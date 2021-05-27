package changchun;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;

import changchun.util.CCutil;


/**
 * Take the alignment output. Format it for viewing
 * 
 * @author Johan Henriksson
 *
 */
public class FormatAlignment {
	
	public static Connection conn;
	

	
	/**
	 * Connect to the database
	 */
	private static Connection connect(String db) throws SQLException {
		// SQLite connection string
		String url = "jdbc:sqlite:"+db;
		Connection conn = DriverManager.getConnection(url);
		return conn;
	}

	/**
	 * Insert an alignment into the SQL table
	 */
	public static void insertSQL(Connection conn, String strain, OneEdit edit) throws IOException {
		try {
			String sql = "INSERT INTO alignment(strain, wbid, alignment, countread, countins, countdel,   fracins, fracdel, fracfine) VALUES(?,?,?,?,?,?,?,?,?)";

			PreparedStatement pstmt = conn.prepareStatement(sql);

			pstmt.setString(1, strain);
			pstmt.setString(2, edit.geneid);
			pstmt.setString(3, edit.alignment);
			
			pstmt.setInt(4, edit.numread);
			pstmt.setInt(5, edit.sumins);
			pstmt.setInt(6, edit.sumdel);
			
			double numread=edit.numread;
			if(numread==0)
				numread=1;
			pstmt.setDouble(7, edit.sumins/numread);
			pstmt.setDouble(8, edit.sumdel/numread);
			pstmt.setDouble(9, edit.fine/numread);
			
			pstmt.executeUpdate();
		} catch (SQLException ex) {
			throw new IOException(ex);
		}
	}
	    
	
	public static void createTable(Connection conn) throws SQLException {
        String sql_create_rawreads_table = "CREATE TABLE alignment ("+
        		"strain    text NOT NULL,"+
                "wbid      text NOT NULL,"+
                "alignment text NOT NULL,"+
                "countread number NOT NULL,"+
                "countins  number NOT NULL,"+
                "countdel  number NOT NULL,"+
                "fracins   number NOT NULL,"+
                "fracdel   number NOT NULL,"+
                "fracfine   number NOT NULL"+
                ");";

        Statement stmt = conn.createStatement();
        stmt.execute(sql_create_rawreads_table);
        
	}
	
	
	public static void processFile(
			File fSerial, String strain) throws IOException {

		@SuppressWarnings("unchecked")
		ArrayList<OneEdit> allEdits=(ArrayList<OneEdit>)CCutil.readObjectFile(fSerial);
	
		int iii=0;
	    for(OneEdit e:allEdits) {
	    	iii++;
	    	if(!e.alignment.contentEquals("")) {
	    		//System.out.println(e.geneid);
	    		if(iii%100==0)
	    			System.out.println("   "+iii);
	    		
		    	BufferedReader br=new BufferedReader(new StringReader(e.alignment));
		    	
		    	//Remove the first lines
		    	for(int i=0;i<3;i++)
		    		br.readLine();
		    	
		    	//Read the first block of alignments
		    	ArrayList<StringBuilder> lines=new ArrayList<StringBuilder>(100);
		    	for(;;) {
		    		String line=br.readLine();
		    		if(line==null || line.contentEquals("") || line.startsWith(" "))
		    			break;
		    		StringBuilder sb=new StringBuilder();
		    		sb.append(line.substring(8));
		    		lines.add(sb);
		    	}
		    	br.readLine();
		    	
		    	//Read the other blocks
		    	for(;;) {
		    		String line=br.readLine();
		    		if(line==null) {
		    			break;
		    		} else {
		    			for(int i=0;i<lines.size();i++) {
				    		lines.get(i).append(line.substring(8));
				    		line=br.readLine();
		    			}
		    			//Skip one line
			    		line=br.readLine();
		    		}
		    	}
		    	
		    	//Concatenate the new alignment
		    	StringBuilder sb=new StringBuilder();
				for(int i=0;i<lines.size();i++) {
					sb.append(lines.get(i));
					sb.append("\n");
				}
		    	
				e.alignment=sb.toString();
				
		    	br.close();
		    	
    			insertSQL(conn, strain, e);
	    	} else {
	    		e.alignment="";
	    	}
	    }
	}
	
	
	
	/**
	 * 
	 * Entry point
	 */
	public static void main(String[] args) throws IOException, ClassNotFoundException, SQLException {
		//HashMap<String,String> attemptedGeneEdits=CCutil.readMapWbidEditedseq();
		//CCutil.readGeneFastaWithExtra();
		
		/*
		String name="P19764_154_S54_L001_R1_001.fastq.gz.out.bam.mut.serial";
		name=name.substring(name.indexOf('S')+1);
		name=name.substring(0,name.indexOf('_'));
		int id=Integer.parseInt(name);
		System.out.println("CHS "+(1000+id));
		System.exit(0);*/
				
		File rootdir=new File("/media/mahogny/TOSHIBA/changchun/mut");
		if(args.length!=0)
			rootdir=new File(args[0]);
		
		File fSQL=new File(rootdir, "alignments.sqlite");
		if(fSQL.exists())
			fSQL.delete();
		
		conn=connect(fSQL.getAbsolutePath());
		createTable(conn);

		//101,CHS 1001
		for(File fSerial:rootdir.listFiles()) {
			if(fSerial.getName().endsWith(".serial")) {
				System.out.println(fSerial);
				String name=fSerial.getName();
				name=name.substring(name.indexOf('S')+1);
				name=name.substring(0,name.indexOf('_'));
				int id=Integer.parseInt(name);
				String chsid="CHS "+(1000+id);
				processFile(fSerial, chsid);

			}
		}
		/*
		for(int i=0;i<280;i++) {
			System.out.println(i);
			int ngi_id=101+i;
			int ccid=1001+i;
			
			File fSerial=new File(rootdir, "P19764_"+ngi_id+"_S"+(i+1)+"_L001_R1_001.fastq.gz.out.bam.mut.serial");		
			if(fSerial.exists()) {
				processFile(fSerial, "CHS "+ccid);
			} else {
				System.out.println("Missing file "+fSerial);
			}
		}*/
		
	    
		
	}
}
