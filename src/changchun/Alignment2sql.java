package changchun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.TreeMap;

import changchun.util.CCutil;
import changchun.util.ClustalOmega;
import changchun.util.MatchSeq;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Statement;

/**
 * Take the alignment output. Format it for viewing
 * 
 * @author Johan Henriksson
 *
 */
public class Alignment2sql {
	
	
	
	/**
	 * Connect to the database
	 * @throws SQLException 
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
	public static void insertSQL(Connection conn, String strain, String wbid, String alignment) throws IOException {
		try {
			String sql = "INSERT INTO alignment(strain, wbid,alignment) VALUES(?,?,?)";

			PreparedStatement pstmt = conn.prepareStatement(sql);

			pstmt.setString(1, strain);
			pstmt.setString(2, wbid);
			pstmt.setString(3, alignment);
			pstmt.executeUpdate();
		} catch (SQLException e) {
			throw new IOException(e);
		}
	}
	    
	
	public static void createTable(Connection conn) throws SQLException {
        String sql_create_rawreads_table = "CREATE TABLE alignment ("+
                "wbid text NOT NULL,"+
                "alignment  text NOT NULL);";

        Statement stmt = conn.createStatement();
        stmt.execute(sql_create_rawreads_table);
        
	}
	
	/**
	 * Read what genes should be edited to. wbid -> sequence
	 */
	public static HashMap<String, String> readAttemptedGeneEdits() throws IOException {
		HashMap<String,String> attemptedGeneEdits=new HashMap<String, String>();
		TreeMap<String, String> geneMap=CCutil.getNameMap();

		@SuppressWarnings("resource")
		BufferedReader br=new BufferedReader(new FileReader(new File("/home/mahogny/Desktop/celegans/gene_edits.csv")));
				
		br.readLine();
		String line=null;
		while((line=br.readLine())!=null) {
			String[] col=line.split(",", 0);
			
			String geneName=col[0];
			String wbid=geneMap.get(geneName.toUpperCase());
			
			if(wbid==null) {
				System.out.println(geneName);
				throw new RuntimeException("missing wbid");
			}
			
			String newseq=col[3];
			
			
			newseq=newseq.replaceAll(" ", "");
			newseq=newseq.replaceAll("\\.", "");
			newseq=newseq.toUpperCase();
			
			attemptedGeneEdits.put(wbid, newseq);
		}
		br.close();
		return attemptedGeneEdits;
	}
	
	
	
	public static void main(String[] args) throws IOException, ClassNotFoundException, SQLException {
		
		HashMap<String,String> attemptedGeneEdits=readAttemptedGeneEdits();
		CCutil.readGeneFasta();
		
		String strain="101";
		
		File fSQL=new File("/media/mahogny/TOSHIBA/changchun/mut/P19764_101_S1_L001_R1_001.fastq.gz.out.bam.mut.sqlite");
		if(fSQL.exists())
			fSQL.delete();
		
		Connection conn=connect(fSQL.getAbsolutePath());
		createTable(conn);
		
		File fSerial=new File("/media/mahogny/TOSHIBA/changchun/mut/P19764_101_S1_L001_R1_001.fastq.gz.out.bam.mut.serial");
		//File fSerial=new File("P19764_101_S1_L001_R1_001.fastq.gz.out.bam.mut.serial");
	    ObjectInputStream ois = new ObjectInputStream(new FileInputStream(fSerial));
	    @SuppressWarnings("unchecked")
		ArrayList<OneEdit> allEdits=(ArrayList<OneEdit>)ois.readObject();
	    ois.close();

	    int iii=0;
	    for(OneEdit e:allEdits) {
	    	System.out.println(iii);
	    	iii++;
	    	
	    	if(!e.alignment.contentEquals("")) {
	    		//Add genome fastq
		    	e.reads.add(CCutil.mapGeneOrigfasta.get(e.geneid));
		    	
		    	String fastaSeq=CCutil.mapGeneOrigfasta.get(e.geneid);
		    	String editSeq=attemptedGeneEdits.get(e.geneid);
		    	if(editSeq==null)
		    		throw new RuntimeException("missing "+e.geneid);
				MatchSeq m=new MatchSeq();
				if(!m.match(fastaSeq, editSeq, 30)) {
					m=new MatchSeq();
					editSeq=CCutil.revcomp(e.geneid);
					if(!m.match(fastaSeq, editSeq, 30)) {
						throw new RuntimeException("unable to match "+e.geneid);
					}
				}
		    	
				//Add sequence to substitute
		    	e.reads.add(editSeq);
		    	
		    	//Align
		    	e.alignment=ClustalOmega.call(e.reads);

	    		
	    		System.out.println(e.geneid);
	    		
		    	BufferedReader br=new BufferedReader(new StringReader(e.alignment));
		    	
		    	//Remove the first lines
		    	for(int i=0;i<3;i++)
		    		br.readLine();
		    	
		    	//Read the first block of alignments
		    	ArrayList<String> lines=new ArrayList<String>();
		    	for(;;) {
		    		String line=br.readLine();
		    		if(line.contentEquals("") || line.startsWith(" "))
		    			break;
		    		lines.add(line.substring(8));
		    	}
		    	br.readLine();
		    	
		    	//Read the other blocks
		    	for(;;) {
		    		String line=br.readLine();
		    		if(line==null) {
		    			break;
		    		} else {
		    			for(int i=0;i<lines.size();i++) {
				    		lines.set(i,lines.get(i)+line.substring(8));
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
    	    	
    			insertSQL(conn, strain, e.geneid, sb.toString());
    	    	//String genomeSeqAlign=lines.get(lines.size()-2);
    	    	//String replaceSeqAlign=lines.get(lines.size()-1);
    	    	////////////// Figure out the edit positions!
    	    	
		    	br.close();
		    	
		    	
		    	//System.exit(0);
	    	}

	    	
	    }
	    
	}

}
