package changchun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.io.StringReader;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import changchun.util.CCutil;

/**
 * 
 * look for: TAAGGATCCTAA
 * 
 * 
 * 
 * @author mahogny
 *
 */
public class ExtractSeqFromMSAOld {

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
	
	public static class ReplacedSeqMSA {
		
		String from, to;
		int count;
		private int posOrig, posIntend;
		
		
		public ReplacedSeqMSA(String from, String to, int count, int posOrig, int posIntend) {
			this.from = from;
			this.to = to;
			this.count = count;
			this.posOrig=posOrig;
			this.posIntend=posIntend;
		}
		

		public String toString() {
			return from+"/"+to+"@"+posOrig+"@"+posIntend+"#"+count;
		}
	}
	
	
	public static int countN(String s) {
		int cnt=0;
		for(char c:s.toCharArray())
			if(c=='N')
				cnt++;
		return cnt;
	}
	
	
	public static List<String> readAlignment(String alignment) throws IOException{
    	BufferedReader br=new BufferedReader(new StringReader(alignment));
    	
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
    	
    	//////// After this, the alignment is a list of lines!
    	
    	//Turn to string
    	String[] line=new String[lines.size()];
		for(int i=0;i<lines.size();i++) {
			line[i]=lines.get(i).toString().replaceAll(" ", "");
			//System.out.println(line[i]);
		}
		return Arrays.asList(line);
	}
	
	public static ArrayList<OneEdit> processFile(
			File fSerial, String strain, PrintWriter pw) throws IOException {

		@SuppressWarnings("unchecked")
		ArrayList<OneEdit> allEdits=(ArrayList<OneEdit>)CCutil.readObjectFile(fSerial);

		ArrayList<OneEdit> updatedEdits=new ArrayList<>();
		
		int iii=0;
	    for(OneEdit e:allEdits) {
	    	iii++;
	    	if(!e.alignment.contentEquals("")) {
	    		//System.out.println(e.geneid);
	    		if(iii%100==0)
	    			System.out.println("   "+iii);

	    		//System.out.println("################################# one alignment");

	    		List<String> arrayLines=readAlignment(e.alignment);

	    		
	    		int perfectFitOriginal=0;
	    		int perfectFitEdit=0;
	    		
	    		//Prepare modified edit
	    		OneEdit newOneEdit=new OneEdit();
	    		newOneEdit.geneid=e.geneid;
	    		newOneEdit.alignment="";

				//Algorithm to split up MSA
				ArrayList<ReplacedSeqMSA> replacements=new ArrayList<>();
				int numCroppedToFit=0;
		    	int numNotConsensus=0;
		    	int numFilteredPoorFit=0;
				if(arrayLines.size()!=0) {
					
					//////// Remove reads that have too many errors. These seem to be in the wrong place
					
					//List<String> arrayLines = Arrays.asList(line);
					
			    	numFilteredPoorFit=arrayLines.size();
			    	ArrayList<String> cleanedReads=removePoorFit(arrayLines);
					numFilteredPoorFit-=cleanedReads.size();
					
					//System.out.println("#Poor fit: "+numFilteredPoorFit);
			    	
					//////// Now focus on the region around the expected edit
			    	//System.out.println("num reads before crop: "+lines.size());
			    	ArrayList<String> cropped=cropToIntendedSequence(cleanedReads);

					//Store new alignment
			    	StringBuilder sbNewEdit=new StringBuilder();
			    	for(String s:cropped) {
			    		if(sbNewEdit.length()>0)
			    			sbNewEdit.append("\n");
			    		sbNewEdit.append(s);			    		
			    	}
					newOneEdit.alignment=sbNewEdit.toString();
					updatedEdits.add(newOneEdit);
					
			    	//Last two reads are special so separate them
			    	ArrayList<String> fitseqs=new ArrayList<String>(cropped);
			    	String intendSeq=fitseqs.remove(fitseqs.size()-1);
			    	String origSeq = fitseqs.remove(fitseqs.size()-1);
			    	numCroppedToFit=fitseqs.size();

			    	
			    	//Is there a read that perfect fits the original or edited sequence, over the expected edit region?
			    	//Where do diffs start?
			    	int firstDiff=0;
		    		//System.out.println(origSeq);
		    		//System.out.println(intendSeq);
			    	while(true) {
			    		char charOrig=origSeq.charAt(firstDiff);
			    		char charIntend=intendSeq.charAt(firstDiff);
			    		if(charOrig==charIntend) {
			    			firstDiff++;
			    		} else
			    			break;
			    	}
			    	int lastDiff=firstDiff+1;
			    	for(int i=lastDiff;i<origSeq.length();i++) {
			    		char charOrig=origSeq.charAt(i);
			    		char charIntend=intendSeq.charAt(i);
			    		if(charOrig!=charIntend) {
			    			lastDiff=i;
			    		} 
			    	}
			    	String subOrig=origSeq.substring(firstDiff, lastDiff);
			    	String subIntend=intendSeq.substring(firstDiff, lastDiff);
			    	for(String s:fitseqs) {
			    		if(s.substring(firstDiff, lastDiff).equals(subOrig)) {
			    			//System.out.println(s+"\tto orig");
			    			perfectFitOriginal++;
			    		}
			    	}
			    	for(String s:fitseqs) {
			    		if(s.substring(firstDiff, lastDiff).equals(subIntend)) {
			    			//System.out.println(s+"\tto intend");
			    			perfectFitEdit++;			    			
			    		}
			    	}
			    	
			    	//System.out.println("Diffs are "+firstDiff+"  "+lastDiff);
			    	//System.out.println(""+perfectFitOriginal+"\t"+perfectFitEdit);
			    	
		    		//boolean perfectFitOriginal=false;
		    		//boolean perfectFitEdit=false;

					
			    	//Find consensus sequence. Then extract reads that agree with it. 
			    	//System.out.println("## Find all consensus sequences");
			    	/*
			    	for(String s:fitseqs)
			    		System.out.println(s);*/
			    	
					while(fitseqs.size()>0) {
						ArrayList<String> notfittingseqs=new ArrayList<>();
				    	int numFittingToConsensusOneRound=fitToConcensus(fitseqs, origSeq, intendSeq, notfittingseqs, replacements);
				    	fitseqs=notfittingseqs;
				    	numNotConsensus=fitseqs.size();
				    	if(numFittingToConsensusOneRound==0) {
				    		//Cannot fit anymore; give up
				    		//System.out.println("Left: "+fitseqs.size());
				    		break;
				    	}
					}
					

				}
				
				//Print out results
				pw.print(strain+"\t"+e.geneid+"\t"+numCroppedToFit +"\t"+numNotConsensus+"\t"+numFilteredPoorFit+"\t"+perfectFitOriginal+"\t"+perfectFitEdit+"\t");
				boolean first=true;
				for(ReplacedSeqMSA s:replacements) {
					
					if(s.count<2)
						continue;
						
					if(countN(s.to)>3)
						continue;
					
					if(s.from.length()<=1 && s.to.length()<=1)
						continue;
					
					if(!first)
						pw.print(",");
					pw.print(s);
					first=false;
				}
				pw.println();

	    	} else {
	    		e.alignment="";
	    	}
	    }
	    return updatedEdits;
	}
	
	

	/**
	 * Remove reads that have too many errors. These seem to be in the wrong place
	 */
	private static ArrayList<String> removePoorFit(List<String> line) {
		line=new ArrayList<>(line);
		
    	//Last two reads are special so separate them
    	String intendSeq = line.remove(line.size()-1);
    	String origSeq = line.remove(line.size()-1);

    	ArrayList<String> filtered=new ArrayList<>();
    	for(String sRead:line) {
        	if(numdiffToOriginalOrIntend(sRead, origSeq, intendSeq)<10) {
        		filtered.add(sRead);
        	}
    	}

    	//Readd special seqs
    	filtered.add(origSeq);
    	filtered.add(intendSeq);
    	
    	return filtered;
	}


	/**
	 * Extract sequences that fit to consensus, remove these, and report what is the replaced sequence
	 */
	private static int fitToConcensus(ArrayList<String> fitseqs, String origSeq, String intendSeq, ArrayList<String> remainingSeqs, ArrayList<ReplacedSeqMSA> replacements) {
		int len=fitseqs.get(0).length();
		int countA[]=new int[len];
		int countT[]=new int[len];
		int countG[]=new int[len];
		int countC[]=new int[len];
		
		//Count occurrences
		for(String s:fitseqs) {
			for(int i=0;i<s.length();i++) {
				char c=s.charAt(i);
				if(c=='A')
					countA[i]++;
				else if(c=='T')
					countT[i]++;
				else if(c=='G')
					countG[i]++;
				else if(c=='C')
					countC[i]++;
			}
		}
		
		//Create consensus. N if no seq present
		char[] consensus=new char[len];
		for(int i=0;i<len;i++) {
			char c='N';
			int tot=countA[i]+countT[i]+countG[i]+countC[i];
			if(tot!=0) {
				int max=countA[i];
				c='A';
				if(countT[i]>max) {
					max=countT[i];
					c='T';
				}
				if(countG[i]>max) {
					max=countG[i];
					c='G';
				}
				if(countC[i]>max) {
					max=countC[i];
					c='C';
				}
			}
			consensus[i]=c;
		}
		
		//Check which seqs agree with consensus
		//ArrayList<String> notfittingseqs=new ArrayList<>();
		int numfitting=0;
		for(String s:fitseqs) {
			int numdiff=numdiffToConsensus(s, consensus);
			
			if(numdiff<3)
				numfitting++;
			else
				remainingSeqs.add(s);
		}
		
		if(numfitting>0) {
			//Figure out what this replaced sequence is

			/*
			System.out.println("consensus");
	    	System.out.println(new String(consensus));
	    	System.out.println(origSeq);
	    	*/
			

	    	int firstIndex=0;
	    	while(firstIndex<origSeq.length() && (origSeq.charAt(firstIndex)==consensus[firstIndex] || consensus[firstIndex]=='N'))
	    		firstIndex++;
	    	
	    	if(firstIndex==origSeq.length()) {
	    		//System.out.println("Agrees with original, #fitting: " +numfitting+"\t#notfitting: "+remainingSeqs.size());
	    		//But this is a complicated case actually. more tests needed
	    	} else {
		    	int lastIndex=0;
		    	for(int i=firstIndex;i<origSeq.length();i++)
		    		if(origSeq.charAt(i)!=consensus[i] && consensus[i]!='N')
		    			lastIndex=i+1;  

		    	//Figure out substitution. remove gaps
		    	String origSub=origSeq.substring(firstIndex, lastIndex).replaceAll("-", "");
		    	String newSub=new String(consensus).substring(firstIndex, lastIndex).replaceAll("-", "");
						    	
		    	//To find position, remove gaps up until this place
		    	int posOrig=origSeq.substring(0,firstIndex).replaceAll("-", "").length();
		    	int posIntend=intendSeq.substring(0,firstIndex).replaceAll("-", "").length();
		    	
		    	//System.out.println(pos);
		    	/*
		    	if(posOrig>1000) {
		    		
		    		System.out.println("pos");
		    		System.out.println(origSeq);
		    		System.out.println(intendSeq);
		    		
		    		System.exit(0);
		    	}*/
		    	
		    	//System.out.println(new String(consensus));
		    	//System.out.println("======================================================> "+origSub+"/"+newSub+"#"+numfitting);
		    	replacements.add(new ReplacedSeqMSA(origSub,newSub,numfitting,posOrig,posIntend));
		    	//System.exit(0);

	    	}
		} else {
			//System.out.println("None fitting");
			
			/*
    		System.out.println("############");
    		for(String s:remainingSeqs)
    			System.out.println(s);
    		System.out.println("////////////");
    		System.out.println(new String(consensus));
    		*/
		}
		
		return numfitting;
	}

	
	/**
	 * Get the number of differences for a read vs consensus sequence. 
	 * If the consensus has an N-base then ignore this position. Though this might never happen?
	 */
	private static int numdiffToConsensus(String s, char[] consensus) {
		int numdiff=0;
		for(int i=0;i<consensus.length;i++) {
			char c=consensus[i];
			char r=s.charAt(i);
			if(c!='N' && c!=r && r!='-') {
				numdiff++;
			}
		}
		return numdiff;
	}


	/**
	 * Get the number of differences between read vs original sequence. Ignore whenever there is a gap in either sequence
	 */
	private static int numdiffToOriginalOrIntend(String sRead, String sOrig, String sIntend) {
		int numdiff=0;
		for(int i=0;i<sRead.length();i++) {
			char r=sRead.charAt(i);
			if(r!='-') {
				char o=sOrig.charAt(i);
				char intend=sIntend.charAt(i);
				if(r==o || r==intend) {
					//Fits orig or intended. good!
				} else if(o=='-' && intend=='-') {
					//orig and intended both have gaps. ignore this position
				} else
					numdiff++;
			}
		}
		return numdiff;
	}


	// update future faculty site regarding sequencing
	
	/**
	 * Crop the alignment to only cover region of the intended sequence
	 */
	private static ArrayList<String> cropToIntendedSequence(List<String> line) {
    	//Crop to last line (to replace)
    	String expectedEdit = line.get(line.size()-1);
    	//System.out.println("e "+expectedEdit);
    	
    	int firstIndex=0;
    	while(expectedEdit.charAt(firstIndex)=='-')
    		firstIndex++;

    	int lastIndex=0;
    	for(int i=firstIndex;i<expectedEdit.length();i++)
    		if(expectedEdit.charAt(i)!='-')
    			lastIndex=i+1;

    	//System.out.println(firstIndex);
    	//System.out.println(lastIndex);
    	//System.out.println(expectedEdit.substring(firstIndex, lastIndex));
    	
    	//String[] cropped=new String[line.length];
    	ArrayList<String> cropped=new ArrayList<>();
    	for(String oneline:line) {
    	//for(int i=0;i<line.length;i++) {
    		String sub=oneline.substring(firstIndex, lastIndex).replaceAll(" ",""); //Spaces added into alignment. removed here
    		if(!isBlank(sub)) {
    			cropped.add(sub);
    		}
    		//System.out.println(cropped[i]);
    	}
    	
    	/*
    	for(String s:cropped)
    		System.out.println(s);*/

    	return cropped;
	}


	public static boolean isBlank(String s) {
		for(int i=0;i<s.length();i++)
			if(s.charAt(i)!='-')
				return false;
		return true;
	}

	

	public static void createCroppedTable(Connection conn) throws SQLException {
        String sql_create_rawreads_table = "CREATE TABLE cropped_alignment ("+
        		"strain    text NOT NULL,"+
                "wbid      text NOT NULL,"+
                "alignment text NOT NULL"+
                ");";

        Statement stmt = conn.createStatement();
        stmt.execute(sql_create_rawreads_table);
	}
	
	
	
	/**
	 * Insert an alignment into the SQL table
	 */
	public static void insertCroppedSQL(Connection conn, String strain, OneEdit edit) throws IOException {
		try {
			String sql = "INSERT INTO cropped_alignment(strain, wbid, alignment) VALUES(?,?,?)";

			PreparedStatement pstmt = conn.prepareStatement(sql);

			pstmt.setString(1, strain.replaceAll("_", " "));
			pstmt.setString(2, edit.geneid);
			pstmt.setString(3, edit.alignment);
			
			pstmt.executeUpdate();
		} catch (SQLException ex) {
			throw new IOException(ex);
		}
	}
	
	/**
	 * 
	 * Entry point
	 */
	public static void main(String[] args) throws IOException, ClassNotFoundException, SQLException {
		//File rootdir=new File("/media/mahogny/TOSHIBA/changchun/mut");
		File rootdir=new File("/home/mahogny/temp/serial");
		File outfile=new File("/home/mahogny/temp/serial/diffs.csv");
		//File fSQL=new File("/corgi/websites/changchun_mutant/website/data/cropped_alignments.sqlite");
		File fSQL=new File("/home/mahogny/temp/serial/cropped_alignments.sqlite");
		
		if(args.length!=0) {
			rootdir=new File(args[0]);
			outfile=new File(args[1]);
			fSQL=new File(args[2]);
		}
		

		if(fSQL.exists())
			fSQL.delete();
		
		conn=connect(fSQL.getAbsolutePath());
		createCroppedTable(conn);
		
		
		PrintWriter pw=new PrintWriter(outfile);
		pw.println("ccid\twbid\tcoveringreads\tunexplainedreads\tfilteredpoorfit\tperfectFitOriginal\tperfectFitEdit\tlistmuts");
		
		//For all files
		for(File fSerial:rootdir.listFiles()) {
			if(fSerial.getName().endsWith(".serial")) {
				System.out.println(fSerial);
				String name=fSerial.getName();
				name=name.substring(0, name.indexOf("_R1"));
				String chsid=name;
				ArrayList<OneEdit> newEdits = processFile(fSerial, chsid, pw);
				
				
				for(OneEdit e:newEdits) {
					insertCroppedSQL(conn, name, e);
				}
				
				/*
				System.out.println("serializing new edits");
				File fNewSerial=new File(fSerial.getParentFile(), fSerial.getName()+".edited");
			    FileOutputStream fileOutputStream = new FileOutputStream(fNewSerial);
			    ObjectOutputStream objectOutputStream = new ObjectOutputStream(fileOutputStream);
			    objectOutputStream.writeObject(newEdits);
			    objectOutputStream.flush();
			    objectOutputStream.close();
			    */


			}
		}
		pw.close();
		System.out.println("done");
	}
}
