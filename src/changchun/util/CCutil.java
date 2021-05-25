package changchun.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.HashMap;
import java.util.TreeMap;


/**
 * Common tools
 * 
 * @author Johan Henriksson
 *
 */
public class CCutil {

	
	public static HashMap<String,String> mapGeneOrigfastaNoExtra=new HashMap<String, String>();
	public static TreeMap<String, String> mapGeneOrigfastaWithExtra=new TreeMap<String, String>();
	public static TreeMap<String, String> mapGeneChrom=new TreeMap<String, String>();
	public static TreeMap<String, Integer> mapGeneStart=new TreeMap<String, Integer>();

	
	/**
	 * Get the complement
	 */
	public static char comp(char c) {
		if(c=='A')
			return 'T';
		else if(c=='T')
			return 'A';
		else if(c=='C')
			return 'G';
		else if(c=='G')
			return 'C';
		else if(c=='N')
			return 'N';
		else {
            System.out.println("base: "+c);
            return 'N';
		}
	
	}

	/**
	 * Get the reverse complement
	 */
	public static String revcomp(String s) {
		char[] c=new char[s.length()];
		
		for(int i=0;i<s.length();i++) {
			c[c.length-i-1] = comp(s.charAt(i));
		}
		return new String(c);
	}

	
	
	/**
	 * Read FASTA with the sequence of each gene
	 */
	public static void readGeneFastaNoExtra() throws IOException {
		BufferedReader br=new BufferedReader(new FileReader(new File("genes.fa")));
		
		String line;
		while((line=br.readLine())!=null) {
			String wbid=line.substring(1);
			wbid=wbid.substring(0,wbid.indexOf(":"));
			line=br.readLine();
			mapGeneOrigfastaNoExtra.put(wbid,line);
		}
		
		br.close();
	}

	
	/**
	 * Get FASTA seq for genes that we have previously cut out
	 */
	public static void readGeneFastaWithExtra() throws IOException {
		//BufferedReader br=new BufferedReader(new FileReader(new File("genes.fa")));
		BufferedReader br=new BufferedReader(new FileReader(new File("all.bed.fa")));  //cut out a bit extra here... 300bp?
				
	
		StringBuilder sb=new StringBuilder();
		String thisSeq="";
		String line=null;
		while((line=br.readLine())!=null) {
			if(line.startsWith(">")) {
				if(sb.length()>0) {
					mapGeneOrigfastaWithExtra.put(thisSeq,sb.toString());
				}
				
				//>WBGene00199986::X:15627747-15628189
				thisSeq=line.substring(1);
				thisSeq=thisSeq.substring(0,thisSeq.indexOf(":"));
				sb=new StringBuilder();
				
				//Coord: 
				String coord=line.substring(line.indexOf(":")+2);
				
				String chr=coord.substring(0,coord.indexOf(':'));
				CCutil.mapGeneChrom.put(thisSeq,chr);
				
				String pos=coord.substring(coord.indexOf(':')+1,coord.indexOf('-'));
				CCutil.mapGeneStart.put(thisSeq,Integer.parseInt(pos));
				
			} else {
				sb.append(line);
			}
		}
		
		if(sb.length()>0) {
			mapGeneOrigfastaWithExtra.put(thisSeq,sb.toString());
		}
		
		br.close();
	}

	
	/**
	 * Get map symbol -> wbid
	 */
	public static TreeMap<String, String> getNameMap() throws IOException {
		BufferedReader br=new BufferedReader(new FileReader(new File("map_wbid_transcid.csv")));
		
		TreeMap<String, String> map=new TreeMap<String, String>();
		
		br.readLine(); //header
		String line=null;
		while((line=br.readLine())!=null) {
			String[] col=line.split(",", 0);
			
			String wbid=col[0];
			
			String tr_name=col[1];
			String gene_name=col[2];
			
			//System.out.println(gene_name);
			
			map.put(tr_name.toUpperCase(),wbid);
			map.put(gene_name.toUpperCase(),wbid);
	
		}
		
		
		br.close();
		return map;
	}


	public static Object readObjectFile(File fSerial) throws IOException {
	    Object allEdits;
		try {
			ObjectInputStream ois = new ObjectInputStream(new FileInputStream(fSerial));
			allEdits = ois.readObject();
			ois.close();
		} catch (ClassNotFoundException e) {
			throw new IOException(e);
		}
	    return allEdits;
	}

	public static void writeObjectFile(File fSerialOut, Object allEdits) throws IOException {
	    FileOutputStream fileOutputStream = new FileOutputStream(fSerialOut);
	    ObjectOutputStream objectOutputStream = new ObjectOutputStream(fileOutputStream);
	    objectOutputStream.writeObject(allEdits);
	    objectOutputStream.flush();
	    objectOutputStream.close();
	}

	/**
	 * Read what genes should be edited to. wbid -> sequence
	 */
	public static HashMap<String, String> readMapWbidEditedseq() throws IOException {
		HashMap<String,String> attemptedGeneEdits=new HashMap<String, String>();
		TreeMap<String, String> geneMap=getNameMap();
	
		@SuppressWarnings("resource")
		BufferedReader br=new BufferedReader(new FileReader(new File("gene_edits.csv")));
				
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


}
