package changchun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.StringTokenizer;
import java.util.TreeMap;

/**
 * 
 * Get coordinates of each gene and store in BED file - for later FASTA extraction
 * 
 * @author Johan Henriksson
 *
 */
public class GetGeneCoord {

	public static void main(String[] args) throws IOException {
		PrintWriter pw=new PrintWriter(new File("/home/mahogny/Desktop/celegans/all.bed"));
		
		BufferedReader br=new BufferedReader(new FileReader(new File("/home/mahogny/Desktop/celegans/all.gff")));
		
		String line=null;
		while((line=br.readLine())!=null) {
			
			if(!line.startsWith("#")) {
				
				StringTokenizer stok=new StringTokenizer(line,"\t");
				
				String chrom=stok.nextToken();
				stok.nextToken();
				String t=stok.nextToken();
				
				if(t.contentEquals("gene")) {
					
					String sfrom=stok.nextToken();
					String sto=stok.nextToken();
					/*String sdot=*/stok.nextToken();
					String sstrand=stok.nextToken();
					/*String sdot2=*/stok.nextToken();
					String sattr=stok.nextToken();

					int from=Integer.parseInt(sfrom);
					int to=Integer.parseInt(sto);
					//from=Math.max(1,from-150);
					//to=to+150; //pray...
					
					TreeMap<String,String> attr=parseAttr(sattr);
					
					String gene_id=attr.get("gene_id");
					pw.println(
							chrom+"\t"+from+"\t"+to+"\t"+
									gene_id+"\t"+
									"0\t"+sstrand);
					//https://en.wikipedia.org/wiki/BED_(file_format)
					
					//	System.out.println(line);
				}
				
				
			}
			
			
		}
		
		
		br.close();
		pw.close();
	}
	
	

	/**
	 * Parse an attribute string into a map
	 */
	private static TreeMap<String,String> parseAttr(String s){
		TreeMap<String,String> m=new TreeMap<String, String>();
		StringTokenizer stok=new StringTokenizer(s,";");
		while(stok.hasMoreTokens()) {
			String tok=stok.nextToken();
			tok=tok.trim();
			int ind=tok.indexOf("=");
			if(ind==-1) {
				//KEY value
				ind=tok.indexOf(" ");
				if(ind==-1)
					throw new RuntimeException("Unknown format of attributes; example: "+s);
				String key=tok.substring(0,ind);
				String value=tok.substring(ind+1);
				value=parseOutCitation(value);
				m.put(key, value);
			} else {
				//KEY=value
				String key=tok.substring(0,ind);
				//key=key.substring(0,key.length()-1);
				String value=tok.substring(ind+1);
				//value=value.substring(0,value.length()-1);
				m.put(key, value);
			}
		}
		return m;
	}
	
	/**
	 * Remove citation signs around "foo" if detected
	 */
	private static String parseOutCitation(String value) {
		if(value.startsWith("\""))
			value=value.substring(1,value.length()-1);
		return value;
	}

	
}


/**
 * 
 *map transcript_id -> gene_id

take all lines "gene"

extract geneseq +-300bp
	getfasta

replace region

new gff file - used for feature counting


map the reads
remove multimapping reads

perform feature counting





 
 * 
 * 
*/