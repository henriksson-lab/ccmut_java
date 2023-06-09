package changchun.crap;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.TreeMap;

import changchun.util.CCutil;
import changchun.util.MatchSeq;

public class ReadGeneEdit {

	/**
	 * Get map symbol -> wbid
	 */
	/*
	public static TreeMap<String, String> getNameMap() throws IOException {
		BufferedReader br=new BufferedReader(new FileReader(new File("/home/mahogny/Desktop/celegans/map_wbid_transcid.csv")));
		
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
	*/
	
	/*
	public static TreeMap<String, String> geneChrom=new TreeMap<String, String>();
	public static TreeMap<String, Integer> geneStart=new TreeMap<String, Integer>();

	public static TreeMap<String, String> getCutFasta() throws IOException {
		BufferedReader br=new BufferedReader(new FileReader(new File("/home/mahogny/Desktop/celegans/all.bed.fa")));
				
		TreeMap<String, String> map=new TreeMap<String, String>();

		StringBuilder sb=new StringBuilder();
		String thisSeq="";
		String line=null;
		while((line=br.readLine())!=null) {
			if(line.startsWith(">")) {
				if(sb.length()>0) {
					map.put(thisSeq,sb.toString());
				}
				thisSeq=line.substring(1); //not quite
				thisSeq=thisSeq.substring(0,thisSeq.indexOf(":"));
				//>WBGene00199986::X:15627747-15628189
				sb=new StringBuilder();
				
				//Coord: 
				String coord=line.substring(line.indexOf(":")+2);
				
				String chr=coord.substring(0,coord.indexOf(':'));
				geneChrom.put(thisSeq,chr);
				
				String pos=coord.substring(coord.indexOf(':')+1,coord.indexOf('-'));
				geneStart.put(thisSeq,Integer.parseInt(pos));
				
				
				//V:30-479
				
			} else {
				sb.append(line);
			}
		}
		br.close();
		return map;
	}
	
	*/

	
	
	/////////////////// TODO: calculate substitute regions
	/*
	
	public static int moveUntilEndOfMatch(String fasta, String seq, int offset) {
		int i=0;
		while(seq.charAt(i)==fasta.charAt(offset+i))
			i++;
		return i;
	}
	
	public static int moveUntilStartOfMatch(String fasta, String seq, int offset) {
		int i=0;
		while(seq.charAt(i)==fasta.charAt(offset-i))
			i++;
		return i;
	}
	
	*/
	
	
	/**
	 * 
	 * @param pwGTF
	 * @param pwFasta
	 * @param pwListGene
	 * @param fasta
	 * @param newseq
	 * @param wbid
	 * @return
	 */
	public static boolean out(
			PrintWriter pwGTF, PrintWriter pwFasta, PrintWriter pwListGene, PrintWriter pwBed,
			
			String fasta, String newseq, String wbid) {
	
		MatchSeq m=new MatchSeq();
		if(m.match(fasta, newseq, 30)) {
			String newid="rep_"+wbid;
			pwFasta.println(">"+newid);
			pwFasta.println(m.fastaSub);
			
			//V	WormBase	gene	180	329	.	+	.	gene_id "WBGene00197333"; gene_name "cTel3X.2"; gene_source "WormBase"; gene_biotype "ncRNA";
			pwGTF.println("rep_"+wbid+"\t"+
					"mutant\t"+"gene\t"+
					"1\t"+m.fastaSub.length()+"\t"+
					".\t"+"+\t"+".\t"+
					"gene_id \""+newid+"\";");

			pwGTF.println("rep_"+wbid+"\t"+
					"mutant\t"+"exon\t"+
					"1\t"+m.fastaSub.length()+"\t"+
					".\t"+"+\t"+".\t"+
					"gene_id \""+newid+"\";");

			
			//V       WormBase        transcript      180     329     .       +       .       gene_id "WBGene00197333"; transcript_id "cTel3X.2"; gene_name "cTel3X.2"; gene_source "WormBase"; gene_biotype "ncRNA"; transcript_source "WormBase"; transcript_biotype "ncRNA";
			//V       WormBase        exon    180     329     .       +       .       gene_id "WBGene00197333"; transcript_id "cTel3X.2"; exon_number "1"; gene_name "cTel3X.2"; gene_source "WormBase"; gene_biotype "ncRNA"; transcript_source "WormBase"; transcript_biotype "ncRNA"; exon_id "cTel3X.2.e1";

			
			
			pwListGene.println(wbid);
			
			//Bed region for original fasta
			int gs=CCutil.mapGeneStart.get(wbid); //fair chance of some off-by-1 error here
			pwBed.println(CCutil.mapGeneChrom.get(wbid)+"\t"+(gs+m.fasta_i1+1)+"\t"+(gs+m.fasta_i2)+"\torig_"+wbid);  //store location of gene!
			
			//Bed region for substitute
			int dx=m.newseq_i2-m.newseq_i1;
			pwBed.println(newid+"\t"+(m.fasta_i1+1)+"\t"+(m.fasta_i1+1+dx)+"\t"+newid);

			
			///wrong!!!
			
			return true;
		} else
			return false;
		
	}
	
	/**
	 * Entry point
	 */
	public static void main(String[] args) throws IOException {
		TreeMap<String, String> geneMap=CCutil.getNameMap();
		CCutil.readGeneFastaWithExtra();
		BufferedReader br=new BufferedReader(new FileReader(new File("/home/mahogny/Desktop/celegans/gene_edits.csv")));
		
		PrintWriter pwGTF=new PrintWriter(new File("/home/mahogny/Desktop/celegans/mutant.gtf"));
		PrintWriter pwFasta=new PrintWriter(new File("/home/mahogny/Desktop/celegans/mutant.fa"));
		PrintWriter pwListGene=new PrintWriter(new File("/home/mahogny/Desktop/celegans/mutant.csv"));
		PrintWriter pwBed=new PrintWriter(new File("/home/mahogny/Desktop/celegans/mutant.bed"));
		
		br.readLine();
		String line=null;
		while((line=br.readLine())!=null) {
			String[] col=line.split(",", 0);
			
			String geneName=col[0];
			String wbid=geneMap.get(geneName.toUpperCase());
			
			if(wbid==null) {
				System.out.println(geneName);
				System.exit(0);
			}
			
			String newseq=col[3];
			
			
			newseq=newseq.replaceAll(" ", "");
			newseq=newseq.replaceAll("\\.", "");
			newseq=newseq.toUpperCase();
			
			//System.out.println(wbid+","+newseq);
			String fasta=CCutil.mapGeneOrigfastaWithExtra.get(wbid);
			if(fasta!=null) {
				//Replace the content. Return new seq
				if(!out(pwGTF, pwFasta, pwListGene, pwBed, fasta, newseq, wbid)) {
					if(!out(pwGTF, pwFasta, pwListGene, pwBed, fasta, CCutil.revcomp(newseq), wbid)) {
						System.out.println("failing with "+wbid);
						System.out.println(fasta);
						System.out.println(newseq);
						System.out.println(CCutil.revcomp(newseq));
						System.out.println("===");
						
						//failing: WBGene00012181   --- actually does not have the correct end, first newseq first part matches
					}
				}
			} else {
				System.out.println("No such wbid "+wbid);
				
				
				//no such: WBGene00305048 . in neither all.fa or all.gff
				
			}
			
			
		}
		
		
		br.close();
		
		pwListGene.close();
		pwGTF.close();
		pwFasta.close();
		
	}
	
	
	
}
