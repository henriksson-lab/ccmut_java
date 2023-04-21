package changchun;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.TreeMap;

import changchun.util.CCutil;
import changchun.util.MatchSeq;


/**
 * 
 * Write a BED file that point to regions that have been gene edited. This is precisely the edit position, not the entire gene
 * 
 * 
 * TODO need a new mapping wbid - symbol to run this
 * 
 * @author Johan Henriksson
 *
 */
public class WriteGeneEditPositions {

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
	public static boolean attemptFitGeneEdit(
			PrintWriter pwListGene, PrintWriter pwBed,String fasta, String newseq, String wbid) {
	
		MatchSeq m=new MatchSeq();
		if(m.match(fasta, newseq, 30)) {
			
			pwListGene.println(wbid);
			
			//Bed region for original fasta
			int gs=CCutil.mapGeneStart.get(wbid); //fair chance of some off-by-1 error here
			pwBed.println(
				CCutil.mapGeneChrom.get(wbid)+"\t"+
				(gs+m.fasta_i1+1)+"\t"+
				(gs+m.fasta_i2)+"\t"+
				wbid);  //store location of gene!
			
			return true;
		} else
			return false;
		
	}
	
	/**
	 * Entry point
	 */
	public static void main(String[] args) throws IOException {
		TreeMap<String, String> geneMap=CCutil.getNameMap();
		//uses map_wbid_transcid.csv   contains zzz-1! with WBGene00009333
		
		
		CCutil.readGeneFastaWithExtra();
		BufferedReader br=new BufferedReader(new FileReader(new File("/home/mahogny/Desktop/celegans/gene_edits.csv")));
		//now: /home/mahogny/Dropbox/umea/project/changchun/gpcr/mutantanalysis/gene_edits.csv
		//has: 1692 (and a header)
		//this is a list of genesym, from, to, sequences
		
		PrintWriter pwListGene=new PrintWriter(new File("/home/mahogny/Desktop/celegans/mutant.csv"));
		//now: /home/mahogny/oldhome/Desktop/alldesktop/celegans/mutant.csv
		//has: 1690
		//this is a list of WBIDs ... why is 1 id missing?
		
		
		PrintWriter pwBed=new PrintWriter(new File("/home/mahogny/Desktop/celegans/edit_positions.bed"));
		//now: /home/mahogny/oldhome/Desktop/alldesktop/celegans/edit_positions.bed
		//has: 1441
		
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
				if(!attemptFitGeneEdit(pwListGene, pwBed, fasta, newseq, wbid)) {
					if(!attemptFitGeneEdit(pwListGene, pwBed, fasta, CCutil.revcomp(newseq), wbid)) {
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
				
				break;//throw new IOException("No wbid!");
				
				//TODO these must all be handled!
				
				//no such: WBGene00305048 . in neither all.fa or all.gff
				
			}
			
			
		}
		
		
		br.close();
		
		pwListGene.close();
		
	}
	
	
	
}
