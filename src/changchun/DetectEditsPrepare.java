package changchun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.zip.GZIPOutputStream;

import htsjdk.samtools.util.Interval;

/**
 * 
 * This file extract the number of reads, deletions etc stats for every gene in the list.
 * It also aligns the reads to the region
 * 
 * 
 * Pre-filter to get speed: samtools view P19764_101_S1_L001_R1_001.fastq.gz.out.bam -L test.bed
 * samtools view P19764_101_S1_L001_R1_001.fastq.gz.out.bam -F 4 -L edit_positions.bed -o subset.bam
 * 
 * 
 * Also extract reads covering the region
 * 
 * @author Johan Henriksson
 *
 */
public class DetectEditsPrepare {
	
	
	

	
	/**
	 * Read FASTA with the sequence of each gene -- no extra to aid alignment
	 */
	public static HashMap<String,String> readGeneFastaNoExtra(File fFASTA) throws IOException {
		HashMap<String,String> mapGeneOrigfastaNoExtra=new HashMap<String, String>();
		
		BufferedReader br=new BufferedReader(new FileReader(fFASTA));
		
		String line;
		while((line=br.readLine())!=null) {
			String wbid=line.substring(1);
			if(wbid.indexOf(":")!=-1)
				wbid=wbid.substring(0,wbid.indexOf(":"));
			line=br.readLine();
			mapGeneOrigfastaNoExtra.put(wbid,line);
		}
		
		br.close();
		
		return mapGeneOrigfastaNoExtra;
	}
	
	
	/**
	 * Read what genes should be edited to. wbid -> new sequence ---- NEW, ONLY RELIES ON WBID
	 */
	public static HashMap<String, String> readMapWbidEditedseqNEW(File fCSV) throws IOException {
		HashMap<String,String> attemptedGeneEdits=new HashMap<String, String>();
	
		//@SuppressWarnings("resource")
		BufferedReader br=new BufferedReader(new FileReader(fCSV));
				
		br.readLine(); //Fine should have a header
		String line=null;
		while((line=br.readLine())!=null) {
			String[] col=line.split(",", 0);
			
			String wbid=col[0];  /// Should have WBID in first column
			
			String newseq=col[1];  // template in second column
			
			
			newseq=newseq.replaceAll(" ", "");
			newseq=newseq.replaceAll("\\.", "");
			newseq=newseq.toUpperCase();
			
			attemptedGeneEdits.put(wbid, newseq);
		}
		br.close();
		return attemptedGeneEdits;
	}
	
	
	public static void main(String[] args) throws IOException {
		
		//Default files
		File fPos=new File("/home/mahogny/Desktop/celegans/edit_positions.bed");
		File fFastaNoExtra=new File("/home/mahogny/mystore/dataset/changchun/gpcr/refgenome/list_geneedit_pos_NEW.fa");
		File fTemplates=new File("gene_edits_20230331.csv");
		File fPrep=new File("/home/mahogny/Desktop/celegans/prepped_detection");
		
		//Get name of files
		if(args.length!=0) {
			fPos=new File(args[0]);
			fFastaNoExtra=new File(args[1]);
			fTemplates=new File(args[2]);
			fPrep=new File(args[3]);
		}
		
		System.out.println("Reading fasta, no extra "+fFastaNoExtra);
		HashMap<String,String> mapGeneOrigfastaNoExtra=readGeneFastaNoExtra(fFastaNoExtra);
		
		
		System.out.println("Read map wbid to template");
		HashMap<String,String> mapWbidTemplate=readMapWbidEditedseqNEW(fTemplates);



		/**
		 * Read info about positions that could have been edited
		 */
		System.out.println("Reading gene edit positions");
		ArrayList<OneEdit> allEdits=new ArrayList<OneEdit>();
		String line;
		BufferedReader br=new BufferedReader(new FileReader(fPos));
		while((line=br.readLine())!=null) {

			OneEdit e=new OneEdit();
			allEdits.add(e);
			
			String[] toks=line.split("\t", 4);
			
			
			Interval interval=new Interval(
					toks[0],
					Integer.parseInt(toks[1]),
					Integer.parseInt(toks[2]));
			e.chrom=interval.getContig();
			e.from=interval.getStart();
			e.to=interval.getEnd();
			e.geneid=toks[3];

			System.out.println("Preparing "+e.geneid);

			///// Figure out interval more precisely
			GenomeEdit ge=new GenomeEdit(
					mapGeneOrigfastaNoExtra.get(e.geneid),
					mapWbidTemplate.get(e.geneid)
			);
			e.reducedFrom=e.from+ge.from;
			e.reducedTo=e.from+ge.to;
			e.alignedGenome=ge.alignedGenome;
			e.alignedTemplate=ge.alignedTemplate;
			
			System.out.println(e.geneid+"\t"+e.chrom+"\t"+e.from+":"+e.to+"\t"+e.reducedFrom+":"+e.reducedTo);
		}
		br.close();

		
		
		System.out.println("Write prep object");
	    ObjectOutputStream objectOutputStream = new ObjectOutputStream(new GZIPOutputStream(new FileOutputStream(fPrep)));
	    objectOutputStream.writeObject(allEdits);
	    objectOutputStream.close();
		
	    System.out.println("done");
		
	}
	


}
