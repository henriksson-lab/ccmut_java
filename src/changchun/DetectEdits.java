package changchun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Interval;

/**
 * 
 * Pre-filter to get speed: samtools view P19764_101_S1_L001_R1_001.fastq.gz.out.bam -L test.bed
 * samtools view P19764_101_S1_L001_R1_001.fastq.gz.out.bam -F 4 -L edit_positions.bed -o subset.bam
 * 
 * This file extract the number of reads, deletions etc stats for every gene in the list
 * 
 * Also extract reads covering the region
 * 
 * @author Johan Henriksson
 *
 */
public class DetectEdits {
	
	public static class OneEdit {
		String geneid;
		Interval interval;
		public int sumdel;
		public int sumins;
		public int numread;
		public int fine;
		
		LinkedList<String> reads=new LinkedList<String>();		
	}
	
	public static boolean showFound=false;
	
	
	public static void main(String[] args) throws IOException {
		
		File fPos=new File("/home/mahogny/Desktop/celegans/edit_positions.bed");
		File fBAM=new File("/home/mahogny/Desktop/celegans/P19764_101_S1_L001_R1_001.fastq.gz.out.bam");
		File fOut=new File("/home/mahogny/Desktop/celegans/mutant_summary");
		
		
		if(args.length!=0) {
			fPos=new File(args[0]);
			fBAM=new File(args[1]);
			fOut=new File(args[2]);
		}
		
		BufferedReader br=new BufferedReader(new FileReader(fPos));


		/**
		 * Read positions that could have been edited
		 */
		ArrayList<OneEdit> allEdits=new ArrayList<OneEdit>();
		HashMap<String, ArrayList<OneEdit>> mapChromEdits=new HashMap<String, ArrayList<OneEdit>>();
		String line;
		while((line=br.readLine())!=null) {
			OneEdit e=new OneEdit();
			//System.out.println(line);
			String[] toks=line.split("\t", 4);
			e.interval=new Interval(
					toks[0],
					Integer.parseInt(toks[1]),
					Integer.parseInt(toks[2]));
			e.geneid=toks[3];
			allEdits.add(e);
			
			ArrayList<OneEdit> elist=mapChromEdits.get(e.interval.getContig());
			if(elist==null)
				mapChromEdits.put(e.interval.getContig(),elist=new ArrayList<OneEdit>());
			elist.add(e);	
		}
		
		
		
		
		/**
		 * Loop through BAM file
		 */
		final SamReader reader = SamReaderFactory.makeDefault().open(fBAM);
		int readRecords=0;
		for (final SAMRecord samRecord : reader) {
			readRecords++;
			if(readRecords%1000000 == 0){
				System.out.println(readRecords);
			}

			if(readRecords%3000000 == 0){
				//System.exit(0);
			}

			
			//check that read mapped in proper pair (0x2)
			if((samRecord.getFlags() & 0x2)!=0) {
				
				//with bowtie2, can remove the second read like this. but some argue that this is overkill
				//-F "[XS] == null and not unmapped  and not duplicate"
								
				boolean found=false;
				ArrayList<OneEdit> elist=mapChromEdits.get(samRecord.getContig());
				if(elist!=null) {
					editloop: for(OneEdit edit:elist) {
						if(samRecord.overlaps(edit.interval)) {
							

							/**
							 * Gather the reads. reverse complement if on the other strand
							 * 
							 * read reverse strand (0x10) 
							 */
							if((samRecord.getFlags() & 0x10)!=0) { 
								edit.reads.add(revcomp(samRecord.getReadString()));
							} else {
								edit.reads.add(samRecord.getReadString());
							}
							
							/**
							 * Count how many deletions in the sequence based on the cigar
							 */
							int numdel=0, numins=0;
							Cigar cigar=samRecord.getCigar();
							//System.out.println(samRecord.getCigarString());
							for(CigarElement e:cigar.getCigarElements()) {
								CigarOperator op=e.getOperator();
								if(op==CigarOperator.D) {
									numdel+=e.getLength();
								}
								if(op==CigarOperator.I) {
									numins+=e.getLength();
								}
							}
							
							if(showFound) {
								System.out.println("wbid "+edit.geneid+"   "+numdel+"  "+numins);
							}
						
							edit.sumdel+=numdel;
							edit.sumins+=numins;
							edit.numread++;
							if(numdel+numins==0) {
								edit.fine++;
							}
							
							found=true;
							break editloop;
						}
					}
				}
				if(!found && showFound) {
					System.out.println("no overlap");
					System.out.println(samRecord);
				}
				
			}
		}
		br.close();
		
		
		/**
		 * Write statistics
		 */
		PrintWriter pw=new PrintWriter(fOut);
		for(OneEdit e:allEdits) {
			pw.println("reads\t"+e.geneid+"\t"+e.numread);
		}
		for(OneEdit e:allEdits) {
			pw.println("del\t"+e.geneid+"\t"+e.sumdel);
		}
		for(OneEdit e:allEdits) {
			pw.println("ins\t"+e.geneid+"\t"+e.sumins);
		}
		for(OneEdit e:allEdits) {
			pw.println("fine\t"+e.geneid+"\t"+e.fine);
		}
		pw.close();
		
		
		/**
		 * Write the reads for each gene
		 */
		for(OneEdit edit:allEdits) {
			File fFastq=new File(fOut.getParentFile(), fOut.getName()+".gene."+edit.geneid);
			pw=new PrintWriter(fFastq);
			int i=0;
			for(String seq:edit.reads) {
				pw.println(">"+i);
				pw.println(seq);
				i++;
			}
			pw.close();
		}
		
	}
	

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
		else throw new RuntimeException(":(");
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

}
