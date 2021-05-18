package changchun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

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
 * @author Johan Henriksson
 *
 */
public class DetectEdits {
	
	public static class OneEdit {
		String geneid;
		Interval interval;
		public int sumdel;
		public int numread;
		
		
		
		
	}
	
	public static boolean showFound=false;
	
	
	@SuppressWarnings("unused")
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
		ArrayList<OneEdit> allEdits=new ArrayList<DetectEdits.OneEdit>();
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
				mapChromEdits.put(e.interval.getContig(),elist=new ArrayList<DetectEdits.OneEdit>());
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

			if(readRecords%100000000 == 0){
				System.out.println(readRecords);
			}

			boolean found=false;
			ArrayList<OneEdit> elist=mapChromEdits.get(samRecord.getContig());
			if(elist!=null) {
				editloop: for(OneEdit edit:elist) {
					if(samRecord.overlaps(edit.interval)) {
						

						/**
						 * Count how many deletions in the sequence based on the cigar
						 */
						int numdel=0;
						Cigar cigar=samRecord.getCigar();
						for(CigarElement e:cigar.getCigarElements()) {
							CigarOperator op=e.getOperator();
							if(op==CigarOperator.D) {
								numdel+=e.getLength();
							}
						}
						
						System.out.println("wbid "+edit.geneid+"   "+numdel);
					
						edit.sumdel+=numdel;
						edit.numread++;
						
						
						
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
		br.close();
		
		
		PrintWriter pw=new PrintWriter(fOut);
		for(OneEdit e:allEdits) {
			pw.println("reads\t"+e.geneid+"\t"+e.numread);
		}
		for(OneEdit e:allEdits) {
			pw.println("del\t"+e.geneid+"\t"+e.sumdel);
		}
		pw.close();
	}

}
