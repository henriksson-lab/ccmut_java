package changchun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

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
		
		
		
		
	}
	
	
	
	
	public static void main(String[] args) throws IOException {
		BufferedReader br=new BufferedReader(new FileReader("/home/mahogny/Desktop/celegans/edit_positions.bed"));
		File fBAM=new File("/home/mahogny/Desktop/celegans/P19764_101_S1_L001_R1_001.fastq.gz.out.bam");
		
		
		ArrayList<OneEdit> edits=new ArrayList<DetectEdits.OneEdit>();
		
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
			
			edits.add(e);
		}
		
		
		final SamReader reader = SamReaderFactory.makeDefault().open(fBAM);
		int readRecords=0;
		for (final SAMRecord samRecord : reader) {
			readRecords++;
			if(readRecords%1000000 == 0){
				System.out.println(readRecords);
			}
			
			
			
			editloop: for(OneEdit edit:edits) {
				boolean found=false;
				if(samRecord.overlaps(edit.interval)) {
					

					int numdel=0;
					Cigar cigar=samRecord.getCigar();
					for(CigarElement e:cigar.getCigarElements()) {
						CigarOperator op=e.getOperator();
						if(op==CigarOperator.D) {
							numdel+=e.getLength();
						}
					}
					
					System.out.println("numdel "+numdel);
				
					found=true;
					break editloop;
				}
				if(!found)
					System.out.println("no overlap");
				
			}
			
			
		}
		
		
		br.close();
	}

}
