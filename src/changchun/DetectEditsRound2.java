package changchun;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.TreeSet;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import changchun.util.CCutil;
import changchun.util.ClustalOmega;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
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
public class DetectEditsRound2 {
	
	public static boolean showFound=false;



	
	
	
	
	
	public static boolean useProf=false;
	
	/**
	 * Align all the sequences
	 */
	public static void align(OneEdit e) throws IOException {

		//Deduplicate the reads
		TreeSet<String> uniqueReads = new TreeSet<String>(e.reads);
		e.reads.clear();
		e.reads.addAll(uniqueReads);


		if(useProf) {
			//Add original genome sequence 
			//Add sequence to substitute
			LinkedList<String> alignedReads = new LinkedList<>();
			alignedReads.add(e.alignedGenome);
			alignedReads.add(e.alignedTemplate);

			//edit.alignment=Kalign.call(edit.reads);
			e.alignment=ClustalOmega.call(uniqueReads, alignedReads);
		} else {
			//Add original genome sequence 
			//Add sequence to substitute
			e.reads.add(e.alignedGenome);
			e.reads.add(e.alignedTemplate);

			//edit.alignment=Kalign.call(edit.reads);
			e.alignment=ClustalOmega.call(e.reads);

		}
	}
	
	

	
	
	
	public static void main(String[] args) throws IOException, ClassNotFoundException {
		File fPrep=new File("/home/mahogny/Desktop/celegans/prepped_detection");
		File fBAM=new File("/home/mahogny/Desktop/celegans/P19764_101_S1_L001_R1_001.fastq.gz.out.bam");
		File fOut=new File("/home/mahogny/Desktop/celegans/mutant_summary");

		//Get name of files
		if(args.length!=0) {
			fPrep=new File(args[0]);
			fBAM=new File(args[1]);
			fOut=new File(args[2]);
		}
		
		
		System.out.println("=== reading prep object =========");
		
		ObjectInputStream ois = new ObjectInputStream(new GZIPInputStream(new FileInputStream(fPrep)));
		ArrayList<OneEdit> allEdits=(ArrayList<OneEdit>)ois.readObject();
		ois.close();
		
		HashMap<String, ArrayList<OneEdit>> mapChromEdits=new HashMap<String, ArrayList<OneEdit>>();
		for(OneEdit e:allEdits) {
			
			e.interval=new Interval(
					e.chrom, 
					e.reducedFrom - 50, 
					e.reducedTo + 50);
			
			ArrayList<OneEdit> elist=mapChromEdits.get(e.interval.getContig());
			if(elist==null)
				mapChromEdits.put(e.interval.getContig(),elist=new ArrayList<OneEdit>());
			elist.add(e);				
		}
		
		
		/**
		 * Loop through BAM file and extract overlapping reads
		 */
		System.out.println("============= Reading BAM");
		final SamReader reader = SamReaderFactory.makeDefault().open(fBAM);
		int readRecords=0;
		int usedread=0;
		for (final SAMRecord samRecord : reader) {
			readRecords++;
			if(readRecords%1000000 == 0){
				System.out.println(readRecords);
			}

			if(usedread> 5000){
				//break;
			}

			
			//check that read mapped in proper pair (0x2)
			if((samRecord.getFlags() & 0x2)!=0) {
				
				//with bowtie2, can remove the second read like this. but some argue that this is overkill
				//-F "[XS] == null and not unmapped  and not duplicate"
								
				boolean found=false;
				ArrayList<OneEdit> elist=mapChromEdits.get(samRecord.getContig());
				if(elist!=null) {
					editloop: for(OneEdit edit:elist) {
						

						/**
						 * Count how many deletions etc in the sequence based on the cigar
						 */
						int numdel=0, numins=0, numgap=0;
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
							if(op==CigarOperator.N) {
								numgap+=e.getLength();
							}
						}
						
						//Discard reads with too large gaps
						if(numgap>20) {
							break editloop;
						}
						
						
						if(samRecord.overlaps(edit.interval)) {
							usedread++;

							/**
							 * Gather the reads. reverse complement if on the other strand
							 * 
							 * read reverse strand (0x10) 
							 */
							String readString=CCutil.revcomp(samRecord.getReadString());
							if((samRecord.getFlags() & 0x10)!=0) { 
								readString = CCutil.revcomp(samRecord.getReadString());
								//System.out.println("rc");
																
								//Is this truly needed?
								
							} else {
								//System.out.println("not_rc");
								readString = samRecord.getReadString();
							}

							
							edit.reads.add(readString);

							/*
							System.out.println(edit.geneid);
							System.out.println("CIGAR "+samRecord.getCigarString());
							System.out.println(readString);
							System.out.println();*/
							
							
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
		
		
		/**
		 * Write statistics
		 */
		PrintWriter pw=new PrintWriter(fOut);
		for(OneEdit e:allEdits) {
			pw.println("reads\t"+e.geneid+"\t"+e.numread);
			pw.println("del\t"+e.geneid+"\t"+e.sumdel);
			pw.println("ins\t"+e.geneid+"\t"+e.sumins);
			pw.println("fine\t"+e.geneid+"\t"+e.fine);
			pw.println("totc\t"+e.geneid+"\t"+readRecords);
		}
		pw.close();
				
		
		/**
		 * Write all the objects; new method that conserves memory!
		 */
		System.out.println("serializing and aligning");
		File fSerial=new File(fOut.getParentFile(), fOut.getName()+".serial");
	    FileOutputStream fileOutputStream = new FileOutputStream(fSerial);
	    GZIPOutputStream gzipOutputStream = new GZIPOutputStream(fileOutputStream);
	    ObjectOutputStream objectOutputStream = new ObjectOutputStream(gzipOutputStream);

	    objectOutputStream.writeInt(allEdits.size());

	    //Transfer to a list so we can kick alignments ASAP
	    LinkedList<OneEdit> listAllEdits=new LinkedList<>(allEdits);
	    allEdits.clear();
	    mapChromEdits.clear();
	    
	    while(!listAllEdits.isEmpty()) {
	    	OneEdit e = listAllEdits.removeFirst();
			System.out.println("align "+e.geneid);
			
			align(e);
			
			//Remove interval variable - not serializable
			e.interval=null; 
			
		    objectOutputStream.writeObject(e);

	    }
	    
		
	    objectOutputStream.close();
	    gzipOutputStream.close();
	    fileOutputStream.close();
	    
	    System.out.println("done");
		
	}
	


}
