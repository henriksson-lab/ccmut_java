package changchun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.TreeMap;
import java.util.TreeSet;
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
public class DetectEditsRound2_20230406 {
	
	public static boolean showFound=false;

	public static HashMap<String,String> mapWbidTemplate;

	public static HashMap<String,String> mapGeneOrigfastaNoExtra=new HashMap<String, String>();
	public static TreeMap<String, String> mapGeneOrigfastaWithExtra=new TreeMap<String, String>();


	/**
	 * Dirty check if proper direction
	 */
	public static boolean match(String fasta, String editSeq, int len) {
		String s1=editSeq.substring(0,len);
		int start=editSeq.length()-len;
		String s2=editSeq.substring(start,start+len);
		return fasta.indexOf(s1) != -1 || fasta.indexOf(s2) != -1;
	}
	
	/**
	 * Align all the sequences
	 */
	public static void align(OneEdit e) throws IOException {
		
    	String fastaSeqWithExtra=mapGeneOrigfastaWithExtra.get(e.geneid);
    	if(fastaSeqWithExtra==null)
    		throw new RuntimeException("missing fastaseq extra "+e.geneid);

    	String editSeq=mapWbidTemplate.get(e.geneid);
    	if(editSeq==null)
    		throw new RuntimeException("missing editseq "+e.geneid);
    	
    	
    	//Figure out if the edit sequence should be reverse complemented or not
		//MatchSeq m=new MatchSeq();
		if(!match(fastaSeqWithExtra, editSeq, 20)) {
			editSeq=CCutil.revcomp(editSeq);
			if(!match(fastaSeqWithExtra, editSeq, 20)) {
				editSeq=CCutil.revcomp(editSeq);
				if(!match(fastaSeqWithExtra, editSeq, 10)) {  // Pray that shorter sequence will do
					editSeq=CCutil.revcomp(editSeq);
					if(!match(fastaSeqWithExtra, editSeq, 10)) {
						System.out.println("unable to match "+e.geneid+"  "
								+editSeq
								+CCutil.revcomp(editSeq)
								+"\n\n"+fastaSeqWithExtra);
						editSeq=fastaSeqWithExtra;  //At least it aligns!
						//throw new RuntimeException(
						//		);
					
					}
				}
			}
		}


		//Deduplicate the reads
		TreeSet<String> uniqueReads = new TreeSet<String>(e.reads);
		e.reads.clear();
		e.reads.addAll(uniqueReads);
		
		//Add original genome sequence 
    	e.reads.add(mapGeneOrigfastaWithExtra.get(e.geneid));
		//Add sequence to substitute
    	e.reads.add(editSeq);

		//edit.alignment=Kalign.call(edit.reads);
		e.alignment=ClustalOmega.call(e.reads);
	}
	
	
	
	
	/**
	 * Read FASTA with the sequence of each gene -- some extra to aid alignment
	 */
	public static void readGeneFastaWithExtra(File fFASTA) throws IOException {
		BufferedReader br=new BufferedReader(new FileReader(fFASTA));
		
		String line;
		while((line=br.readLine())!=null) {
			String wbid=line.substring(1);
			wbid=wbid.substring(0,wbid.indexOf(":"));
			line=br.readLine();
			mapGeneOrigfastaWithExtra.put(wbid,line);
		}
		
		br.close();
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

		
		// /home/mahogny/mystore/dataset/changchun/gpcr/refgenome/list_geneedit_pos_NEW.bed
		
		//Default files
		File fPos=new File("/home/mahogny/Desktop/celegans/edit_positions.bed");
		File fFastaWithExtra=new File("/home/mahogny/mystore/dataset/changchun/gpcr/refgenome/list_geneedit_pos_wextra_NEW.fa");
		File fTemplates=new File("gene_edits_20230331.csv");
		File fBAM=new File("/home/mahogny/Desktop/celegans/P19764_101_S1_L001_R1_001.fastq.gz.out.bam");
		File fOut=new File("/home/mahogny/Desktop/celegans/mutant_summary");
		
		//Get name of files
		if(args.length!=0) {
			fPos=new File(args[0]);
			fFastaWithExtra=new File(args[1]);
			fTemplates=new File(args[2]);
			fBAM=new File(args[3]);
			fOut=new File(args[4]);
		}
		
		readGeneFastaWithExtra(fFastaWithExtra);
		mapWbidTemplate=readMapWbidEditedseqNEW(fTemplates);



		/**
		 * Read info about positions that could have been edited
		 */
		ArrayList<OneEdit> allEdits=new ArrayList<OneEdit>();
		HashMap<String, ArrayList<OneEdit>> mapChromEdits=new HashMap<String, ArrayList<OneEdit>>();
		String line;
		BufferedReader br=new BufferedReader(new FileReader(fPos));
		while((line=br.readLine())!=null) {
			OneEdit e=new OneEdit();
			//System.out.println(line);
			String[] toks=line.split("\t", 4);
			e.interval=new Interval(
					toks[0],
					Integer.parseInt(toks[1]),
					Integer.parseInt(toks[2]));
			e.chrom=e.interval.getContig();
			e.from=e.interval.getStart();
			e.to=e.interval.getEnd();
			e.geneid=toks[3];
			allEdits.add(e);
			
			ArrayList<OneEdit> elist=mapChromEdits.get(e.interval.getContig());
			if(elist==null)
				mapChromEdits.put(e.interval.getContig(),elist=new ArrayList<OneEdit>());
			elist.add(e);	
		}
		
		
		
		
		/**
		 * Loop through BAM file and extract overlapping reads
		 */
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
						if(samRecord.overlaps(edit.interval)) {
							usedread++;

							/**
							 * Gather the reads. reverse complement if on the other strand
							 * 
							 * read reverse strand (0x10) 
							 */
							if((samRecord.getFlags() & 0x10)!=0) { 
								edit.reads.add(CCutil.revcomp(samRecord.getReadString()));
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

	    //Transfer to a list so we can kick them out soon enough
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
