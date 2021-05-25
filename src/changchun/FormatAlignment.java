package changchun;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.HashMap;

import changchun.util.CCutil;
import changchun.util.ClustalOmega;
import changchun.util.MatchSeq;


/**
 * Take the alignment output. Format it for viewing
 * 
 * @author Johan Henriksson
 *
 */
public class FormatAlignment {
	
	
	public static void main(String[] args) throws IOException, ClassNotFoundException {
		
		HashMap<String,String> attemptedGeneEdits=CCutil.readMapWbidEditedseq();
		CCutil.readGeneFasta();
		
		
		
		
		File fSerial=new File("/media/mahogny/TOSHIBA/changchun/mut/P19764_101_S1_L001_R1_001.fastq.gz.out.bam.mut.serial");
		File fSerialOut=new File("/media/mahogny/TOSHIBA/changchun/mut/P19764_101_S1_L001_R1_001.fastq.gz.out.bam.mut.serial2");
		
		
		if(args.length!=0) {
			fSerial=new File(args[0]);
			fSerialOut=new File(args[0]+".out");
		}
		
		@SuppressWarnings("unchecked")
		ArrayList<OneEdit> allEdits=(ArrayList<OneEdit>)CCutil.readObjectFile(fSerial);
	
	    int iii=0;
	    for(OneEdit e:allEdits) {
	    	System.out.println(iii);
	    	iii++;
	    	
	    	if(!e.alignment.contentEquals("")) {
		    	
		    	String fastaSeq=CCutil.mapGeneOrigfasta.get(e.geneid);
		    	String editSeq=attemptedGeneEdits.get(e.geneid);
		    	if(editSeq==null)
		    		throw new RuntimeException("missing "+e.geneid);
				MatchSeq m=new MatchSeq();
				if(!m.match(fastaSeq, editSeq, 30)) {
					m=new MatchSeq();
					editSeq=CCutil.revcomp(e.geneid);
					if(!m.match(fastaSeq, editSeq, 30)) {
						throw new RuntimeException("unable to match "+e.geneid);
					}
				}
		    	
	    		//Add genome fastq
		    	e.reads.add(CCutil.mapGeneOrigfasta.get(e.geneid));
				//Add sequence to substitute
		    	e.reads.add(editSeq);
		    	
		    	//Align
		    	e.alignment=ClustalOmega.call(e.reads);
	
	    		
	    		System.out.println(e.geneid);
	    		
		    	BufferedReader br=new BufferedReader(new StringReader(e.alignment));
		    	
		    	//Remove the first lines
		    	for(int i=0;i<3;i++)
		    		br.readLine();
		    	
		    	//Read the first block of alignments
		    	ArrayList<String> lines=new ArrayList<String>();
		    	for(;;) {
		    		String line=br.readLine();
		    		if(line==null || line.contentEquals("") || line.startsWith(" "))
		    			break;
		    		lines.add(line.substring(8));
		    	}
		    	br.readLine();
		    	
		    	//Read the other blocks
		    	for(;;) {
		    		String line=br.readLine();
		    		if(line==null) {
		    			break;
		    		} else {
		    			for(int i=0;i<lines.size();i++) {
				    		lines.set(i,lines.get(i)+line.substring(8));
				    		line=br.readLine();
		    			}
		    			//Skip one line
			    		line=br.readLine();
		    		}
		    	}
		    	
		    	//Concatenate the new alignment
		    	StringBuilder sb=new StringBuilder();
				for(int i=0;i<lines.size();i++) {
					sb.append(lines.get(i));
					sb.append("\n");
				}
		    	
				e.alignment=sb.toString();
				
		    	br.close();
	    	}
	
	    	
	    }
	    
	    
		System.out.println("serializing");
		CCutil.writeObjectFile(fSerialOut, allEdits);
		
	}
}
