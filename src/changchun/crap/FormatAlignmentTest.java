package changchun.crap;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.HashMap;

import changchun.OneEdit;
import changchun.util.CCutil;
import changchun.util.ClustalOmega;

public class FormatAlignmentTest {
	
	public static void main(String[] args) throws IOException, ClassNotFoundException {
		
		
		CCutil.readGeneFastaWithExtra();
		
		
		File fSerial=new File("P19764_101_S1_L001_R1_001.fastq.gz.out.bam.mut.serial");
	    ObjectInputStream ois = new ObjectInputStream(new FileInputStream(fSerial));
	    @SuppressWarnings("unchecked")
		ArrayList<OneEdit> allEdits=(ArrayList<OneEdit>)ois.readObject();
	    ois.close();

	    for(OneEdit e:allEdits) {
	    	
	    	//WBGene00005319
	    	
	    	//e.reads.add("ATGATATCATCTATAGATAACTATTATTCTACAAATTATTCTAAATGTAATCTCAGTGAGAGCTTCTTAGCTTCTTGGAAAGGAGTTGCGTACCCTACTGATTTTATTCAAATATTTTCCCTACCACTTCAAATCCTCGCATTTTACATTATCCTAACAAAAACTCCAGTTCAAATGAAAAGTATGAAGTGGCCATTGTTCTATAATCACCTATTTTGTTCCATTTTCGACGTGATTCTATGCACTTTTTCAACTATTTACATAATTCTTCCGATGCTGGGGGTCTTCACTGTCGGAGTTTTTAGTTGGTTGGGGATTCCTATTATCCTGGAGTTGATTCTCCTTACGTGCTCTTTGCTGTGTGAGTATTTAAATTTGGGGAGCCAGGAAAGTGTCGGCTGCTTCAGTAAAACTGCCGCCAGGAAGTGTCAGTTGCTGAGAAGCCTCATCAATATTGAAAAAATTGAAGCACCAGGGATATGTCGGCTGCTTCTGTAAACTTTAAAAACAAAAACTACACTTGACGGGTGGCCCGTATTGTGTCGGCTGCTTAGTCCACGATCGCATATTGATATGCCAACCCGGAAACTGTCGGCTGCACTTCATGATCTTCATTTAAAAACACAGTATACTATTTTCAGCACTTGCCTTCTCCTACATCTACCTCTTCGAGAGTCGCAGCCGTGCGGTTTCACAGAATCCGTTTAAAATGAGCAGGAAAAGCACGAGAATCAAATATTATTCCTATTTAATATTATCTTACTCCACTATTTTCATATTTCTTATAATTATCCCTTCCGACCAGGAAACTGCAAAACTTCAAGCGTTGCAAGCTTACCCCTGCCCGACACAAGAGTTTTTCACGTTTCCAATTCTGATTCTAAACTCTGATTCAACAACTTCCACATTTATAGTGATTATTTTTATGCCTATATTCATAGCTCATTCAGTTGGTCATGGCGTTTTTCATGTTACATGGACAATATGGTATTTGTATGTTGCACCTTCAAATCAAGTTTCGATAGAGACTCAGAAGAAGCAAAAAACGTTCCTCAAAAATGTTATATTGCAATTTTCGATTCCATCAGTATTCATTTTATTCTCAGTTGCAATTATTTTTACATCTAGTTTTTATTCGCAAGAAATGATGAATCTAGGTGTGGACGTTGCAGGTAAGATTGATGAAAAGCAGCCGACATATCCCGGCTGTCACAATACTTTAGATTGATACGACACTTACCCGGGAAGTGTCGGACGCATCACATAACCTGTAAGTTGTCGGCTGCTTAGATTCATAATTATTCTTCAGGACTTCATGGCATCGGTGAAAGTATTGCAGTTATTTTCGTTCATTCACCTTACCGGAAAGCCGCTGGGCAGCTAATTTTCAGATGGAACGGTGAGTTGGAAGTTAACTCTTGCAAGTGATATAGTATACATTTTAAGATCCACAATCTAGCACACGAGTTTCTGCTGTACGGGGGAATTCCAGTGCGTTCACTATTTG");
	    	e.reads.add(CCutil.mapGeneOrigfastaWithExtra.get("WBGene00005319"));
	    	
	    	e.alignment=ClustalOmega.call(e.reads);
	    	
	    	
	    	if(!e.alignment.contentEquals("")) {
		    	//System.out.println(e.alignment);

	    		System.out.println(e.geneid);
	    		System.out.println();
	    		
		    	BufferedReader br=new BufferedReader(new StringReader(e.alignment));
		    	
		    	//Remove the first lines
		    	for(int i=0;i<3;i++)
		    		br.readLine();
		    	
		    	ArrayList<String> lines=new ArrayList<String>();
		    	for(;;) {
		    		String line=br.readLine();
		    		if(line.contentEquals("") || line.startsWith(" "))
		    			break;
		    		lines.add(line.substring(0));
		    	}
		    	br.readLine();
		    	
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
		    	
		    	
		    	//e.alignment.split("\n", arg1)
		    	
    			for(int i=0;i<lines.size();i++) {
		    		System.out.println(lines.get(i));
    			}
    	    	System.out.println();
    	    	System.out.println();
    	    	System.out.println();
		    	
		    	br.close();
		    	
		    	
		    	System.exit(0);
	    	}

	    	
	    }
	    
	}

}
