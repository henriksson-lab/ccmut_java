package changchun.util;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.LinkedList;
import java.util.concurrent.TimeUnit;


/**
 * Wrapper around Kalign
 * 
 * @author Johan Henriksson
 *
 */
public class ClustalOmega {
	
	
	public static String call(Collection<String> reads) throws IOException {
		
		File tmpIn=File.createTempFile("foo", "align");
		File tmpOut=File.createTempFile("foo", "align");
		
		//clustalo -i my-in-seqs.fa -o my-out-seqs.fa
		//Kalign.java
		
        PrintWriter pw=new PrintWriter(new FileWriter(tmpIn));
        int i=0;
		for(String seq:reads) {
			pw.println(">"+i);
			pw.println(seq);
			i++;
		}
		pw.close();

		
		
		try {
	        ProcessBuilder pb = new ProcessBuilder(
	        		"clustalo",
	        		"--force",
	        		"--outfmt=clu",
	        		"-i",tmpIn.getAbsolutePath(),
	        		"-o",tmpOut.getAbsolutePath());
	        //System.out.println(pb.command()+"   #reads "+reads.size());
	        
	        
	        if(reads.size()==0) {
	        	return "";
	        } else {
		        Process process = pb.start();
		        
		        
				boolean ok=process.waitFor(10,TimeUnit.SECONDS);
				
				String out;
				if(ok) {
					out=new String(Files.readAllBytes(Paths.get(tmpOut.getAbsolutePath())));
				} else {
					out="err";
					System.out.println("error, moving on");
				}
				
				tmpOut.delete();
				tmpIn.delete();

				return out;
	        }
			
		} catch (InterruptedException e) {
			throw new IOException(e);
		}
	}
	
	
	public static void main(String[] args) throws IOException {
		
		LinkedList<String> list=new LinkedList<String>();
		list.add("ATCGATCGATCGATCGATCGATCGATCGATCGATCG");
		list.add("ATCGATCGATCG");

		String out=call(list);
		System.out.println("--------------------");
		System.out.println(out);
		System.out.println("done");
	}

}
