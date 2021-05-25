package changchun.util;

import java.io.File;
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
public class Kalign {
	
	
	public static String call(Collection<String> reads) throws IOException {
		
		File tmpOut=File.createTempFile("foo", "align");
		
		try {
	        ProcessBuilder pb = new ProcessBuilder(
	        		"kalign","-out",tmpOut.getAbsolutePath(),
	        		"-f","aln",
	        		"-t","0",
					"-sort","tree");
	        System.out.println(pb.command()+"   #reads "+reads.size());
	        
	        
	        if(reads.size()==0) {
	        	return "";
	        } else {
		        Process process = pb.start();
		        
		        PrintWriter pw=new PrintWriter(process.getOutputStream());
		        
		        int i=0;
				for(String seq:reads) {
					pw.println(">"+i);
					pw.println(seq);
					i++;
				}
				pw.close();
		        

				//BufferedReader stdInput = new BufferedReader(new InputStreamReader(process.getInputStream()));
				//BufferedReader stdError = new BufferedReader(new InputStreamReader(process.getErrorStream()));

				// Read the output from the command
				//while ((stdInput.readLine()) != null);
				//while ((stdError.readLine()) != null);
				boolean ok=process.waitFor(10,TimeUnit.SECONDS);
				
				
				String out;
				if(ok) {
					out=new String(Files.readAllBytes(Paths.get(tmpOut.getAbsolutePath())));
					tmpOut.delete();
				} else {
					out="err";
					System.out.println("error, moving on");
				}
				
				tmpOut.delete();

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
