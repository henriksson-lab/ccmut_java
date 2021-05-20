package changchun;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.LinkedList;

public class Kalign {
	
	
	public static String call(Collection<String> reads) throws IOException {
		
		File tmpFastq=File.createTempFile("foo", "fasta");
		File tmpOut=File.createTempFile("foo", "align");
		
		
		System.out.println(tmpFastq);
		System.out.println(tmpOut);
		
		PrintWriter pw=new PrintWriter(tmpFastq);
		int i=0;
		for(String seq:reads) {
			pw.println(">"+i);
			pw.println(seq);
			i++;
		}
		pw.close();
		
		try {
			//Runtime runTime = Runtime.getRuntime();
			/*Process process = runTime.exec(new String[] {
					"kalign",
					"-i",tmpFastq.getAbsolutePath(),
					"-o",tmpOut.getAbsolutePath(),
					});
			*/

	        ProcessBuilder pb = new ProcessBuilder(
	        		"kalign","-out",tmpOut.getAbsolutePath(),"-f","aln",
					"-sort","tree");
					//tmpFastq.getAbsolutePath(),
					//);
//					"-f","aln",
	//				"-sort","tree");
	        //pb.redirectErrorStream(true);
	        //pb.inheritIO();
	        System.out.println(pb.command());
	        Process process = pb.start();
	        
	        pw=new PrintWriter(process.getOutputStream());
	        
	        i=0;
			for(String seq:reads) {
				pw.println(">"+i);
				pw.println(seq);
				i++;
			}
			pw.close();
	        

			BufferedReader stdInput = new BufferedReader(new InputStreamReader(process.getInputStream()));
			BufferedReader stdError = new BufferedReader(new InputStreamReader(process.getErrorStream()));

			// Read the output from the command
			System.out.println("Here is the standard output of the command:\n");
			String s = null;
			while ((s = stdInput.readLine()) != null) {
			    //System.out.println(s);
			}
			while ((s = stdError.readLine()) != null) {
			    System.err.println(s);
			}

			process.waitFor();
			//kalign -i $CURRENT_FILE -o alignment/$CURRENT_FILE -f aln -sort tree
			
			String out=new String(Files.readAllBytes(Paths.get(tmpOut.getAbsolutePath())));
			
			
			tmpFastq.delete();
			tmpOut.delete();
			
			return out;
			
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
