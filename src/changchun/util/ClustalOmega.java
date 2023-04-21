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
 * Wrapper around Clustal omega
 * 
 * @author Johan Henriksson
 *
 */
public class ClustalOmega {
	
	// 2 mins per alignment
	
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
	        		"--threads","3",
	        		"--outfmt=clu",
	        		"-i",tmpIn.getAbsolutePath(),
	        		"-o",tmpOut.getAbsolutePath());
	        //System.out.println(pb.command()+"   #reads "+reads.size());
	        
	        System.out.println("in "+tmpIn.getAbsolutePath()+"\t"+reads.size());
	        System.out.println("out "+tmpOut.getAbsolutePath());
	        
	        
	        if(reads.size()==0) {
	        	return "";
	        } else {
		        Process process = pb.start();
		        
		        
				boolean ok=process.waitFor(2000,TimeUnit.SECONDS);
				
				
				String out;
				if(ok) {
					out=new String(Files.readAllBytes(Paths.get(tmpOut.getAbsolutePath())));
				} else {
					out="err";
					System.out.println("clustalo timeout error, moving on");
					//System.exit(1);
					
					try {
						process.destroy();
					} catch (Exception e) {
					}

				}
				
				tmpOut.delete();
				tmpIn.delete();
				
		        System.out.println();

				return out;
	        }
			
		} catch (InterruptedException e) {
			throw new IOException(e);
		}
	}
	
	
	
	
	
	//clustalo -i temp_3.fa --profile1 temp_2.fa > out

	public static String call(Collection<String> reads, Collection<String> readsAligned) throws IOException {
		
		File tmpIn=File.createTempFile("foo", "align");
		File tmpInAligned=File.createTempFile("foo", "align");
		File tmpOut=File.createTempFile("foo", "align");
		
		//clustalo -i my-in-seqs.fa -o my-out-seqs.fa
		//Kalign.java
		System.out.println();
		System.out.println();
		System.out.println();
		System.out.println("in1==============");

		
        PrintWriter pw=new PrintWriter(new FileWriter(tmpIn));
        int i=0;
		for(String seq:reads) {
			pw.println(">"+i);
			pw.println(seq);
			i++;

			//System.out.println(">"+i);
			//System.out.println(seq);

		}
		pw.close();

		System.out.println();
		System.out.println();
		System.out.println();
		System.out.println("in2==============");
		
        pw=new PrintWriter(new FileWriter(tmpInAligned));
        i=0;
		for(String seq:readsAligned) {
			pw.println(">a"+i);
			pw.println(seq);
			i++;
			
			//System.out.println(">"+i);
			//System.out.println(seq);

		}
		pw.close();
		
		
		try {
	        ProcessBuilder pb = new ProcessBuilder(
	        		"clustalo",
	        		"--force",
	        		"--threads","10",
	        		"--outfmt=clu",
	        		"-i",tmpIn.getAbsolutePath(),
	        		"--profile1", tmpInAligned.getAbsolutePath(),
	        		"-o",tmpOut.getAbsolutePath());
	        //System.out.println(pb.command()+"   #reads "+reads.size());
	        
	        System.out.println("Align "+reads.size()+"  "+readsAligned.size());
	        //System.out.println("in "+tmpIn.getAbsolutePath()+"\t"+reads.size());
	        //System.out.println("out "+tmpOut.getAbsolutePath());
	        
	        
	        if(reads.size()==0) {
	        	return "";
	        } else {
		        Process process = pb.start();
		        
		        
				boolean ok=process.waitFor(2000,TimeUnit.SECONDS);
				
				
				String out;
				if(ok) {
					out=new String(Files.readAllBytes(Paths.get(tmpOut.getAbsolutePath())));
				} else {
					out="err";
					System.out.println("clustalo timeout error, moving on");
					//System.exit(1);
					
					try {
						process.destroy();
					} catch (Exception e) {
					}

				}
				
				tmpOut.delete();
				tmpIn.delete();
				tmpInAligned.delete();

				/*
				System.out.println();
				System.out.println();
				System.out.println();
				System.out.println("out=============");

				System.out.println(out);
		        System.out.println();
		        System.out.println();
		        System.out.println();*/

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
