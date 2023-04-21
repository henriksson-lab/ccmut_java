package changchun;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.StringTokenizer;
import java.util.TreeMap;

public class FindCommonBCF {
	
	
	/**
	 * 
	 * Count BCFs to find which are common over many files
	 * 
	 */
	public static void main(String[] args) throws IOException {
		
		File root=new File("/home/mahogny/bcf");
		
		TreeMap<String, TreeMap<Integer, Integer>> mapChrom=new TreeMap<>();
		
		int countFile=0;
		for(File f:root.listFiles()) {
			if(f.getName().endsWith(".bcf")) {
				System.out.println(f);
				System.out.println(countFile);
				
				Runtime rt = Runtime.getRuntime();
				String[] commands = { "bcftools", "view", f.getAbsolutePath() };
				Process proc = rt.exec(commands);

				BufferedReader stdInput = new BufferedReader(new InputStreamReader(proc.getInputStream()));

				String s = null;
				while ((s = stdInput.readLine()) != null) {
				    //System.out.println(s);
				    if(!s.startsWith("#")) {

					    StringTokenizer stok=new StringTokenizer(s, "\t");
					    String contig=stok.nextToken();
					    int start=Integer.parseInt(stok.nextToken());
					    
						TreeMap<Integer, Integer> mapLoc=mapChrom.get(contig);
						if(mapLoc==null)
							mapChrom.put(contig,mapLoc=new TreeMap<>());

						Integer cnt=mapLoc.get(start);
						if(cnt==null)
							cnt=1;
						else
							cnt++;
						
						mapLoc.put(start, cnt);
						
						//System.out.println("666");

				    }
				}				
			}
			countFile++;
			/*
			if(countFile>10)
				break;
				*/
			//System.out.println(mapChrom.get("I").size());
		}
	
		//TreeMap<String, TreeMap<Integer, Integer>> mapChrom=new TreeMap<>();

		//Write to output
		PrintWriter pw=new PrintWriter(new File("/home/mahogny/bcf/commonlist.csv"));
		for(String contig:mapChrom.keySet()) {
			TreeMap<Integer, Integer> mapLoc=mapChrom.get(contig);
			for(int pos:mapLoc.keySet()) {
				int cnt=mapLoc.get(pos);
				if(cnt>1) {
					pw.println(""+contig+"\t"+pos+"\t"+cnt);
				}
			}
		}
		pw.close();
		
		
		
	}

	
	


}
