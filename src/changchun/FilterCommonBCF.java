package changchun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.StringTokenizer;
import java.util.TreeMap;

public class FilterCommonBCF {

	
	public static void main(String[] args) throws IOException {

		
		File root=new File("/home/mahogny/bcf");

		//Read list of common SNPs
		System.out.println("Reading list");
		TreeMap<String, TreeMap<Integer, Integer>> mapChrom=new TreeMap<>();
		BufferedReader br=new BufferedReader(new FileReader(new File("/home/mahogny/bcf/commonlist.csv")));
		String line;
		while((line=br.readLine())!=null) {
			StringTokenizer stok=new StringTokenizer(line,"\t");
		    String contig=stok.nextToken();
		    int start=Integer.parseInt(stok.nextToken());
		    int cnt=Integer.parseInt(stok.nextToken());
			TreeMap<Integer, Integer> mapLoc=mapChrom.get(contig);
			if(mapLoc==null)
				mapChrom.put(contig,mapLoc=new TreeMap<>());
			if(cnt>2) {
				mapLoc.put(start, cnt);
			}
		}
		br.close();
		
		
		
		int countFile=0;
		for(File f:root.listFiles()) {
			if(f.getName().endsWith(".bcf")) {
				System.out.println(f);
				System.out.println(countFile);
				
				File outfile=new File(f.getParentFile(), f.getName().replaceAll("bcf", "vcf"));
				PrintWriter pw=new PrintWriter(outfile);
				
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
					    
					    if(mapChrom.containsKey(contig) && mapChrom.get(contig).containsKey(start)) {
					    	//System.out.println("Exists");
					    } else {
					    	pw.println(s);
					    }
					    
						//System.out.println("666");
				    }
				}
				pw.close();
				countFile++;
			}
			/*
			if(countFile>10)
				break;
				*/
			//System.out.println(mapChrom.get("I").size());
		}
		
		
		
	}
	
}
