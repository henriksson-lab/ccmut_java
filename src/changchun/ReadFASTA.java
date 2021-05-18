package changchun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.TreeMap;

public class ReadFASTA {
	
	public static void main(String[] args) throws IOException {
		TreeMap<String, String> map=new TreeMap<String, String>();
		BufferedReader br=new BufferedReader(new FileReader(new File("/home/mahogny/Desktop/celegans/all.bed.fa")));
		String toGene=null;
		StringBuilder sb=new StringBuilder();
		
		String line=null;
		while((line=br.readLine())!=null) {
			
			if(line.startsWith(">")) {
				if(toGene!=null) {
					map.put(toGene,sb.toString());
				}
				toGene=line.substring(1).split(":", 0)[0];
				sb=new StringBuilder();
			} else {
				sb.append(line);
			}
		}
		map.put(toGene,sb.toString());
		br.close();
	}
	
	

}
