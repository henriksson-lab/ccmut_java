package changchun;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.TreeMap;

public class ReadGeneMapping {

	public static void main(String[] args) throws IOException {
		BufferedReader br=new BufferedReader(new FileReader(new File("/home/mahogny/Desktop/celegans/map_wbid_transcid.csv")));
		
		System.out.println(br.readLine());
		
		TreeMap<String, String> map=new TreeMap<String, String>();
		
		String line=null;
		while((line=br.readLine())!=null) {
			String[] col=line.split(",", 0);
			
			String wbid=col[0];
			
			String tr_name=col[1];
			String gene_name=col[2];
			
			//System.out.println(gene_name);
			
			map.put(tr_name,wbid);
			map.put(gene_name,wbid);
			
		}
		
		
		br.close();

	}
}
