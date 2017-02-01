package cz1.test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;

public class ContigSplit {

	public static void main(String[] args) {
		split(args[0], args[1], args[2]);
	}
	
	private static void split(String breakageFile, String contigFile, String out) {
		
		try {
			BufferedReader br = getBufferedReader(breakageFile);
			String line;
			String[] s;
			Map<String, int[]> breakages = new HashMap<String, int[]>();
			while ( (line=br.readLine())!=null ) {
				s = line.split("\\s+");
				int[] b = new int[s.length-1];
				for(int i=1; i<s.length; i++)
					b[i-1] = Integer.parseInt(s[i]);
				breakages.put(s[0], b);
			}
			br.close();
			
			br = getBufferedReader(contigFile);
			line = br.readLine();
			while( line!=null ) {
				if( !line.startsWith(">") ||
						!breakages.containsKey(line.substring(1)) ) {
					line = br.readLine();
					continue;
				}

				String original_id = line.substring(1);

				StringBuilder contig_sb = new StringBuilder(); 
				while( (line=br.readLine())!=null && !line.startsWith(">"))
					contig_sb.append(line);
				String contig_str = contig_sb.toString();
				
				String id = original_id.replaceAll("^Itr_sc0{0,6}", "").
						replaceAll(".1$", "");
				int[] breakage = breakages.get(original_id);
				int[] interval = new int[breakage.length+2];
				for(int i=1; i<interval.length-1; i++) interval[i] = breakage[i-1]-1;
				interval[interval.length-1] = contig_str.length();
				
				
				for(int i=0; i<interval.length-1; i++) {
					BufferedWriter bw = new BufferedWriter(new FileWriter(
							out+"/"+id+"_"+(i+1)+".fa"));
					bw.write(">"+id+"_"+(i+1)+"\n");
					bw.write(contig_str.substring(interval[i], interval[i+1]));
					bw.write("\n");
					bw.close();
				}
			}
			br.close();
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}
	
	public static BufferedReader getBufferedReader(String path) throws IOException {
		BufferedReader br = null;
		if(path.endsWith(".gz") || path.endsWith(".zip")){
			br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new 
					FileInputStream(path))), 65536);
		}else{
			br = new BufferedReader(new FileReader(path), 65536);
		}
		return br;
	}
}
