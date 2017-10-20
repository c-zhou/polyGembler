package cz1.hmm.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.TreeBidiMap;

import cz1.util.Executor;
import cz1.util.Utils;

public class DataRetrievor extends Executor {

	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		
	}
	
	public static void main(String[] args) {
		String in_vcf = "C:\\Users\\chenxi.zhou\\Desktop\\new_alg\\____2.vcf";
		String ncsu_sample = "C:\\Users\\chenxi.zhou\\Desktop\\new_alg\\ncsu_samples.txt";
		String out_txt = "C:\\Users\\chenxi.zhou\\Desktop\\new_alg\\gbs_nsp306_02.txt";
		String[] parent_sample = new String[]{"440132", "440166"};
		int M = 10;
		
		BidiMap<String, Integer> target_sample = new TreeBidiMap<String, Integer>();
		target_sample.put(parent_sample[0], 0);
		target_sample.put(parent_sample[1], 1);
		int idx = 2;
		String line;
		try {
			BufferedReader br_s = Utils.getBufferedReader(ncsu_sample);
			while( (line=br_s.readLine())!=null ) {
				if(!target_sample.containsKey(line)) 
					target_sample.put(line, idx++);
			}
			br_s.close();
			
			// C:\\Users\\chenxi.zhou\\Desktop\\new_alg\\gbs_nsp306_01.txt
			//int[] retrieve_line = new int[]{
			//		134, 160, 184, 186, 189, 191, 192, 201, 214, 215,
			//		217, 219, 221, 222, 227, 229, 230, 233, 256, 273
			//};
			int[] retrieve_line = new int[]{
					160, 200, 213, 215, 216, 227, 229, 304, 348, 462
			};
			int N = target_sample.size();
			BufferedReader br_v = Utils.getBufferedReader(in_vcf);
			int skip = 0;
			final Map<Integer, Integer> sample_index = new HashMap<Integer, Integer>();
			String[] s;
			while( (line=br_v.readLine())!=null ) {
				skip++;
				if(line.startsWith("#")
						&& !line.startsWith("##")) {
					s = line.split("\\s+");
					for(int i=0; i!=s.length; i++)
						if(target_sample.containsKey(s[i]))
							sample_index.put(target_sample.get(s[i]), i);
					break;
				}
			}
			BufferedWriter bw_txt = Utils.getBufferedWriter(out_txt);
			bw_txt.write(N+"\n");
			bw_txt.write(M+"\n");
			bw_txt.write(parent_sample[0]);
			for(int i=1; i!=N; i++) 
				bw_txt.write("\t"+target_sample.getKey(i));
			bw_txt.write("\n");
			bw_txt.write(parent_sample[0]+"\t"+parent_sample[0]+"\n");
			bw_txt.write("M1");
			for(int i=1; i!=M; i++) bw_txt.write("\tM"+(i+1));
			bw_txt.write("\n");
			
			for(int i=0; i!=M; i++) {
				if(i==0) skip = retrieve_line[0]-skip-1;
				else skip = retrieve_line[i]-
						retrieve_line[i-1]-1;
				for(int j=0; j!=skip; j++) br_v.readLine();
				line = br_v.readLine();
				s = line.split("\\s+");
				bw_txt.write("M"+(i+1)+"\t"+s[1]+"\tA");
				int m = s[3].split(";").length;
				for(int j=1; j!=m; j++) 
					bw_txt.write(","+(char)('A'+j));
				for(int j=0; j!=N; j++) {
					idx = sample_index.get(j);
					if(s[idx].equals(".")) {
						bw_txt.write("\t0");
						for(int k=1; k!=m; k++)
							bw_txt.write(",0");
					} else bw_txt.write("\t"+s[idx].replaceAll(";", ","));
				}
				bw_txt.write("\n");
			}
			br_v.close();
			bw_txt.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
