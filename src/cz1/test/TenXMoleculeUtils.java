package cz1.test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.Map;

public class TenXMoleculeUtils {
	
	private static void variantParser(String stats_file, String out_file) {
		// TODO Auto-generated method stub
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(new File(stats_file))));
			String line, chr;
			String[] s;
			String position;
			final Map<String, Map<String, int[]>> depthStats = new HashMap<String, Map<String, int[]>>();
			Map<String, int[]> stats;
			int[] stat;
			int a;
			while( (line=br.readLine())!=null ) {
				s = line.split("\\|");
				chr = s[4].trim().split(":")[0];
				if(!depthStats.containsKey(chr)) depthStats.put(chr, new HashMap<String, int[]>());
				stats = depthStats.get(chr);
				position = s[6].trim();
				if(!stats.containsKey(position)) stats.put(position, new int[2]);
				stat = stats.get(position);
				a = Integer.parseInt(s[9].trim());
				if(a<2) ++stat[a];
			}
			br.close();
			
			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(out_file))));
			for(String key : depthStats.keySet()) {
				stats = depthStats.get(key);
				for(String pos : stats.keySet()) {
					stat = stats.get(pos);
					bw.write(key+"\t"+pos+"\t"+stat[0]+"\t"+stat[1]+"\n");
				}
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	
	}
	
	public static void main(String[] args) {
		switch(args[0].toLowerCase()) {
		case "depth":
			// args: variants statistics
			// args: output file
			variantParser(args[1], args[2]);	
			break;
		default:
			throw new RuntimeException("!!!");	
		}
	}
}
