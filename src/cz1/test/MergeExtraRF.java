package cz1.test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class MergeExtraRF {
	
	public static void main(String[] args) {
		new MergeExtraRF().merge(args[0], 
				args[1], args[2], args[3]);
	}
	
	private final Set<String> exclude;
	
	public MergeExtraRF() {
		exclude = new HashSet<String>();
		exclude.add("1");
		exclude.add("3");
		exclude.add("5");
		exclude.add("47");
	}
	
	private void merge(String contigFile, 
			String originRF, 
			String extraRF,
			String newRF) {
		Map<String, Integer> lineN = 
				new HashMap<String, Integer>();
		try {
			BufferedReader br = new BufferedReader(
					new FileReader(contigFile));
			String line;
			String[] s;
			int l=0;
			while( (line = br.readLine()) !=null) 
				lineN.put(line.split("\\s+")[0], 
						++l);
			br.close();
			BufferedWriter bw = new BufferedWriter(
					new FileWriter(newRF));
			br = new BufferedReader(
					new FileReader(originRF));
			StringBuilder os = new StringBuilder();
			while( (line=br.readLine())!=null ) {
				s = line.split("\\s+");
				if( exclude.contains(s[6]) ||
						exclude.contains(s[7]))
					continue;
				os.delete(0, os.length());
				for(int i=0; i<8; i++) {
					os.append(s[i]);
					os.append("\t");
				}
				os.append( lineN.get(s[8]) );
				os.append( "\t" );
				os.append( lineN.get(s[9]) );
				os.append("\n");
				bw.write(os.toString());
			}
			br.close();
			br = new BufferedReader(
					new FileReader(extraRF));
			while( (line=br.readLine())!=null ) {
				s = line.split("\\s+");
				if( exclude.contains(s[6]) ||
						exclude.contains(s[7]))
					continue;
				os.delete(0, os.length());
				for(int i=0; i<8; i++) {
					os.append(s[i]);
					os.append("\t");
				}
				os.append( lineN.get(s[8]) );
				os.append( "\t" );
				os.append( lineN.get(s[9]) );
				os.append("\n");
				bw.write(os.toString());
			}
			br.close();
			bw.close();
		
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
