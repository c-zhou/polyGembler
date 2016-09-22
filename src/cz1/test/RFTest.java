package cz1.test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class RFTest {
	public static void main(String[] args) {
		
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(args[0]));
			String line;
			String[] s;
			Map<String, Integer> idxMap = new HashMap<String, Integer>();
			int i=0;
			while( (line=br.readLine())!=null ){
				s = line.split("\\s+");
				idxMap.put(s[0], ++i);
			}
			br.close();
			br = new BufferedReader(new FileReader(args[1]));
			BufferedWriter bw = new BufferedWriter(new FileWriter(args[1]));
			while( (line=br.readLine())!=null ) {
				s = line.split("\\s+");
				for(int k=0; k<8; k++)
					bw.write(s[k]+"\t");
				bw.write(idxMap.get(s[6])+"\t");
				bw.write(idxMap.get(s[7])+"\n");
			}
			bw.close();
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}
}
