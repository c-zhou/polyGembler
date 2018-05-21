package cz1.test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;

public class MUMmerUtils {

	public static void variantParser(String snp_file, String out_file) {
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(new File(snp_file))));
			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(out_file))));
			for(int i=0; i<5; i++) br.readLine();
			int P1, P2, tmp_P1, tmp_P2;
			String line = br.readLine();
			String[] s = line.trim().split("\\s+");
			final int sLen = s.length-1;
			String allele1, allele2, sequence;
			while( line!=null ) {
				P1 = Integer.parseInt(s[0]);
				P2 = Integer.parseInt(s[3]);
				sequence = s[sLen];
				allele1 = "";
				allele2 = "";
				allele1 += s[1];
				allele2 += s[2];
				while( (line=br.readLine())!=null ) {
					s = line.trim().split("\\s+");
					tmp_P1 = Integer.parseInt(s[0]);
					tmp_P2 = Integer.parseInt(s[3]);
					if(tmp_P1==P1&&s[1].equals(".")||
							tmp_P2==P2&&s[2].equals(".")) {
						if(!s[1].equals(".")) allele1 += s[1];
						if(!s[2].equals(".")) allele2 += s[2];
					}
					else break;
				}
				if(!allele1.contains("N")&&!allele2.contains("N")) 
					bw.write(sequence+"\t"+(P1+1)+"\t"+(P2+1)+"\t"+allele1+"\t"+allele2+"\n");
			}
			bw.close();
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		variantParser(args[0], args[1]);
	}
}
