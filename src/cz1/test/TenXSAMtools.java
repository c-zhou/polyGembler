package cz1.test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;

import cz1.util.Utils;

public class TenXSAMtools {

	private final static int chunk  = 40000;

	public static void main(String[] args) {
		String fastq_in = args[0];
		String out_prefix = args[1];

		long seq_id = 0;
		int vol = 1;
		int bc = 0;
		String line;
		String[] s;
		String bc_str=null;
		try {
			BufferedReader br = Utils.getBufferedReader(fastq_in);
			BufferedWriter bw = Utils.getBufferedWriter(out_prefix+".volume"+String.format("%04d", vol)+".fastq");
			while( (line=br.readLine())!=null ) {
				s = line.trim().split("\\s+");
				if(s.length<2) {
					for(int i=0; i<7; i++) br.readLine();
					continue;
				}
				if(bc_str==null) {
					bc_str = s[1];
					++bc;
				}
				if(!s[1].equals(bc_str)) {
					bc_str = s[1];
					++bc;
					if(bc>chunk) {
						++vol;
						bw.close();
						bw = Utils.getBufferedWriter(out_prefix+".volume"+String.format("%04d", vol)+".fastq");
						bc = 1;
					}
				}
				++seq_id;
				bw.write("@"+String.format("%012d", seq_id)+" "+s[1]+"\n");
				bw.write(br.readLine()+"\n");
				bw.write(br.readLine()+"\n");
				bw.write(br.readLine()+"\n");
				bw.write("@"+String.format("%012d", seq_id)+" "+s[1]+"\n");
				br.readLine();
				bw.write(br.readLine()+"\n");
				bw.write(br.readLine()+"\n");
				bw.write(br.readLine()+"\n");
			}
			bw.close();
			br.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
	}
}
