package cz1.test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;

import cz1.util.IO;

public class PGSC {
	
	
	public static void main(String[] args) {
		/**
		String[] clFile = new String[]{
				"C:\\Users\\chenxi.zhou\\Desktop\\genetic mapping writing-up\\meta_final\\potatoPGSC.cl",
				"C:\\Users\\chenxi.zhou\\Desktop\\genetic mapping writing-up\\meta_final\\potatoPGSC_Stirlingx12601ab1.cl"};
		
		String line;
		String[] s;
		Map<String, int[]> scaffs = new HashMap<String, int[]>();
		try {	
			for(int i=0; i<clFile.length; i++) {
				BufferedReader br = IO.getBufferedReader(clFile[i]);
				line = br.readLine();
				while( line!=null ) {
					if(line.startsWith("#")) {
						String chr = "CH"+padding(line.replaceAll("#", ""), 2)+"PGSC0003DMB";
						while( (line=br.readLine())!=null && !line.startsWith("#")) {
							s = line.split("\\s+");
							int a = Integer.parseInt(s[1])-1;
							int b = Integer.parseInt(s[2])-a;
							scaffs.put(chr+padding(s[0],9), new int[]{a,b});
						}
					}
				}
				br.close();
			}
			**/
		//String[] clFile = new String[]{
		//		"C:\\Users\\chenxi.zhou\\Desktop\\genetic mapping writing-up\\meta_final\\PGSC_DM_v4.03_pseudomolecules.agp2",
		//		"C:\\Users\\chenxi.zhou\\Desktop\\genetic mapping writing-up\\meta_final\\PGSC_DM_v4.03_unanchored_regions_chr00.agp2"};
		
		//String[] faFile = new String[] {
		//		"C:\\Users\\chenxi.zhou\\Desktop\\genetic mapping writing-up\\meta_final\\PGSC_DM_v4.03_pseudomolecules.fasta",
		//		"C:\\Users\\chenxi.zhou\\Desktop\\genetic mapping writing-up\\meta_final\\PGSC_DM_v4.03_unanchored_regions_chr00.fasta"};
		String[] clFile = new String[]{args[0], args[1]};
		String[] faFile = new String[]{args[2], args[3]};
		
		String line;
		String[] s;
		Map<String, int[]> scaffs = new HashMap<String, int[]>();
		Map<String, String> pseudoM = new HashMap<String, String>();
		try {	
			for(int i=0; i<clFile.length; i++) {
				int idx = i==0 ? 5 : 4;
				BufferedReader br = IO.getBufferedReader(clFile[i]);
				while( (line=br.readLine())!=null ) {
					if(line.startsWith("#") || 
							line.length()==0 || 
							line.contains("clone")) 
						continue;
					s = line.split("\\s+");
					String chr = "CH"+s[0].substring(8)+s[idx];
					int a = Integer.parseInt(s[1])-1;
					int b = Integer.parseInt(s[2]);
					scaffs.put(chr, new int[]{a,b});
					
				}
				br.close();
			}
			
			System.err.println("Loading agp done...");
			
			for(int i=0; i<faFile.length; i++) {
				BufferedReader br = IO.getBufferedReader(faFile[i]);
				line = br.readLine();
				while( line!=null ) {
					if(line.startsWith(">")) {
						String key = line.substring(7).toUpperCase();
						StringBuilder pm = new StringBuilder();
						while( (line=br.readLine())!=null &&
								!line.startsWith(">"))
							pm.append(line);
						pseudoM.put(key, pm.toString());
					}
				}
				br.close();
			}
			System.err.println("Loading fasta done...");
			
			//BufferedWriter bw = IO.getGZIPBufferedWriter("C:\\Users\\chenxi.zhou\\Desktop\\genetic mapping writing-up\\meta_final\\PGSC_DM_v4.03.fasta.gz");
			BufferedWriter bw = IO.getGZIPBufferedWriter(args[4]);
			for(String scaff : scaffs.keySet()) {
				System.err.println(scaff);
				bw.write(">");
				bw.write(scaff);
				bw.write("\n");
				int[] x = scaffs.get(scaff);
				bw.write(pseudoM.get(scaff.substring(0, 4)).substring(x[0], x[1]));
				bw.write("\n");
			}
			bw.close();
			System.err.println("Writing scaffolds done...");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static String padding(String s, int l) {
		int i = s.indexOf('_');
		if(i<0) return StringUtils.leftPad(s, l, "0");
		return StringUtils.leftPad(s, l+s.length()-i, "0");
	}
}
