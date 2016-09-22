package cz1.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.ArrayUtils;

public class BioString {

	public static void main(String[] args) throws IOException {
		BioString bs = new BioString();
		//bs.scafflod2contig("C:\\Users\\chenxi.zhou\\Desktop\\"
		//		+ "putty\\v2_collapsed.fa.gz", 
		//		"C:\\Users\\chenxi.zhou\\Desktop\\"
		//		+ "v2_collapsed.ctg.gz");
		
		//bs.scafflod2contig("C:\\Users\\chenxi.zhou\\Desktop\\"
		//		+ "putty\\NSP306_trifida_v2.scf.gz", 
		//		"C:\\Users\\chenxi.zhou\\Desktop\\"
		//		+ "NSP306_trifida_v2.ctg.gz");
		//bs.contigstats("C:\\Users\\chenxi.zhou\\Desktop\\"
		//		+ "NSP306_trifida_v2.ctg.gz", 
		//		"C:\\Users\\chenxi.zhou\\Desktop\\"
		//s		+ "NSP306_trifida_v2.stat.gz");
		/***
		bs.rename("C:\\Users\\chenxi.zhou\\Desktop\\"
				+ "putty\\snp_trifida_1.txt", 
				"C:\\Users\\chenxi.zhou\\Desktop\\"
				+ "putty\\map_trifida_1.txt",
				"C:\\Users\\chenxi.zhou\\Desktop\\"
				+ "putty\\trifida_final_contig_all.txt",
				"C:\\Users\\chenxi.zhou\\Desktop\\"
				+ "putty\\snp_trifida_2.txt", 
				"C:\\Users\\chenxi.zhou\\Desktop\\"
				+ "putty\\map_trifida_2.txt");
		**/
		
		//bs.contigdistance("C:\\Users\\chenxi.zhou\\Desktop\\"
		//		+ "NSP306_trifida_v2.stat.gz",
		//		"C:\\Users\\chenxi.zhou\\Desktop\\"
		//		+ "putty\\snp_trifida.txt",
		//		"C:\\Users\\chenxi.zhou\\Desktop\\"
		//		+ "putty\\map_trifida.txt",
		//		"C:\\Users\\chenxi.zhou\\Desktop\\"
		//		+ "contig_distance.txt");
		
		bs.contigdistance("C:\\Users\\chenxi.zhou\\Desktop\\"
				+ "NSP306_trifida_v2.stat.gz",
				"C:\\Users\\chenxi.zhou\\Desktop\\"
				+ "new_tmp\\dlpbwa2m.snp",
				"C:\\Users\\chenxi.zhou\\Desktop\\"
				+ "new_tmp\\rf_for_r_dlp_final_cz1_mt_1.map",
				"C:\\Users\\chenxi.zhou\\Desktop\\"
				+ "new_tmp\\contig_distance_dlp_mt_1.txt");
		
		/**
		BufferedReader br = IO.getBufferedReader(
				"C:\\Users\\chenxi.zhou\\Desktop\\v2_collapsed.fa.gz");
		BufferedWriter bw = IO.getBufferedWriter(
				"C:\\Users\\chenxi.zhou\\Desktop\\v2_collapsed.fa");
		bw.write(br.readLine()+"\n");
		String line;
		StringBuilder sb = new StringBuilder(10000000);
		int L = 60;
		while( (line=br.readLine())!=null ) {
			sb.append(line);
			if(sb.length()>=9000000) { 
				for(int i=0; i<150000; i++) 
					bw.write(sb.substring(i*L, L*(i+1))+"\n");
				sb.delete(0, 9000000);
			}
		}
		for(int i=0; i<sb.length()/L; i++) 
			bw.write(sb.substring(i*L, L*(i+1))+"\n");
		bw.write(sb.substring(sb.length()/L*L)+"\n");
		br.close();
		bw.close();
		*/
		
		/**
		System.out.println("Trifida_004	Trifida_0042".replaceAll("^"+"Trifida_0{0,}", ""));
		BufferedReader br = IO.getBufferedReader(
				"C:\\Users\\chenxi.zhou\\Desktop\\v2_collapsed.ctg.gz");
		int n=5511917, l=64;
		while(--n>0) br.readLine();
		
		String ref = "";
		String line;
		for(int i=0; i<l; i++) {
			line = br.readLine();
			ref += line.split("\\s+")[2];
		}
		String seq="TGTAAAGGCTGTAGATATCAATCGCCTTGTACTGTTTCAGCAGTTCATTCACAGCTTGGCTGCA";
		br.close();
		System.out.println("ref "+ref);
		System.out.println("seq "+seq);
		int sc = 0;
		for(int i=0; i<ref.length(); i++)
			sc += seq.charAt(i)==ref.charAt(i) ? 0 :1;
		System.out.println(sc);
		*/
	}

	private void contigdistance(String stat,
			String snp,
			String map,
			String output) throws IOException {
		// TODO Auto-generated method stub
		BufferedReader br = IO.getBufferedReader(snp);
		Map<String, int[]> snpMap = new HashMap<String, int[]>();
		String[] s;
		String line = br.readLine();
		String scaffold = null;
		List<Integer> position = new ArrayList<Integer>();
		while( line!=null ) {
			s = line.split("_");
			if( !s[0].equals(scaffold) ) {
				if(scaffold!=null) {
					int l = scaffold.length();
					for(int i=0; i<5-l; i++)
						scaffold = "0"+scaffold;
					snpMap.put(scaffold, ArrayUtils.toPrimitive(
							position.toArray(
									new Integer[position.size()])));
				}
				scaffold = s[0];
				position = new ArrayList<Integer>();
			}
			position.add(Integer.parseInt(s[1]));
			line = br.readLine();
		}
		br.close();
		
		br = IO.getBufferedReader(map);
		Map<String, double[]> mapMap = new HashMap<String, double[]>();
		while( (line=br.readLine())!=null ) {
			if(line.startsWith("#"))
				continue;
			s = line.split("\\s+");
			scaffold = s[0].replace("*", "");
			int l = scaffold.length();
			for(int i=0; i<5-l; i++)
				scaffold = "0"+scaffold;
			String[] s5 = s[5].split(",");
			double[] rf = new double[s5.length];
			for(int i=0; i<rf.length; i++)
				rf[i] = Double.parseDouble(s5[i]);
			mapMap.put(scaffold, rf);
		}
		br.close();
		
		br = IO.getBufferedReader(stat);
		Map<String, Contig> contigMap = new HashMap<String, Contig>();
		line = br.readLine();
		scaffold = null;
		List<Contig> contigs = new ArrayList<Contig>();
		while( line!=null ) {
			s = line.split("_|\\s+");
			if( !s[0].equals(scaffold) ) {
				if(scaffold!=null) {
					int l = scaffold.length();
					for(int i=0; i<5-l; i++)
						scaffold = "0"+scaffold;
					contigMap.put(scaffold, contigBTree(
							contigs.toArray(new Contig[contigs.size()])));
				}
				scaffold = s[0];
				contigs = new ArrayList<Contig>();
			}
			contigs.add(new Contig(s[0]+"_"+s[1],
					new int[]{Integer.parseInt(s[2]),Integer.parseInt(s[3])},
					null,
					null));
			line = br.readLine();
		}
		br.close();
		
		BufferedWriter bw = IO.getBufferedWriter(output);
		for(String scaf : snpMap.keySet()) {
			int[] posi = snpMap.get(scaf);
			double[] rf = mapMap.get(scaf);
			if(rf==null) continue;
			for(int i=0; i<posi.length; i++) {
				//System.out.println(posi[i]);
				String contig = 
						search(contigMap.get(scaf), posi[i]);
				if(i<rf.length)
					bw.write(contig+"\t"+rf[i]+
							"\t"+posi[i]+"\t"+scaf+"\n");
				else
					bw.write(contig+"\tNA\t"+posi[i]+
							"\t"+scaf+"\n");
			}
		}
		bw.close();
	}
	
	private String search(Contig contig, int i) {
		// TODO Auto-generated method stub
		if(contig==null) return "NULL";
		int[] range = contig.range;
		if(i>=range[0] && i<=range[1])
			return contig.id;
		if(i<range[0]) 
			return search(contig.left, i);
		if(i>range[1]) 
			return search(contig.right, i);
		return null;
	}

	private Contig contigBTree(Contig[] array) {
		// TODO Auto-generated method stub
		if(array.length==1) {
			array[0].left = null;
			array[0].right = null;
			return array[0];
		} 
		int b = array.length/2;
		array[b].left = b-1<0 ? null : 
			contigBTree(ArrayUtils.subarray(
					array, 0, b));
		array[b].right = b+1>array.length-1 ? null : 
			contigBTree(ArrayUtils.subarray(
					array, b+1, array.length));
		return array[b];
	}

	private void contigstats(String input, String output) 
			throws IOException {
		// TODO Auto-generated method stub
		BufferedReader br = IO.getBufferedReader(input);
		BufferedWriter bw = IO.getGZIPBufferedWriter(output);

		String scaffold=null, contig = null;
		String[] s;
		String line = br.readLine();
		String start = null, end = null;
		while( line!=null ) {

			s = line.split("\\s+|_");
			if( s[1].equals("0")) {
				line = br.readLine();
				continue;
			}
			if(  !s[0].equals(scaffold)
					|| !s[1].equals(contig) ) {
				
				if( contig!=null ) {
					bw.write(scaffold+"_"+contig+
							"\t"+start+"\t"+end+"\n");
					System.out.println(scaffold+"_"+contig+
							"\t"+start+"\t"+end);
				}
				scaffold = s[0];
				contig = s[1];
				start = s[2];
				end = s[2];
			} else end = s[2];
			line = br.readLine();
		}
		br.close();
		bw.close();
	}

	public void scafflod2contig(String input, String output) 
			throws IOException {
		BufferedReader br = IO.getBufferedReader(input);
		BufferedWriter bw = IO.getGZIPBufferedWriter(output);

		String scaffold = null;
		int contig = -1, position = -1;
		boolean next = false;
		String line = br.readLine();
		while( line!=null ) {
			if(line.startsWith(">")) {
				scaffold = line.replace(">Trifida_", "");
				System.err.println(scaffold);
				contig = 1;
				position = 1;
				next = false;
			}
			while( (line=br.readLine())!=null ) {
				if(line.startsWith(">"))
					break;
				char[] nuls = line.toCharArray();
				for(char nul : nuls) {
					if(nul=='N') {
						bw.write(scaffold+
								"_0\t"+position+"\tN\n");
						next = true;
					} else {
						if(next) {
							contig++;
							next = false;
						}
						bw.write(scaffold+"_"+contig+
								"\t"+position+"\t"+nul+"\n");
					}
					position++;
				}
			}
		}
		br.close();
		bw.close();
	}
	
	private class Contig {
		private Contig left;
		private Contig right;
		private final String id;
		private final int[] range;
		
		public Contig(String id,
				int[] range,
				Contig left,
				Contig right) {
			this.id = id;
			this.range = range;
			this.left = left;
			this.right = right;
		}
		
		public boolean leaf() {
			return this.left==null &&
					this.right==null;
		}
	}
	
	private void rename(String vcf, String contig,
			String vcf_new) throws IOException {
		BufferedReader br = IO.getBufferedReader(contig);
		List<String> contigs = new ArrayList<String>();
		contigs.add("DUMMY");
		String line;
		while( (line=br.readLine())!=null ) 
			contigs.add(line.replaceAll("^Trifida_0{0,}", ""));
		br.close();
		br = IO.getBufferedReader(vcf);
		BufferedWriter bw = IO.getBufferedWriter(vcf_new);
		String[] s;
		Map<String, String> all_snps = new HashMap<String, String>();
		while( (line=br.readLine())!=null ) {
			if(line.startsWith("#")) { 
				bw.write(line+"\n");
				continue;
			}
			s = line.split("\\s+");
			String new_contig = contigs.get(Integer.parseInt(s[0]));
			line = line.replaceAll("^"+s[0], new_contig);
			line = line.replace("S"+s[0]+"_", "S"+new_contig+"_");
			line = line+"\n";
			if(all_snps.get(new_contig)==null)
				all_snps.put(new_contig, line);
			else
				all_snps.put(new_contig, all_snps.get(new_contig)+line);
		}
		br.close();
		for(int i=1; i<contigs.size(); i++) 
			if( (line=all_snps.get(contigs.get(i))) !=null )
					bw.write(line);
		bw.close();
		return;
	}
	
	private void rename(String snp, String map, String contig,
			String snp_new, String map_new) 
			throws IOException {
		BufferedReader br = IO.getBufferedReader(contig);
		List<String> contigs = new ArrayList<String>();
		contigs.add("DUMMY");
		String line;
		while( (line=br.readLine())!=null ) 
			contigs.add(line.replaceAll("^Trifida_0{0,}", ""));
		br.close();
		br = IO.getBufferedReader(snp);
		BufferedWriter bw = IO.getBufferedWriter(snp_new);
		String[] s;
		while( (line=br.readLine())!=null ) {
			s = line.split("_");
			bw.write(contigs.get(
					Integer.parseInt(s[0]))+"_"+s[1]+"\n");
		}
		bw.close();
		br.close();
		br = IO.getBufferedReader(map);
		bw = IO.getBufferedWriter(map_new);
		while( (line=br.readLine())!=null ) {
			if(line.startsWith("#")) { 
				bw.write(line+"\n");
				continue;
			}
			s = line.split("\\s+");
			bw.write(line.replace(s[0],
					"*"+contigs.get(
							Integer.parseInt(s[0].substring(1))))
							+"\n");
		}
		bw.close();
		br.close();
		return;
	}
}












