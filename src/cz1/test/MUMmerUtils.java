package cz1.test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import cz1.ngs.model.Sequence;

public class MUMmerUtils {

	private static void variantParser(String snp_file, String out_file) {
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(new File(snp_file))));
			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(out_file))));
			for(int i=0; i<5; i++) br.readLine();
			int P01, P02, P1, P2, frm, tmp_P1, tmp_P2;
			String line = br.readLine();
			String[] s = line.trim().split("\\s+");
			final int sLen = s.length-1;
			String allele1, allele2, sequence;
			boolean d1, d2;
			while( line!=null ) {
				P01 = Integer.parseInt(s[0]);
				P02 = Integer.parseInt(s[3]);
				frm = Integer.parseInt(s[sLen-2]);
				sequence = s[sLen];
				P1  = P01;
				P2  = P02;
				allele1 = s[1];
				allele2 = s[2];
				d1 = allele1.equals(".");
				d2 = allele2.equals(".");
				while( (line=br.readLine())!=null ) {
					s = line.trim().split("\\s+");
					tmp_P1 = Integer.parseInt(s[0]);
					tmp_P2 = Integer.parseInt(s[3]);
					if(d1&&tmp_P1==P1&&s[1].equals(".")&&Math.abs(tmp_P2-P2)==1) {
						allele2 += s[2];
						P2 = tmp_P2;
					} else if(d2&&tmp_P2==P2&&s[2].equals(".")&&Math.abs(tmp_P1-P1)==1) {
						allele1 += s[1];
						P1 = tmp_P1;
					}
					else break;
				}
				
				if(!allele1.contains("N")&&!allele2.contains("N")) {
					if(P02>P2) {
						P02 = P2;
						allele2 = new StringBuilder(allele2).reverse().toString();
					}
					if(frm<0) allele2 = cmp(allele2);
					bw.write(sequence+"\t"+P01+"\t"+P02+"\t"+allele1+"\t"+allele2+"\t"+frm+"\n");
				}
			}
			bw.close();
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private final static Map<Character, Character> cmp = new HashMap<Character, Character>();
	static {
		cmp.put('A', 'T');
		cmp.put('C', 'G');
		cmp.put('G', 'C');
		cmp.put('T', 'A');
		cmp.put('N', 'N');
		cmp.put('a', 't');
		cmp.put('c', 'g');
		cmp.put('g', 'c');
		cmp.put('t', 'a');
		cmp.put('n', 'n');
	}
	
	private static String cmp(String seq) {
		if(seq.equals(".")) return seq;
		StringBuilder seq_buf = new StringBuilder(seq);
		int n = seq.length();
		for(int i=0; i!=n; i++) 
			seq_buf.setCharAt(i, cmp.get(seq_buf.charAt(i)));
		return seq_buf.toString();
	}
	
	private final static int min_gap = 1;
	
	private static void variantMerger(String var_in, String ref_in1, String ref_in2, String var_out) {
		// TODO Auto-generated method stub
		final Map<String, Sequence> refSequence1 = Sequence.parseFastaFileAsMap(ref_in1);
		final Map<String, Sequence> refSequence2 = Sequence.parseFastaFileAsMap(ref_in2);
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(new File(var_in))));
			BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(var_out))));
			
			String line = br.readLine();
			String[] s = line.trim().split("\\s+");
			String seq_sn = s[0];
			String seq1 = refSequence1.get(seq_sn).seq_str();
			String seq2 = refSequence2.get(seq_sn).seq_str();
			int tmp_pos1, tmp_pos2, frm;
			final List<Integer> pos1 = new ArrayList<Integer>();
			final List<Integer> pos2 = new ArrayList<Integer>();
			final List<String> var1 = new ArrayList<String>();
			final List<String> var2 = new ArrayList<String>();
			String var, allele1, allele2;
			int pos, p10, p20, p11, p21, n;
			
			while(line!=null) {
				pos1.clear();
				pos2.clear();
				var1.clear();
				var2.clear();
				pos1.add(Integer.parseInt(s[1]));
				pos2.add(Integer.parseInt(s[2]));
				var1.add(s[3]);
				var2.add(s[4]);
				frm = Integer.parseInt(s[5]);
				
				while( (line=br.readLine())!=null ) {
					s = line.trim().split("\\s+");
					tmp_pos1 = Integer.parseInt(s[1]);
					tmp_pos2 = Integer.parseInt(s[2]);
					
					if( Math.abs(pos1.get(pos1.size()-1)-tmp_pos1)-s[3].length()<=min_gap && 
							Math.abs(pos2.get(pos2.size()-1)-tmp_pos2)-s[4].length()<=min_gap &&
							Integer.parseInt(s[5])==frm)  {
						pos1.add(tmp_pos1);
						pos2.add(tmp_pos2);
						var1.add(s[3]);
						var2.add(s[4]);
					} else
						break;
				}
				
				n = pos1.size();
				
				if(pos1.get(0)>pos1.get(n-1)) {
					Collections.reverse(pos1);
					Collections.reverse(var1);
				}
				
				if(pos2.get(0)>pos2.get(n-1)) {
					Collections.reverse(pos2);
					Collections.reverse(var2);
				}
				
				for(int i=0; i<n; i++) {
					pos = pos1.get(i);
					var = var1.get(i);
					if(var.equals(".")) continue;
					if(!seq1.substring(pos-1, pos+var.length()-1).toUpperCase().equals(var)) 
						throw new RuntimeException("1. "+pos+","+var);
				}
				
				for(int i=0; i<n; i++) {
					pos = pos2.get(i);
					var = var2.get(i);
					if(var.equals(".")) continue;
					if(!seq2.substring(pos-1, pos+var.length()-1).toUpperCase().equals(var)) 
						throw new RuntimeException("2. "+pos+","+var);
				}
				
				p10 = pos1.get(0);
				p11 = pos1.get(n-1)+
						(var1.get(n-1).equals(".")?0:var1.get(n-1).length());
				
				p20 = pos2.get(0);
				p21 = pos2.get(n-1)+
						(var2.get(n-1).equals(".")?0:var2.get(n-1).length());
				
				if(p10>p11||p20>p21) throw new RuntimeException("!!!");
				
				allele1 = p10==p11?".":seq1.substring(p10-1,p11-1).toUpperCase();
				allele2 = p20==p21?".":seq2.substring(p20-1,p21-1).toUpperCase();
			
				bw.write(seq_sn+"\t"+p10+"\t"+p20+"\t"+allele1+"\t"+allele2+"\t"+frm+"\n");
			}
			
			bw.close();
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
	}
	
	public static void main(String[] args) {
		switch(args[0].toLowerCase()) {
		case "parser":
			// args: mummer snps
			// args: output file
			variantParser(args[1], args[2]);	
			break;
		case "merger":
			// args: input variants
			// args: reference seq1
			// args: reference seq2
			// args: output file
			variantMerger(args[1], args[2], args[3], args[4]);
			break;
		default:
			throw new RuntimeException("!!!");	
		}
	}
}
