package cz1.appl;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;

import cz1.util.ArgsEngine;
import cz1.util.Utils;

public class VcfToolPlugin {
	
	private final static Logger myLogger = 
			Logger.getLogger(VcfToolPlugin.class);
	private static ArgsEngine myArgsEngine = null;
	static {
		BasicConfigurator.configure();
	}
	
	/**
	public static void main(String[] args) {
		//recode(args[0], args[1]);
		//recode(args[0], args[1], args[2], 
		//		Double.parseDouble(args[3]), 
		//		Integer.parseInt(args[4]));
		//statistic("C:\\Users\\chenxi.zhou\\Desktop\\genetic mapping writing-up\\metafile\\tetrasim.fb.m50.f192.vcf.gz",
		//		"C:\\Users\\chenxi.zhou\\Desktop\\genetic mapping writing-up\\metafile\\Tetra_Trifida_out_alleledose.dat.gz",
		//		"C:\\Users\\chenxi.zhou\\Desktop\\genetic mapping writing-up\\metafile\\tetrasim.f192.log.txt",
		//		4);
		//statistic("C:\\Users\\chenxi.zhou\\Desktop\\genetic mapping writing-up\\metafile\\disim.fb.m50.vcf.gz",
		//		"C:\\Users\\chenxi.zhou\\Desktop\\genetic mapping writing-up\\metafile\\Trifida_D_out_alleledose.dat.gz",
		//		"C:\\Users\\chenxi.zhou\\Desktop\\genetic mapping writing-up\\metafile\\disim.log.txt",
		//		2);
		recode(args[0], args[1]);
	}
	*/
	private static void printUsage() {
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i  Input vcf file\n"
						+ " -o  Output recoded vcf file.\n\n");
	}

	public static void main(String[] args) {
		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}

		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-i", "--input-vcf", true);
			myArgsEngine.add("-o", "--output-vcf", true);
			myArgsEngine.parse(args);
		}

		String in=null, out=null;
		if (myArgsEngine.getBoolean("-i")) {
			in = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify input vcf file.");
		}

		if (myArgsEngine.getBoolean("-o")) {
			out = myArgsEngine.getString("-o");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify output vcf file.");
		}
		
		recode(in, out);
	}
	
	
	private static void statistic(String vcf, String dose, 
			String out, int ploidy) {
		
		BufferedReader br_vcf = Utils.getBufferedReader(vcf);
		BufferedReader br_dose = Utils.getBufferedReader(dose);
		BufferedWriter wr_log = Utils.getBufferedWriter(out);
		String line_snp, line_dosa;
		String[] s, snp, dosa;
		int all=0, fp=0, tp=0, error=0, miss=0;
		int[][] tr = new int[ploidy+1][ploidy+1];;
		try {
			line_dosa = br_dose.readLine();
			s = line_dosa.split("\\s+");
			List<String> header = new ArrayList<String>();
			for(int i=0; i<s.length; i++) 
				header.add(s[i]);
			int[] indices = null;
			boolean shift = true;
			while( (line_snp=br_vcf.readLine())!=null ) {
				//if(line_snp.startsWith("scaffold_chrom2_1535_244830"))
				//	System.out.println(line_snp);
				if( line_snp.startsWith("##"))
					continue;
				if( line_snp.startsWith("#")) {
					s = line_snp.split("\\s+");
					indices =  new int[s.length];
					for(int i=9; i<s.length; i++)
						indices[i] = header.indexOf(s[i]);
					continue;
				}
				all++;
				if(all%1000==0) 
					System.out.println(all+" done."+
							tp+" true positive."+
							error+" errors."+
							miss+" missings.");
				snp = line_snp.split("\\s+");
				s = snp[0].split("_");
				int chrom = Integer.parseInt(
						s[1].replace("chrom", ""));
				int position = Integer.parseInt(s[2])+
						Integer.parseInt(snp[1])-1;
				boolean target = false;
				if(shift) line_dosa = br_dose.readLine();
				while( line_dosa!=null ) {
					//if(line_dosa.startsWith("CHROM2"))
					//	System.out.println("CHROM2");
					int c = compare(line_dosa, 
							chrom, 
							position);
					if(c==0) {
						target = true;
						shift = true;
						break;
					}
					if(c>0) {
						shift = false;
						break;
					}
					line_dosa=br_dose.readLine();
				}
				if(!target) {
					fp++;
				} else {
					tp++;
					//System.out.println(line_dosa);
					//System.out.println(line_snp);
					
					dosa = line_dosa.split("\\s+");
					StringBuilder os = new StringBuilder(); 
					os.append(snp[0]);
					os.append("\t");
					os.append(snp[1]);
					os.append("\t");
					os.append(dosa[0]);
					int e0=0, e1=0;
					int[][] tr0 = new int[ploidy+1][ploidy+1],
							tr1 = new int[ploidy+1][ploidy+1];
					for(int i=9; i<snp.length; i++) {
						s = snp[i].split(":");
						miss += miss(s[0]);
						int[] e = error(s[0], 
								dosa[indices[i]], '0');
						if(e!=null) {
							if(e[0]!=e[1]) e0++;
							tr0[e[0]][e[1]]++;
						}
						e = error(s[0],
								dosa[indices[i]], '1');
						if(e!=null) {
							if(e[0]!=e[1]) e1++;
							tr1[e[0]][e[1]]++;
						}
						os.append("\t");
						os.append(s[0]);
						os.append("(");
						os.append(dosa[indices[i]]);
						os.append(")");
					}
					os.append("\n");
					error += Math.min(e0, e1);
					if(e0<e1) {
						for(int i=0; i<tr0.length; i++)
							for(int j=0; j<tr0[i].length; j++)
								tr[i][j] += tr0[i][j];
					} else {
						for(int i=0; i<tr1.length; i++)
							for(int j=0; j<tr1[i].length; j++)
								tr[i][j] += tr1[i][j];
					}
					wr_log.write(os.toString());
				}
			}
			wr_log.write("##SNPs  "+all+"\n");
			wr_log.write("##FP 	  "+fp+"\n");
			wr_log.write("##TP 	  "+tp+"\n");
			wr_log.write("##MISS  "+miss+"\n");
			wr_log.write("##ERROR "+error+"\n");
			wr_log.write("##TRANS "+"\n");
			wr_log.write("#\t");
			for(int i=0; i<=ploidy; i++)
				wr_log.write("\t\t"+i);
			wr_log.write("\n");
			for(int i=0; i<tr.length; i++) {
				wr_log.write("#"+i+"\t");
				for(int j=0; j<tr[i].length; j++)
					wr_log.write("\t\t"+tr[i][j]);
				wr_log.write("\n");
			}
			br_vcf.close();
			br_dose.close();
			wr_log.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private static int miss(String snp) {
		// TODO Auto-generated method stub
		return snp.startsWith(".") ? 1 : 0;
	}

	private static int[] error(String snp, 
			String dose, char ref) {
		// TODO Auto-generated method stub
		if(snp.startsWith("."))
			return null;
		int d = 0;
		for(int i=0; i<snp.length(); i+=2)
			if(snp.charAt(i)==ref)
				d++;
		return new int[]{Integer.parseInt(dose), d};
	}
	
	private static int compare(String line, 
			int chrom, 
			int position) {
		// TODO Auto-generated method stub
		String s[] = line.split("\\s+")[0].split("\\.");
		int chr = Integer.parseInt(s[0].replace("CHROM", ""));
		int pos = Integer.parseInt(s[1]);
		if(chr<chrom) return -1;
		if(chr==chrom) return pos-position;
		if(chrom>chr) return 1;
		return 0;
	}

	private static void recode(String in, String out, String log, double thres, int nF1) {
		BufferedReader br = Utils.getBufferedReader(in);
		BufferedWriter bw = Utils.getBufferedWriter(out);
		BufferedWriter log_f = Utils.getBufferedWriter(log);
		String snp;
		String[] s, info;
		StringBuilder os = new StringBuilder();
	
		Map<String, String> pl = new HashMap<String, String>();
		pl.put("0/0/0/0", "0,255,255,255,255");
		pl.put("0/0/0/1", "255,0,255,255,255");
		pl.put("0/0/1/1", "255,255,0,255,255");
		pl.put("0/1/1/1", "255,255,255,0,255");
		pl.put("1/1/1/1", "255,255,255,255,0");
		
		try {

			Set<String> f1 = new HashSet<String>();
			for(int i=1; i<=nF1-2; i++)
				f1.add("F"+i);
			f1.add("P1");
			f1.add("P2");
			int[] index = new int[nF1];
			while( (snp=br.readLine())!=null ) {
				if(snp.startsWith("##")) {
					bw.write(snp+"\n");
					continue;
				}
				if(snp.startsWith("#")) {
					os.setLength(0);
					s = snp.split("\\s+");
					for(int i=0; i<8; i++) {
						os.append(s[i]);
						os.append("\t");
					}
					os.append(s[8]);
					int k=0;
					for(int i=9; i<s.length; i++)
						if(f1.contains(s[i])) {
							index[k++] = i;
							os.append("\t");
							os.append(s[i]);
						}
					os.append("\n");
					bw.write(os.toString());
					continue;
				}
				s = snp.split("\\s+");
				if(Double.parseDouble(s[5])<100) {
					log_f.write("QUAL| "+snp+"\n");
					continue;
				}
				double af = Double.parseDouble(s[7].
						split(";")[3].
						split("=|,")[1]);
				if(af>.9 || af<.1) {
					log_f.write("MAF | "+snp+"\n");
					continue;
				}
				
				os.setLength(0);
				os.append(s[0]);
				os.append("\t");
				os.append(s[1]);
				os.append("\t");
				os.append(s[2]);
				os.append("\t");
				os.append('A');
				os.append("\t");
				os.append('B');
				info = s[4].split(",");
				if(info.length>2) {
					log_f.write(snp+"\n");
					continue;
				}
				boolean MNP = info.length>1,
						MN = false;
				for(int i=5; i<8; i++) {
					os.append("\t");
					os.append(s[i]);
				}
				os.append("\tGT:AD:DP:GQ:PL");
				double m=0;
				for(int i : index) {
					os.append("\t");
					info = s[i].split(":");
					if(s[i].startsWith(".") || 
							Double.parseDouble(info[1])<20) {
						os.append("./.:0,0:0:0:0,0,0,0,0");
						m+=1.0;
						continue;
					}
					if(MNP && info[0].indexOf("0")>-1) {
						log_f.write(snp+"\n");
						MN = true;
						break;
					}
										
					String gt = null;
					gt = MNP ? info[0]
							.replaceAll("1", "0")
							.replaceAll("2", "1") : info[0];
					os.append(gt);
					os.append(":");
					os.append(info[3]);
					os.append(",");
					os.append(info[5]);
					os.append(":");
					os.append(info[2]);
					os.append(":");
					os.append((int)Math.round(Double.parseDouble(info[1])));
					os.append(pl.get(gt));
				}
				os.append("\n");
				if(m/nF1>thres && !MN) {
					log_f.write(snp+"\n");
					continue;
				} 
				if(!MN)	bw.write(os.toString());
			}
			br.close();
			bw.close();
			log_f.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static void recode(String in, String out) {
		BufferedReader br = Utils.getBufferedReader(in);
		BufferedWriter bw = Utils.getBufferedWriter(out);
		String snp;
		String[] s, info, qs;
		double[] q;
		Map<String, Integer> fieldMap = new HashMap<String, Integer>();
		StringBuilder os = new StringBuilder();
		try {
			while( (snp=br.readLine())!=null ) {
				if(snp.startsWith("#")) {
					bw.write(snp+"\n");
					continue;
				}
				s = snp.split("\\s+");
				os.setLength(0);
				os.append(s[0].replace(".", "_"));
				os.append("\t");
				os.append(s[1]);
				os.append("\t");
				os.append(s[2]);
				os.append("\t");
				os.append('A');
				os.append("\t");
				os.append('B');
				for(int i=5; i<8; i++) {
					os.append("\t");
					os.append(s[i]);
				}
				if(fieldMap.isEmpty()) {
					info = s[8].split(":");
					for(int i=0; i<info.length; i++) {
						fieldMap.put(info[i], i);
					}
				}
				os.append("\tGT:AD:DP:GQ:PL");
				for(int i=9; i<s.length; i++) {
					os.append("\t");
					if(s[i].startsWith(".")) {
						os.append("./.:0,0:0:33:0,0,0");
						continue;
					}
					info = s[i].split(":");
					os.append(info[fieldMap.get("GT")]);
					os.append(":");
					os.append(info[fieldMap.get("RO")]);
					os.append(",");
					os.append(info[fieldMap.get("AO")]);
					os.append(":");
					os.append(info[fieldMap.get("DP")]);
					os.append(":");
					os.append(info[fieldMap.get("GQ")]);
					os.append(":");
					qs = info[fieldMap.get("GL")].split(",");
					q = new double[qs.length];
					for(int k=0; k<q.length; k++)
						q[k] = -10*Double.parseDouble(qs[k]);
					/***
					int gq = Integer.MAX_VALUE;
					boolean b = false;
					for(int k=0; k<q.length; k++) 
						if(q[k]<gq) {
							if(q[k]!=0 || b) 
								gq = q[k];
							if(q[k]==0) 
								b = true;
						}
					os.append(gq);
					os.append(":");
					*/
					os.append(String.format("%.3f", q[0]));
					for(int k=1; k<q.length; k++) {
						os.append(",");
						os.append(String.format("%.3f", q[k]));
					}
				}
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
