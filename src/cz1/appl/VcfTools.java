package cz1.appl;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;

import cz1.util.IO;

public class VcfTools {

	private final static int ploidy = 4;
	
	public static void main(String[] args) {
		//homozygotes("C:\\Users\\chenxi.zhou\\desktop\\potato\\"
		//		+ "PGSC_DM_v403_8303SNPs_Stirlingx12601ab1.vcf", 
		//		"Stirling:12601ab1");
		//statistics("C:\\Users\\chenxi.zhou\\Desktop\\genetic mapping writing-up\\metafile\\tetrasim.fb.m50.f192.vcf.gz",
		//		"C:\\Users\\chenxi.zhou\\Desktop\\genetic mapping writing-up\\metafile\\Tetra_Trifida_out_alleledose.dat.gz",
		//		"C:\\Users\\chenxi.zhou\\Desktop\\genetic mapping writing-up\\metafile\\tetrasim.f192.log.txt",
		//		4);
		//statistics("C:\\Users\\chenxi.zhou\\Desktop\\putty\\vd_comp\\out.recode.vcf",
		//		"C:\\Users\\chenxi.zhou\\Desktop\\putty\\vd_comp\\Trifida_D_out_alleledose.dat.gz",
		//		"C:\\Users\\chenxi.zhou\\Desktop\\putty\\vd_comp\\out.recode.txt",
		//		2);
		vcf2Onemap(args[0], args[1].split(":"), args[2]);
		
	}
	
	private static void vcf2Onemap(String vcfFile, String[] parents_str, String out) {
		try {
			
			int[] parents = new int[2];
			final BufferedReader br = IO.getBufferedReader(vcfFile);
			String[] s;
			String line;
			int[] f1 = null;
			String scaffold = null;
			List<String[]> scaff = null;
			int nF1 = -1;
			while( (line=br.readLine())!=null ) {
				if(line.startsWith("##")) continue;
				if(line.startsWith("#")) {
					s = line.split("\\s+");
					List<Integer> f1_list = new ArrayList<Integer>();
					for(int i=9; i<s.length; i++) {
						if(s[i].equals(parents_str[0]))
							parents[0] = i;
						else if(s[i].equals(parents_str[1]))
							parents[1] = i;
						else f1_list.add(i);
					}
					Integer[] f1I = new Integer[f1_list.size()];
					f1_list.toArray(f1I);
					f1 = ArrayUtils.toPrimitive(f1I);
					nF1 = f1.length;
					continue;
				}
				s = line.split("\\s+");
				if(scaffold==null) {
					scaff = new ArrayList<String[]>();
					scaffold = s[0];
					scaff.add(s);
					continue;
				} else if(scaffold.equals(s[0])) {
					scaff.add(s);
					continue;
				} else {
					System.out.println(scaffold);
					BufferedWriter bw = IO.getBufferedWriter(out+"/"+scaffold+".dat");
					StringBuilder oos = new StringBuilder();
					int n=0;
					String segregation = null;
					for(String[] snp : scaff) {
						
						if(!snp[parents[0]].startsWith("./.") &&
								!snp[parents[1]].startsWith("./.")) {
							String p1 = snp[parents[0]].startsWith("0/0") ? "a" : 
								((snp[parents[0]].startsWith("0/1")||
										snp[parents[0]].startsWith("1/0")) ? "ab" : "b");
							String p2 = snp[parents[1]].startsWith("0/0") ? "a" : 
								((snp[parents[1]].startsWith("0/1")||
										snp[parents[1]].startsWith("1/0")) ? "ab" : "b");
							if(p1.equals("a")&&
									p2.equals("ab")||
									p1.equals("b")&&
									p2.equals("ab")) {
								oos.append("*"+snp[0]+"_"+snp[1]+" ");
								oos.append("D2.15\t");
								segregation = "D2.15";
							} else if(p1.equals("ab")&&
									p2.equals("a")||
									p1.equals("ab")&&
									p2.equals("b")) {
								oos.append("*"+snp[0]+"_"+snp[1]+" ");
								oos.append("D1.10\t");
								segregation = "D1.10";
							} else if(p1.equals("ab")&&
									p2.equals("ab")) {
								oos.append("*"+snp[0]+"_"+snp[1]+" ");
								oos.append("B3.7\t");
								segregation = "B3.7";
							} else {
								continue;
							}
							n++;
							String[] f1g = new String[nF1];
							String g;
							for(int i=0; i<f1.length; i++) {
								g = snp[f1[i]].substring(0, 3);
								if(g.equals("./.")) { f1g[i] = "-"; continue; }
								if(p1.equals("a")&&p2.equals("ab") || 
										p1.equals("ab")&&p2.equals("a")) {
									if(g.equals("b")) { 
										f1g[i] = "-"; 
										continue; 
									}
								}
								if(p1.equals("b")&&p2.equals("ab") || 
										p1.equals("ab")&&p2.equals("b")) {
									if(g.equals("a")) { 
										f1g[i] = "-"; 
										continue; 
									}
								}
								
								switch(segregation) {
								case "D2.15":
								case "D1.10":
									switch(g) {
									case "0/0":
									case "1/1":
										f1g[i] = "a";
										break;
									case "0/1":
									case "1/0":
										f1g[i] = "ab";
										break;
									default:
										throw new RuntimeException("!!!");
									}
									break;
								case "B3.7":
									switch(g) {
									case "0/0":
										f1g[i] = "a";
										break;
									case "0/1":
									case "1/0":
										f1g[i] = "ab";
										break;
									case "1/1":
										f1g[i] = "b";
										break;
									default:
										throw new RuntimeException("!!!");
									}
									break;
								}
							}
							oos.append(cat(f1g,","));
							oos.append("\n");
						}
					}
					bw.write(nF1+" "+n+"\n");
					bw.write(oos.toString());
					bw.close();
					
					scaff = new ArrayList<String[]>();
					scaffold = s[0];
					scaff.add(s);					
				}
			}
			
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static Object cat(String[] f1g, String delimeter) {
		// TODO Auto-generated method stub
		StringBuilder oos = new StringBuilder();
		oos.append(f1g[0]);
		for(int i=1; i<f1g.length; i++) {
			oos.append(",");
			oos.append(f1g[i]);
		}
		return oos.toString();
	}

	private static void homozygotes(String vcfFile, String parents) {
		
		try {
			final BufferedReader br = IO.getBufferedReader(vcfFile);
			int[] index_p = new int[2];
			String[] s;
			String line = br.readLine(); 
			while( line!=null ) {
				if(line.startsWith("##")) {
					line = br.readLine();
					continue;
				}
				if(line.startsWith("#")) {
					String[] p = parents.split(":");
					s = line.split("\\s+");
					for(int i=0; i<s.length; i++) {
						if(s[i].equals(p[0]))
							index_p[0] = i;
						if(s[i].equals(p[1]))
							index_p[1] = i;
					}
					line = br.readLine();
					continue;
				}
				s = line.split("\\s+");
				String scaff = s[0];
				int n = 1, h = homo(s, index_p);
				while( (line=br.readLine())!=null ) {
					s = line.split("\\s+");
					if( !s[0].equals(scaff) )
						break;
					n++;
					h+=homo(s,index_p);
				}
				//if(((double)h)/n<=.5 && n>=10) 
				//	System.out.println(scaff+"\t"+(n-h));
				System.out.println(scaff+"\t"+h+"/"+n);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	private static int homo(String[] s, int[] index_p) {
		// TODO Auto-generated method stub
		int dosage_p0 = dosage(s[index_p[0]].split(":")[0]),
				dosage_p1 = dosage(s[index_p[1]].split(":")[0]);
		if(dosage_p0==-1 && dosage_p1==-1) 
			return 0;
		if(dosage_p0!=-1 && dosage_p1!=-1) {
			if(dosage_p0!=ploidy && dosage_p0!=0 &&
					dosage_p1!=ploidy && dosage_p1!=0)
				return 0;
			if(dosage_p0==ploidy || dosage_p0==0 ||
					dosage_p1==ploidy || dosage_p1==0)
				return 1;
		}
		if(dosage_p0==-1) return guess(s, index_p, dosage_p1);
		else return guess(s, index_p, dosage_p0);
	}

	private static int guess(String[] s, int[] index_p, 
			int dosage_known) {
		// TODO Auto-generated method stub
		int H = 0;
		int dosage_sum = 0;
		for(int i=0; i<s.length; i++) {
			if(i!=index_p[0]&&i!=index_p[1]) {
				int dosage = dosage(s[i].split(":")[0]);
				if(dosage==-1) continue;
				dosage_sum += dosage;
				H += ploidy;
			}
		}
		double frac = ((double) dosage_sum)/H;
		if( Math.abs(frac-((double)dosage_known)/2/ploidy)<0.01 || 
				Math.abs(frac-((double)dosage_known)/2/ploidy-.5)<0.01)
			return 1;
		else return 0;
	}

	private static int dosage(String genotype) {
		int dosage = 0;
		if(genotype.startsWith("."))
			return -1;
		for(int i=0; i<genotype.length(); i+=2)
			dosage += genotype.charAt(i)=='0' ? 1 : 0;
		return dosage;
	}
	
	private static void statistics(String vcf, String dose, 
			String out, int ploidy) {
		
		
		String line_snp, line_dosa;
		String[] s, snp, dosa;
		int all=0, fp=0, tp=0, error=0, miss=0;
		int[][] tr = new int[ploidy+1][ploidy+1];;
		try {
			BufferedReader br_vcf = IO.getBufferedReader(vcf);
			BufferedReader br_dose = IO.getBufferedReader(dose);
			BufferedWriter wr_log = IO.getBufferedWriter(out);
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
				snp = line_snp.split("\\s+");
				if( snp[3].length()>1 ||
						snp[4].length()>1 ) continue;
				
				all++;
				if(all%1000==0) 
					System.out.println(all+" done."+
							tp+" true positive."+
							error+" errors."+
							miss+" missings.");
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
}
