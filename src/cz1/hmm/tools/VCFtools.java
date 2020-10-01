package cz1.hmm.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;


public class VCFtools extends Executor {
	private final static String[] vcf_header = new String[]{"CHROM","POS","ID",
			"REF","ALT","QUAL","FILTER","INFO","FORMAT"};

	private final static double err = 0.01;

	private int min_depth = 0;
	private int max_depth = Integer.MAX_VALUE;
	private int min_qual = 0;
	private double min_maf = 0;
	private double max_missing = 1.0;
	private String vcf_in = null;
	private String vcf_out = null;

	public VCFtools(int min_depth, 
			int max_depth, int min_qual, 
			double min_maf, double max_missing, 
			String vcf_in, String vcf_out) {
		this.min_depth = min_depth;
		this.max_depth = max_depth;
		this.min_qual = min_qual;
		this.min_maf = min_maf;
		this.max_missing = max_missing;
		this.vcf_in = vcf_in;
		this.vcf_out = vcf_out;
	}

	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--vcf            Input VCF file.\n"
						+ " -l/--min-depth      Minimum depth to keep a SNP (0).\n"
						+ " -u/--max-depth      Maximum depth to keep a SNP ("+Integer.MAX_VALUE+").\n"
						+ " -q/--min-qual       Minimum quality to keep a SNP (0).\n"
						+ " -f/--min-maf        Minimum minor allele frequency to keep a SNP (default 0).\n"
						+ " -m/--max-missing    Maximum proportion of missing data to keep a SNP (default 1.0).\n"
						+ " -o/--prefix         Prefix for output VCF file.\n\n");
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub

		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}

		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-i", "--vcf", true);
			myArgsEngine.add("-l", "--min-depth", true);
			myArgsEngine.add("-u", "--max-depth", true);
			myArgsEngine.add("-q", "--min-qual", true);
			myArgsEngine.add("-f", "--min-maf", true);
			myArgsEngine.add("-m", "--max-missing", true);
			myArgsEngine.add("-o", "--prefix", true);
		}
		myArgsEngine.parse(args);

		if (myArgsEngine.getBoolean("-i")) {
			vcf_in = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your VCF file.");
		}

		if (myArgsEngine.getBoolean("-o")) {
			vcf_out = myArgsEngine.getString("-o");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the output file.");
		}

		if (myArgsEngine.getBoolean("-l")) {
			min_depth = Integer.parseInt(myArgsEngine.getString("-l"));
		}

		if (myArgsEngine.getBoolean("-u")) {
			max_depth = Integer.parseInt(myArgsEngine.getString("-u"));
		}

		if (myArgsEngine.getBoolean("-q")) {
			min_qual = Integer.parseInt(myArgsEngine.getString("-q"));
		}

		if (myArgsEngine.getBoolean("-f")) {
			min_maf = Double.parseDouble(myArgsEngine.getString("-f"));
		}

		if (myArgsEngine.getBoolean("-m")) {
			max_missing = Double.parseDouble(myArgsEngine.getString("-m"));
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		// recoding and filtering

		try {
			BufferedReader br = Utils.getBufferedReader(vcf_in);
			BufferedWriter bw = Utils.getBufferedWriter(vcf_out);
			int af_i = -1;
			int dp_i = -1;
			int nS = -1;
			final int header_field = vcf_header.length;
			boolean first_line = true;
			final Map<String, Integer> format_set = new HashMap<String, Integer>();
			List<String> format_arr = new ArrayList<>();
			String format_str = null;
			String snp;
			String[] s, info;
			StringBuilder os = new StringBuilder();
			int filteredByMNP = 0;
			int filteredByQual = 0;
			int filteredByMAF = 0;
			int filteredByDP = 0;
			int filteredByMIS = 0;

			while( (snp=br.readLine())!=null ) {

				if(snp.startsWith("##")) {
					bw.write(snp+"\n");
					continue;
				}
				
				if(snp.startsWith("#")) {
					s = snp.split("\\s+");
					String[] out_str = new String[s.length];
					for(int i=0; i<header_field; i++) {
						if(!s[i].toUpperCase().
								contains(vcf_header[i])) {
							String err_msg = "Corrupted VCF file. The first "+
									header_field+" columns should be,\n";
							for(int j=0; j<header_field; j++)
								err_msg += vcf_header[j]+"\n";
							br.close();
							bw.close();
							throw new RuntimeException(err_msg);
						}
						out_str[i] = s[i];
					}
					
					nS = s.length-header_field;
					
					Set<String> allSamples = new HashSet<>();
					String sample;
					for(int i=header_field; i<s.length; i++) {
						sample = s[i];
						if(sample.contains(":")) {
							sample = sample.split(":")[0];
							myLogger.warn("sample renamed: "+s[i]+" -> "+sample);
						}
						if(allSamples.contains(sample)) {
							br.close();
							bw.close();
							throw new RuntimeException("Duplicate samples "+sample+"!!!");
						}
						allSamples.add(sample);
						out_str[i] = sample;
					}
					
					bw.write(cat(out_str, "\t"));
					bw.write("\n");
					continue;
				}

				s = snp.split("\\s+");

				if(first_line) {
					first_line = false;

					info = s[7].split(";");
					
					for(int i=0; i<info.length; i++) {
						if(info[i].startsWith("AF="))
							af_i = i;
						if(info[i].startsWith("DP="))
							dp_i = i;
					}
					if(af_i==-1)  myLogger.warn("No AF field in VCF file. "
							+ "Filtering by SNP allele frequency disabled.");
					if(dp_i==-1)  myLogger.warn("No DP field in VCF file. "
							+ "Filtering by SNP allele depth disabled.");
					
					info = s[8].split(":");
					
					for(int i=0; i<info.length; i++) {
						format_set.put(info[i], i);
					}

					if(format_set.containsKey("GT")) 
						format_arr.add("GT");

					if(format_set.containsKey("AD")
							|| format_set.containsKey("RO") && 
							format_set.containsKey("AO")) {
						format_arr.add("AD");
					}

					/***
					if(format_set.containsKey("DP") 
							||format_set.containsKey("AD")
							|| format_set.containsKey("RO") && 
							format_set.containsKey("AO")) {
						format_arr.add("DP");
					}

					if(format_set.containsKey("GQ")) {
						format_arr.add("GQ");
					}

					if(format_set.containsKey("PL")
							|| format_set.containsKey("GL")) {
						format_arr.add("PL");
					}
					**/
					
					if(format_arr.isEmpty()) {
						br.close();
						bw.close();
						throw new RuntimeException("NO GT or AD field in the VCF file!!!");
					}
					
					format_str = cat(format_arr, ":");
				}

				// filtered by QUAL
				if(!s[5].equals(".") && Double.parseDouble(s[5])<min_qual) {
					filteredByQual++;
					continue;
				}

				// filtered by MAF
				if(af_i>0) {
					double maf = Double.parseDouble(s[7].
							split(";")[af_i].
							replaceAll("^AF=","").split(",")[0]);
					if(maf<min_maf || maf>1-min_maf) {
						filteredByMAF++;
						continue;
					}
				}

				// filtered by DP
				if(dp_i>0) {
					double d = Double.parseDouble(s[7].
							split(";")[dp_i].
							replaceAll("^DP=", ""));
					if(d>max_depth || d<min_depth) {
						filteredByDP++;
						continue;
					}
				}

				os.setLength(0);
				os.append(s[0].replace('|', '_').replace('.', '_'));
				os.append("\t");
				os.append(s[1]);
				os.append("\t");
				os.append(s[2]);
				info = s[4].split(",");
				// filtered by MNP
				if(info.length>2) {
					filteredByMNP++;
					continue;
				}
				boolean isPotentialMNP = info.length>1;
				
				if(isPotentialMNP) {
					os.append("\t");
					os.append(info[0]);
					os.append("\t");
					os.append(info[1]);
				} else {
					os.append("\t");
					os.append(s[3]);
					os.append("\t");
					os.append(s[4]);
				}
				
				for(int i=5; i<header_field-1; i++) {
					os.append("\t");
					os.append(s[i]);
				}
				os.append("\t");
				os.append(format_str);

				int[] ad = new int[2];
				//int dp = -1;
				int m = 0;
				boolean isMNP = false;
				for(int i=header_field; i<s.length; i++) {
					os.append("\t");
					info = s[i].split(":");
					if(s[i].startsWith(".")) {
						os.append(".");
						m++;
						continue;
					}
					
					if(isPotentialMNP && 
							format_set.containsKey("GT") && 
							info[format_set.get("GT")].contains("0")) {
						isMNP = true;
                        break;
                    }
					
					if(format_set.containsKey("AD")) {
						String[] ad_str = info[format_set.get("AD")].split(",");
						if(isPotentialMNP && 
								!ad_str[0].equals(".") &&
								Integer.parseInt(ad_str[0])>0) {
							isMNP = true;
	                        break;
						}
						int shift = isPotentialMNP ? 1 : 0;
						if(ad_str.length<2) {
                            ad[0] = 0;
                            ad[1] = 0;
                        } else {
                    	    ad[0] = ad_str[0+shift].equals(".")?0:Integer.parseInt(ad_str[0+shift]);
                		    ad[1] = ad_str[1+shift].equals(".")?0:Integer.parseInt(ad_str[1+shift]);
                        }
						//dp = ad[0]+ad[1];
					} else if(format_set.containsKey("RO") &&
							format_set.containsKey("AO")) {
						String ro_str = info[format_set.get("RO")];
						String ao_str = info[format_set.get("AO")];
						
						if(isPotentialMNP) {
							if(!ro_str.equals(".") && Integer.parseInt(ro_str)>0) {
								isMNP = true;
		                        break;
							}
                    		String[] ad_str = ao_str.split(",");
                    		ad[0] = ad_str[0].equals(".")?0:Integer.parseInt(ad_str[0]);
                    		ad[1] = ad_str[1].equals(".")?0:Integer.parseInt(ad_str[1]);
                    	} else {
                    		ad[0] = ro_str.equals(".")?0:Integer.parseInt(ro_str);
                    		ad[1] = ao_str.equals(".")?0:Integer.parseInt(ao_str);
                    	}
						//dp = ad[0]+ad[1];
					}
					//else if(format_set.containsKey("DP")) {
					//	dp = Integer.parseInt(info[format_set.get("DP")]);
					//}

					String[] out_arr = new String[format_arr.size()];
					for(int j=0; j<format_arr.size(); j++) {
						String target_str = format_arr.get(j);
						switch(target_str) {
						case "GT":
							out_arr[j] = info[format_set.get("GT")];
							if(isPotentialMNP) {
								out_arr[j] = out_arr[j].replace('1', '0').replace('2', '1');
							}
							break;
						case "AD":
							out_arr[j] = cat(ad,",");
							break;
						/***
						case "DP":
							out_arr[j] = dp+"";
							break;
						case "GQ":
							out_arr[j] = info[format_set.get("GQ")];
							break;
						case "PL":
						    if(isPotentialMNP) {
						    	br.close();
								bw.close();
								throw new RuntimeException("Cannot rescue MNP with PL/GL field!!!");
						    }
							if(format_set.containsKey("PL")) {
								out_arr[j] = info[format_set.get("PL")];
							} else {
								String[] ss = info[format_set.get("GL")].split(",");
								long[] pl = new long[ss.length];
								for(int k=0; k<pl.length; k++) {
									pl[k] = Math.round(Double.parseDouble(ss[j])*-10);
								}
								out_arr[j] = cat(pl,",");
							}
							break;
						**/
						default:
							br.close();
							bw.close();
							throw new RuntimeException("Undefined VCF field!!!");
						}
					}
					os.append(cat(out_arr, ":"));
				}
				os.append("\n");
				
				// filtered by MNP
				if(isMNP) {
					filteredByMNP++;
					continue;
				}
				
				// filtered by missing data
				if((double)m/nS>max_missing) {
					filteredByMIS++;
					continue;
				}
				
				bw.write(os.toString());
			}
			bw.close();
			br.close();
			
			int filteredTotal = filteredByMNP+filteredByQual+filteredByMAF+filteredByDP+filteredByMIS;
			myLogger.info("#Filtered by Multi-allelic: "+filteredByMNP);
			myLogger.info("#Filtered by Quality      : "+filteredByQual);
			myLogger.info("#Filtered by MAF          : "+filteredByMAF);
			myLogger.info("#Filtered by Allele Depth : "+filteredByDP);
			myLogger.info("#Filtered by Missing      : "+filteredByMIS);
			myLogger.info("--------------------------- ");
			myLogger.info("#Filtered Total           : "+filteredTotal);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static String misSiteStr(final String[] format_arr, final int ploidy) {
		// TODO Auto-generated method stub
		StringBuilder os = new StringBuilder();
		for(String target_str : format_arr) {
			switch(target_str) {
			case "GT":
				os.append(".");
				for(int i=1; i<ploidy; i++)
					os.append("/.");
				break;
			case "AD":
				os.append(":0,0");
				break;
			case "DP":
				os.append(":0");
				break;
			case "GQ":
				os.append(":0");
				break;
			case "PL":
				os.append(":0");
				for(int i=0; i<ploidy; i++)
					os.append(",0");
				break;
			default:
				throw new RuntimeException("Undefined VCF field!!!");
			}
		}
		return os.toString();
	}

	public static double fit(double[] ll, int[] depth, int ploidy) {
		int d = depth[0]+depth[1];
		double maxLL = Double.NEGATIVE_INFINITY;
		for(int i=0; i<ll.length; i++) {
			double pa = ((ploidy-i)*(1-err)+i*err)/ploidy;
			//double nk = Math.log10(nchoosek(d, depth[0]));
			double nk = nchoosekL(d, depth[0]);
			ll[i] = nk+depth[0]*Math.log10(pa)+
					depth[1]*Math.log10(1-pa);
			if(ll[i]>maxLL) maxLL = ll[i];
		}
		for(int i=0; i<ll.length; i++) ll[i] = 10*(maxLL-ll[i]);
		double gq = Double.POSITIVE_INFINITY;
		double pseudo_zero = Double.POSITIVE_INFINITY;
		for (int i = 0; i < ll.length; i++) {
			if (ll[i] < pseudo_zero) {
				gq = pseudo_zero;
				pseudo_zero = ll[i];
			} else if (ll[i] < gq) {
				gq = ll[i];
			}
		}
		return gq;
	}

	public static double[] fit(int[] depth, int ploidy) {
		double[] ll = new double[ploidy+1];
		fit(ll, depth, ploidy);
		return ll;
	}

	public static double[] fitGL(int[] depth, int ploidy) {
		// TODO Auto-generated method stub
		return PL2GL(fit(depth, ploidy));
	}

	public static void fit(double[] ll, String[] genotype, String ref_allele) {
		Arrays.fill(ll, 255);
		if(genotype[0].equals(".")) {
			Arrays.fill(ll, 0);
			return;
		}
		int k = ll.length-1;
		for(String g : genotype)
			if(g.equals(ref_allele))
				k--;
		ll[k] = 0;
	}

	public static double[] fit(String[] genotype, String ref_allele) {
		double[] ll = new double[genotype.length+1];
		fit(ll, genotype, ref_allele);
		return ll;
	}

	public static double[] fitGL(String[] genotype, String ref_allele) {
		// TODO Auto-generated method stub
		return PL2GL(fit(genotype, ref_allele));
	}

	private static double nchoosek(int n, int k) {
		if (k<0||k>n) return 0;
		if (k>n/2) k=n-k;
		double choose = 1.0;
		for (int i=1; i<=k; i++) {
			choose *= (n+1-i);
			choose /= i;
		}
		return choose;
	}

	private static double nchoosekL(int n, int k) {
		// TODO Auto-generated method stub
		if (k<0||k>n) throw new RuntimeException("!!!");
		if (k>n/2) k=n-k;
		double choose = 0;
		for (int i=1; i<=k; i++) {
			choose += Math.log10(n+1-i);
			choose -= Math.log10(i);
		}
		return choose;
	}

	private String cat(List<String> ss, String delimeter) {
		// TODO Auto-generated method stub
		StringBuilder os = new StringBuilder();
		os.append(ss.get(0));
		for(int i=1; i<ss.size(); i++) {
			os.append(delimeter);
			os.append(ss.get(i));
		}
		return os.toString();
	}

	private static String cat(double[] ds, String delimeter) {
		// TODO Auto-generated method stub
		StringBuilder os = new StringBuilder();
		os.append(ds[0]);
		for(int i=1; i<ds.length; i++) {
			os.append(delimeter);
			os.append(ds[i]);
		}
		return os.toString();
	}

	private static String cat(int[] ds, String delimeter) {
		// TODO Auto-generated method stub
		StringBuilder os = new StringBuilder();
		os.append(ds[0]);
		for(int i=1; i<ds.length; i++) {
			os.append(delimeter);
			os.append(ds[i]);
		}
		return os.toString();
	}
	
	private static String cat(long[] ds, String delimeter) {
		// TODO Auto-generated method stub
		StringBuilder os = new StringBuilder();
		os.append(ds[0]);
		for(int i=1; i<ds.length; i++) {
			os.append(delimeter);
			os.append(ds[i]);
		}
		return os.toString();
	}
	
	private static String cat(String[] ds, String delimeter) {
		// TODO Auto-generated method stub
		StringBuilder os = new StringBuilder();
		os.append(ds[0]);
		for(int i=1; i<ds.length; i++) {
			os.append(delimeter);
			os.append(ds[i]);
		}
		return os.toString();
	}
	
	private static String uniGT(double[] ll) {
		// TODO Auto-generated method stub
		int l = ll.length-1;
		int[] g = new int[l];
		for(int i=0; i<l; i++) 
			if(ll[i]!=0)
				g[l-1-i]=1;
			else break;
		return cat(g, "/");
	}

	public static String uniPL(String gt_str, int ploidy) {
		// TODO Auto-generated method stub
		int[] g = new int[ploidy+1];
		Arrays.fill(g, 255);
		int k = 0;
		for(int i=0; i<gt_str.length(); i+=2) {
			k += gt_str.charAt(i)=='0' ? 1 : 0;
		}
		g[k] = 0;
		return cat(g, ",");
	}

	public static double[] PL2GL(double[] pl) {
		// TODO Auto-generated method stub
		double[] ll = new double[pl.length];
		for(int i=0; i!=ll.length; i++)
			ll[i] = Math.pow(10,-pl[i]/10);
		return ll;
	}

}
