package cz1.hmm.tools;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.HashMap;
import java.util.List;
import java.util.Set;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;

import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Executor;
import cz1.util.Utils;


public class VCFtools extends Executor {
	
	private int ploidy = -1;
	private int min_depth = 0;
	private int max_depth = Integer.MAX_VALUE;
	private int min_qual = 0;
	private double min_maf = 0.1;
	private double max_missing = 0.5;
	private String vcf_in = null;
	private String vcf_out = null;
	
	public VCFtools(int ploidy, int min_depth, 
			int max_depth, int min_qual, 
			double min_maf, double max_missing, 
			String vcf_in, String vcf_out) {
		this.ploidy = ploidy;
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
						+ " -i/--vcf			Input VCF file.\n"
						+ " -p/--ploidy			Real ploidy of genome. "
						+ "						NOTE: You may call variant as diploid and the program will "
						+ "							  fit a binomial model to call genotypes and genotype "
						+ "							  qualities from allele depth.\n"
						+ " -l/--min-depth		Minimum depth to keep a SNP (DP).\n"
						+ " -u/--max-depth		Maximum depth to keep a SNP (DP).\n"
						+ " -q/--min-qual  		Minimum quality to keep a SNP (QUAL).\n"
						+ " -f/--min-maf		Minimum minor allele frequency to keep a SNP (default 0.1).\n"
						+ " -m/--max-missing	Maximum proportion of missing data to keep a SNP (default 0.5).\n"
						+ " -o/--prefix			Prefix for output VCF file.\n\n");
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
			myArgsEngine.add("-p", "--ploidy", true);
			myArgsEngine.add("-u", "--max-depth", true);
			myArgsEngine.add("-q", "--min-qual", true);
			myArgsEngine.add("-f", "--min-maf", true);
			myArgsEngine.add("-m", "--max-missing", true);
			myArgsEngine.add("-o", "--prefix", true);
			myArgsEngine.parse(args);
		}

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
		
		if (myArgsEngine.getBoolean("-p")) {
			ploidy = Integer.parseInt(myArgsEngine.getString("-p"));
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
			int ploidy_obs = -1;
			final int header_field = Constants._vcf_header.length;
			boolean switchoff = true;
			final Map<String, Integer> format_set = new HashMap<String, Integer>();
			String[] format_arr = null;
			String format_str = "";
			String snp;
	        String[] s, info;
	        StringBuilder os = new StringBuilder();
	        boolean contain_ad_info = false;
			while( (snp=br.readLine())!=null ) {
				
				if(snp.startsWith("##")) {
					bw.write(snp+"\n");
					continue;
				}
				if(snp.startsWith("#")) {
					s = snp.split("\\s+");
					for(int i=0; i<header_field; i++) {
						if(!s[Constants._vcf_header_i[i]].toUpperCase().
								contains(Constants._vcf_header[i])) {
							String err_msg = "Corrupted VCF file. The first "
						+header_field+" columns should be,\n";
							for(int j=0; j<header_field; j++)
								err_msg += Constants._vcf_header[j]+"\n";
							throw new RuntimeException(err_msg);
						}
					}
					nS = s.length-header_field;
                    bw.write(snp+"\n");
                    continue;
                }
				
				s = snp.split("\\s+");
				
				if(switchoff) {
					switchoff = false;
					info = s[Constants._vcf_header_i[7]].split(";");
					for(int i=0; i<info.length; i++) {
						if(info[i].startsWith("AF="))
							af_i = i;
						if(info[i].startsWith("DP="))
							dp_i = i;
					}
					if(af_i==-1)  myLogger.warn("No AF field in VCF file. "
                    		+ "Filtering by SNP allele frequency disabled.");
					if(dp_i==-1)  myLogger.warn("No AF field in VCF file. "
                    		+ "Filtering by SNP read depth disabled.");
					info = s[Constants._vcf_header_i[8]].split(":");
					
					for(int i=0; i<info.length; i++) {
						format_set.put(info[i], i);
					}
					if(!format_set.keySet().contains("GT"))
						throw new RuntimeException("No GT field in VCF file. Progam exit.");
					info = s[Constants._vcf_header_i[8]+1].split(":");
					ploidy_obs = info[format_set.get("GT")].split("\\||/").length;
					
					if(ploidy<0) ploidy = ploidy_obs;
					
					List<String> tmp_format = new ArrayList<String>();
					if(ploidy_obs!=ploidy) {
						if(!format_set.containsKey("AD")
								&& !(format_set.containsKey("RO") && 
										format_set.containsKey("AO")) ) {
							throw new RuntimeException("No AD field in VCF file. AD field "
									+ "is required to change the ploidy.");
						}
						format_str = Constants._vcf_format_str;
						format_arr = Constants._vcf_format;
					} else {
						for(String target_str : Constants._vcf_format) {
							switch(target_str) {
							case "GT":
								tmp_format.add("GT");
								format_str += "GT";
								break;
							case "AD":
								if(format_set.containsKey("AD")
										|| format_set.containsKey("RO") && 
										format_set.containsKey("AO")) {
									tmp_format.add("AD");
									format_str += ":AD";
									contain_ad_info = true;
								}
								break;
							case "DP":
								if(format_set.containsKey("DP") 
										||format_set.containsKey("AD")
										|| format_set.containsKey("RO") && 
										format_set.containsKey("AO")) {
									tmp_format.add("DP");
									format_str += ":DP";
								}
								break;
							case "GQ":
								if(format_set.containsKey("GQ")) {
									tmp_format.add("GQ");
									format_str += ":GQ";
								}
								break;
							case "PL":
								if(format_set.containsKey("PL")
										|| format_set.containsKey("GL")) {
									tmp_format.add("PL");
									format_str += ":PL";
								}
								break;
							default:
								throw new RuntimeException("Undefined VCF field!!!");
							}
						}
						format_arr = format_str.split(":");
					}
				}
				
				if(Double.parseDouble(s[Constants._vcf_header_i[5]])<min_qual) {
                    continue;
                }
				
				if(af_i>0) {
					double maf = Double.parseDouble(s[Constants._vcf_header_i[7]].
							split(";")[af_i].
							replaceAll("^AF=","").split(",")[0]);
					if(maf<min_maf || maf>1-min_maf)
						continue;
				}
				
				if(dp_i>0) {
					double d = Double.parseDouble(s[Constants._vcf_header_i[7]].
							split(";")[dp_i].
							replaceAll("^DP=", ""));
					if(d>max_depth || d<min_depth)
	                    continue;
				}
				
				os.setLength(0);
				os.append(s[Constants._vcf_header_i[0]].replace('|', '_'));
				os.append("\t");
				os.append(s[Constants._vcf_header_i[1]]);
				os.append("\t");
				os.append(s[Constants._vcf_header_i[2]]);
				os.append("\t");
				os.append('A');
				os.append("\t");
				os.append('B');
				info = s[Constants._vcf_header_i[4]].split(",");
				if(info.length>2) continue;
				boolean MNP = info.length>1, MN = false;
				for(int i=5; i<header_field-1; i++) {
					os.append("\t");
					os.append(s[Constants._vcf_header_i[i]]);
				}
				os.append("\t");
				os.append(format_str);
				
				double m=0;
				int[] ad = new int[2];
				int dp = -1;
                for(int i=header_field; i<s.length; i++) {
                    os.append("\t");
                    info = s[i].split(":");
                    if(s[i].startsWith(".")) {
                        os.append(misSiteStr(format_arr, ploidy));
                        m += 1.0;
                        continue;
                    }
                    if(MNP && info[0].indexOf("0")>-1) {
                        MN = true;
                        break;
                    }
                    if(format_set.containsKey("RO") &&
                    		format_set.containsKey("AO")) {
                    	if(MNP) {
                    		String[] ad_str = info[format_set.get("AO")].split(",");
                    		ad[0] = Integer.parseInt(ad_str[0]);
                    		ad[1] = Integer.parseInt(ad_str[1]);
                    	} else {
                    		ad[0] = Integer.parseInt(info[format_set.get("RO")]);
                    		ad[1] = Integer.parseInt(info[format_set.get("AO")]);
                    	}
                    	dp = ad[0]+ad[1];
                    } else if(format_set.containsKey("AD")) {
                    	String[] ad_str = info[format_set.get("AD")].split(",");
                    	int shift = 0;
                    	if(MNP) shift = 1;
                    	ad[0] = Integer.parseInt(ad_str[0+shift]);
                		ad[1] = Integer.parseInt(ad_str[1+shift]);
                    	dp = ad[0]+ad[1];
                    } else if(format_set.containsKey("DP")) {
                    	dp = Integer.parseInt(info[format_set.get("DP")]);
                    }
                    
                    if(ploidy_obs != ploidy) {
                    	double[] ll = new double[ploidy+1];
                    	double gq = fit(ll, ad, ploidy);
                    	for(String target_str : format_arr) {
							switch(target_str) {
							case "GT":
								os.append(uniGT(ll));
								break;
							case "AD":
								os.append(":");
								os.append(cat(ad,","));
		                    	break;
							case "DP":
								os.append(":");
								os.append(dp);
								break;
							case "GQ":
								os.append(":");
								os.append(gq);
								break;
							case "PL":
								os.append(":");
								os.append(cat(ll,","));
								break;
							default:
								throw new RuntimeException("Undefined VCF field!!!");
							}
						}
                    } else {
                    	
                    	String gt_str = info[format_set.get("GT")];
                    	for(String target_str : format_arr) {
							switch(target_str) {
							case "GT":
								if(!MNP)
									os.append(gt_str);
								else {
									gt_str = gt_str.replace('1', '0').
											replace('2', '1');
									os.append(gt_str);
								}
								break;
							case "AD":
								os.append(":");
								os.append(cat(ad,","));
		                    	break;
							case "DP":
								os.append(":");
								os.append(dp);
								break;
							case "GQ":
								os.append(":");
								os.append(info[format_set.get("GQ")]);
								break;
							case "PL":
								os.append(":");
								if(!MNP) {
									if(format_set.containsKey("PL"))
										os.append(info[format_set.get("PL")]);
									else {
										double[] ll = new double[ploidy+1];
										String[] tmp = info[format_set.get("GL")].split(",");
										for(int j=0; j<ploidy+1; j++)
											ll[j] = Double.parseDouble(tmp[j])*-10.0;
										os.append(cat(ll,","));
									}
								} else {
									if(ploidy_obs==2) {
										String[] qual_str = format_set.containsKey("GL") ? 
												info[format_set.get("GL")].split(","):
												info[format_set.get("PL")].split(",");
										double[] ll = new double[3];
										ll[0] = Double.parseDouble(qual_str[2]);
										ll[1] = Double.parseDouble(qual_str[4]);
										ll[2] = Double.parseDouble(qual_str[5]);
										if(format_set.containsKey("GL"))
											for(int j=0; j<3; j++)
												ll[j] *= -10.0;
										os.append(cat(ll,","));
									} else {
										if(contain_ad_info) {
											double[] ll = new double[ploidy+1];
					                    	fit(ll, ad, ploidy);
											os.append(cat(ll,","));
										} else {
											os.append(uniPL(gt_str));
										}
									}
								}
								break;
							default:
								throw new RuntimeException("Undefined VCF field!!!");
							}
						}
                    }
                }
                os.append("\n");
                if(m/nS>max_missing && !MN) continue;
                if(!MN) bw.write(os.toString());
            }
			bw.close();
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private Object misSiteStr(final String[] format_arr, final int ploidy) {
		// TODO Auto-generated method stub
		StringBuilder os = new StringBuilder();
		for(String target_str : format_arr) {
			switch(target_str) {
			case "GT":
				os.append("./.");
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
    
	private final static double err = Constants.seq_err;
	
    public static double fit(double[] ll, int[] depth, int ploidy) {
    	int d = depth[0]+depth[1];
    	double maxLL = Double.NEGATIVE_INFINITY;
        for(int i=0; i<ll.length; i++) {
    		double pa = ((ploidy-i)*(1-err)+i*err)/ploidy;
    		double nk = Math.log10(nchoosek(d, depth[0]));
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
    
    public static void fit(double[] ll, char[] genotype) {
    	Arrays.fill(ll, 255);
    	int k = ll.length-1;
    	for(char g : genotype)
    		if(g==Constants._universal_A_allele)
    			k--;
    	ll[k] = 0;
    }
    
    public static double[] fit(char[] genotype) {
    	double[] ll = new double[genotype.length+1];
    	fit(ll, genotype);
    	return ll;
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
    
    private Object uniPL(String gt_str) {
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
}
