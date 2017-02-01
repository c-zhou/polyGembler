package cz1.tools;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.HashMap;
import java.util.Set;
import java.util.Map;


public class VcfToolPlugin {
	
	public static void main(String[] args) {
		/**
		System.out.println(nchoosek(60,30));
		double[] ll = fit(new int[]{20,40},6);
		for(int i=0; i<ll.length; i++)
			System.out.println(ll[i]);
		**/
		if(args.length==2)
            recode(args[0], args[1]);
        else if(args.length==7)
            if(Integer.parseInt(args[5])==4)
	    	    recode4(args[0], args[1], args[2], Double.parseDouble(args[3]), Integer.parseInt(args[4]), 
	    	    		Integer.parseInt(args[5]), Double.parseDouble(args[6]));
	        else if(Integer.parseInt(args[5])==2)
                recode2(args[0], args[1], args[2], Double.parseDouble(args[3]), Integer.parseInt(args[4]), 
                		Integer.parseInt(args[5]), Double.parseDouble(args[6]));
            else if(Integer.parseInt(args[5])==6)
                recode6(args[0], args[1], args[2], Double.parseDouble(args[3]), Integer.parseInt(args[4]),
                        Integer.parseInt(args[5]), Double.parseDouble(args[6]));
        else
            throw new RuntimeException("!!!");
    }

    private static void recode6(String in, String out, String log, double thres, int nF1, int ploidy, double avg_dp) {
        String snp;
        String[] s, info;
        double max_dp = avg_dp*nF1*1.5;
        StringBuilder os = new StringBuilder();
        
        try {
            BufferedReader br = new BufferedReader(new FileReader(in));
            BufferedWriter bw = new BufferedWriter(new FileWriter(out));
            BufferedWriter log_f = new BufferedWriter(new FileWriter(log));
            int[] dp = new int[2];
            Set<String> f1 = new HashSet<String>();
            for(int i=1; i<=nF1-2; i++)
                f1.add("Mx23Hm6F"+i);
            f1.add("Mx23Hm6P1");
            f1.add("Mx23Hm6P2");
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
                    for(int i : index)
                        System.out.print(i+"\t");
                    System.out.println();
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
                
                double d = Double.parseDouble(s[7].
                        split(";")[7].
                        split("=")[1]);
                if(d>max_dp) {
                    log_f.write("DUP | "+snp+"\n");
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
                    log_f.write("MNP | "+snp+"\n");
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
                    if(s[i].startsWith(".")) {
                        os.append("./.:0,0:0:0:0,0,0,0,0,0,0");
                        m+=1.0;
                        continue;
                    }
                    if(MNP && info[0].indexOf("0")>-1) {
                        log_f.write("MNP | "+snp+"\n");
                        MN = true;
                        break;
                    }
                    if(MNP) {
                        dp[0] = Integer.parseInt(info[6].split(",")[0]);
                        dp[1] = Integer.parseInt(info[6].split(",")[1]);
                        double[] ll = new double[ploidy+1];
                        double gq = fit(ll, dp, ploidy);
                        os.append(uniGT(ll));
                        os.append(":");
                        os.append(info[6]);
                        os.append(":");
                        os.append(dp[0]+dp[1]);
                        os.append(":");
                        os.append(gq);
                        os.append(":");
                        os.append(cat(ll,","));
                    } else {
                        dp[0] = Integer.parseInt(info[4]);
                        dp[1] = Integer.parseInt(info[6]);
                        double[] ll = new double[ploidy+1];
                        double gq = fit(ll, dp, ploidy);
                        os.append(uniGT(ll));
                        os.append(":");
                        os.append(info[4]);
                        os.append(",");
                        os.append(info[6]);
                        os.append(":");
                        os.append(dp[0]+dp[1]);
                        os.append(":");
                        os.append(gq);
                        os.append(":");
                        os.append(cat(ll,","));
                    }
                }
                os.append("\n");
                if(m/nF1>thres && !MN) {
                    log_f.write("MIS | "+snp+"\n");
                    continue;
                }
                if(!MN) bw.write(os.toString());
            }
            br.close();
            bw.close();
            log_f.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    private static void recode2(String in, String out, String log, double thres, int nF1, int ploidy, double avg_dp) {
        String snp;
        String[] s, s0, info;
        StringBuilder os = new StringBuilder();
        double max_dp = avg_dp*nF1*1.5;
        try {
            BufferedReader br = new BufferedReader(new FileReader(in));
            BufferedWriter bw = new BufferedWriter(new FileWriter(out));
            BufferedWriter log_f = new BufferedWriter(new FileWriter(log));
            int[] dp = new int[2];
            
            while( (snp=br.readLine())!=null ) {
                if(snp.startsWith("#")) {
                    bw.write(snp+"\n");
                    continue;
                }                                                                                    
                s = snp.split("\\s+");
                
                if(Double.parseDouble(s[5])<100)
                    continue;
                
                double af = Double.parseDouble(s[7].
                        split(";")[3].
                        split("=|,")[1]);
                if(af>.9 || af<.1) {
                    log_f.write("MAF | "+snp+"\n");
                    continue;
                }
                
                double d = Double.parseDouble(s[7].
                        split(";")[7].
                        split("=")[1]);
                if(d>max_dp) {
                    log_f.write("DUP | "+snp+"\n");
                    continue;
                }

                os.setLength(0);
                os.append(s[0].replaceAll("^Itr_sc0{0,9}","").replaceAll("\\.1$",""));
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
                    log_f.write("MNP | "+snp+"\n");
                    continue;
                }
                if(info.length>2) continue;
                boolean MNP = info.length>1,
                        MN = false;
                for(int i=5; i<8; i++) {
                    os.append("\t");
                    os.append(s[i]);
                }
                
                os.append("\tGT:AD:DP:GQ:PL");
                double m=0;
                
                for(int i=9; i<s.length; i++) {
                    os.append("\t");
                    info = s[i].split(":");
                    if(s[i].startsWith(".")) {
                        os.append("./.:0,0:0:0:0,0,0");
                        m+=1.0;
                        continue;
                    }
                    if(MNP && info[0].indexOf("0")>-1) {
                        log_f.write("MNP | "+snp+"\n");
                        MN = true;
                        break;
                    }
                    
                    if(MNP) {
                    	dp[0] = Integer.parseInt(info[6].split(",")[0]);
                        dp[1] = Integer.parseInt(info[6].split(",")[1]);
                        double[] ll = new double[ploidy+1];
                        double gq = fit(ll, dp, ploidy);
                        os.append(uniGT(ll));
                        os.append(":");
                        os.append(info[6]);
                        os.append(":");    
                        os.append(dp[0]+dp[1]);
                        os.append(":");
                        os.append(gq);
                        os.append(":");
                        os.append(cat(ll,","));
                        /**
                        s0 = info[8].split(",");
                        os.append(-10.0*Double.parseDouble(s0[2]));
                        os.append(",");
                        os.append(-10.0*Double.parseDouble(s0[4]));
                        os.append(",");
                        os.append(-10.0*Double.parseDouble(s0[5]));
                    	**/
                    } else {
                    	dp[0] = Integer.parseInt(info[4]);
                        dp[1] = Integer.parseInt(info[6]);
                        double[] ll = new double[ploidy+1];
                        double gq = fit(ll, dp, ploidy);
                        os.append(uniGT(ll));
                        os.append(":");
                        os.append(info[4]);
                        os.append(",");
                        os.append(info[6]);
                        os.append(":");
                        os.append(dp[0]+dp[1]);
                        os.append(":");
                        os.append(gq);
                        os.append(":");
                        os.append(cat(ll,","));

                        /**
                        s0 = info[8].split(",");
                        os.append(-10.0*Double.parseDouble(s0[0]));
                        os.append(",");
                        os.append(-10.0*Double.parseDouble(s0[1]));
                        os.append(",");
                        os.append(-10.0*Double.parseDouble(s0[2]));
                    	**/
                    }
                }
        
                os.append("\n");
                if(m/nF1>thres && !MN) {
                    log_f.write("MIS | "+snp+"\n");
                    continue;
                }
                if(!MN) bw.write(os.toString());
            }
                  
            br.close();
            bw.close();
            log_f.close();
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
        }

    }
    
    private final static double err = 0.01;
    private static double fit(double[] ll, int[] depth, int ploidy) {
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
    
	private static void recode4(String in, String out, String log, double thres, int nF1, int ploidy, double avg_dp) {
		String snp;
		String[] s, info;
		double max_dp = avg_dp*nF1*1.5;
        StringBuilder os = new StringBuilder();
        try {
            BufferedReader br = new BufferedReader(new FileReader(in));
            BufferedWriter bw = new BufferedWriter(new FileWriter(out));
            BufferedWriter log_f = new BufferedWriter(new FileWriter(log));
			
            int[] dp = new int[2];
            Set<String> f1 = new HashSet<String>();
			for(int i=1; i<=nF1-2; i++)
				f1.add("Mx23Hm4F"+i);
			f1.add("Mx23Hm4P1");
			f1.add("Mx23Hm4P2");
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
                    for(int i : index)
                        System.out.print(i+"\t");
                    System.out.println();
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

                double d = Double.parseDouble(s[7].
                        split(";")[7].
                        split("=")[1]);
                if(d>max_dp) {
                    log_f.write("DUP | "+snp+"\n");
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
                    log_f.write("MNP | "+snp+"\n");
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
					if(s[i].startsWith(".")) {
						os.append("./.:0,0:0:0:0,0,0,0,0");
						m+=1.0;
						continue;
					}
					if(MNP && info[0].indexOf("0")>-1) {
                        log_f.write("MNP | "+snp+"\n");
						MN = true;
                        break;
					}
					if(MNP) {
						dp[0] = Integer.parseInt(info[6].split(",")[0]);
                        dp[1] = Integer.parseInt(info[6].split(",")[1]);
                        double[] ll = new double[ploidy+1];
                        double gq = fit(ll, dp, ploidy);
                        os.append(uniGT(ll));
                        os.append(":");
                        os.append(info[6]);
                        os.append(":");    
                        os.append(dp[0]+dp[1]);
                        os.append(":");
                        os.append(gq);
                        os.append(":");
                        os.append(cat(ll,","));
                    } else {
						dp[0] = Integer.parseInt(info[4]);
                        dp[1] = Integer.parseInt(info[6]);
                        double[] ll = new double[ploidy+1];
                        double gq = fit(ll, dp, ploidy);
                        os.append(uniGT(ll));
                        os.append(":");
                        os.append(info[4]);
                        os.append(",");
                        os.append(info[6]);
                        os.append(":");
                        os.append(dp[0]+dp[1]);
                        os.append(":");
                        os.append(gq);
                        os.append(":");
                        os.append(cat(ll,","));
                    }
                }
				os.append("\n");
				if(m/nF1>thres && !MN) {
					log_f.write("MIS | "+snp+"\n");
                    continue;
				}
                if(!MN) bw.write(os.toString());
			}
			br.close();
			bw.close();
            log_f.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
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

	private static void recode(String in, String out) {
		String snp;
		String[] s, info, qs;
	    double[] q;
		StringBuilder os = new StringBuilder();
		try {
		        
            BufferedReader br = new BufferedReader(new FileReader(in));
            BufferedWriter bw = new BufferedWriter(new FileWriter(out));

            while( (snp=br.readLine())!=null ) {
				if(snp.startsWith("#")) {
					bw.write(snp+"\n");
					continue;
				}
				s = snp.split("\\s+");
				os.setLength(0);
				os.append(s[0].replaceAll("^Itr_sc0{0,9}", "").
						replace(".1", ""));
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
				os.append("\tGT:AD:DP:GQ:PL");
				for(int i=9; i<s.length; i++) {
					os.append("\t");
					if(s[i].startsWith(".")) {
						os.append("./.:0,0:0:33:0,0,0");
						continue;
					}
					info = s[i].split(":");
					/**
                    os.append(info[0]);
					os.append(":");
					os.append(info[2]);
					os.append(",");
					os.append(info[4]);
					os.append(":");
					os.append(info[1]);
					os.append(":");
					qs = info[6].split(",");
					q = new int[qs.length];
					for(int k=0; k<q.length; k++)
						q[k] = (int) Math.round(
								-10*Double.parseDouble(qs[k]));
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
					os.append(q[0]);
					for(int k=1; k<q.length; k++) {
						os.append(",");
						os.append(q[k]);
					}
                    **/

                    os.append(info[0]);
                    os.append(":");
                    os.append(info[4]);
                    os.append(",");
                    os.append(info[6]);
                    os.append(":");
                    os.append(info[2]);
                    os.append(":");
                    os.append(info[1]);
                    qs = info[8].split(",");
                    q = new double[qs.length];
                    for(int k=0; k<q.length; k++)
                        q[k] = -10.0*Double.parseDouble(qs[k]);
                    os.append(":");
                    os.append(q[0]);
                    for(int k=1; k<q.length; k++) {
                        os.append(",");
                        os.append(q[k]);
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
