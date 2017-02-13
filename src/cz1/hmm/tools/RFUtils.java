package cz1.hmm.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Arrays;

import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.apache.commons.math3.stat.inference.GTest;

import cz1.util.Algebra;
import cz1.util.Constants;
import cz1.util.Executor;
import cz1.util.Utils;

public abstract class RFUtils extends Executor {
	
	protected static final double kb = 1024;
	protected static final double mb = 1024*1024;
	protected static final double gb = 1024*1024*1024;
	
	protected String in_haps;
	protected String out_prefix;
	protected String[] founder_haps;
	protected String expr_id = null;
	protected int drop_thres = 1;
	protected double skew_phi = 2;
	protected int best_n = 10;
	protected final String goodness_of_fit = "fraction";
	
	protected double[] probs_uniform = new double[]{.5,.5,.5,.5};

	protected static NumberFormat formatter = new DecimalFormat("#0.000");
	
	protected static String[][] dc = null;
	protected static int nF1;
	
	protected class FileLoader implements Runnable {
		private final String id;
		private final File[] files;
		private final int i;

		public FileLoader(String id, File[] files, int i) {
			this.id = id;
			this.files = files;
			this.i = i;
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub
			double[] ll = new double[this.files.length];
			for(int k=0; k<ll.length; k++) 
				ll[k]=Double.NEGATIVE_INFINITY;

			for(int k=0; k<this.files.length; k++) {

				File file = this.files[k];
				if(!new File(file.getAbsolutePath()+
						"/phasedStates/"+expr_id+".txt").exists()) {
					System.err.println("warning: "+
							file.getName()+
							" exsits, but phased states do not.");
					continue;
				}
				try {
					BufferedReader br = 
							Utils.getBufferedReader(file.getAbsolutePath()+
									"/phasedStates/"+expr_id+".txt");
					br.readLine();
					String mak = br.readLine();
					br.close();
					if( mak==null ) {
						System.err.println("warning: "+
								file.getName()+
								" exists, but phased states are NULL.");
						continue;
					}
					if( Integer.parseInt(mak)<2 ) {
						System.err.println("warning: "+
								file.getName()+
								" exists, but #marker is less than 2.");
						continue;
					}

					String res = file.getAbsolutePath()+"/stderr_true";
					//IO.println(res);
					File _res_hmm = new File(res);
					String line;
					if( _res_hmm.exists() ) {
						BufferedReader br2 = Utils.getBufferedReader(res);
						String lprob=null;
						while( (line=br2.readLine()) !=null) {
							if(line.startsWith("log"))
								lprob = line;
						}
						br2.close();
						if( lprob!=null ) {
							String[] s0 = lprob.split("\\s+");
							ll[k] = Double.parseDouble(s0[3]);
						}
					} else {
						BufferedReader br2 = 
								Utils.getBufferedReader(file.getAbsolutePath()+
										"/phasedStates/"+
										expr_id+".txt");
						String lprob=br2.readLine();
						br2.close();
						if( lprob!=null ) 
							ll[k] = Double.parseDouble(lprob);
						System.out.println(ll[k]);
					}
				} catch (IOException e) {
					e.printStackTrace();
					System.exit(1);
				}
			}

			int[] maxN = maxN(ll);
			StringBuilder oos = new StringBuilder();
			boolean[] drop = new boolean[maxN.length];
			int dropped = 0;
			for(int k=0; k<maxN.length; k++) {

				int[] haps_observed = readHaplotypes(maxN[k]);
				
				long[] observed;
				double p;
				switch(goodness_of_fit) {
				case "fraction":
					double[] phases = new double[Constants._haplotype_z];
					for(int z=0; z<phases.length; z++) 
						phases[z] = (double) haps_observed[z];
					double expected = StatUtils.sum(phases)/Constants._ploidy_H/2;
					double maf = StatUtils.max(phases)/expected, 
							mif = StatUtils.min(phases)/expected;
					if( maf>1/skew_phi || mif<skew_phi) {
						System.out.println(this.files[maxN[k]].getName()+
								" was dropped due to large haploptype frequency variance. (" +
								Utils.cat(phases, ",") +")");
						drop[k] = true;
					}
					if(drop[k]) dropped++;
					oos.append("["+(drop[k]?"drop](maf,":"keep](maf,")+maf+";mif,"+mif+") "+
							Utils.cat(haps_observed,",")+"\t"+this.files[maxN[k]].getName()+"\n");
					break;
				case "chisq":
					observed = new long[Constants._haplotype_z];
					for(int z=0; z<observed.length; z++) 
						observed[z] = (long) haps_observed[z];
					p = new ChiSquareTest().chiSquareTest(probs_uniform, observed);
					if(p<skew_phi) drop[k] = true;
					if(drop[k]) dropped++;
					oos.append("["+(drop[k]?"drop](p,":"keep](p,")+formatter.format(p)+") "+
							Utils.cat(haps_observed,",")+"\t"+this.files[maxN[k]].getName()+"\n");
					break;
				case "gtest":
					observed = new long[Constants._haplotype_z];
					for(int z=0; z<observed.length; z++) 
						observed[z] = (long) haps_observed[z];
					p = new GTest().gTest(probs_uniform, observed);
					if(p<skew_phi) drop[k] = true;
					if(drop[k]) dropped++;
					oos.append("["+(drop[k]?"drop](p,":"keep](p,")+formatter.format(p)+") "+
							Utils.cat(haps_observed,",")+"\t"+this.files[maxN[k]].getName()+"\n");
					break;
				default:
					System.err.println("Goodness-of-fit test should be fraction, chisq or gTest.");
					System.exit(1);
				}

			}
			System.out.print(oos.toString());
			System.err.println(this.id+" - dropped "+dropped);
			if( drop.length-dropped<drop_thres ) {
				System.err.println("Scaffold "+this.id+" dropped.");
			} else {
				int kk=0;
				for(int k=0; k<drop.length; k++) {
					if(!drop[k]) {
						dc[i][kk] = this.files[maxN[k]]
								.getName();
						
						kk++;
					}
					if(kk>=best_n) break;
				}
			}
		}
		
		private int[] readHaplotypes(final int i) {
			// TODO Auto-generated method stub
			try {
				BufferedReader br_states = Utils.getBufferedReader(this.files[i]+
								"/phasedStates/"+expr_id+".txt");;
				String line, stateStr;
				String[] s;
				int[] haps_observed = new int[Constants._haplotype_z];
				while( (line=br_states.readLine())!=null ) {
					if(!line.startsWith("#")) continue;
					s = line.split("\\s+|:");
					stateStr = s[s.length-1];
					for(char h : stateStr.toCharArray())
						haps_observed[h>'9'?(h-'a'+9):(h-'1')]++;
				}
				br_states.close();
				return haps_observed;
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				System.exit(1);
			}
			return null;
		}

		private int[] maxN(double[] ll) {
			double[] ll0 = Arrays.copyOf(ll, ll.length);
			int n = ll.length;//best_n_phases[0].length;
			int[] maxN = new int[n];
			Arrays.fill(maxN, -1);
			for(int k=0; k<n; k++) {
				if(k>=ll0.length) return maxN;
				int p = 0;
				double e = Double.NEGATIVE_INFINITY;
				for(int s=0; s<ll0.length; s++)
					if(ll0[s]>e) {
						e = ll0[s];
						p = s;
					}
				maxN[k] = p;
				ll0[p] = Double.NEGATIVE_INFINITY;
			}
			return maxN;
		}
	}
	
	public void calcGDsAll(String phasedStates, int ploidy, 
			String[] parents, int nF1, int nSNP, 
			double[][] rfAll, int s) {
		// TODO Auto-generated method stub
		char[][] h = readHaplotypes(phasedStates, 
				ploidy, parents, nF1);
		int c = 0;
		for(int i=0; i<nSNP; i++) { 
			for(int j=i+1; j<nSNP; j++) {
				double r = 0;
				for(int k=0; k<h.length; k++) 
					r += h[k][i]==h[k][j] ? 0 : 1;
				rfAll[c++][s] = r/h.length;
			}
		}
	}

	public static double[] calcGDs(String phasedStates,
			String[] parents, int nF1) {
		// TODO Auto-generated method stub
		char[][] h = readHaplotypes(phasedStates, Constants._ploidy_H, parents, nF1);
		double[] d = new double[h[0].length-1];
		for(int i=0; i<d.length; i++) {
			double c = 0;
			for(int j=0; j<h.length; j++) 
				c += h[j][i]==h[j][i+1] ? 0 : 1;
			d[i] = c/h.length;
		}
		return d;
	}

	protected static double calcGD(String phasedStates,
			String[] parents, int nF1) {
		char[][] h = readHaplotypes(phasedStates, Constants._ploidy_H, parents, nF1);
		double c = 0;
		for(int i=0; i<h.length; i++)
			c += h[i][0]==h[i][h[i].length-1] ? 0 : 1;
		//return geneticDistance( c/h.length, mapFunc);
		return c/h.length;	
	}

	protected static char[][] readHaplotypes(String phasedStates, int ploidy,
			String[] parents, int nF1) {
		// TODO Auto-generated method stub
		try {
			BufferedReader br = Utils.getBufferedReader(phasedStates);
			br.readLine();
			int m = Integer.parseInt(br.readLine());
			char[][] h = new char[nF1*ploidy][m];
			String line, stateStr;
			String[] s;
			int c = 0;
			while( (line=br.readLine())!=null ) {
				if(!line.startsWith("#")) continue;
				//if(skip++<2) continue;
				s = line.split("\\s+|:");
				if(Arrays.asList(parents).contains(s[2])) continue;
				stateStr = s[s.length-1];
				h[c++] = stateStr.toCharArray();
			}
			br.close();
			return h;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.exit(1);
		}

		return null;
	}
	
	protected static BufferedWriter mapWriter;
	
	protected class mapCalculator implements Runnable {
		private final int i;
		public mapCalculator(int i) {
			this.i = i;
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub
			double[][] rf_all = new double[dc[this.i].length][];

			for(int k=0; k<dc[this.i].length; k++) 
				if( dc[i][k]!=null ) 
					rf_all[k] = calcGDs(in_haps+"/"+dc[i][k]+
							"/phasedStates/"+expr_id+".txt",
							founder_haps,
							nF1);

			String contig = dc[i][0].replace(expr_id,"experiment")
					.split("\\.")[1];
			rf_all = Algebra.transpose(rf_all);
			
			try {
				StringBuilder os = new StringBuilder();
				os.append("*"+contig+"\n");
				for(int k=0; k<rf_all.length; k++)
					if(rf_all[k]!=null)
						os.append(Utils.cat(rf_all[k], ",")+"\n");
				mapWriter.write(os.toString());
			} catch (MathIllegalArgumentException | IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	protected static double geneticDistance(double r, String mapFunc) {
		// TODO Auto-generated method stub
		switch(mapFunc.toUpperCase()) {
		case "KOSAMBI":
			return .25*Math.log((1+2*r)/(1-2*r));
		case "HALDANE":
			return -.5*Math.log(1-2*r);	
		default:
			System.err.println("Error - Undefined genetic mapping function.");
			System.exit(1);
		}
		return -1;
	}
}
