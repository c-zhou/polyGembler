package cz1.hmm.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.log4j.Logger;

import cz1.hmm.data.DataEntry;
import cz1.math.Combination;
import cz1.math.Permutation;
import cz1.math.SaddlePointExpansion;
import cz1.util.Constants;
import cz1.util.Constants.Field;

public abstract class HiddenMarkovModel {
	protected final static Logger myLogger = Logger.getLogger(HiddenMarkovModel.class);

	protected final Field field;
	protected final DataEntry de;

	protected int M; // #markers
	protected int N; // #individuals
	protected int H; // #founder haplotypes
	protected int K; // #compound hidden states
	protected String[] samples;
	protected String[] parents;
	protected int[] parents_i;
	protected int[] progeny_i;
	protected double[] distance;
	
	protected List<Integer[]> sspace; // state space for each sample
	
	protected StateUnit state;
	protected ObUnit[][] obs;
	protected EmissionUnit[] emission;
	protected TransitionUnit[] transition;
	protected ViterbiUnit[] vbs;
	
	public HiddenMarkovModel(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse,
			Field field) {
		this.field = field;
		this.de = this.catDE(de, seperation, reverse);
		this.initialise();
	}
	
	public abstract void train();
	public abstract double loglik();
	public abstract void write(String output, 
			String experiment, 
			String contig);
	public abstract void write(String output, 
			String experiment, 
			String contig,
			double loglik_diff);
	public abstract void print(boolean details);
	
	public void print() {this.print(false);}

	protected DataEntry catDE(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse) {
		// TODO Auto-generated method stub
		for(int i=0; i<de.length; i++)
			if(reverse[i]) de[i].reverse();
		for(int i=1; i<de.length; i++) {
			de[0].addAll(de[i], seperation[i-1]);
		}
		return de[0];
	}

	private void initialise() {
		// TODO Auto-generated method stub
		this.samples = de.getSample();
		this.M = this.de.getAllele().size();
		this.N = this.samples.length;
		this.H = Constants._ploidy_H;
		this.state = new StateUnit(H);
		this.K = state.hsc.length;
		
		this.parents = Constants._founder_haps.split(":");
		List<Integer> pars = new ArrayList<Integer>();
		List<Integer> pros = new ArrayList<Integer>();
		for(int i=0; i<samples.length; i++) {
			boolean isprogeny = true;
			for(int j=0; j<parents.length; j++) {
				if(samples[i].equals(parents[j])) {
					pars.add(i);
					isprogeny = false;
					break;
				}
			}
			if(isprogeny) pros.add(i);
		}
		this.parents_i = ArrayUtils.toPrimitive(pars.toArray(new Integer[pars.size()]));
		this.progeny_i = ArrayUtils.toPrimitive(pros.toArray(new Integer[pros.size()]));
		double[] position = de.getPosition();
		this.distance = new double[this.M-1];
		for(int i=0; i<distance.length; i++)
			distance[i] = position[i+1]-position[i];
		this.sspace = new ArrayList<>(N);
		for(int i=0; i<N; i++) sspace.add(null);
		sspace.set(parents_i[0], new Integer[]{0});
		sspace.set(parents_i[1], new Integer[]{1});
		Integer[] progeny_s = new Integer[K-2];
		for(int i=0; i<K-2; i++) progeny_s[i] = i+2;
		for(int i : progeny_i) sspace.set(i, progeny_s);
		
		this.makeObUnits();
		this.makeEmissionUnits();
		this.makeTransitionUnits();
		this.makeViterbiUnits();
	}

	private void makeViterbiUnits() {
		// TODO Auto-generated method stub
		vbs = new ViterbiUnit[N];
		for(int i=0; i<N; i++) vbs[i] = new ViterbiUnit();
	}

	private final int allele_d = 10; 
	private void makeObUnits() {
		// TODO Auto-generated method stub
		this.obs = new ObUnit[N][M];
		switch(this.field) {
		case AD:
			List<List<int[]>> ad = this.de.getAlleleDepth();
			if(ad==null) throw new RuntimeException("AD feild not available!!! Try GT (-G/--genotype) options.");
			for(int i=0; i<N; i++) {
				for(int j=0; j<M; j++) {
					int[] dp = ad.get(j).get(i);
					obs[i][j] = new ObUnit(dp[0]+dp[1], dp[0], K);
				}
			}
			break;
		case GT:
			List<List<String[]>> gt = this.de.getGenotype();
			if(gt==null) throw new RuntimeException("GT field not available!!! Try AD (-D/--allele-depth) option.");
			List<String[]> allele = this.de.getAllele();
			int acnt, bcnt;
			for(int i=0; i<this.M; i++) {
				String[] a = allele.get(i);
				for(int j=0; j<this.N; j++) {
					String[] g = gt.get(i).get(j);
					acnt = 0;
					bcnt = 0;
					for(int k=0; k<H; k++) {
						acnt += (g[k].equals(a[0]) ? 1 : 0);
						bcnt += (g[k].equals(a[1]) ? 1 : 0);
					}
					acnt *= allele_d;
					bcnt *= allele_d;
					obs[j][i] = new ObUnit(acnt+bcnt, acnt, K);
				}
			}
			break;
		default:
			throw new RuntimeException("!!!");
		}
	}

	private void makeEmissionUnits() {
		// TODO Auto-generated method stub
		this.emission = new EmissionUnit[M];
		for(int i=0; i<M; i++) 
			emission[i] = new EmissionUnit(de.getAllele().get(i),
					H*2,
					K,
					bfrac(i));
	}
	
	private double bfrac(final int i) {
		// TODO Auto-generated method stub
		double baf = 0.5, acnt, bcnt;
		switch(this.field) {
		case AD:
			if(this.de.getAlleleDepth()==null)
				throw new RuntimeException("AD feild not available!!! Try GT (-G/--genotype) options.");
			List<int[]> ad = this.de.getAlleleDepth().get(i);
			acnt = 0;
			bcnt = 0;
			for(int j=0; j<N; j++) {
				int[] aa = ad.get(j);
				acnt += aa[0];
				bcnt += aa[1];
			}
			if(acnt!=0&&bcnt!=0) baf = bcnt/(acnt+bcnt);
			break;
		case GT:
			if(this.de.getGenotype()==null)
				throw new RuntimeException("GT field not available!!! Try AD (-D/--allele-depth) option.");
			List<String[]> gt = this.de.getGenotype().get(i);
			String[] allele = this.de.getAllele().get(i);
			acnt = 0;
			bcnt = 0;
			for(int j=0; j<N; j++) {
				String[] g = gt.get(j);
				for(int k=0; k<H; k++) {
					acnt += (g[k].equals(allele[0]) ? 1 : 0);
					bcnt += (g[k].equals(allele[1]) ? 1 : 0);
				}
			}
			if(acnt!=0&&bcnt!=0) baf = bcnt/(acnt+bcnt);
			break;
		default:
			throw new RuntimeException("!!!");
		}
		return baf;
	}

	private void makeTransitionUnits() {
		// TODO Auto-generated method stub
		this.transition = new TransitionUnit[M-1];
		for(int i=0; i<M-1; i++)
			transition[i] = new TransitionUnit(distance[i], state.confs);
	}
	
	protected void refresh() {
		// refresh emission probabilities for obs
		double[] emissc;
		for(int i=0; i<M; i++) {
			emissc = emission[i].emissc;
			for(int j=0; j<N; j++) {
				obs[j][i].updateEmiss(emissc);
			}
		}
	}
	
	protected double findViterbiPath() {
		// TODO Auto-generated method stub
		// not sure if this makes much difference
		double probability = 0.0;
		for(int i=0; i<N; i++) {
			Integer[] ss = this.sspace.get(i);
			ObUnit[] ob = obs[i];
			double pi = 1.0/ss.length; // init probability
			
			this.vbs[i].clear();
			double[][] v = this.vbs[i].v;
			int[][] trace = this.vbs[i].trace;
			double[] logscale = this.vbs[i].logscale;
			
			double[] emiss = ob[0].emiss;;
			TransitionUnit t;
			
			for(int k : ss) v[0][k] = pi*emiss[k];
			logscale[0] = ob[0].getLogScale(i);
			for(int j=1; j<M; j++) {
				emiss = ob[j].emiss;
				t = transition[j-1];

				for(int k : ss) {
					double a, b = 0, c = 0;
					int s = k;
					for(int z : ss) {
						a = v[j-1][z]*t.trans(z, k);
						if(a > c) {
							c = a;
							s = z;
						}
						if(z==k) b = a;
					}

					if(b==c) s = k;
					a = emiss[k]*c;
					trace[j-1][k] = s;
					v[j][k] = a;
				}
				logscale[j] = ob[j].getLogScale(i);
				vbs[i].scale(j);
			}

			this.vbs[i].finalise();
			probability += this.vbs[i].probability();
		}

		return probability;
	}
	
	protected double findViterbiPath1() {
		// TODO Auto-generated method stub
		// not sure yet if this makes much difference
		// cause we need to calculate a lot of logs
		double probability = 0.0;
		for(int i=0; i<N; i++) {
			Integer[] ss = this.sspace.get(i);
			ObUnit[] ob = obs[i];
			double pi = 1.0/ss.length; // init probability
			
			// pre-calculation for a lower bound for Viterbi path
			// used to get rid of of some paths later
			double lower_bound = Double.NEGATIVE_INFINITY;
			for(int k : ss) {
				double x = pi*ob[0].emiss[k];
				double logscale = ob[0].getLogScale(i);
				
				for(int j=1; j<M; j++) {
					x *= ob[j].emiss[k]*transition[j-1].trans(k,k);
					logscale += ob[j].getLogScale(i);
					if(x<Constants.threshMin) {
						logscale += Constants.logThreshMax;
						x /= Constants.threshMax;
					}
				}
				x = Math.log(x)+logscale;
				if(lower_bound < x) lower_bound = x;
			}
			lower_bound -= 0.04139268515;

			this.vbs[i].clear();
			double[][] v = this.vbs[i].v;
			int[][] trace = this.vbs[i].trace;
			double[] logscale = this.vbs[i].logscale;
			
			Set<Integer> ss_copy   = new HashSet<Integer>(Arrays.asList(ss));
			Set<Integer> ss_from = new HashSet<Integer>();
			double[] emiss;
			TransitionUnit t;
			
			for(int k : ss) v[0][k] = pi*ob[0].emiss[k];
			logscale[0] = ob[0].getLogScale(i);
			for(int j=1; j<M; j++) {
				ss_from.clear();
				ss_from.addAll(ss_copy);
				ss_copy.clear();
				emiss = ob[j].emiss;
				t = transition[j-1];

				for(int k : ss) {
					double a, b = 0, c = 0;
					int s = k;
					for(int z : ss_from) {
						a = v[j-1][z]*t.trans(z, k);
						if(a > c) {
							c = a;
							s = z;
						}
						if(z==k) b = a;
					}

					if(b==c) s = k;
					a = emiss[k]*c;
					
					if(Math.log(a)+logscale[j-1]+ob[j].getLogScale(i)>=lower_bound) {
						trace[j-1][k] = s;
						v[j][k] = a;
						ss_copy.add(k);
					}
				}
				logscale[j] = ob[j].getLogScale(i);
				vbs[i].scale(j);
			}

			this.vbs[i].finalise();
			probability += this.vbs[i].probability();
		}

		return probability;
	}
	
	protected static void clear(double[][][] matrix) {
		// TODO Auto-generated method stub
		int r = matrix.length, 
				c = matrix[0].length;
		for(int i=0; i<r; i++)
			for(int j=0; j<c; j++)
				Arrays.fill(matrix[i][j], 0);

	}

	protected static void clear(double[][] matrix) {
		// TODO Auto-generated method stub
		int r = matrix.length;
		for(int i=0; i<r; i++)
			Arrays.fill(matrix[i], 0);
	}

	protected static void clear(int[][] matrix) {
		// TODO Auto-generated method stub
		int r = matrix.length;
		for(int i=0; i<r; i++)
			Arrays.fill(matrix[i], 0);
	}
	
	public static void clear(double[] array) {
		// TODO Auto-generated method stub
		Arrays.fill(array, 0);
	}
	
	protected class ObUnit {
		private final int cov;	// depth of coverage
		private final int aa;	// A-allele depth
		protected final double[] emiss; // place holder for emission probs
		protected final double[] logscale;
		
		public ObUnit(final int cov,
				final int aa, 
				final int k) {
			this.cov = cov;
			this.aa = aa;
			this.emiss = new double[k];
			this.logscale = new double[3];
		}

		public double getLogScale(final int i) {
			// TODO Auto-generated method stub
			if(i==parents_i[0])
				return logscale[0];
			else if(i==parents_i[1])
				return logscale[1];
			else
				return logscale[2];
		}

		public void updateEmiss(double[] emissA) {
			if(emiss.length!=emissA.length) 
				throw new RuntimeException("!!!");
			if(cov==0) {
				Arrays.fill(emiss, 1.0/emiss.length);
				return;
			}
			for(int i=0; i<emiss.length; i++)
				emiss[i] = SaddlePointExpansion.logBinomialProbability(aa, cov, emissA[i]);
			
			logscale[0] = emiss[0];
			emiss[0] = 1.0;
			logscale[1] = emiss[1];
			emiss[1] = 1.0;
			logscale[2] = StatUtils.max(emiss, 2, K-2);
			for(int i=2; i<K; i++) {
				emiss[i] = Math.exp(emiss[i]-logscale[2]);
			}
		}

		public int getCov() {
			return this.cov;
		}

		public int getAa() {
			return this.aa;
		}

		public double[] getEmiss() {
			return this.emiss;
		}
	}

	public class StateUnit {
		private final char[] hs; // hidden states
		private final int[][] hsc; // compound hidden states
		private final String[] hsc_str; // compound hidden states str
		private final int confs; // #confs
		private final int[][] confs_hsc; // confs of compound hidden states
		
		private final int[][] confs_tab; // summary of confs
		private final int[][] confs_cnt; // confs for probability calculation
		/***
		 * p is the probability of jumps
		 * 
		 * 			(1-p)^0*p^k		(1-p)^1*p^(k-1)		(1-p)^2*p^(k-2)		...		(1-p)^k*p^0
		 * 	1			
		 * 	2
		 * 	3
		 * 	...
		 */
		
		public StateUnit(final int h) {
			hs = makeHs(h*2);
			hsc = makeHsc();
			hsc_str = makeHscStr();
			confs = (h/2+1)*(h/2+2)/2+1;
			confs_tab = new int[confs][2];
			confs_hsc = new int[hsc.length][hsc.length];
			confs_cnt = new int[confs][h+1];
			
			this.configure();
		}

		private void configure() {
			// TODO Auto-generated method stub
			final int h = hs.length/4;
			final int k = hsc.length;
			// make confs look up table
			final Map<String, Integer> confs_str = new HashMap<>();
			int c = 0;
			for(int i=0; i<=h; i++) {
				for(int j=i; j<=h; j++) {
					++c;
					confs_str.put(i+"_"+j, c);
					confs_str.put(j+"_"+i, c);
					confs_tab[c][0] = i;
					confs_tab[c][1] = j;
				}
			}
			// make compound states look up table
			int comn0, comn1;
			int[] h0i = new int[h], 
					h1i = new int[h],
					h0j = new int[h], 
					h1j = new int[h];
			for(int i=0; i<k; i++) {
				if(i<2) {
					Arrays.fill(confs_hsc[i], -1);
					confs_hsc[i][i] = 0;
					continue;
				}
				System.arraycopy(hsc[i], 0, h0i, 0, h);
				System.arraycopy(hsc[i], h, h1i, 0, h);
				for(int j=0; j<k; j++) {
					if(j<2) {
						confs_hsc[i][j] = -1;
						continue;
					}
					System.arraycopy(hsc[j], 0, h0j, 0, h);
					System.arraycopy(hsc[j], h, h1j, 0, h);
					comn0 = hsComn(h0i, h0j);
					comn1 = hsComn(h1i, h1j);
					confs_hsc[i][j] = confs_str.get(comn0+"_"+comn1);
				}
			}
			// make probability calculation look up table
			int[][] comns = new int[h+1][h+1];
			for(int i=0; i<=h; i++) {
				for(int j=0; j<=i; j++) {
					comns[i][j] = f(h, i, j);
				}
			}
			
			int[] probs;
			int[] comns0, comns1;
			for(int i=1; i<confs; i++) {
				comn0 = confs_tab[i][0];
				comn1 = confs_tab[i][1];
				comns0 = comns[comn0];
				comns1 = comns[comn1];
				probs = confs_cnt[i];
				for(int u=0; u<=comn0; u++)
					for(int v=0; v<=comn1; v++) 
						probs[u+v] += comns0[u]*comns1[v];
			}
		}
		
		private int f(int n, int a, int b) {
			// TODO Auto-generated method stub
			// inefficient but straightforward
			/***
			 * f(n,a,b) = c(a,b)*f(n-b,a-b,0)
			 * f(n,a,0) = n!-f(n,a,1)-f(n,a,2)-...-f(n,a,a) 
			 */
			if(n<a||n<b||a<b) return 0;
			if(b==0) {
				long f = Permutation.factorial(n);
				if(f>Integer.MAX_VALUE) throw new RuntimeException("!!!");
				for(int i=1; i<=a; i++) f -= f(n,a,i);
				return (int) f;
			}
			return Combination.nchoosek(a, b)*f(n-b,a-b,0);
		}

		public void calcProbs(double[] trans, double[][] cnts, double jump) {
			final int h = hs.length/2;
			trans[0] = 1.0;
			cnts[0][0] = 0.0;
			cnts[0][1] = h;
			int[] cnt;
			final double stay = 1-jump;
			for(int i=1; i<confs; i++) {
				cnt = confs_cnt[i];
				double[] probs = new double[h+1];
				for(int j=0; j<=h; j++) 
					probs[j] = cnt[j]*Math.pow(jump, h-j)*Math.pow(stay, j);
				trans[i] = StatUtils.sum(probs);
				for(int j=0; j<=h; j++) {
					cnts[i][0] += (h-j)*probs[j]/trans[i];
					cnts[i][1] += j*probs[j]/trans[i];
				}
			}
		}
		
		private char[] makeHs(int h) {
			// TODO Auto-generated method stub
			char[] hs = new char[h];
			for(int k=0; k<h; k++) {
				hs[k] = k<10 ? (char)('1'+k) : (char)('a'+k-10);
			}
			return hs;
		}

		private int[][] makeHsc() {
			// TODO Auto-generated method stub
			int mid1 = this.hs.length/2;
			int[] p1 = new int[mid1], p2 = new int[mid1];
			for(int i=0; i<mid1; i++) {
				p1[i] = i;
				p2[i] = mid1+i;
			}

			List<List<Integer>> com1 = Combination.combination(p1, mid1/2),
					com2 = Combination.combination(p2, mid1/2);
			int k = com1.size();
			int k2 = k*k+2;
			int[][] hsc = new int[k2][mid1];
			hsc[0] = p1;
			hsc[1] = p2;
			for(int i=0; i<k; i++) {
				for(int j=0; j<k; j++) {
					List<Integer> com = new ArrayList<Integer>(com1.get(i));
					com.addAll(com2.get(j));
					hsc[i*k+j+2] = ArrayUtils.toPrimitive(com.toArray(new Integer[mid1]));
				}
			}
			return hsc;
		}

		private String[] makeHscStr() {
			// TODO Auto-generated method stub
			int K = this.hsc.length;
			int h = this.hs.length/2;
			String[] hsc_str = new String[K];
			char[] st = new char[h];
			for(int i=0; i<K; i++) {
				for(int j=0; j<h; j++) {
					st[j] = hs[hsc[i][j]];
				}
				hsc_str[i] = StringUtils.join(st, '_');
			}
			return hsc_str;
		}

		private int hsComn(int[] h1, int[] h2) {
			// TODO Auto-generated method stub
			int h = h1.length;
			int u = 0, v = 0, c = 0;
			while (u<h && v<h) {
				if(h1[u]==h2[v]) {
					++c;
					++u;
					++v;
				} else if (h1[u] < h2[v]) {
					++u;
				} else {
					++v;
				}
			}
			return c;
		}
		
		protected char[] getHs() {
			return this.hs;
		}
		
		protected int[][] getHsc() {
			return this.hsc;
		}
		
		protected String[] getHscStr() {
			return this.hsc_str;
		}

		protected int getConfs() {
			return this.confs;
		}
		
		protected int[][] getConfsHsc() {
			return this.confs_hsc;
		}
		
		protected int hsc(int i, int j) {
			return confs_hsc[i][j];
		}
	}

	protected class EmissionUnit {
		private final String[] allele;
		private final double[] emiss;
		private final double[] emissc;
		private final double pseudo;
		private final double[][] count;
		private final double[][][] cnts_prior; // priors for counting
		
		public EmissionUnit(final String[] allele,
				int H,
				int K,
				final double pseudo) {
			this.allele = allele;
			this.emiss = new double[H];
			this.emissc = new double[K];
			this.pseudo = pseudo;
			this.count = new double[H][allele.length];
			this.cnts_prior = new double[K][H/2][2];
			this.prior();
			this.updatec();
		}

		protected void prior() {
			// TODO Auto-generated method stub
			BetaDistribution beta = new BetaDistribution(Constants.rg, 
					(1-pseudo)*Constants._mu_A_e, pseudo*Constants._mu_A_e);
			for(int i=0; i<emiss.length; i++) {
				emiss[i] = beta.sample();
				if(emiss[i]==0) emiss[i] = 0.001;
				if(emiss[i]==1) emiss[i] = 0.999;
			}
		}

		public void addCount(int s, double acnt, double bcnt) {
			// TODO Auto-generated method stub
			int[] hs = state.hsc[s];
			double[][] prior = cnts_prior[s];
			for(int i=0; i<H; i++) {
				count[hs[i]][0] += acnt*prior[i][0];
				count[hs[i]][1] += bcnt*prior[i][1];
			}
		}
		
		protected void update() {
			// TODO Auto-generated method stub
			for(int i=0; i<emiss.length; i++) {
				emiss[i] = count[i][0]/
				(count[i][0]+count[i][1]);
			}
			this.updatec();
		}

		protected void updatec() {
			// TODO Auto-generated method stub
			int[][] hsc = state.hsc;
			for(int i=0; i<emissc.length; i++) {
				double p = 0, p1;
				for(int j=0; j<H; j++)
					p += emiss[hsc[i][j]];
				p1 = H-p;
				for(int j=0; j<H; j++) {
					cnts_prior[i][j][0] = emiss[hsc[i][j]]/p;
					cnts_prior[i][j][1] = (1-emiss[hsc[i][j]])/p1;
				}
				emissc[i] = p/H;
			}
		}
		
		protected String[] getAllele() {
			return this.allele;
		}
		
		protected double[] getEmiss() {
			return this.emiss;
		}
		
		protected double[] getEmissc() {
			return this.emissc;
		}
		
		protected double getPseudoCount() {
			return this.pseudo;
		}
		
		protected double[][] getCount() {
			return this.count;
		}
		
		protected void pseudo() {
			for(int i=0; i<count.length; i++) {
				count[i][0] = (1-pseudo)*Constants._mu_A_m;
				count[i][1] = pseudo*Constants._mu_A_m;
			}
		}
	}

	protected class TransitionUnit {
		private final double[] trans;
		private final double distance;
		private final double pseudo;
		private double jump;
		protected final double[] count;
		private final double[][] cnts_prior; // priors for counting jumps
		
		
		public TransitionUnit(final double distance,
				final int K) {
			this.distance = distance;
			this.pseudo = Math.exp(-distance*Constants._con_base_r);
			this.jump = prior();
			this.trans = new double[K];
			this.count = new double[2];
			this.cnts_prior = new double[K][2];
			this.updatec();
		}

		private void updatec() {
			// TODO Auto-generated method stub
			state.calcProbs(trans, cnts_prior, jump);
		}

		protected void update(double jump) {
			// TODO Auto-generated method stub
			this.jump = jump;
			this.updatec();
		}
		
		protected void update() {
			// TODO Auto-generated method stub
			double c = count[0]+count[1];
			if(c==0) jump = 0.5;
			else jump = count[0]/c;
			this.updatec();
		}
		
		private double prior() {
			// TODO Auto-generated method stub
			double p = new BetaDistribution(Constants.rg, 
					(1-pseudo)*Constants._mu_J_e, pseudo*Constants._mu_J_e).sample();
			if(p==0) p = 1e-12;
			if(p==1) p = 1-1e-12;
			return p;
		}
		
		protected void addCount(int from, int to, double n) {
			int hsc = state.hsc(from, to);
			count[0] += cnts_prior[hsc][0]*n;
			count[1] += cnts_prior[hsc][1]*n;
		}
		
		protected double trans(int from, int to) {
			return trans[state.hsc(from, to)];
		}
		
		protected double[] getTrans() {
			return trans;
		}
		
		protected double[] getCount() {
			return count;
		}
		
		protected double getDistance() {
			return distance;
		}
		
		protected double getJump() {
			return jump;
		}
		
		protected void pseudo() {
			count[0] = (1-pseudo)*Constants._mu_J_m;
			count[1] = pseudo*Constants._mu_J_m;
		}
	}
	
	protected class ViterbiUnit {
		protected double[][] v;
		protected int[][] trace;
		protected double[] logscale;
		protected String[] path_str;
		protected int[] path;
		protected int ends = -1; // end state
		protected double probability = 0;
		
		public ViterbiUnit() {
			// TODO Auto-generated constructor stub
			this.v = new double[M][K];
			this.trace = new int[M-1][K];
			this.logscale = new double[M];
			this.path_str = new String[M];
			this.path = new int[M];
		}

		public double probability() {
			// TODO Auto-generated method stub
			return probability;
		}

		public void finalise() {
			for(int i=0; i<K; i++) {
				if(v[M-1][i]>probability) {
					ends = i;
					probability = v[M-1][i];
				}
			}
			probability = Math.log(probability)+logscale[M-1];
			this.trace();
		}
		
		protected void trace() {
			int tr = ends;
			this.path[M-1] = tr;
			this.path_str[M-1] = state.hsc_str[tr];
			for(int i=M-2; i>=0; i--) {
				tr = trace[i][tr];
				this.path[i] = tr;
				this.path_str[i] = maxMatch(state.hsc_str[tr], 
						this.path_str[i+1]);
			}
			return;
		}
		
		protected String maxMatch(String str1, String str2) {
			// TODO Auto-generated method stub
			final String[] s1 = str1.split("_");
			final String[] s2 = str2.split("_");
			final Map<String, Integer> a = new HashMap<String, Integer>();
			for(int i=0; i<s2.length; i++) a.put(s2[i], i);
			final String[] s = new String[s1.length];
			for(int i=0; i<s.length; i++) {
				if(a.containsKey(s1[i])) {
					s[a.get(s1[i])] = s1[i];
					s1[i] = null;
				}
			}
			int i = 0, j = 0;
			while(i<s.length) {
				if(s1[i]!=null) {
					while(j<s.length && s[j]!=null)
						++j;
					s[j] = s1[i];
				}
				++i;
			}
			return StringUtils.join(s, '_');
		}

		protected void scale(final int i) {
			// TODO Auto-generated method stub
			double[] probs = this.v[i];
			double min = Double.POSITIVE_INFINITY,
					max = Double.NEGATIVE_INFINITY;
			for(int k=0; k<probs.length; k++) {
				if(probs[k]>0) {
					min = probs[k]<min ? probs[k] : min;
					max = probs[k]>max ? probs[k] : max;
				}
			}

			this.logscale[i] += this.logscale[i-1];
			if(min<Constants.threshMin &&
					max<Constants.threshMax) {
				this.logscale[i] += Constants.logThreshMax;
				for(int k=0; k<probs.length; k++)
					probs[k] /= Constants.threshMax;
			}
		}
		
		protected void clear() {
			// TODO Auto-generated method stub
			HiddenMarkovModel.clear(v);
			HiddenMarkovModel.clear(trace);
			HiddenMarkovModel.clear(logscale);
		}
	}
	
	protected class FBUnit { /** forward/backward unit */
		protected double[][] probsMat;
		protected double[] logscale;
		protected double probability;
		protected final boolean backward;

		public FBUnit(boolean backward) {
			this.backward = backward;
			this.probability = 0;
			this.probsMat = new double[M][K];
			this.logscale = new double[M];
		}

		public double probability() {
			// TODO Auto-generated method stub
			return this.probability;
		}
		
		public void probability(double p) {
			// TODO Auto-generated method stub
			if(this.backward)
				this.probability = Math.log(p)+this.logscale[0];
			else
				this.probability = Math.log(p)+
				this.logscale[this.logscale.length-1];
		}
		
		public void probability(double p, double logscale) {
			// TODO Auto-generated method stub
			if(this.backward)
				this.probability = Math.log(p)+this.logscale[0];
			else
				this.probability = Math.log(p)+
				this.logscale[this.logscale.length-1];
			this.probability += logscale;
		}


		protected void scale() {
			// TODO Auto-generated method stub
			this.logscale = new double[this.probsMat.length];
			if(this.backward)
				for(int i=this.logscale.length-1; i>=0; i++)
					this.scale(i);
			else
				for(int i=0; i<this.logscale.length; i++)
					this.scale(i);
		}

		protected void scale(final int i) {
			// TODO Auto-generated method stub
			double[] probs = this.probsMat[i];
			double min = Double.POSITIVE_INFINITY,
					max = Double.NEGATIVE_INFINITY;
			for(int k=0; k<probs.length; k++) {
				if(probs[k]>0) {
					min = probs[k]<min ? probs[k] : min;
					max = probs[k]>max ? probs[k] : max;
				}
			}

			int dv = this.backward ? 1 : -1;
			this.logscale[i] += this.logscale[i+dv];
			if(min<Constants.threshMin &&
					max<Constants.threshMax) {
				this.logscale[i] += 
						Constants.logThreshMax;
				for(int k=0; k<probs.length; k++)
					probs[k] /= Constants.threshMax;
			}
		}
	}
}
