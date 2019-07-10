package cz1.hmm.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.log4j.Logger;

import cz1.hmm.data.DataEntry;
import cz1.math.Combination;
import cz1.math.SaddlePointExpansion;
import cz1.util.Constants;
import cz1.util.Constants.Field;

public class EmissionModel {
	protected final static Logger myLogger = Logger.getLogger(EmissionModel.class);
	protected static int iteration = 0;
	
	protected Field field;
	protected DataEntry de;

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
	
	protected boolean logspace;
	
	public EmissionModel(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse,
			Field field,
			boolean logspace) {
		this.field = field;
		this.de = this.catDE(de, seperation, reverse);
		this.logspace = logspace;
		this.initialise();
	}
	
	public EmissionModel(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse,
			Field field) {
		this(de, seperation, reverse, field, true);
	}
	
	public EmissionModel() {
		// TODO Auto-generated constructor stub
		super();
	}

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

	protected void initialise() {
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
	
	protected void switchNumericalSpace() {
		// TODO Auto-generated method stub
		switchNumericalSpace(!logspace);
		logspace = !logspace;
	}
	
	protected void switchNumericalSpace(boolean logspace) {
		// TODO Auto-generated method stub
		if(logspace==this.logspace)
			return;
		for(int i=0; i<N; i++)
			for(int j=0; j<M; j++) obs[i][j].switchNumericSpace();
		this.logspace = logspace;
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
	
	protected static void clear(double[] array) {
		// TODO Auto-generated method stub
		Arrays.fill(array, 0);
	}
	
	public int iteration() {
		// TODO Auto-generated method stub
		return iteration;
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
			
			if(!logspace) switchToNormalSpace();
		}

		public void switchNumericSpace() {
			// TODO Auto-generated method stub
			if(logspace)
				switchToNormalSpace();
			else
				switchToLogSpace();
		}

		private void switchToNormalSpace() {
			// TODO Auto-generated method stub
			// if not in log space
			logscale[0] = emiss[0];
			emiss[0] = 1.0;
			logscale[1] = emiss[1];
			emiss[1] = 1.0;
			logscale[2] = StatUtils.max(emiss, 2, K-2);
			for(int i=2; i<K; i++) {
				emiss[i] = Math.exp(emiss[i]-logscale[2]);
			}
		}

		private void switchToLogSpace() {
			// TODO Auto-generated method stub
			emiss[0] = logscale[0];
			emiss[1] = logscale[1];
			for(int i=2; i<K; i++) {
				emiss[i] = Math.log(emiss[i])+logscale[2];
			}
			Arrays.fill(logscale, 0);
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

	protected class StateUnit {
		protected final char[] hs; // hidden states
		protected final int[][] hsc; // compound hidden states
		protected final String[] hsc_str; // compound hidden states str
		
		public StateUnit(final int h) {
			hs = makeHs(h*2);
			hsc = makeHsc();
			hsc_str = makeHscStr();
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
		
		protected char[] getHs() {
			return this.hs;
		}
		
		protected int[][] getHsc() {
			return this.hsc;
		}
		
		protected String[] getHscStr() {
			return this.hsc_str;
		}
	}

	protected class EmissionUnit {
		protected final String[] allele;
		protected final double[] emiss;
		protected final double[] emissc;
		protected final double bfrac;
		protected final double[][] count;
		protected final double[][][] cnts_prior; // priors for counting
		
		public EmissionUnit(final String[] allele,
				int H,
				int K,
				final double bfrac) {
			this.allele = allele;
			this.emiss = new double[H];
			this.emissc = new double[K];
			this.bfrac = bfrac;
			this.count = new double[H][allele.length];
			this.cnts_prior = new double[K][H/2][2];
			this.prior();
			this.updatec();
		}

		protected void prior() {
			// TODO Auto-generated method stub
			BetaDistribution beta = new BetaDistribution(Constants.rg, 
					(1-bfrac)*Constants._mu_A_e, bfrac*Constants._mu_A_e);
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
		
		protected double getBfrac() {
			return this.bfrac;
		}
		
		protected double[][] getCount() {
			return this.count;
		}
		
		protected void pseudo() {
			for(int i=0; i<count.length; i++) {
				count[i][0] = (1-bfrac)*Constants._mu_A_m;
				count[i][1] = bfrac*Constants._mu_A_m;
			}
		}
	}
}
