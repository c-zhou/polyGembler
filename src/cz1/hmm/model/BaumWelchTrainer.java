package cz1.hmm.model;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.zip.ZipEntry;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import cz1.hmm.data.DataEntry;
import cz1.math.Combination;
import cz1.math.Permutation;
import cz1.util.Constants;
import cz1.util.Constants.Field;

public class BaumWelchTrainer extends EmissionModel implements ForwardBackwardTrainer {
	private final static Logger myLogger = LogManager.getLogger();
	
	protected final static double mu_J_e = 1e5;
	protected final static double mu_J_m = 0.1;
	protected final static double mu_J_p = 1e-8; //precision
	protected final static double con_base_r = 1e-8;

	private static int bwt_iter = 0;
	private int trans_alter = Integer.MAX_VALUE;
	
	protected StateUnit1 state1;
	protected TransitionUnit[] transition;
	protected ViterbiUnit[] vbs;
	private FBUnit[] forward, backward;
	
	public BaumWelchTrainer(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse,
			Field field,
			int ploidy,
			String[] parents) {
		super(de, seperation, reverse, field, ploidy, parents, false);
		initialise1();
	}
	
	public BaumWelchTrainer() {
		// TODO Auto-generated constructor stub
		super();
	}

	protected void initialise1() {
		// TODO Auto-generated method stub
		this.state1 = new StateUnit1(H);
		this.makeTransitionUnits();
		this.makeViterbiUnits();
		this.makeNaiveTrainer();
	}

	public static BaumWelchTrainer copyOf(EmissionModel model) {
		BaumWelchTrainer hmm = new BaumWelchTrainer();
		hmm.field = model.field;
		hmm.de = model.de;
		hmm.M = model.M; // #markers
		hmm.N = model.N; // #individuals
		hmm.H = model.H; // #founder haplotypes
		hmm.K = model.K; // #compound hidden states
		hmm.Nf1 = model.Nf1;
		hmm.state = model.state;
		hmm.samples = model.samples;
		hmm.parents = model.parents;
		hmm.fi1ter = model.fi1ter;
		hmm.parents_i = model.parents_i;
		hmm.progeny_i = model.progeny_i;
		hmm.weights = model.weights;
		hmm.distance = model.distance;
		hmm.sspace = model.sspace; // state space for each sample
		if(iteration>0) 
			model.switchNumericalSpace(false);
		else
			model.logspace = false;
		hmm.obs = model.obs;
		hmm.pas = model.pas;
		hmm.emission = model.emission;
		hmm.logspace = model.logspace;
		hmm.conjs.addAll(model.conjs);
		hmm.true_id = model.true_id;
		hmm.true_pos = model.true_pos;
		hmm.chrs = model.chrs;
		hmm.chrs_rev = model.chrs_rev;
		hmm.Ms = model.Ms;
		hmm.initialise1();
		hmm.makeNaiveTrainer();
		return hmm;
	}
	
	@Override
	public void makeNaiveTrainer() {
		// TODO Auto-generated method stub
		this.forward = new FBUnit[N];
		for(int i=0; i<N; i++) 
			this.forward[i] = new FBUnit(false);
		this.backward = new FBUnit[N];
		for(int i=0; i<N; i++) 
			this.backward[i] = new FBUnit(true);
		return;
	}

	private void makeViterbiUnits() {
		// TODO Auto-generated method stub
		vbs = new ViterbiUnit[N];
		for(int i=0; i<N; i++) vbs[i] = new ViterbiUnit();
	}

	private void makeTransitionUnits() {
		// TODO Auto-generated method stub
		this.transition = new TransitionUnit[M-1];
		for(int i=0; i<M-1; i++)
			transition[i] = new TransitionUnit(distance[i], state1.confs);
	}
	
	@Override
	public double findPath() {
		// TODO Auto-generated method stub
		// not sure if this makes much difference
		double probability = 0.0;
		for(int i=0; i<N; i++) {
			if(fi1ter[i]) continue;
			
			Integer[] ss = this.sspace.get(i);
			ObUnit[] ob = obs[i];
			double pi = 1.0/ss.length; // init probability
			
			this.vbs[i].clear();
			double[][] v = this.vbs[i].v;
			int[][] trace = this.vbs[i].trace;
			double[] logscale = this.vbs[i].logscale;
			
			double[] emiss = ob[0].emiss;
			TransitionUnit t;
			
			for(int k : ss) v[0][k] = pi*emiss[k];
			logscale[0] = ob[0].getLogScale();
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
				logscale[j] = ob[j].getLogScale();
				vbs[i].scale(j);
			}

			this.vbs[i].finalise();
			this.vbs[i].trace(pas[i].path, pas[i].path_str);
			probability += this.vbs[i].probability();
		}

		return probability;
	}
	
	public double findPath1() {
		// TODO Auto-generated method stub
		// not sure yet if this makes much difference
		// cause we need to calculate a lot of logs
		double probability = 0.0;
		for(int i=0; i<N; i++) {
			if(fi1ter[i]) continue;
			
			Integer[] ss = this.sspace.get(i);
			ObUnit[] ob = obs[i];
			double pi = 1.0/ss.length; // init probability
			
			// pre-calculation for a lower bound for Viterbi path
			// used to get rid of of some paths later
			double lower_bound = Double.NEGATIVE_INFINITY;
			for(int k : ss) {
				double x = pi*ob[0].emiss[k];
				double logscale = ob[0].getLogScale();
				
				for(int j=1; j<M; j++) {
					x *= ob[j].emiss[k]*transition[j-1].trans(k,k);
					logscale += ob[j].getLogScale();
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
			logscale[0] = ob[0].getLogScale();
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
					
					if(Math.log(a)+logscale[j-1]+ob[j].getLogScale()>=lower_bound) {
						trace[j-1][k] = s;
						v[j][k] = a;
						ss_copy.add(k);
					}
				}
				logscale[j] = ob[j].getLogScale();
				vbs[i].scale(j);
			}

			this.vbs[i].finalise();
			this.vbs[i].trace(pas[i].path, pas[i].path_str);
			probability += this.vbs[i].probability();
		}

		return probability;
	}

	@Override
	public void train() {
		// TODO Auto-generated method stub
		++iteration;
		++bwt_iter;
		
		refresh();
		forward();
		backward();
		check();
		em();
	}
	
	@Override
	public void em() {
		// TODO Auto-generated method stub
		for(int i=0; i<M; i++) updateEmiss(i);
		
		if(bwt_iter%trans_alter==0) {
			for(int i=0; i<M-1; i++) updateTrans(i);
			myLogger.info("jump probabilities updated.");
		} else {
			double jump1, jump2;
			for(int i : conjs) {
				jump1 = transition[i].jump;
				updateTrans(i);
				jump2 = transition[i].jump;
				myLogger.info("jump probability at conjunction #"+i+" updated: "+jump1+"->"+jump2+";");
			}
		}
	}

	private void updateTrans(final int i) {
		// TODO Auto-generated method stub
		FBUnit fw1, bw1;
		ObUnit ob1;
		double exp_c, exp, count;
		Integer[] ss;

		TransitionUnit t1;
		t1 = transition[i];
		t1.pseudoCount();
		for(int j=0;j<N; j++) {
			if(fi1ter[j]) continue;
			
			ss = sspace.get(j);
			fw1 = forward[j];
			bw1 = backward[j];
			ob1 = obs[j][i+1];
			exp_c = fw1.logscale[i]+
					bw1.logscale[i+1]+
					ob1.getLogScale()-
					fw1.probability;

			if(exp_c>Constants.MAX_EXP_DOUBLE) { 
				for(int a : ss) {
					for(int b : ss) { 
						count = Math.exp(Math.log(
								fw1.probsMat[i][a]*
								t1.trans(a, b)*
								ob1.emiss[b]*
								bw1.probsMat[i+1][b])+
								exp_c);
						t1.addCount(a, b, count);
					}
				}
			} else {
				exp = Math.exp(exp_c);
				for(int a : ss) {
					for(int b : ss) { 
						count = fw1.probsMat[i][a]*
								t1.trans(a, b)*
								ob1.emiss[b]*
								bw1.probsMat[i+1][b]*
								exp;
						t1.addCount(a, b, count);
					}
				}
			}
		}
		t1.update();
	}

	private void updateEmiss(final int i) {
		// TODO Auto-generated method stub
		FBUnit fw1, bw1;
		ObUnit ob1;
		double exp_c, exp, count, coeff;
		Integer[] ss;
		int acnt, bcnt;
	
		EmissionUnit e1;
		e1 = emission[i];
		e1.pseudoCount();
		for(int j=0;j<N; j++) {
			if(fi1ter[j]) continue;
			
			ss = sspace.get(j);
			fw1 = forward[j];
			bw1 = backward[j];
			ob1 = obs[j][i];
			acnt = ob1.getAa();
			bcnt = ob1.getCov()-acnt;
			coeff = weights[j==parents_i[0]||j==parents_i[1]?0:1];
			exp_c = fw1.logscale[i]+
					bw1.logscale[i]-
					fw1.probability;

			if(exp_c>Constants.MAX_EXP_DOUBLE) {
				for(int a : ss) {
					count = coeff*
							Math.exp(Math.log(
							fw1.probsMat[i][a]*
							bw1.probsMat[i][a])+
							exp_c);
					e1.addCount(a, acnt, bcnt, count);
				}
			} else {
				exp = Math.exp(exp_c);
				for(int a : ss) {
					count = coeff*
							fw1.probsMat[i][a]*
							bw1.probsMat[i][a]*
							exp;
					e1.addCount(a, acnt, bcnt, count);
				}
			}
		}
		e1.update();
	}

	@Override
	public void backward() {
		// TODO Auto-generated method stub
		for(int i=0; i<N; i++) {
			if(fi1ter[i]) continue;
			
			Integer[] ss = sspace.get(i);
			double[][] probsMat = backward[i].probsMat;
			double[] logscale = backward[i].logscale;
			ObUnit[] ob = obs[i];
			
			double[] emiss;
			TransitionUnit t;
			
			for(int k : ss) probsMat[M-1][k] = 1.0;
			logscale[M-1] = 0;
			double tmp; 
			
			for(int j=M-2; j>=0; j--) {	
				emiss = ob[j+1].emiss;
				t = transition[j];
				for(int k : ss) {
					tmp = 0;
					for(int z : ss) 
						tmp += t.trans(k, z)*emiss[z]*probsMat[j+1][z];
					probsMat[j][k] = tmp;
				}
				
				logscale[j] = ob[j+1].getLogScale();
				backward[i].scale(j);
			}
			
			double pi = 1.0/ss.length;
			double p = 0.0;
			emiss = ob[0].emiss;
			for(int z : ss)
				p += pi*emiss[z]*probsMat[0][z];
			backward[i].probability(p, ob[0].getLogScale());
		}
		return;
	}

	@Override
	public void forward() {
		// TODO Auto-generated method stub
		for(int i=0; i<N; i++) {
			if(fi1ter[i]) continue;
			
			Integer[] ss = sspace.get(i);
			double pi = 1.0/ss.length;
			
			double[][] probsMat = forward[i].probsMat;
			double[] logscale = forward[i].logscale;
			ObUnit[] ob = obs[i];
			
			double[] emiss = ob[0].emiss;
			TransitionUnit t;
			
			for(int k : ss) probsMat[0][k] = pi*emiss[k];
			logscale[0] = ob[0].getLogScale();
			double tmp; 
			
			for(int j=1; j<M; j++) {
				
				emiss = ob[j].emiss;
				t = transition[j-1];
				
				for(int k : ss) {
					tmp = 0;
					for(int z : ss)
						tmp += probsMat[j-1][z]
								*t.trans(z, k);
					probsMat[j][k] = emiss[k]*tmp;
				}
				
				logscale[j] = ob[j].getLogScale();
				forward[i].scale(j);
			}
			forward[i].probability(StatUtils.sum(probsMat[M-1]));
		}
		return;
	}

	@Override
	protected double loglik(int fromIndex, int toIndex) {
		// TODO Auto-generated method stub
		if(fromIndex<0||fromIndex>=toIndex||toIndex>M)
			throw new RuntimeException("!!!");
		
		int m = toIndex - fromIndex;
		double probability = 0;

		for(int i=0; i<N; i++) {
			if(fi1ter[i]) continue;
			
			Integer[] ss = sspace.get(i);
			double pi = 1.0/ss.length;
			
			double[][] probsMat = new double[m][K];
			double[] logscale = new double[m];
			ObUnit[] ob = obs[i];
			
			double[] emiss = ob[fromIndex].emiss;
			TransitionUnit t;
			
			for(int k : ss) probsMat[0][k] = pi*emiss[k];
			logscale[0] = ob[fromIndex].getLogScale();
			double tmp; 
			
			for(int j=1; j<m; j++) {
				
				emiss = ob[fromIndex+j].emiss;
				t = transition[fromIndex+j-1];
				
				for(int k : ss) {
					tmp = 0;
					for(int z : ss)
						tmp += probsMat[j-1][z]
								*t.trans(z, k);
					probsMat[j][k] = emiss[k]*tmp;
				}
				
				logscale[j] = ob[fromIndex+j].getLogScale();
				scale(logscale, probsMat, j);
			}
			probability += Math.log(StatUtils.sum(probsMat[m-1]))+logscale[m-1];
		}

		return probability;	
	}
	
	protected void scale(final double[] logscale, final double[][] probsMat, final int i) {
		// TODO Auto-generated method stub
		double[] probs = probsMat[i];
		double min = Double.POSITIVE_INFINITY,
				max = Double.NEGATIVE_INFINITY;
		for(int k=0; k<probs.length; k++) {
			if(probs[k]>0) {
				min = probs[k]<min ? probs[k] : min;
				max = probs[k]>max ? probs[k] : max;
			}
		}

		logscale[i] += logscale[i-1];
		if(min<Constants.threshMin &&
				max<Constants.threshMax) {
			logscale[i] += 
					Constants.logThreshMax;
			for(int k=0; k<probs.length; k++)
				probs[k] /= Constants.threshMax;
		}
	}
	
	@Override
	public void check() {
		// TODO Auto-generated method stub
		if(iteration==0) return;
		for(int i=0; i<this.forward.length; i++) {
			double r = Math.abs(forward[i].probability-backward[i].probability);
			if(r>1e-6) 
				throw new RuntimeException("Different likelihood by forward and backward algorithm: forward, "+
						forward[i].probability+"; backward, "+backward[i].probability);
		}
	}

	public double loglik1() {
		// TODO Auto-generated method stub
		if(iteration==0)
			return Double.NEGATIVE_INFINITY;
		else {
			double probability = 0;
			for(int i=0; i<N; i++) 
				probability += weights[i==parents_i[0]||i==parents_i[1]?0:1]*this.backward[i].probability;
			return probability;
		}
	}

	@Override
	public double loglik() {
		// TODO Auto-generated method stub
		if(iteration==0)
			return Double.NEGATIVE_INFINITY;
		else {
			double probability = 0;
			for(int i=0; i<N; i++)
				probability += weights[i==parents_i[0]||i==parents_i[1]?0:1]*this.forward[i].probability;
			return probability;
		}
	}

	@Override
	public void write(String output, String experiment, String scaff) {
		// TODO Auto-generated method stub
		this.findPath();

		String root = experiment+"."+scaff+"."+System.nanoTime();
		ModelWriter1 writer1 = new ModelWriter1(output+"/"+root+".zip");
		writer1.writeHaplotype();
		writer1.writeDosage();
		writer1.writeGenotype();
		writer1.writeEmissionModel();
		writer1.writeTransitionModel();
		writer1.writeSNP();
		writer1.writeRunInfo();
		writer1.close();		
	}
	
	protected class ModelWriter1 extends ModelWriter {

		public ModelWriter1(String file) {
			// TODO Auto-generated constructor stub
			super(file);
		}
		
		private void writeTransitionModel() {
			try {
				out.putNextEntry(new ZipEntry("transition.txt"));
				for(int i=0; i<M-1; i++)
					out.write((de.getId()+"_"+de.getPosition()[i]+"\t\t\t"+transition[i].getJump()+"\n").getBytes());
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	protected class StateUnit1 extends StateUnit {
		private final int confs; // #confs
		private final int[][] confs_hsc; // confs of compound hidden states
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
		
		public StateUnit1(final int h) {
			super(h);
			confs = (h/2+1)*(h/2+2)/2+1;
			confs_hsc = new int[hsc.length][hsc.length];
			confs_cnt = new int[confs][h+1];
			
			this.configure();
		}

		private void configure() {
			// TODO Auto-generated method stub
			final int h = hs.length/4;
			final int k = hsc.length;
			// make confs look up table
			int[][] confs_tab = new int[confs][2];
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
			
			// hidden states for parents 
			int[] probs = confs_cnt[0];
			int h2 = h*2;
			for(int i=0; i<=h2; i++)
				probs[i] = f(h2, h2, i);
			
			// hidden states for f1 progeny
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
			// two sets of /n elements with /a elements in common
			// calculates #permutation with /b elements matches 
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

		public void calcProbs(double[] trans, double[][] cnts, double p) {
			clear(cnts);
			
			final int h = hs.length/2;
			int[] cnt;
			final double stay = 1-p;
			final double jump = p/(h-1);
			double p1, pA;
			for(int i=0; i<confs; i++) {
				cnt = confs_cnt[i];
				pA = 0;
				for(int j=0; j<=h; j++) {
					p1 = cnt[j]*Math.pow(jump, h-j)*Math.pow(stay, j);
					cnts[i][0] += (h-j)*p1;
					cnts[i][1] += j*p1;
					pA += p1;
				}
				
				trans[i] = pA;
				cnts[i][0] /= pA;
				cnts[i][1] /= pA;
			}
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

	protected class TransitionUnit {
		private final double[] trans;
		private final double distance;
		private final double base_r;
		private double jump;
		protected final double[] count;
		private final double[][] cnts_prior; // priors for counting jumps
		
		
		public TransitionUnit(final double distance,
				final int K) {
			this.distance = distance;
			this.base_r = Math.exp(-2*distance*con_base_r);
			this.jump = prior();
			this.trans = new double[K];
			this.count = new double[2];
			this.cnts_prior = new double[K][2];
			this.updatec();
		}

		private void updatec() {
			// TODO Auto-generated method stub
			state1.calcProbs(trans, cnts_prior, jump);
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
			jump = Math.max(jump, mu_J_p);
			jump = Math.min(jump, 0.5-mu_J_p);
			this.updatec();
		}
		
		private double prior() {
			// TODO Auto-generated method stub
			double p = new BetaDistribution(Constants.rg, 
					(1-base_r)*mu_J_e, base_r*mu_J_e).sample()*0.5;
			p = Math.max(p, mu_J_p);
			p = Math.min(p, 0.5-mu_J_p);
			return p;
		}
		
		protected void addCount(int from, int to, double n) {
			int hsc = state1.hsc(from, to);
			count[0] += cnts_prior[hsc][0]*n;
			count[1] += cnts_prior[hsc][1]*n;
		}
		
		protected double trans(int from, int to) {
			return trans[state1.hsc(from, to)];
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
		
		protected void pseudoCount() {
			count[0] = (1-base_r)*mu_J_m;
			count[1] = base_r*mu_J_m;
		}
	}
	
	protected class ViterbiUnit {
		protected double[][] v;
		protected int[][] trace;
		protected double[] logscale;
		protected int ends = -1; // end state
		protected double probability = 0;
		
		public ViterbiUnit() {
			// TODO Auto-generated constructor stub
			super();
			this.v = new double[M][K];
			this.trace = new int[M-1][K];
			this.logscale = new double[M];
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
		}
		
		protected void trace(int[] path, String[] path_str) {
			int tr = ends;
			path[M-1] = tr;
			path_str[M-1] = state1.hsc_str[tr];
			for(int i=M-2; i>=0; i--) {
				tr = trace[i][tr];
				path[i] = tr;
				path_str[i] = maxMatch(state1.hsc_str[tr], 
						path_str[i+1]);
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
			BaumWelchTrainer.clear(v);
			BaumWelchTrainer.clear(trace);
			BaumWelchTrainer.clear(logscale);
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

	public void modifyTransAlter(int trans_alter) {
		// TODO Auto-generated method stub
		this.trans_alter = trans_alter;
	}
}
