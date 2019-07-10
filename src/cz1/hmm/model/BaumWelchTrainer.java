package cz1.hmm.model;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.log4j.Logger;

import cz1.hmm.data.DataEntry;
import cz1.math.Combination;
import cz1.math.Permutation;
import cz1.util.Constants;
import cz1.util.Constants.Field;

public class BaumWelchTrainer extends EmissionModel implements ForwardBackwardTrainer {
	protected final static Logger myLogger = Logger.getLogger(BaumWelchTrainer.class);

	final boolean updateEmiss;
	final boolean updateTrans;
	
	protected StateUnit1 state;
	protected TransitionUnit[] transition;
	protected ViterbiUnit[] vbs;
	private FBUnit[] forward, backward;
	
	public BaumWelchTrainer(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse,
			Field field) {
		super(de, seperation, reverse, field, false);
		this.updateEmiss = true;
		this.updateTrans = false;
		initialise1();
	}
	
	public BaumWelchTrainer(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse,
			Field field,
			boolean updateEmiss, 
			boolean updateTrans) {
		super(de, seperation, reverse, field, false);
		this.updateEmiss = updateEmiss;
		this.updateTrans = updateTrans;
		initialise1();
	}
	
	public BaumWelchTrainer() {
		// TODO Auto-generated constructor stub
		super();
		this.updateEmiss = true;
		this.updateTrans = false;
	}

	public BaumWelchTrainer(boolean updateEmiss, boolean updateTrans) {
		// TODO Auto-generated constructor stub
		super();
		this.updateEmiss = updateEmiss;
		this.updateTrans = updateTrans;
	}
	
	protected void initialise1() {
		// TODO Auto-generated method stub
		this.state = new StateUnit1(H);
		this.makeTransitionUnits();
		this.makeViterbiUnits();
		this.makeNaiveTrainer();
	}

	public static BaumWelchTrainer copyOf(EmissionModel model) {
		return copyOf(model, true, false);
	}
	
	public static BaumWelchTrainer copyOf(EmissionModel model, 
			boolean updateEmiss, boolean updateTrans) {
		BaumWelchTrainer hmm = new BaumWelchTrainer(updateEmiss, updateTrans);
		hmm.field = model.field;
		hmm.de = model.de;
		hmm.M = model.M; // #markers
		hmm.N = model.N; // #individuals
		hmm.H = model.H; // #founder haplotypes
		hmm.K = model.K; // #compound hidden states
		hmm.samples = model.samples;
		hmm.parents = model.parents;
		hmm.parents_i = model.parents_i;
		hmm.progeny_i = model.progeny_i;
		hmm.distance = model.distance;
		hmm.sspace = model.sspace; // state space for each sample
		if(iteration>0) 
			model.switchNumericalSpace(false);
		else
			model.logspace = false;
		hmm.obs = model.obs;
		hmm.emission = model.emission;
		hmm.logspace = model.logspace;
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
			transition[i] = new TransitionUnit(distance[i], state.confs);
	}
	
	@Override
	public double findPath() {
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
	
	public double findPath1() {
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

	@Override
	public void train() {
		// TODO Auto-generated method stub
		refresh();
		forward();
		backward();
		check();
		em();
		++iteration;
	}
	
	@Override
	public void em() {
		// TODO Auto-generated method stub
		
		FBUnit fw1, bw1;
		ObUnit ob1;
		double exp_c, exp, count;
		Integer[] ss;
		
		if(updateEmiss) {
			
			int acnt, bcnt;
			EmissionUnit e1;
			for(int i=0; i<M; i++) {
				e1 = emission[i];
				e1.pseudo();

				for(int j=0;j<N; j++) {
					ss = sspace.get(j);
					fw1 = forward[j];
					bw1 = backward[j];
					ob1 = obs[j][i];
					acnt = ob1.getAa();
					bcnt = ob1.getCov()-acnt;
					exp_c = fw1.logscale[i]+
							bw1.logscale[i]-
							fw1.probability;

					if(exp_c>Constants.MAX_EXP_DOUBLE) {
						for(int a : ss) {
							count = Math.exp(Math.log(
									fw1.probsMat[i][a]*
									bw1.probsMat[i][a])+
									exp_c);
							e1.addCount(a, count*acnt, count*bcnt);
						}
					} else {
						exp = Math.exp(exp_c);
						for(int a : ss) {
							count = fw1.probsMat[i][a]*
									bw1.probsMat[i][a]*
									exp;
							e1.addCount(a, count*acnt, count*bcnt);
						}
					}
				}
				e1.update();
			}
		}
		
		if(updateTrans) {

			TransitionUnit t1;
			for(int i=0; i<M-1; i++) {
				// transitions
				t1 = transition[i];
				t1.pseudo();
				for(int j=0;j<N; j++) {
					ss = sspace.get(j);
					fw1 = forward[j];
					bw1 = backward[j];
					ob1 = obs[j][i+1];
					exp_c = fw1.logscale[i]+
							bw1.logscale[i+1]+
							ob1.getLogScale(j)-
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
		}
	}

	@Override
	public void backward() {
		// TODO Auto-generated method stub
		for(int i=0; i<N; i++) {
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
				
				logscale[j] = ob[j+1].getLogScale(i);
				backward[i].scale(j);
			}
			
			double pi = 1.0/ss.length;
			double p = 0.0;
			emiss = ob[0].emiss;
			for(int z : ss)
				p += pi*emiss[z]*probsMat[0][z];
			backward[i].probability(p, ob[0].getLogScale(i));
		}
		return;
	}

	@Override
	public void forward() {
		// TODO Auto-generated method stub
		for(int i=0; i<N; i++) {
			
			Integer[] ss = sspace.get(i);
			double pi = 1.0/ss.length;
			
			double[][] probsMat = forward[i].probsMat;
			double[] logscale = forward[i].logscale;
			ObUnit[] ob = obs[i];
			
			double[] emiss = ob[0].emiss;
			TransitionUnit t;
			
			for(int k : ss) probsMat[0][k] = pi*emiss[k];
			logscale[0] = ob[0].getLogScale(i);
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
				
				logscale[j] = ob[j].getLogScale(i);
				forward[i].scale(j);
			}
			forward[i].probability(StatUtils.sum(probsMat[M-1]));
		}
		return;
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
			for(FBUnit bw : this.backward) probability += bw.probability;
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
			for(FBUnit fw : this.forward) probability += fw.probability;
			return probability;
		}
	}

	@Override
	public void write(String output, String experiment, String scaff) {
		// TODO Auto-generated method stub
		this.findPath();

		String root = experiment+"."+scaff+"."+System.nanoTime();

		try {
			ZipOutputStream out = new ZipOutputStream(new BufferedOutputStream(new 
					FileOutputStream(output+"/"+root+".zip"), 65536));

			out.putNextEntry(new ZipEntry("phasedStates/"+experiment+".txt"));
			out.write((""+this.loglik()+"\n").getBytes());
			out.write((""+M+"\n").getBytes());
			for(int i=0; i<N; i++) {
				String[] path = this.vbs[i].path_str;
				List<String[]> path_s = new ArrayList<String[]>();
				for(int k=0; k<path.length; k++)
					path_s.add(path[k].split("_"));
				for(int k=0; k<Constants._ploidy_H; k++) {
					out.write(("# id "+this.samples[i]+":"+(k+1)+"\t\t\t").getBytes());
					for(int s=0; s<path_s.size(); s++)
						out.write(path_s.get(s)[k].getBytes());
					out.write("\n".getBytes());
				}
			}
			
			out.putNextEntry(new ZipEntry("phased_genotypes"+experiment+".txt"));
			for(int i=0; i<N; i++) {
				int[] path = this.vbs[i].path;
				for(int j=0; j<M; j++) {
					int[] states = this.state.getHsc()[path[j]];
					double[] emiss = this.emission[j].getEmiss();
					int dosa = 0;
					for(int k=0; k<H; k++)
						if(emiss[states[k]]<0.5)
							++dosa;
					out.write((""+dosa+"\t").getBytes());
				}
				out.write("\n".getBytes());
			}

			out.putNextEntry(new ZipEntry("results_hmm/emissionModel.txt"));
			for(int i=0; i<M; i++) {
				double[] emiss = this.emission[i].getEmiss();
				String[] allele = this.de.getAllele().get(i);
				out.write((this.de.getId()+"_"+this.de.getPosition()[i]+"\t\t\t").getBytes());
				for(int j=0; j<emiss.length; j++) {
					out.write((this.state.getHs()[j]+"-> {").getBytes());
					out.write((allele[0]+","+emiss[j]+";").getBytes());
					out.write((allele[1]+","+(1-emiss[j])+";").getBytes());
					out.write("} ".getBytes());
				}
				out.write("\n".getBytes());
			}

			out.putNextEntry(new ZipEntry("results_hmm/transitionModel.txt"));
			for(int i=0; i<M-1; i++)
				out.write((this.de.getId()+"_"+this.de.getPosition()[i]+"\t\t\t"+this.transition[i].getJump()+"\n").getBytes());
			
			out.putNextEntry(new ZipEntry("stderr_true"));

			out.write("cz1.model.HiddenMarkovModel:\n".getBytes());
			out.write("cz1.model.HiidenMarkovModel$EM:\n".getBytes());
			out.write(("log prob is "+this.loglik()+" at "+iteration+"\n").getBytes());
			out.write(("random seed is "+Constants.seed).getBytes());

			out.putNextEntry(new ZipEntry("snp_"+experiment+".txt"));

			StringBuilder os = new StringBuilder();
			List<String[]> allele = de.getAllele();
			double[] position = de.getPosition();
			String[] markers_id = de.getMarker();
			String id = de.getId();
			for(int i=0; i<position.length; i++) {
				os.setLength(0);
				os.append(id);
				os.append("\t");
				int p = (int) position[i];
				os.append(p);
				os.append("\t");
				os.append(p);
				os.append("\t");
				os.append(markers_id[i]);
				os.append("\t");
				for(int j=0; j<allele.get(i).length; j++) {
					os.append(allele.get(i)[j]);
					os.append("\t");
				}
				os.append(experiment);
				os.append("\t");
				os.append("false\n");
				out.write(os.toString().getBytes());
			}			
			out.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
	}
	
	protected class StateUnit1 extends StateUnit {
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
		
		public StateUnit1(final int h) {
			super(h);
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

		public void calcProbs(double[] trans, double[][] cnts, double p) {
			clear(cnts);
			
			final int h = hs.length/2;
			trans[0] = 1.0;
			cnts[0][0] = 0.0;
			cnts[0][1] = h;
			int[] cnt;
			final double stay = 1-p;
			final double jump = p/(h-1);
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
			//double p = new BetaDistribution(Constants.rg, 
			//		(1-pseudo)*Constants._mu_J_e, pseudo*Constants._mu_J_e).sample();
			//if(p==0) p = 1e-16;
			//if(p==1) p = 1-1e-16;
			//return p;
			return 1-pseudo;
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
}
