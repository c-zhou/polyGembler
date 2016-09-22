package cz1.model;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import cz1.data.DataEntry;
import cz1.util.Algebra;
import cz1.util.Combination;
import cz1.util.Constants;
import cz1.util.IO;
import cz1.util.Permutation;
import cz1.util.Dirichlet;

public class HiddenMarkovModelDupSHEXA extends HiddenMarkovModel {
	private final static Logger logger = LogManager.getLogger(HiddenMarkovModelDupSHEXA.class.getName());
	private final static Runtime runtime = Runtime.getRuntime();

	private final DataEntry de;
	private final DP[][] dp;
	private final ST[] statespace;
	private final OB[] obspace;
	private final String[] hs; // all hidden states 1,2,3,4,5,6,7,8
	private final String[] str_statespace;
	private final double[] bfrac;
	private final boolean train_exp;
	private final Set<Integer> exp_b = new HashSet<Integer>();

	private double[] distance;
	private TP[] transProbs;
	private EP[] emissProbs;

	private String[] sample;
	private String[] parent;
	private int N;
	private int M;
	private Viterbi vb[];

	private TP[] compoundTransProbs;
	private EP[] compoundEmissProbs;

	private FB[] forward, backward;
	
	private final int[] pedigree; //0 and 1 represent two parents, and -1s represent F1 offspring
	private final List<Integer[]> validStateSpace; //valid statespace for each sample

	private final Map<String, int[][]> spMap;

	public HiddenMarkovModelDupSHEXA(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse,
			boolean trainExp) {
		this.train_exp = trainExp;
		this.de = this.catDE(de, seperation, reverse);
		this.initialise();
		this.hs = this.makeHS();
		this.bfrac = this.bfrac();
		this.obspace = this.makeOB();
		this.statespace = this.makeST();
		this.str_statespace = this.makeSS();
		this.dp = this.makeDP();
		this.spMap = this.makeSpMap();
		this.pedigree = this.pedigree();
		this.validStateSpace = this.validStateSpace();
		this.makeBWT();
		//this.makeViterbi();
		//this.print();
		//this.train();
	}

	private DataEntry catDE(DataEntry[] de2, 
			double[] seperation, 
			boolean[] reverse) {
		// TODO Auto-generated method stub
		for(int i=0; i<de2.length; i++)
			if(reverse[i]) de2[i].reverse();
		for(int i=1; i<de2.length; i++) {
			if(this.train_exp)
				this.exp_b.add(de2[0].modelLength());
			de2[0].addAll(de2[i], seperation[i-1]);
		}
		return de2[0];
	}

	private double[] bfrac() {
		// TODO Auto-generated method stub
		String[] sample = this.de.getSample();
		Set<Integer> p = new HashSet<Integer>();
		for(int i=0; i<sample.length; i++)
			if(isParent(this.sample[i]))
				p.add(i);
		List<List<char[]>> gt = this.de.getGenotype();
		double[] baf = new double[gt.size()];
		for(int i=0; i<baf.length; i++) {
			double cnt = 0, bcnt = 0;
			String[] allele = this.de.getAllele().get(i);
			for(int j : p) {
				char[] g = gt.get(i).get(j);
				for(int k=0; k<g.length; k++) {
					cnt += 1;
					bcnt += ((""+g[k]).equals(allele[1]) ? 1 : 0);
				}
			}
			if(bcnt==cnt || bcnt==0) baf[i] = .5;
			else baf[i] = bcnt/cnt;
		}
		return baf;
	}

	private String[] makeSS() {
		// TODO Auto-generated method stub
		String[] str_states = new String[this.statespace.length];
		for(int i=0; i<str_states.length; i++) 
			str_states[i] = StringUtils.join(
					this.statespace[i].state_str,"_");
		return str_states;
	}

	private Map<String, int[][]> makeSpMap() {
		// TODO Auto-generated method stub
		Map<String, int[][]> spMap = 
				new HashMap<String, int[][]>();
				for(ST s : this.statespace) {
					ST[] perms = s.stateperm();
					int[][] p = new int[perms.length] 
							[s.state.length];
					for(int i=0; i<p.length; i++)
						p[i] = perms[i].state;
					spMap.put(StringUtils.join(
							s.state_str,"_"),p);
				}
				return spMap;
	}

	private List<Integer[]> validStateSpace() {
		// TODO Auto-generated method stub
		Integer[] validP0, validP1, validF1;
		validP0 = new Integer[]{0,1};
		validP1 = new Integer[]{0,2};
		validF1 = new Integer[this.statespace.length-2];
		validF1[0] = 0;
		for(int i=1; i<validF1.length; i++)
			validF1[i] = i+2;
		List<Integer[]> validSS = new ArrayList<Integer[]>();
		for(int i=0; i<this.pedigree.length; i++)
			if(this.pedigree[i]==0)
				validSS.add(validP0);
			else if(this.pedigree[i]==1)
				validSS.add(validP1);
			else
				validSS.add(validF1);
		return validSS;
	}

	private int[] pedigree() {
		// TODO Auto-generated method stub
		int pedigree[] = new int[this.sample.length];
		Arrays.fill(pedigree, -1);
		for(int i=0; i<pedigree.length; i++)
			for(int j=0; j<this.parent.length; j++)
				if(this.sample[i].equals(this.parent[j]))
					pedigree[i] = j;
		return pedigree;
	}

	private void initialise() {
		// TODO Auto-generated method stub
		this.sample = de.getSample();
		this.parent = Constants._founder_haps.split(":");
		this.N = this.sample.length;
		this.M = this.de.getAllele().size();
		double[] position = de.getPosition();
		this.distance = new double[this.M-1];
		for(int i=0; i<distance.length; i++)
			distance[i] = position[i+1]-position[i];
	}

	public void train() {
		Constants.iteration++;
		logger.info("###################");
		logger.info("train: "+Constants.iteration);
		long[] tic = new long[10];
		int k=0;
		tic[k++] = System.nanoTime();
		this.makeForward();
		tic[k++] = System.nanoTime();
		logger.info("forward done "+(tic[k-1]-tic[k-2])+"ns");
		this.makeBackward();
		tic[k++] = System.nanoTime();
		logger.info("backward done "+(tic[k-1]-tic[k-2])+"ns");
		this.checkFW();
		tic[k++] = System.nanoTime();
		this.EM();
		tic[k++] = System.nanoTime();
		logger.info("EM algorithm "+(tic[k-1]-tic[k-2])+"ns");
		this.makeCompoundEP();
		tic[k++] = System.nanoTime();
		logger.info("make compound ep "+(tic[k-1]-tic[k-2])+"ns");
		this.makeCompoundTP();
		tic[k++] = System.nanoTime();
		logger.info("make compound tp "+(tic[k-1]-tic[k-2])+"ns");
		//this.makeViterbi();
		this.updateDP();
		tic[k++] = System.nanoTime();
		logger.info("update dp "+(tic[k-1]-tic[k-2])+"ns");
		this.printExpTrain();
		return;
	}

	private void printExpTrain() {
		// TODO Auto-generated method stub
		for(final Integer i : this.exp_b) {
			double cumR = 0.0;
			double[][] tp_i = this.transProbs[i].probsMat;
			for(int j=0; j<tp_i.length; j++)
				for(int k=0; k<tp_i[j].length; k++)
					if(j!=k) cumR+=tp_i[j][k];
			logger.info("********** sep "+i+" "+cumR/(tp_i.length-1));
		}
	}

	private void updateDP() {
		// TODO Auto-generated method stub
		for(int i=0; i<M; i++) {
			String[] dosaS = this.obspace[i].dosaS;
			int _d_ = dosaS.length;
			Map<String, Integer[]> dgMap = this.obspace[i].dgMap;
			EP ep1 = this.compoundEmissProbs[i];
			for(int j=0; j<N; j++) {
				Integer[] vst = this.validStateSpace.get(j);
				
				DP dp1 = this.dp[i][j];
				for(int k : vst) {
					double ss = 0;
					for(int a=0; a<_d_; a++) {
						Integer[] ix = dgMap.get(dosaS[a]);
						for(int c : ix)
							ss += dp1.likelihood[a]
									*ep1.probsMat[k][c];
					}
					dp1.emiss[k] = ss;
					//this.dp[i][j].emiss[k] = 
					//		StatUtils.sum(this.dp[i][j].emissG[k]);
					//for(int b=0; b<this.dp[i][j].emissG[k].length; b++)
					//	this.dp[i][j].emissG[k][b] /= this.dp[i][j].emiss[k];
				}
			}
		}
	}

	private double[][] normalise(double[][] mat, boolean byrow,
			boolean logspace) {
		// TODO Auto-generated method stub
		if(!byrow) mat = transpose(mat);
		for(int i=0; i<mat.length; i++) {
			double s = Algebra.sum(mat[i]);
			if(s==0) continue;
			for(int j=0; j<mat[i].length; j++)
				mat[i][j] /= s;
		}
		if(!byrow) mat = transpose(mat);
		if(logspace) mat = logspace(mat);
		return mat;
	}

	private double[][] logspace(double[][] mat) {
		// TODO Auto-generated method stub
		for(int i=0; i<mat.length; i++)
			mat[i] = logspace(mat[i]);
		return mat;
	}

	private double[] logspace(double[] array) {
		// TODO Auto-generated method stub
		for(int i=0; i<array.length; i++)
			array[i] = Math.log(array[i]);
		return array;
	}

	private double[][] normalspace(double[][] mat) {
		// TODO Auto-generated method stub
		for(int i=0; i<mat.length; i++)
			mat[i] = normalspace(mat[i]);
		return mat;
	}

	private double[] normalspace(double[] array) {
		// TODO Auto-generated method stub
		for(int i=0; i<array.length; i++)
			array[i] = Math.exp(array[i]);
		return array;
	}

	private static double[][] transpose(double[][] mat) {
		// TODO Auto-generated method stub
		double[][] tMat = new double[mat[0].length][mat.length];
		for (int i = 0; i < mat.length; i++)
			for (int j = 0; j < mat[0].length; j++)
				tMat[j][i] = mat[i][j];
		return tMat;
	}

	private void clear(double[][][] matrix) {
		// TODO Auto-generated method stub
		int r = matrix.length, 
				c = matrix[0].length;
		for(int i=0; i<r; i++)
			for(int j=0; j<c; j++)
				Arrays.fill(matrix[i][j], 0);

	}

	private void clear(double[][] matrix) {
		// TODO Auto-generated method stub
		int r = matrix.length;
		for(int i=0; i<r; i++)
			Arrays.fill(matrix[i], 0);
	}

	private void clear(int[][] matrix) {
		// TODO Auto-generated method stub
		int r = matrix.length;
		for(int i=0; i<r; i++)
			Arrays.fill(matrix[i], 0);
	}

	private void EM() {
		// TODO Auto-generated method stub
		int _n_ = this.forward.length;
		int _m_ = this.dp.length;
		int[] s_a, s_b; 
		String[] s_l;
		ST[] st_b;
		double tmp, A;
		int _k_ = this.statespace.length;
		int _p_ = Constants._ploidy_H;
		int _g_, _d_, _t_;
		//long elapsed_emiss = 0, 
		//		elapsed_trans = 0, 
		//		tic, toc;

		for(int i=0; i<_m_; i++) {
			_g_ = this.obspace[i].genotype.length;
			OB ob = this.obspace[i];
			Map<String, Integer[]> dgMap = ob.dgMap;
			String[] dosaS = ob.dosaS;
			_d_ = dosaS.length;
			String[][] genotype = ob.genotype;
			String[] allele = ob.allele;
			Map<String, Integer> alMap = new HashMap<String, Integer>();
			for(int j=0; j<allele.length; j++) alMap.put(allele[j], j);
			
			//tic = System.nanoTime();
			if(i==0 || this.exp_b.contains(i)) {
				double[][] trans_count = this.transProbs[i].pseudo();
				TP tp1 = this.compoundTransProbs[i];
				for(int j=0;j<_n_; j++) {
					Integer[] vst = this.validStateSpace.get(j);
					
					DP dp1 = this.dp[i][j];
					FB fw1 = this.forward[j];
					FB bw1 = this.backward[j];
					double exp_c = fw1.logscale[i]+
							bw1.logscale[i]-
							fw1.probability;
					Integer[] vst2 = (i==0 ? 
							new Integer[]{0} : 
								this.validStateSpace.get(j));
					if(exp_c>Constants.MAX_EXP_DOUBLE) { 
						// cannot calculate exponential directly
						// logarithm and then exponential 
						// time consuming
						for(int a : vst2) {
							s_a = this.statespace[a].state;
							for(int b : vst) { 
								st_b = this.statespace[b].stateperm();
								_t_ = st_b.length;
								double[] prior = this.transProbs[i].prior[a][b];
								A = Math.exp(Math.log(
										fw1.probsMat[i][a]*
										tp1.probsMat[a][b]*
										dp1.emiss[b]*
										bw1.probsMat[i][b])+
										exp_c);
								for(int l=0; l<_t_; l++) {
									tmp = prior[l]*A;
									s_b = st_b[l].state;
									for(int m=0; m<_p_; m++) {
										trans_count[s_a[m]][s_b[m]] += tmp;
									}
								}
							}
						}
					} else { // exponential is safe
						// ideal way but dangerous!!!
						double exp = Math.exp(exp_c);
						for(int a : vst2) {
							s_a = this.statespace[a].state;
							for(int b : vst) { 
								st_b = this.statespace[b].stateperm();
								_t_ = st_b.length;
								double[] prior = this.transProbs[i].prior[a][b];
								A = fw1.probsMat[i][a]*
										tp1.probsMat[a][b]*
										dp1.emiss[b]*
										bw1.probsMat[i][b]*
										exp;
								for(int l=0; l<_t_; l++) {
									tmp = prior[l]*A;
									s_b = st_b[l].state;
									for(int m=0; m<_p_; m++) {
										trans_count[s_a[m]][s_b[m]] += tmp;
									}
								}
							}
						}
					}
				}
				this.transProbs[i].posterior();
			}
			//toc = System.nanoTime();
			//elapsed_trans += toc-tic;

			//tic = System.nanoTime();
			EP ep1 = this.compoundEmissProbs[i];
			double[][] emiss_count = this.emissProbs[i].pseudo();
			double[][] emissG = new double[_k_][_g_];
			for(int j=0;j<_n_; j++) {
				
				Integer[] vst = this.validStateSpace.get(j);
				DP dp1 = this.dp[i][j];
				FB fw1 = this.forward[j];
				FB bw1 = this.backward[j];
				double exp_c = fw1.logscale[i+1]+
						bw1.logscale[i]-
						fw1.probability;

				for(int k : vst) {
					Arrays.fill(emissG[k], 0);
					for(int a=0; a<_d_; a++) {
						Integer[] idx = dgMap.get(dosaS[a]);
						for(int c : idx)
							emissG[k][c] = dp1.likelihood[a]*
								ep1.probsMat[k][c]/dp1.emiss[k];
					}
				}
				
				if(exp_c>Constants.MAX_EXP_DOUBLE) {
					// cannot calculate exponential directly
					// logarithm and then exponential 
					// time consuming
					for(int a : vst) {
						if(a==0) continue;
						s_a = this.statespace[a].state;
						for(int c=0; c<_g_; c++) {
							A = Math.exp(Math.log(
									emissG[a][c]*
									fw1.probsMat[i+1][a]*
									bw1.probsMat[i][a])+
									exp_c);
							s_l = genotype[c];
							for(int l=0; l<_p_; l++) {
								emiss_count[s_a[l]]
										[alMap.get(s_l[l])] += A;
							}
						}
					}
				} else { // exponential is safe
					// ideal way but dangerous!!!
					double exp = Math.exp(exp_c);
					for(int a : vst) {
						if(a==0) continue;
						s_a = this.statespace[a].state;
						for(int c=0; c<_g_; c++) {
							A = emissG[a][c]*
									fw1.probsMat[i+1][a]*
									bw1.probsMat[i][a]*
									exp;
							s_l = genotype[c];
							for(int l=0; l<_p_; l++) {
								emiss_count[s_a[l]]
										[alMap.get(s_l[l])] += A;
							}
						}
					}
				}
			}
			this.emissProbs[i].posterior();

			//toc = System.nanoTime();
			//elapsed_emiss += toc-tic;
		}
		
		//System.out.println("elapsed emiss "+elapsed_emiss);
		//System.out.println("elapsed trans "+elapsed_trans);
	}

	private void makeBackward() {
		// TODO Auto-generated method stub
		int _m_ = this.dp.length;
		int _n_ = this.dp[0].length;
		int _k_ = this.statespace.length;
		double tmp;
		
		for(int i=0; i<_n_; i++) {
			//logger.info(i);
			double[][] probsMat = this.backward[i].probsMat;
			Arrays.fill(probsMat[_m_-1], 1.0);
			Integer[] vst = this.validStateSpace.get(i);
			
			for(int j=_m_-2; j>=0; j--) {	
				Arrays.fill(probsMat[j], 0);
				DP dp1 = this.dp[j+1][i];
				TP tp1 = this.compoundTransProbs[j+1];

				for(int k : vst) {
					tmp = 0;
					for(int a : vst) 
						tmp += dp1.emiss[a]
								*probsMat[j+1][a]*
								tp1.probsMat[k][a];
					probsMat[j][k] = tmp;
				}
				backward[i].scale(j);
			}
			DP dp1 = this.dp[0][i];
			int _d_ = dp1.likelihood.length;
			EP ep1 = this.compoundEmissProbs[0];
			TP tp1 = this.compoundTransProbs[0];
			double s = 0.0;
			for(int j=0; j<_k_; j++)
				for(int b=0; b<_d_; b++)
					s += dp1.likelihood[b]
							*ep1.probsDosaMat[j][b]*
							tp1.probsMat[0][j]*
							probsMat[0][j];
			backward[i].probability(s);
		}
		return;
	}

	private void makeForward() {
		// TODO Auto-generated method stub
		int _m_ = this.dp.length+1;
		int _n_ = this.dp[0].length;
		int _k_ = this.statespace.length;
		double tmp;
		
		for(int i=0; i<_n_; i++) {
			double[][] probsMat = this.forward[i].probsMat;
			probsMat[0][0] = 1.0;
			Arrays.fill(probsMat[0], 1, _k_, 0);
			Integer[] vst = this.validStateSpace.get(i);
			
			for(int j=1; j<_m_; j++) {
				Arrays.fill(probsMat[j], 0);
				TP tp1 = this.compoundTransProbs[j-1];
				for(int k : vst) {
					tmp = 0;
					for(int a : vst)
						tmp += probsMat[j-1][a]
								*tp1.probsMat[a][k];
					probsMat[j][k] = this.dp[j-1][i].emiss[k]*tmp;
				}

				this.forward[i].scale(j);
			}
			this.forward[i].probability(
					StatUtils.sum(probsMat[_m_-1]));
		}
		return;
	}

	private void checkFW() {
		// TODO Auto-generated method stub
		if(Constants.iteration==0) return;

		System.out.println(this.loglik()+"---"+this.loglik1());
		for(int i=0; i<this.forward.length; i++) {
			double r = Math.abs(this.forward[i].probability-
					this.backward[i].probability);
			if(r>1e-6)
				System.out.println(i+" | "+r+
						" --- FORWARD-BACKWARD PRECISION NOT RIGHT!!!");
		}
	}

	private double makeViterbi() {
		// TODO Auto-generated method stub
		int _m_ = this.dp.length+1;
		int _n_ = this.dp[0].length;
		int _k_ = this.statespace.length;

		double probability = 0.0;
		for(int i=0; i<_n_; i++) {
			Integer[] vst = this.validStateSpace.get(i);
			double[][] v = this.vb[i].v;
			int[][] trace = this.vb[i].trace;
			clear(v);
			clear(trace);
			v[0][0] = 1.0;
			
			double[] ss = new double[_k_];
			for(int j=1; j<_m_; j++) {
				DP dp1 = this.dp[j-1][i];
				TP tp1 = this.compoundTransProbs[j-1];

				for(int k : vst) {
					Arrays.fill(ss, 0);;
					for(int a : vst)
						ss[a] = v[j-1][a]*tp1.probsMat[a][k];
					final int ii = maxIndex(trace, j, k, ss);
					v[j][k] = dp1.emiss[k]*ss[ii];
				}
				this.vb[i].scale(j);
			}
			this.vb[i].trace();
			probability += this.vb[i].probability;
		}
		//System.out.println("############# "+probability);
		return probability;
	}

	private int maxIndex(int[][] trace, 
			int j, 
			int k, 
			double[] ss) {
		// TODO Auto-generated method stub
		int i = Algebra.maxIndex(ss);
		if(j>1 && ss[trace[j-2][k]]==ss[i])
			i = trace[j-2][k];
		trace[j-1][k] = i;
		return i;
	}

	private String[] makeHS() {
		// TODO Auto-generated method stub
		String[] hs = new String[Constants._haplotype_z+1];
		hs[0] = ".";
		for(int k=1; k<hs.length; k++)
			hs[k] = k<10 ? ""+k : ""+(char)('a'+k-10);
		return hs;
	}

	private void makeBWT() {
		// TODO Auto-generated method stub
		this.transProbs = new TP[this.M];
		transProbs[0] = new TP(this.statespace,
				this.hs, -1, true, false);
		for(int i=1; i<this.M; i++)
			transProbs[i] = exp_b.contains(i) ? 
					new TP(this.statespace, 
							this.hs, this.distance[i-1], 
							false, true) : 
								new TP(this.statespace, 
										this.hs, this.distance[i-1], 
										false, false);
					this.emissProbs = new EP[this.M];
					for(int i=0; i<this.M; i++)
						emissProbs[i] = new EP(
								this.hs, 
								this.obspace[i], 
								this.bfrac[i]);

					this.compoundTransProbs = new TP[this.M];
					int _n_ = this.statespace.length;
					for(int i=0; i<this.M; i++)
						this.compoundTransProbs[i] = new TP(
								str_statespace, 
								false, 
								new double[_n_][_n_]);
					this.makeCompoundTP();
					//logger.info("Free MEM: "+this.memory("used"));
					//this.makeCompoundAJTP();
					//logger.info("Free MEM: "+this.memory("used"));
					this.compoundEmissProbs = new EP[this.M];
					for(int i=0; i<this.M; i++) {
						OB ob = this.obspace[i];
						this.compoundEmissProbs[i] = new EP(
								str_statespace, 
								ob.genoS, 
								new double[_n_][ob.genotype.length], 
								new double[_n_][ob.dosage.length], 
								true);
					}
					this.makeCompoundEP();
					this.updateDP();
					//logger.info("Free MEM: "+this.memory("used"));
					//this.vb = new Viterbi(this);
					int _m_ = this.dp.length;
					int _k_ = this.dp[0].length;
					this.forward = new FB[_k_];
					for(int i=0; i<_k_; i++) 
						this.forward[i] = new FB(false, _m_+1, _n_);
					this.backward = new FB[_k_];
					for(int i=0; i<_k_; i++) 
						this.backward[i] = new FB(true, _m_+1, _n_);

					this.vb = new Viterbi[_k_];
					for(int i=0; i<_k_; i++)
						this.vb[i] = new Viterbi(_m_+1, _n_, this.str_statespace);
					return;
	}

	private void makeCompoundEP() {
		// TODO Auto-generated method stub
		int _n_ = this.statespace.length;
		for(int i=0; i<this.M; i++) 
			for(int j=0; j<_n_; j++) 
				this.compoundEP(i, j);
	}

	private void compoundEP(int i, int j) {
		// TODO Auto-generated method stub
		EP ep = this.emissProbs[i];
		OB ob = this.obspace[i];
		int[] state = this.statespace[j].state;
		String[][] geno = ob.genotype;
		int _m_ = geno.length;
		int _n_ = geno[0].length;
		Map<String, Integer> esMap = ep.esMap;
		for(int a=0; a<_m_; a++) {
			String[] g = geno[a];
			double ss=1.0;
			for(int b=0; b<_n_; b++) 
				ss *= ep.probsMat[state[b]]
						[esMap.get(g[b])];
			this.compoundEmissProbs[i].probsMat[j][a] = ss;
		}

		Map<String, Integer[]> dgMap = ob.dgMap;
		String[] dosa = ob.dosaS;
		int _k_ = dosa.length;
		int _i_;
		for(int k=0; k<_k_; k++) {
			Integer[] idx = dgMap.get(dosa[k]);
			_i_ = idx.length;
			double ss = 0.0;
			for(int l=0; l<_i_; l++)
				ss += this.compoundEmissProbs[i].probsMat[j][idx[l]];
			this.compoundEmissProbs[i].probsDosaMat[j][k] = ss;
		}
	}

	private void makeCompoundTP() {
		// TODO Auto-generated method stub
		int _n_ = this.statespace.length;
		for(int i=0; i<this.M; i++)
			for(int j=0; j<_n_; j++)
				for(int k=0; k<_n_; k++)
					this.compoundTP(j, k, i, true);
	}

	private void compoundTP(int j, int k, int i, boolean permutate) {
		// TODO Auto-generated method stub

		int[] states_0 = this.statespace[j].state;
		int[] states_1 = this.statespace[k].state;
		double[][] probs = this.transProbs[i].probsMat;

		int[][] perms;
		if(permutate) 
			//perms = Permutation.permutation(states_1);
			perms = this.spMap.get(StringUtils.
					join(this.statespace[k].state_str,
							"_"));
		else {
			perms = new int[1][states_1.length];
			perms[0] = states_1;
		}

		int _a_ = perms.length, _b_ = states_0.length;
		double ss = 0.0;
		for(int a=0; a<_a_; a++) {
			double s = 1.0;
			for(int b=0; b<_b_; b++)
				s *= probs[states_0[b]]
						[perms[a][b]];
			ss += s;
		}
		this.compoundTransProbs[i].probsMat[j][k] = ss;
	}

	private OB[] makeOB() {
		// TODO Auto-generated method stub
		List<String[]> allele = this.de.getAllele();
		OB[] obs = new OB[allele.size()];
		for(int i=0; i<obs.length; i++)
			obs[i] = new OB(allele.get(i));
		return obs;
	}

	private ST[] makeST() {
		// TODO Auto-generated method stub
		int mid1 = this.hs.length/2;
		int[] dot = new int[mid1];
		String[] dot_str = new String[mid1];
		Arrays.fill(dot_str, ".");
		Arrays.fill(dot, 0);
		int[] p1 = new int[mid1], p2 = new int[mid1];
		for(int i=0; i<mid1; i++) {
			p1[i] = i+1;
			p2[i] = mid1+i+1;
		}
		String[] p1_str = new String[mid1],
				p2_str = new String[mid1];
		System.arraycopy(this.hs, 1, p1_str, 0, mid1);
		System.arraycopy(this.hs, mid1+1, p2_str, 0, mid1);
		
		//List<List<String>> com1 = Permutation.multiPermutation(p1, mid1/2),
		//		com2 = Permutation.multiPermutation(p2, mid1/2);
		List<List<Integer>> com1 = Combination.combination(p1, mid1/2),
				com2 = Combination.combination(p2, mid1/2);
		List<Integer> com;
		int _n_ = com1.size();
		_n_ = _n_*_n_+3;
		ST[] st = new ST[_n_];
		int k=0;
		st[k++] = new ST(dot, dot_str, false, false);
		st[k++] = new ST(p1, p1_str, true, true);
		st[k++] = new ST(p2, p2_str, true, true);
		for(int i=0; i<com1.size(); i++) {
			for(int j=0; j<com2.size(); j++) {
				com = new ArrayList<Integer>(com1.get(i));
				com.addAll(com2.get(j));
				Integer[] tmp = new Integer[mid1];
				com.toArray(tmp);
				String[] tmp_str = new String[mid1];
				for(int m=0; m<mid1; m++)
					tmp_str[m] = this.hs[tmp[m]];
				st[k++] = new ST(ArrayUtils.toPrimitive(tmp),
						tmp_str,
						false,
						true);
			}
		}
		return st;
	}

	private DP[][] makeDP() {
		// TODO Auto-generated method stub
		DP[][] dp = new DP[this.M][this.N];
		switch(Constants._use_field_) {
		case "PL":
			List<List<double[]>> pl = this.de.getPhredScaledLikelihood();
			for(int i=0; i<pl.size(); i++)
				for(int j=0; j<pl.get(i).size(); j++)
					dp[i][j] = new DP(pl.get(i).get(j),
							this.statespace.length,
							this.obspace[i].genotype.length,
							isParent(this.sample[j]),
							false);
			break;
		case "AD":
			//TO-DO: need a model to determine the 
			// genotype likelihood from Allele Depth
			throw new RuntimeException(
					"check DataCollection line 406, 407!!!\n"
					+ "need to convert AD field!!!");
			//break;
		case "GT":
			//TO-DO: need a model to cast genotype
			// to likelihood
			throw new RuntimeException("need to convert GT field!!!");
			//break;
		default:
			logger.error("Invalid field provided. "
					+ "Should be \"PL\", \"AD\", or \"GT\". "
					+ "Program halted.");
			System.exit(1);
		}
		return dp;
	}

	private boolean isParent(String sample) {
		// TODO Auto-generated method stub
		return new HashSet<String>(Arrays.asList(this.parent)).contains(sample);
	}

	private long memory(String type) {
		switch(type) {
		case "total":
			return runtime.totalMemory();
		case "free":
			return runtime.freeMemory();
		case "used":
			return runtime.totalMemory()-
					runtime.freeMemory();
		case "max":
			return runtime.maxMemory();
		default:
			System.err.println("Error.");
			System.exit(1);
		}
		return -1;
	}

	private class DP {

		private final boolean isparent; // is parent
		private final double[] likelihood; // likelihood of different dosage
		private final double[] hweDist1;
		private boolean allowtrans; // if transition is allowed
		private boolean logspace; // if is in log space
		private double[] emiss; // product of likelihood and emiss for dosa
		//private double[][] emissG; // product of likelihood and emiss for genotype

		public DP(double[] likelihood, 
				int es,
				int esg,
				boolean isparent, 
				boolean logspace) {
			// TODO Auto-generated constructor stub
			this.isparent = isparent;
			this.allowtrans = !isparent;
			this.emiss = new double[es];
			//this.emissG = new double[es][esg];
			this.hweDist1 = new double[likelihood.length];
			Arrays.fill(hweDist1, 1.0/likelihood.length);
			if(logspace) likelihood = normalspace(likelihood);
			if(StatUtils.max(likelihood)<=0.0) 
				Arrays.fill(likelihood, 1.0/likelihood.length);
			likelihood = Algebra.normalize(likelihood);
			for(int i=0; i<likelihood.length; i++)
				likelihood[i] = hweDist1[i]*Constants._soften_ + 
				likelihood[i]*(1-Constants._soften_);
			this.likelihood = likelihood;
			if(logspace) this.setNormalspace();
		}

		private void setLogspace() {
			for(int i=0; i<likelihood.length; i++) 
				likelihood[i] = Math.log(likelihood[i]);
			this.logspace = true;
		}

		private void setNormalspace() {
			for(int i=0; i<likelihood.length; i++) 
				likelihood[i] = Math.exp(likelihood[i]);
			this.logspace = false;
		}
	}

	private class OB {
		private final String[] allele;
		private final String[][] genotype;
		private final String[][] dosage;
		private final String[] dosaS;
		private final String[] genoS;

		private Map<String, Integer[]> dgMap; // dosage and genotype map

		public OB(String[] allele) {
			// TODO Auto-generated constructor stub
			this.allele = allele;
			this.dosage = this.dosage();
			this.genotype = this.genotype();
			this.dosaS = this.makeDosaS();
			this.genoS = this.makeGenoS();
		}

		private String[] makeGenoS() {
			// TODO Auto-generated method stub
			String[] genoS = new String[this.genotype.length];	
			for(int i=0; i<genoS.length; i++) 
				genoS[i] = StringUtils.join(this.genotype[i],"_");

			return genoS;
		}

		private String[] makeDosaS() {
			// TODO Auto-generated method stub
			String[] dosaS = new String[this.dosage.length];
			for(int i=0; i<dosaS.length; i++) 
				dosaS[i] = StringUtils.join(this.dosage[i],"_");
			return dosaS;
		}

		private String[][] genotype() {
			// TODO Auto-generated method stub
			String[][] genotype = new String[(int)Math.pow(
					this.allele.length, 
					Constants._ploidy_H)][Constants._ploidy_H];
			this.dgMap = new HashMap<String, Integer[]>();
			int k = 0;
			for(int a=0; a<this.dosage.length; a++) {
				String[] dos = dosage[a];
				List<List<String>> perms = 
						Permutation.permutation(dos);
				Integer[] idx = new Integer[perms.size()];
				int i = 0;
				for(List<String> perm : perms) {
					perm.toArray(genotype[k]);
					idx[i++] = k++;
				}
				this.dgMap.put(StringUtils.join(dos, "_"), idx);
			}
			return genotype;
		}

		private String[][] dosage() {
			// TODO Auto-generated method stub
			List<List<String>> mc = Combination.multiCombination(this.allele,
					Constants._ploidy_H);
			String[][] dosage = new String[mc.size()][Constants._ploidy_H];
			for(int i=0; i<mc.size(); i++) 
				mc.get(i).toArray(dosage[i]);
			return dosage;
		}
	}

	private class ST {
		private final boolean isparent; // is valid for parents only
		private final int[] state;
		private final String[] state_str;
		private final boolean isdupstate;
		private final ST[] stateperm;

		public ST(int[] state,
				String[] state_str, 
				boolean isparent, 
				boolean isdupstate) {
			// TODO Auto-generated constructor stub
			this.isparent = isparent;
			this.state = state;
			this.state_str = state_str;
			this.isdupstate = isdupstate;
			this.stateperm = this.stateperm();
		}

		public ST(int[] state,
				String[] state_str, 
				boolean isparent) {
			// TODO Auto-generated constructor stub
			this.isparent = isparent;
			this.state_str = state_str;
			this.state = state;
			this.isdupstate = false;
			this.stateperm = null;
		}

		private ST[] stateperm() {
			// TODO Auto-generated method stub
			if(!isdupstate) return new ST[]{
					new ST(this.state,
							this.state_str,
							this.isparent)};
			if(isparent) {
				int n = this.state.length;
				List<List<Integer>> perm = 
						Permutation.permutation(this.state);
				ST[] stateperm = new ST[perm.size()];
				for(int i=0; i<stateperm.length; i++) {
					Integer[] st = new Integer[n];
					perm.get(i).toArray(st);
					String[] ss = new String[n];
					for(int j=0; j<n; j++)
						ss[j] = hs[st[j]];
					stateperm[i] = new ST(
							ArrayUtils.toPrimitive(st),
							ss,
							this.isparent,
							false);
				}
				return stateperm;
			} else {
				int mid1 = Constants._ploidy_H/2;
				int n = this.state.length;
				int[] p = new int[mid1];
				System.arraycopy(this.state, 0, p, 0, mid1);
				List<List<Integer>> perm1 = Permutation.permutation(p);
				System.arraycopy(this.state, mid1, p, 0, mid1);
				List<List<Integer>> perm2 = Permutation.permutation(p);
				ST[] stateperm = new ST[perm1.size()*perm2.size()];
				int i = 0;
				for(List<Integer> p1 : perm1) {
					for(List<Integer> p2 : perm2) {
						List<Integer> s = new ArrayList<Integer>(p1);
						s.addAll(p2);
						Integer[] st = new Integer[n];
						s.toArray(st);
						String[] ss = new String[n];
						for(int j=0; j<n; j++)
							ss[j] = hs[st[j]];
						stateperm[i++] = new ST(
								ArrayUtils.toPrimitive(st),
								ss,
								this.isparent,
								false);
					}
				}
				return stateperm;
			}
		}
	}

	private class TP {
		private final boolean isDotState;
		private final double distance;
		private final ST[] statespace;
		private final String[] str_statespace;
		private final boolean trainExp;

		private double[][][] prior;
		private boolean logspace;
		private double[][] probsMat;
		private double[][] count;
		private double exp;
		private double[][] alpha;

		public TP(String[] str_statespace, 
				boolean logspace, 
				double[][] probsMat) {
			this.statespace = null;
			this.str_statespace = str_statespace;
			this.isDotState = false;
			this.distance = -1;
			this.logspace = logspace;
			this.probsMat = probsMat;
			if( logspace ) this.setNormalspace();
			this.trainExp = false;
		}

		public TP(ST[] statespace,
				String[] str_statespace, 
				double distance, 
				boolean isDotState,
				boolean trainExp) {
			this.statespace = statespace;
			this.str_statespace = str_statespace;
			this.distance = distance;
			this.isDotState = isDotState;
			this.prior();
			this.trainExp = trainExp;
			int _n_=this.str_statespace.length;
			this.probsMat = new double[_n_][_n_];
			this.count = new double[_n_][_n_];
			this.logspace = false;
			int _m_=this.statespace.length;
			this.prior = new double[_m_][_m_][];
			this.update();
		}

		private void update() {
			// TODO Auto-generated method stub
			int z = Constants._haplotype_z/2;
			if(this.isDotState) {
				for(int i=0; i<2; i++) {
					for(int j=1; j<=z; j++)
						this.probsMat[0][i*z+j] = alpha[i][j-1]/2;
				}
			} else {
				for(int i=0; i<2; i++) {
					int start = 1, end = z;
					if(i==1) {
						start += z;
						end += z;
					}
					for(int j=start; j<=end; j++) {
						for(int k=start; k<=end; k++)
							if(j==k) this.probsMat[j][k] = (1-exp)+exp*alpha[i][k-start];
							else this.probsMat[j][k] = exp*alpha[i][k-start];
					}
				}
			}

			int _a_ = this.statespace.length;
			for(int a=0; a<_a_; a++) {
				int[] s_a = this.statespace[a].state;
				for(int b=0; b<_a_; b++) { 
					ST[] st_b = this.statespace[b].stateperm();
					prior[a][b] = prior(s_a, st_b);
				}
			}
		}

		private double[] prior(int[] s, ST[] st) {
			// TODO Auto-generated method stub
			// TODO Auto-generated method stub
			int _i_ = st.length;
			int _j_ = st[0].state.length;
			double[] prior = new double[_i_];
			Arrays.fill(prior, 1.0);
			for(int i=0; i<_i_; i++) {
				int[] s_i = st[i].state;
				for(int j=0; j<_j_; j++)
					prior[i] *= probsMat[s[j]][s_i[j]];
			}
			return Algebra.normalize(prior);
		}

		private void prior() {
			// TODO Auto-generated method stub	
			double e = Math.exp(-this.distance*Constants._con_base_r);
			double alpha = Math.max(Constants._mu_J_e*(1-e), Constants.eps);
			double beta = Math.max(Constants._mu_J_e*e, Constants.eps);
			this.exp = new BetaDistribution(Constants.rg, alpha, beta).sample();
			int z = Constants._haplotype_z/2;
			this.alpha = new double[2][z];
			//Dirichlet diri = new Dirichlet(distz(z), Constants._mu_alpha_e);
			for(int i=0; i<2; i++)
				//this.alpha[i] = ArrayUtils.toPrimitive(diri.sample());
				Arrays.fill(this.alpha[i], 1.0/z);
			return;
		}

		private void posterior() {
			// TODO Auto-generated method stub				
			int z = Constants._haplotype_z/2;
			if(this.isDotState) {
				double[] c = new double[z];
				for(int i=0; i<2; i++) {
					System.arraycopy(count[0], i*z+1, c, 0, z);
					alpha(i, z, c);
				}
				this.update();
			} else if(this.trainExp) {
				double s = 0, exp_c = 0;
				for(int i=0; i<2; i++) {
					int start = 1, end = z;
					if(i==1) {
						start += z;
						end += z;
					}
					double[] alpha_c = new double[z];
					for(int j=start; j<=end; j++) {
						for(int k=start; k<=end; k++) {
							s += count[j][k];
							if(k==j) {
								exp_c += (1-exp)/probsMat[j][k]*count[j][k];
								alpha_c[k-i*z-1] += exp*alpha[i][k-start]/
										probsMat[j][k]*count[j][k];
							} else {
								alpha_c[k-i*z-1] += count[j][k];
							}
						}
					}
					alpha(i, z, alpha_c);
				}
				exp = 1-exp_c/s;
				this.update();
			}
		}

		public double[][] pseudo() {
			// TODO Auto-generated method stub
			
			//TODO
			//for pair haplotype phasing
			//need to update distance!!!
			
			double e = Math.exp(-this.distance*Constants._con_base_r);
			double alpha_m = Math.max(Constants._pseudo_[2]*(1-e), Constants.eps);
			double beta_m = Math.max(Constants._pseudo_[2]*e, Constants.eps);
			double ab = alpha_m/(alpha_m+beta_m);

			int z = Constants._haplotype_z/2;
			if(this.isDotState) {
				for(int i=0; i<2; i++) {
					for(int j=1; j<=z; j++)
						this.count[0][i*z+j] = ab/z;
				}
			} else {
				for(int i=0; i<2; i++) {
					int start = 1, end = z;
					if(i==1) {
						start += z;
						end += z;
					}
					for(int j=start; j<=end; j++) {
						for(int k=start; k<=end; k++)
							if(j==k) this.count[j][k] = (1-ab)+ab/z;
							else this.count[j][k] = ab/z;
					}
				}
			}

			return count;
		}

		private void alpha(final int i, 
				final int z, 
				final double[] c) {
			// TODO Auto-generated method stub
			for(int j=0; j<c.length; j++)
				c[j] += Constants._pseudo_[0]/z;
			double s = StatUtils.sum(c);
			for(int j=0; j<alpha[i].length; j++)
				alpha[i][j] = c[j]/s;
		}

		private void setLogspace() {
			// TODO Auto-generated method stub
			for(int i=0; i<this.probsMat.length; i++)
				for(int j=0; j<this.probsMat[i].length; j++)
					this.probsMat[i][j] = Math.log(this.probsMat[i][j]);
			this.logspace = true;
		}

		private void setNormalspace() {
			// TODO Auto-generated method stub
			for(int i=0; i<this.probsMat.length; i++)
				for(int j=0; j<this.probsMat[i].length; j++)
					this.probsMat[i][j] = Math.exp(this.probsMat[i][j]);
			this.logspace = false;
		}

		private void setProbsMat(double[][] probsMat) {
			this.probsMat = probsMat;
		}
	}

	private class EP {

		private final String[] statespace;
		private final String[] allele;
		private final double bfrac;
		private double[][] probsMat;
		private boolean logspace;
		private double[][] probsDosaMat;
		private double[][] count;

		private final Map<String, Integer> ssMap;
		private final Map<String, Integer> esMap;

		public EP(String[] statespace, String[] allele, 
				double[][] probsMat,
				double[][] probsDosaMat,
				boolean logspace) {
			// TODO Auto-generated constructor stub
			this.statespace = statespace;
			this.allele = allele;
			this.probsMat = probsMat;
			this.probsDosaMat = probsDosaMat;
			this.logspace = logspace;
			this.bfrac = 0;
			if( logspace ) this.setNormalspace();
			this.ssMap = new HashMap<String, Integer>();
			for(int i=0; i<statespace.length; i++)
				this.ssMap.put(statespace[i], i);
			this.esMap = new HashMap<String, Integer>();
			for(int i=0; i<allele.length; i++)
				this.esMap.put(allele[i], i);
		}

		public EP(String[] statespace, OB obspace, double bfrac) {
			// TODO Auto-generated constructor stub
			this.statespace = statespace;
			this.allele = obspace.allele;
			this.bfrac = bfrac;
			this.prior();
			this.ssMap = new HashMap<String, Integer>();
			for(int i=0; i<statespace.length; i++)
				this.ssMap.put(statespace[i], i);
			this.esMap = new HashMap<String, Integer>();
			for(int i=0; i<allele.length; i++)
				this.esMap.put(allele[i], i);
			this.probsDosaMat = null;
			this.logspace = false;
			this.count = new double[this.statespace.length][this.allele.length];
		}

		private void prior() {
			// TODO Auto-generated method stub
			this.probsMat = new double[this.statespace.length][this.allele.length];
			Arrays.fill(this.probsMat[0], 0);
			for(int i=1; i<this.probsMat.length; i++) {
				//Dirichlet diri = new Dirichlet(distz(this.allele.length), 
				//		Constants._mu_theta_e);
				Dirichlet diri = new Dirichlet(new double[]{1-bfrac, bfrac}, 
						Constants._mu_theta_e);
				this.probsMat[i] = ArrayUtils.toPrimitive(diri.sample());
			}
		}

		private void posterior() {
			// TODO Auto-generated method stub
			int _i_ = count.length,
					_j_ = count[0].length;
			for(int i=1; i<_i_; i++) {
				double s = StatUtils.sum(count[i]);
				if(s==0)
					for(int j=0; j<_j_; j++)
						probsMat[i][j] = 1.0/count[i].length;
				else
					for(int j=0; j<_j_; j++)
						probsMat[i][j] = count[i][j]/s;
			}
		}

		public double[][] pseudo() {
			// TODO Auto-generated method stub
			int _i_ = count.length;
			for(int i=1; i<_i_; i++) {
				//Arrays.fill(count[i], 
				//		1.0/this.allele.length*Constants._mu_theta_m);
				count[i][0] = (1-bfrac)*Constants._pseudo_[1];
				count[i][1] = bfrac*Constants._pseudo_[1];
			}
			return count;
		}

		private void setLogspace() {
			// TODO Auto-generated method stub
			for(int i=0; i<this.probsMat.length; i++)
				for(int j=0; j<this.probsMat[i].length; j++)
					this.probsMat[i][j] = Math.log(this.probsMat[i][j]);
			this.logspace = true;
		}

		private void setNormalspace() {
			// TODO Auto-generated method stub
			for(int i=0; i<this.probsMat.length; i++)
				for(int j=0; j<this.probsMat[i].length; j++)
					this.probsMat[i][j] = Math.exp(this.probsMat[i][j]);
			this.logspace = false;
		}

		private void setProbsMat(double[][] probsMat) {
			this.probsMat = probsMat;
		}
	}

	public DP[][] getDP() {
		return this.dp;
	}

	private static double[] distz() {
		double[] d = new double[Constants._haplotype_z];
		for(int i=0; i<Constants._haplotype_z; i++) 
			d[i] = 1.0/Constants._haplotype_z;
		return d;
	}

	private static double[] distz(int z) {
		double[] d = new double[z];
		for(int i=0; i<z; i++) d[i] = 1.0/z;
		return d;
	}


	private class FB { /** forward/backward algorithm object */
		private double probability;
		private double[][] probsMat;
		private boolean logspace;
		private double[] logscale;
		private final boolean backward;

		public FB(boolean backward,
				int m,
				int s) {
			this.backward = backward;
			this.probsMat = new double[m][s];
			this.probability = 0;
			this.logscale = new double[m];
			Arrays.fill(this.logscale, 0.0);
			this.logspace = false;
		}

		public void probability(double p) {
			// TODO Auto-generated method stub
			if(this.backward)
				this.probability = Math.log(p)+this.logscale[0];
			else
				this.probability = Math.log(p)+
				this.logscale[this.logscale.length-1];
		}

		private void scale() {
			// TODO Auto-generated method stub
			this.logscale = new double[this.probsMat.length];
			if(this.backward)
				for(int i=this.logscale.length-1; i>=0; i++)
					this.scale(i);
			else
				for(int i=0; i<this.logscale.length; i++)
					this.scale(i);
		}

		private void scale(final int i) {
			// TODO Auto-generated method stub
			if(i==this.logscale.length-1 && this.backward) 
				return;
			if(i==0 && !this.backward) return;
			int dv = -1;
			if(this.backward) dv = 1;
			double[] probs = this.probsMat[i];
			double min = Double.POSITIVE_INFINITY,
					max = Double.NEGATIVE_INFINITY;
			for(int k=0; k<probs.length; k++) {
				if(probs[k]>0) {
					min = probs[k]<min ? probs[k] : min;
					max = probs[k]>max ? probs[k] : max;
				}
			}

			this.logscale[i] = this.logscale[i+dv];
			if(min<Constants.threshMin &&
					max<Constants.threshMax) {
				this.logscale[i] += 
						Constants.logThreshMax;
				for(int k=0; k<probs.length; k++)
					probs[k] /= Constants.threshMax;
			}
		}

		private void setLogspace() {
			// TODO Auto-generated method stub
			for(int i=0; i<this.probsMat.length; i++)
				for(int j=0; j<this.probsMat[i].length; j++)
					this.probsMat[i][j] = Math.log(this.probsMat[i][j]);
			this.probability = Math.log(this.probability);
			this.logspace = true;
		}

		private void setNormalspace() {
			// TODO Auto-generated method stub
			for(int i=0; i<this.probsMat.length; i++)
				for(int j=0; j<this.probsMat[i].length; j++)
					this.probsMat[i][j] = Math.exp(this.probsMat[i][j]);
			this.probability = Math.exp(this.probability);
			this.logspace = false;
		}
	}

	private class Viterbi {
		private double[][] v;
		private int m;
		private int[][] trace;
		private String[] statespace;
		private double[] logscale;
		private String[] path;
		private double probability;

		public Viterbi(int _m_, int _n_, String[] statespace) {
			// TODO Auto-generated constructor stub
			this.v = new double[_m_][_n_];
			this.m = _m_;
			this.trace = new int[_m_-1][_n_];
			this.statespace = statespace;
			this.logscale = new double[_m_];
			this.path = new String[_m_-1];
			this.probability = 0.0;
		}

		private void trace() {
			this.probability = Math.log(StatUtils.max(v[m-1]))+
					this.logscale[m-1];
			int tr = Algebra.maxIndex(v[m-1]);
			this.path[m-2] = this.statespace[tr];
			for(int i=m-3; i>=0; i--) {
				tr = trace[i+1][tr];
				this.path[i] = this.statespace[tr];
			}
			return;
		}

		private void scale(final int i) {
			// TODO Auto-generated method stub
			if(i==0) return;
			double[] probs = this.v[i];
			double min = Double.POSITIVE_INFINITY,
					max = Double.NEGATIVE_INFINITY;
			for(int k=0; k<probs.length; k++) {
				if(probs[k]>0) {
					min = probs[k]<min ? probs[k] : min;
					max = probs[k]>max ? probs[k] : max;
				}
			}

			this.logscale[i] = this.logscale[i-1];
			if(min<Constants.threshMin &&
					max<Constants.threshMax) {
				this.logscale[i] += 
						Constants.logThreshMax;
				for(int k=0; k<probs.length; k++)
					probs[k] /= Constants.threshMax;
			}
		}
	}

	public double loglik() {
		if(Constants.iteration==0)
			return Double.NEGATIVE_INFINITY;
		else {
			double probability = 0;
			for(FB fw : this.forward) probability += fw.probability;
			return probability;
		}
	}

	public double loglik1() {
		if(Constants.iteration==0)
			return Double.NEGATIVE_INFINITY;
		else {
			double probability = 0;
			for(FB bw : this.backward) probability += bw.probability;
			return probability;
		}
	}

	public void print(boolean detail) {
		// TODO Auto-generated method stub
		double probability = 0;
		//for(Viterbi vb: this.vb) probability += vb.probability;
		//logger.info(" log-likelihood "+probability);
		if(Constants.iteration==0)
			logger.info(" hmm initialised.");
		else {
			for(FB fw : this.forward) probability += fw.probability;
			logger.info(" log-likelihood "+probability);
		}

		if(detail) {
			this.makeViterbi();
			logger.info("Viterbi path:");
			for(int i=0; i<this.sample.length; i++) {
				String[] path = this.vb[i].path;
				List<String[]> path_s = new ArrayList<String[]>();
				for(int j=0; j<path.length; j++)
					path_s.add(path[j].split("_"));
				for(int j=0; j<Constants._ploidy_H; j++) {
					System.out.print("# id "+this.sample[i]+":"+(j+1)+"\t\t\t");
					for(int k=0; k<path_s.size(); k++)
						System.out.print(path_s.get(k)[j]);
					System.out.println();
				}
			}
			logger.info("Emission probability:");
			for(int i=0; i<this.emissProbs.length; i++) {
				double[][] emissProbs = this.emissProbs[i].probsMat;
				String[] allele = this.emissProbs[i].allele;
				System.out.print("Marker "+(i+1)+"|\t\t");
				for(int j=0; j<emissProbs.length; j++) {
					System.out.print(this.hs[j]+"-> {");
					for(int k=0; k<emissProbs[j].length; k++) 
						System.out.print(allele[k]+","+emissProbs[j][k]+";");
					System.out.print("} ");
				}
				System.out.println();
			}
			logger.info("Transition probability:");
			for(int i=0; i<this.transProbs.length; i++) {
				double[][] transProbs = this.transProbs[i].probsMat;
				System.out.print("Marker "+(i+1)+"|\t\t");
				for(int j=0; j<transProbs.length; j++) {
					System.out.print(this.hs[j]+"-> {");
					for(int k=0; k<transProbs[j].length; k++) 
						System.out.print(this.hs[k]+","+transProbs[j][k]+";");
					System.out.print("} ");
				}
				System.out.println();
			}
		}
	}

	public void print() {
		this.print(false);
	}

	public void write(String output, 
			String experiment, 
			String contig) {
		// TODO Auto-generated method stub
		this.makeViterbi();
		String home = output+
				Constants.file_sep+
				experiment+"."+
				contig+"."+
				Constants.seed+"_1_1_all_0_200mb";
		File fhome = new File(home);
		if(fhome.exists() && fhome.isDirectory()) {
			try {
				FileUtils.deleteDirectory(fhome);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				System.exit(1);
			}
		}
		fhome.mkdir();
		String phasedStates = home+Constants.file_sep+"phasedStates";
		new File(phasedStates).mkdir();
		String results_hmm = home+Constants.file_sep+"results_hmm";
		new File(results_hmm).mkdir();

		try {
			BufferedWriter bw = IO.getBufferedWriter(phasedStates+
					Constants.file_sep+
					experiment+".txt");
			bw.write(""+this.loglik()+"\n");
			bw.write(""+this.dp.length+"\n");
			for(int i=0; i<this.sample.length; i++) {
				String[] path = this.vb[i].path;
				List<String[]> path_s = new ArrayList<String[]>();
				for(int j=0; j<path.length; j++)
					path_s.add(path[j].split("_"));
				for(int j=0; j<Constants._ploidy_H; j++) {
					bw.write("# id "+this.sample[i]+":"+(j+1)+"\t\t\t");
					for(int k=0; k<path_s.size(); k++)
						bw.write(path_s.get(k)[j]);
					bw.write("\n");;
				}
			}
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.exit(1);
		}

		try {
			BufferedWriter bw = IO.getBufferedWriter(results_hmm+
					Constants.file_sep+
					"emissionModel.txt");
			for(int i=0; i<this.emissProbs.length; i++) {
				double[][] emissProbs = this.emissProbs[i].probsMat;
				String[] allele = this.emissProbs[i].allele;
				bw.write(this.de.getId()+"_"+this.de.getPosition()[i]+"\t\t\t");
				for(int j=0; j<emissProbs.length; j++) {
					bw.write(this.hs[j]+"-> {");
					for(int k=0; k<emissProbs[j].length; k++) 
						bw.write(allele[k]+","+emissProbs[j][k]+";");
					bw.write("} ");
				}
				bw.write("\n");
			}
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.exit(1);
		}
		
		try {
			BufferedWriter bw = IO.getBufferedWriter(results_hmm+
					Constants.file_sep+
					"transitionModel.txt");
			for(int i=0; i<this.transProbs.length; i++) {
				double[][] transProbs = this.transProbs[i].probsMat;
				bw.write(this.de.getId()+"_"+this.de.getPosition()[i]+"\t\t\t");
				for(int j=0; j<transProbs.length; j++) {
					bw.write(this.hs[j]+"-> {");
					for(int k=0; k<transProbs[j].length; k++) 
						bw.write(this.hs[k]+","+transProbs[j][k]+";");
					bw.write("} ");
				}
				bw.write("\n");
			}
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.exit(1);
		}
		

		try {
			BufferedWriter bw = IO.getBufferedWriter(fhome+
					Constants.file_sep+
					"stderr_true");
			bw.write("cz1.model.HiddenMarkovModel:\n");
			bw.write("cz1.model.HiidenMarkovModel$EM:\n");
			bw.write("log prob is "+this.loglik()+" at "+Constants.iteration);
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.exit(1);
		}


		try {
			BufferedWriter bw = IO.getBufferedWriter(fhome+
					Constants.file_sep+
					"snp_"+experiment+".txt");
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
				bw.write(os.toString());
			}
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.exit(1);
		}
	}
}
