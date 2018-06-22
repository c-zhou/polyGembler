package cz1.hmm.model;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.log4j.Logger;

import cz1.breeding.data.FullSiblings;
import cz1.hmm.data.DataEntry;
import cz1.hmm.model.HiddenMarkovModel.DP;
import cz1.util.Algebra;
import cz1.util.Combination;
import cz1.util.Constants;
import cz1.util.Utils;
import cz1.util.Permutation;
import cz1.util.Dirichlet;
import cz1.util.Constants.Field;

public class HiddenMarkovModelBWT extends HiddenMarkovModel {
	private FB[] forward, backward;
	
	// ----founder haplotypes coefficients----
	private double founder_hap_coeff = 0.5;
	// higher quality data get more weights for training the model
	private double[][] weights;
	
	public HiddenMarkovModelBWT(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse,
			boolean trainExp,
			Field field,
			double fh_coeff) {
		super(de, seperation, reverse, trainExp, field);
		this.founder_hap_coeff = fh_coeff;
		this.makeWeightingCoeff();
		
		//this.trainAllExp();
		this.makeBWT();
		//this.makeViterbi();
		//this.print();
		//this.train();
	}
	
	public HiddenMarkovModelBWT(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse,
			boolean trainExp,
			Field field) {
		super(de, seperation, reverse, trainExp, field);
		this.makeWeightingCoeff();
		//this.trainAllExp();
		this.makeBWT();
		//this.makeViterbi();
		//this.print();
		//this.train();
	}
	
	public HiddenMarkovModelBWT(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse,
			boolean trainExp,
			Field field,
			double fh_coeff,
			String hmmFile) {
		super(de, seperation, reverse, trainExp, field);
		this.founder_hap_coeff = fh_coeff;
		this.makeWeightingCoeff();
		this.makeBWT(hmmFile);
	}
	
	public HiddenMarkovModelBWT(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse,
			boolean trainExp,
			Field field,
			String hmmFile) {
		super(de, seperation, reverse, trainExp, field);
		this.makeWeightingCoeff();
		this.makeBWT(hmmFile);
	}
	
	public void train() {
		iteration++;
		myLogger.info("###################");
		myLogger.info("train: "+iteration);
		long[] tic = new long[10];
		int k=0;
		tic[k++] = System.nanoTime();
		this.makeForward();
		tic[k++] = System.nanoTime();
		myLogger.info("forward done "+(tic[k-1]-tic[k-2])+"ns");
		this.makeBackward();
		tic[k++] = System.nanoTime();
		myLogger.info("backward done "+(tic[k-1]-tic[k-2])+"ns");
		this.checkFW();
		tic[k++] = System.nanoTime();
		this.EM();
		tic[k++] = System.nanoTime();
		myLogger.info("EM algorithm "+(tic[k-1]-tic[k-2])+"ns");
		this.makeCompoundEP();
		tic[k++] = System.nanoTime();
		myLogger.info("make compound ep "+(tic[k-1]-tic[k-2])+"ns");
		this.makeCompoundTP();
		tic[k++] = System.nanoTime();
		myLogger.info("make compound tp "+(tic[k-1]-tic[k-2])+"ns");
		//this.makeViterbi();
		this.updateDP();
		tic[k++] = System.nanoTime();
		myLogger.info("update dp "+(tic[k-1]-tic[k-2])+"ns");
		this.printExpTrain();
		return;
	}

	private void makeWeightingCoeff() {
		// TODO Auto-generated method stub
		this.weights = new double[M][N];
		
		FullSiblings fullSib = new FullSiblings();
		double[][] pvals = fullSib.calcDosaConfig(this.de, this.parent_i);
		for(int i=0; i<M; i++) {
			if(StatUtils.max(pvals[i])==0) {
				Arrays.fill(weights[i], 1.0);
				continue;
			}
			
			final int p_i = Algebra.maxIndex(pvals[i]);
			double[] weights_i = weights[i];
			
			// parental samples
			if(!miss[i][this.parent_i[0]] && !miss[i][this.parent_i[1]]) {
				int[] parentalDosa = fullSib.getParentalDosaByIndex(p_i);
				double[] ll0 = this.dp[i][this.parent_i[0]].likelihood,
						ll1 = this.dp[i][this.parent_i[1]].likelihood;
				if( ll0[parentalDosa[0]]+ll1[parentalDosa[1]] < 
						ll0[parentalDosa[1]]+ll1[parentalDosa[0]]) {
					int tmp_i = parentalDosa[0];
					parentalDosa[0] = parentalDosa[1];
					parentalDosa[1] = tmp_i;
				}
				weights_i[parent_i[0]] = 1.0+founder_hap_coeff*ll0[parentalDosa[0]]*(N-2);
				weights_i[parent_i[1]] = 1.0+founder_hap_coeff*ll1[parentalDosa[1]]*(N-2);
			} else {
				weights_i[parent_i[0]] = miss[i][this.parent_i[0]]?0.1:1.0;
				weights_i[parent_i[1]] = miss[i][this.parent_i[1]]?0.1:1.0;
			}
			
			boolean[] f1Dosa = fullSib.getF1DosaByIndex(p_i);
			// f1 samples
			for(int j=0; j<N; j++) {
				if(this.pedigree[j]!=-1) continue;
				int dosa = this.getDosage(i, j);
				if(miss[i][j]) weights_i[j] = 0.1;
				weights_i[j] = f1Dosa[dosa] ? 1.0 : 0.01;
			}
			
			Algebra.normalize(weights_i);
			Algebra.multiply(weights_i, N);
		}
		Utils.print(this.weights);
	}

	private boolean isMissing(int i, int j) {
		// TODO Auto-generated method stub
		return false;
	}

	private int getDosage(final int i, final int j) {
		// TODO Auto-generated method stub
		double[] ll = this.dp[i][j].likelihood;
		return Algebra.maxIndex(ll);
	}

	protected void EM() {
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
				
				// ----founder haplotypes coefficients----
				double coeff = weights[i][j];
				
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
							/***
							A = coeff*Math.exp(Math.log(
									emissG[a][c]*
									fw1.probsMat[i+1][a]*
									bw1.probsMat[i][a])+
									exp_c);
							**/
							// ----founder haplotypes coefficients----
							A = coeff*Math.exp(Math.log(
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
							/***
							A = emissG[a][c]*
									fw1.probsMat[i+1][a]*
									bw1.probsMat[i][a]*
									exp;
							 **/
							// ----founder haplotypes coefficients----
							A = coeff*emissG[a][c]*
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
		if(iteration==0) return;

		System.out.println(this.loglik()+"---"+this.loglik1());
		for(int i=0; i<this.forward.length; i++) {
			double r = Math.abs(this.forward[i].probability-
					this.backward[i].probability);
			if(r>1e-6)
				System.out.println(i+" | "+r+
						" --- FORWARD-BACKWARD PRECISION NOT RIGHT!!!");
		}
	}

	private void makeBWT(String hmmFile) {
		// TODO Auto-generated method stub
		
		// load a hmm object from a .zip hmmFile
		try {
			final ZipFile in = new ZipFile(hmmFile);

			this.transProbs = new TP[this.M];
			ZipEntry transEntry = in.getEntry("results_hmm/transitionModel.txt");
			if(transEntry!=null) {
				final BufferedReader br_trans = Utils.getBufferedReader(
						in.getInputStream(transEntry));
				transProbs[0] = new TP(this.statespace,
						this.hs, -1, true, false, br_trans.readLine());
				for(int i=1; i<this.M; i++) {
					transProbs[i] = exp_b.contains(i) ? 
							new TP(this.statespace, 
									this.hs, this.distance[i-1], 
									false, true, br_trans.readLine()) : 
										new TP(this.statespace, 
												this.hs, this.distance[i-1], 
												false, false, br_trans.readLine());	
				} 
				br_trans.close();
			} else {
				transProbs[0] = new TP(this.statespace,
						this.hs, -1, true, false);
				for(int i=1; i<this.M; i++) {
					transProbs[i] = exp_b.contains(i) ? 
							new TP(this.statespace, 
									this.hs, this.distance[i-1], 
									false, true) : 
										new TP(this.statespace, 
												this.hs, this.distance[i-1], 
												false, false);
				}
			}

			this.emissProbs = new EP[this.M];
			ZipEntry emissEntry = in.getEntry("results_hmm/emissionModel.txt");
			if(emissEntry!=null) {
				final BufferedReader br_emiss = Utils.getBufferedReader(
						in.getInputStream(emissEntry));
				for(int i=0; i<this.M; i++) {
					emissProbs[i] = new EP(
							this.hs, 
							this.obspace[i], 
							this.bfrac[i],
							br_emiss.readLine());
				}
				br_emiss.close();
			} else {
				for(int i=0; i<this.M; i++) {
					emissProbs[i] = new EP(
							this.hs, 
							this.obspace[i], 
							this.bfrac[i]);
				}
			}
			in.close();

			this.compoundTransProbs = new TP[this.M];
			int _n_ = this.statespace.length;
			for(int i=0; i<this.M; i++) {
				this.compoundTransProbs[i] = new TP(
						str_statespace, 
						false, 
						new double[_n_][_n_]);
			}
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
		} catch (IOException e) {
			e.printStackTrace();
		}
		return;
	}
	
	private void makeBWT() {
		// TODO Auto-generated method stub
		this.transProbs = new TP[this.M];
		transProbs[0] = new TP(this.statespace,
				this.hs, -1, true, false);
		for(int i=1; i<this.M; i++) {
			transProbs[i] = exp_b.contains(i) ? 
					new TP(this.statespace, 
							this.hs, this.distance[i-1], 
							false, true) : 
								new TP(this.statespace, 
										this.hs, this.distance[i-1], 
										false, false);
		}

		this.emissProbs = new EP[this.M];
		for(int i=0; i<this.M; i++) {
			emissProbs[i] = new EP(
					this.hs, 
					this.obspace[i], 
					this.bfrac[i]);
		}

		this.compoundTransProbs = new TP[this.M];
		int _n_ = this.statespace.length;
		for(int i=0; i<this.M; i++) {
			this.compoundTransProbs[i] = new TP(
					str_statespace, 
					false, 
					new double[_n_][_n_]);
		}
		
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

	public double loglik() {
		if(iteration==0)
			return Double.NEGATIVE_INFINITY;
		else {
			double probability = 0;
			for(FB fw : this.forward) probability += fw.probability;
			return probability;
		}
	}

	public double loglik1() {
		if(iteration==0)
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
		if(iteration==0)
			myLogger.info(" hmm initialised.");
		else {
			for(FB fw : this.forward) probability += fw.probability;
			myLogger.info(" log-likelihood "+probability);
		}

		if(detail) {
			this.makeViterbi();
			myLogger.info("Viterbi path:");
			for(int i=0; i<this.sample.length; i++) {
				String[] path = this.vb[i].path_str;
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
			myLogger.info("Emission probability:");
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
			myLogger.info("Transition probability:");
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
	
	public int[][][] resample(final int n) {
		int[][][] path = new int[this.N][this.M][n];
		for(int i=0; i!=N; i++) path[i] = resample(n, i);
		return path;
	}
	
	public int[][] resample(final int n, final int i) {
		// TODO Auto-generated method stub
		/**
		this.makeViterbi();
		this.printViterbiPath();
		
		Map<Character, Integer> stats = new HashMap<Character, Integer>();
		int w;
		Character c;
		for(int i=0; i<this.sample.length; i++) {
			String[] path = this.vb[i].path_str;
			for(int j=0; j<path.length; j++) {
				w = path[j].length();
				for(int k=0; k<w; k+=2) {
					c = path[j].charAt(k);
					if(!stats.containsKey(c))
						stats.put(c, 0);
					stats.put(c, 1+stats.get(c));
				}
			}
		}
		for(Map.Entry<Character, Integer> entry : stats.entrySet()) 
			System.out.println(entry.getKey()+", "+entry.getValue());
		**/
		
		// this.makeBackward();
		int[][] path = new int[this.M][n];
		int[] path_ij;
		int p;
		double[][] backProbMat_i;
		double[] backProbMat_ij, transProbMat_ij, emissProbMat_ij;
		final int S = this.str_statespace.length;
		final int[] hs_int = new int[S];
		EnumeratedIntegerDistribution dist;
		for(int j=0; j!=S; j++) hs_int[j] = j;
		backProbMat_i = this.backward[i].probsMat;

		backProbMat_ij = backProbMat_i[0];
		emissProbMat_ij = this.dp[0][i].emiss;
		transProbMat_ij = this.compoundTransProbs[0].probsMat[0];

		final double[] enumeratedIntegerProbs = new double[S];
		for(int s=0; s!=S; s++) 
			enumeratedIntegerProbs[s] = transProbMat_ij[s]*emissProbMat_ij[s]*backProbMat_ij[s];
		Algebra.normalize(enumeratedIntegerProbs);
		dist = new EnumeratedIntegerDistribution(hs_int, enumeratedIntegerProbs);
		path[0] = dist.sample(n);

		for(int j=1; j!=M; j++) {
			path_ij = path[j-1];
			final Map<Integer, double[]> enumeratedIntegerProbsBuffer = new HashMap<Integer, double[]>();
			for(int k=0; k!=n; k++) {
				p = path_ij[k];
				if(enumeratedIntegerProbsBuffer.containsKey(p)) {
					dist = new EnumeratedIntegerDistribution(hs_int, enumeratedIntegerProbsBuffer.get(p));
					path[j][k] = dist.sample();
				} else {
					backProbMat_ij = backProbMat_i[j];
					emissProbMat_ij = this.dp[j][i].emiss;
					transProbMat_ij = this.compoundTransProbs[j].probsMat[p];
					final double[] enumeratedIntegerProbs2 = new double[S];
					for(int s=0; s!=S; s++) 
						enumeratedIntegerProbs2[s] = transProbMat_ij[s]*emissProbMat_ij[s]*backProbMat_ij[s];
					Algebra.normalize(enumeratedIntegerProbs2);
					dist = new EnumeratedIntegerDistribution(hs_int, enumeratedIntegerProbs2);
					path[j][k] = dist.sample();
					enumeratedIntegerProbsBuffer.put(p, enumeratedIntegerProbs2);
				}
			}
		}
		
		return path;

		/**
		stats.clear();
		for(int i=0; i!=N; i++) {
			for(int j=0; j!=M; j++) {
				for(int k=0; k!=n; k++) {
					String pa = this.str_statespace[path[i][j][k]];
					w = pa.length();
					for(int l=0; l<w; l+=2) {
						c = pa.charAt(l);
						if(!stats.containsKey(c))
							stats.put(c, 0);
						stats.put(c, 1+stats.get(c));
					}
				}
			}
		}
		for(Map.Entry<Character, Integer> entry : stats.entrySet()) 
			System.out.println(entry.getKey()+", "+entry.getValue());
		
		int[][] head_stats = new int[N][this.hs()];
		for(int i=0; i!=N; i++) {
			for(int k=0; k!=n; k++) {
				String pa = this.str_statespace[path[i][0][k]];
				w = pa.length();
				for(int l=0; l<w; l+=2) {
					int m = pa.charAt(l)-'0';
					if(m>9) m = pa.charAt(l)-'a'+10;
					head_stats[i][m]++;
				}
			}
		}
		Utils.print(head_stats);
		**/
	}
	
	public void write(String output, 
			String experiment, 
			String scaff) {
		// TODO Auto-generated method stub
		write(output, experiment, scaff, 100);
	}
	
	public void write(String output, 
			String experiment, 
			String scaff,
			int resampling) {
		// TODO Auto-generated method stub
		
		this.makeViterbi();
		
		String root = experiment+"."+
				scaff+"."+
				System.nanoTime()+"_1_1_all_0_200mb";
		
		try {
			ZipOutputStream out = new ZipOutputStream(new BufferedOutputStream(new 
					FileOutputStream(output+"/"+root+".zip"), 65536));
			
			out.putNextEntry(new ZipEntry("phasedStates/"+experiment+".txt"));
			out.write((""+this.loglik()+"\n").getBytes());
			out.write((""+this.dp.length+"\n").getBytes());
			for(int i=0; i<this.sample.length; i++) {
				String[] path = this.vb[i].path_str;
				List<String[]> path_s = new ArrayList<String[]>();
				for(int j=0; j<path.length; j++)
					path_s.add(path[j].split("_"));
				for(int j=0; j<Constants._ploidy_H; j++) {
					out.write(("# id "+this.sample[i]+":"+(j+1)+"\t\t\t").getBytes());
					for(int k=0; k<path_s.size(); k++)
						out.write(path_s.get(k)[j].getBytes());
					out.write("\n".getBytes());
				}
			}
			
			if(resampling>0) {
				
				this.makeBackward();
				out.putNextEntry(new ZipEntry("resampling/"+experiment+".txt"));
				out.write((""+resampling+"\n").getBytes());
				out.write((""+this.dp.length+"\n").getBytes());
				final char[][] chars = new char[M][this.hs()];
				int w = Constants._ploidy_H;
				for(int i=0; i<this.sample.length; i++) {
					int[][] path = this.resample(resampling, i);
					for(int j=0; j<resampling; j++) {
						for(int k=0; k<M; k++) {
							String path_s = this.str_statespace[path[k][j]];	
							for(int l=0; l!=w; l++) 
								chars[k][l] = path_s.charAt(l*2);
						}
						for(int k=0; k!=w; k++) {
							out.write(("# id "+this.sample[i]+":"+(k+1)+"\t\t\t").getBytes());
							for(int l=0; l!=M; l++)
								out.write((""+chars[l][k]).getBytes());
							out.write("\n".getBytes());
						}	
					}
				}
			}
			
			out.putNextEntry(new ZipEntry("phased_genotypes"+experiment+".txt"));
			for(int i=0; i<this.N; i++) {
				int[] path = this.vb[i].path;
				for(int j=0; j<this.M; j++) {
					int[] states = this.statespace[path[j]].state;
					int dosa = 0;
					for(int k=0; k<states.length; k++) {
						if(this.emissProbs[j].probsMat[states[k]][0]<
								this.emissProbs[j].probsMat[states[k]][1])
							++dosa;
					}
					out.write((""+dosa+"\t").getBytes());
				}
				out.write("\n".getBytes());
			}
			
			out.putNextEntry(new ZipEntry("results_hmm/emissionModel.txt"));
			for(int i=0; i<this.emissProbs.length; i++) {
				double[][] emissProbs = this.emissProbs[i].probsMat;
				String[] allele = this.emissProbs[i].allele;
				out.write((this.de.getId()+"_"+this.de.getPosition()[i]+"\t\t\t").getBytes());
				for(int j=0; j<emissProbs.length; j++) {
					out.write((this.hs[j]+"-> {").getBytes());
					for(int k=0; k<emissProbs[j].length; k++) 
						out.write((allele[k]+","+emissProbs[j][k]+";").getBytes());
					out.write("} ".getBytes());
				}
				out.write("\n".getBytes());
			}

			out.putNextEntry(new ZipEntry("results_hmm/transitionModel.txt"));
			for(int i=0; i<this.transProbs.length; i++) {
				double[][] transProbs = this.transProbs[i].probsMat;
				out.write((this.de.getId()+"_"+this.de.getPosition()[i]+"\t\t\t").getBytes());
				for(int j=0; j<transProbs.length; j++) {
					out.write((this.hs[j]+"-> {").getBytes());
					for(int k=0; k<transProbs[j].length; k++) 
						out.write((this.hs[k]+","+transProbs[j][k]+";").getBytes());
					out.write("} ".getBytes());
				}
				out.write("\n".getBytes());
			}

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
}
