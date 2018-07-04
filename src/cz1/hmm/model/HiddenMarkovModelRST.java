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

import cz1.hmm.data.DataEntry;
import cz1.hmm.model.HiddenMarkovModel.DP;
import cz1.hmm.model.HiddenMarkovModel.EP;
import cz1.hmm.model.HiddenMarkovModel.Viterbi;
import cz1.util.Algebra;
import cz1.util.Combination;
import cz1.util.Constants;
import cz1.util.Utils;
import cz1.util.Permutation;
import cz1.util.Dirichlet;
import cz1.util.Constants.Field;

public class HiddenMarkovModelRST extends HiddenMarkovModel {
	private FB[] backward;
	private int resample_no = 100;
	
	// ----founder haplotypes coefficients----
	private double founder_hap_coeff = 0.5;
	
	public HiddenMarkovModelRST(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse,
			boolean trainExp,
			Field field,
			double fh_coeff,
			int resample_no) {
		super(de, seperation, reverse, trainExp, field);
		this.founder_hap_coeff = fh_coeff;
		this.resample_no = resample_no;
		//this.trainAllExp();
		this.makeBWT();
		//this.makeViterbi();
		//this.print();
		//this.train();
	}
	
	public HiddenMarkovModelRST(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse,
			boolean trainExp,
			Field field,
			int resample_no) {
		super(de, seperation, reverse, trainExp, field);
		this.resample_no = resample_no;
		//this.trainAllExp();
		this.makeBWT();
		//this.makeViterbi();
		//this.print();
		//this.train();
	}
	
	public HiddenMarkovModelRST(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse,
			boolean trainExp,
			Field field,
			double fh_coeff,
			String hmmFile) {
		super(de, seperation, reverse, trainExp, field);
		this.founder_hap_coeff = fh_coeff;
		this.makeBWT(hmmFile);
	}
	
	public HiddenMarkovModelRST(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse,
			boolean trainExp,
			Field field,
			String hmmFile) {
		super(de, seperation, reverse, trainExp, field);
		this.makeBWT(hmmFile);
	}

	private int[][][] rs_path = null;
			
	public void train() {
		iteration++;
		myLogger.info("###################");
		myLogger.info("train: "+iteration);
		long[] tic = new long[10];
		int k=0;
		tic[k++] = System.nanoTime();
		// this.makeForward();
		// tic[k++] = System.nanoTime();
		// logger.info("forward done "+(tic[k-1]-tic[k-2])+"ns");
		this.makeBackward();
		tic[k++] = System.nanoTime();
		myLogger.info("backward done "+(tic[k-1]-tic[k-2])+"ns");
		// this.checkFW();
		// tic[k++] = System.nanoTime();
		rs_path = this.resample(this.resample_no);
		tic[k++] = System.nanoTime();
		myLogger.info("Resampling "+(tic[k-1]-tic[k-2])+"ns");
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

	protected void EM() {
		// TODO Auto-generated method stub
		
		int[] s_a, s_b; 
		String[] s_l;
		ST[] st_b;
		double tmp, A;
		int _p_ = Constants._ploidy_H;
		int _g_, _d_, _t_;

		for(int i=0; i<M; i++) {
			_g_ = this.obspace[i].genotype.length;
			OB ob = this.obspace[i];
			Map<String, Integer[]> dgMap = ob.dgMap;
			String[] dosaS = ob.dosaS;
			_d_ = dosaS.length;
			String[][] genotype = ob.genotype;
			String[] allele = ob.allele;
			Map<String, Integer> alMap = new HashMap<String, Integer>();
			for(int j=0; j<allele.length; j++) alMap.put(allele[j], j);
			
			EP ep1 = this.compoundEmissProbs[i];
			double[][] emiss_count = this.emissProbs[i].pseudo();
			//clear(emiss_count);
			double[] emissG = new double[_g_];
			
			double[][] trans_count = this.transProbs[i].pseudo();
			for(int j=0;j<N; j++) {
				DP dp1 = this.dp[i][j];
				double coeff = dp1.isparent ? ((N-2)*founder_hap_coeff) : 1.0;
				
				for(int k=0;k<resample_no; k++) {
					// transfer from stats_fr(om) to stats_to
					int stats_fr = i==0 ? 0 : rs_path[j][i-1][k];
					int stats_to = rs_path[j][i][k];
					
					double[] prior = this.transProbs[i].prior[stats_fr][stats_to];
					s_a = this.statespace[stats_fr].state;
					st_b = this.statespace[stats_to].stateperm;
					_t_ = st_b.length;
					for(int l=0; l<_t_; l++) {
						tmp = coeff*prior[l];
						s_b = st_b[l].state;
						for(int m=0; m<_p_; m++)
							trans_count[s_a[m]][s_b[m]] += tmp;
					}
					
					Arrays.fill(emissG, 0);
					for(int l=0; l<_d_; l++) {
						Integer[] idx = dgMap.get(dosaS[l]);
						for(int c : idx)
							emissG[c] = dp1.likelihood[l]*
								ep1.probsMat[stats_to][c]/dp1.emiss[stats_to];
					}
					
					s_a = this.statespace[stats_to].state;
					for(int c=0; c<_g_; c++) {
						A = coeff*emissG[c];
						s_l = genotype[c];
						for(int l=0; l<_p_; l++)
							emiss_count[s_a[l]][alMap.get(s_l[l])] += A;
					}
				}
			}
			
			this.transProbs[i].posterior();
			this.emissProbs[i].posterior();
		}
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

	private void makeBWT(String hmmFile) {
		// TODO Auto-generated method stub
		
		// load a hmm object from a .zip hmmFile
		try {
			final ZipFile in = new ZipFile(hmmFile);

			final BufferedReader br_trans = Utils.getBufferedReader(
					in.getInputStream(in.getEntry("results_hmm/transitionModel.txt")));
			this.transProbs = new TP[this.M];
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

			final BufferedReader br_emiss = Utils.getBufferedReader(
					in.getInputStream(in.getEntry("results_hmm/emissionModel.txt")));
			this.emissProbs = new EP[this.M];
			for(int i=0; i<this.M; i++) {
				emissProbs[i] = new EP(
						this.hs, 
						this.obspace[i], 
						this.bfrac[i],
						br_emiss.readLine());
			}
			br_emiss.close();

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
		this.backward = new FB[_k_];
		for(int i=0; i<_k_; i++) 
			this.backward[i] = new FB(true, _m_+1, _n_);

		this.vb = new Viterbi[_k_];
		for(int i=0; i<_k_; i++)
			this.vb[i] = new Viterbi(_m_+1, _n_, this.str_statespace);
		return;
	}

	public double loglik() {
		return loglik1();
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
			for(FB bw : this.backward) probability += bw.probability;
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
		write(output, experiment, scaff, 100, .0);
	}
	
	public void write(String output, 
			String experiment, 
			String scaff,
			int resampling,
			double loglik_diff) {
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
