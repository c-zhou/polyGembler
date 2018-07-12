package cz1.hmm.model;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.File;
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

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.log4j.Logger;

import cz1.breeding.data.FullSiblings;
import cz1.hmm.data.DataEntry;
import cz1.util.Algebra;
import cz1.util.Combination;
import cz1.util.Constants;
import cz1.util.Constants.Field;
import cz1.util.Utils;
import cz1.util.Permutation;
import cz1.util.Dirichlet;

public class HiddenMarkovModelVBT extends HiddenMarkovModel {
	
	// ----founder haplotypes coefficients----
	private double founder_hap_coeff = 0.5;
	// higher quality data get more weights for training the model
	private double[][] weights;
		
	public HiddenMarkovModelVBT(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse,
			boolean trainExp,
			Field field) {
		super(de, seperation, reverse, trainExp, field);
		
		this.makeWeightingCoeff();
		this.trainAllExp();
		this.trainAllExpExceptDummy();
		this.makeVBT();
		this.makeViterbi();
	}
	
	public HiddenMarkovModelVBT(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse,
			boolean trainExp,
			Field field,
			double fh_coeff) {
		super(de, seperation, reverse, trainExp, field);
		
		this.founder_hap_coeff = fh_coeff;
		this.makeWeightingCoeff();
		this.trainAllExpExceptDummy();
		this.makeVBT();
		this.makeViterbi();
	}

	private void trainAllExpExceptDummy() {
		// TODO Auto-generated method stub
		this.trainAllExp();
		//this.exp_b.remove(0);
	}

	public void train() {
		iteration++;
		myLogger.info("###################");
		myLogger.info("train: "+iteration);
		long[] tic = new long[10];
		int k=0;
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
		
		this.updateDP();
		tic[k++] = System.nanoTime();
		myLogger.info("update dp "+(tic[k-1]-tic[k-2])+"ns");
		
		this.makeViterbi();
		tic[k++] = System.nanoTime();
		myLogger.info("make Viterbi "+(tic[k-1]-tic[k-2])+"ns");
		
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
	
	private int getDosage(final int i, final int j) {
		// TODO Auto-generated method stub
		double[] ll = this.dp[i][j].likelihood;
		return Algebra.maxIndex(ll);
	}
	
	private final double loglik_diff = 0.05129329; // 5% difference of likelihood
	
	protected void EM() {
		// TODO Auto-generated method stub
		int[] s_a, s_b; 
		String[] s_l;
		ST[] st_b;
		double tmp, A;
		int _p_ = Constants._ploidy_H;
		int _g_, _d_, _t_;
		
		final double[][] probs = new double[this.N][];
		for(int i=0; i<this.N; i++) {
			Viterbi vb1 = this.vb[i];
			double ll = vb1.probability(0);
			int z = 0;
			while(ll-vb1.probability(z)<=loglik_diff) ++z;
			final double[] prob = new double[z];
			System.arraycopy(vb1.probability, 0, prob, 0, z);
			for(int k=0; k<z; k++) {
				prob[k] -= ll;
				prob[k] = Math.exp(prob[k]);
			}
			Algebra.normalize(prob);
			probs[i] = prob;
		}
		
		for(int i=0; i<this.M; i++) {
			_g_ = this.obspace[i].genotype.length;
			OB ob = this.obspace[i];
			Map<String, Integer[]> dgMap = ob.dgMap;
			String[] dosaS = ob.dosaS;
			_d_ = dosaS.length;
			String[][] genotype = ob.genotype;
			String[] allele = ob.allele;
			Map<String, Integer> alMap = new HashMap<String, Integer>();
			for(int j=0; j<allele.length; j++) alMap.put(allele[j], j);
			
			double[][] trans_count = this.transProbs[i].pseudo();
			CompoundEP ep1 = this.compoundEmissProbs[i];
			double[][] emiss_count = this.emissProbs[i].pseudo();
			double[] emissG = new double[_g_];
			
			for(int j=0;j<this.N; j++) {
				
				DP dp1 = this.dp[i][j];
				double coeff = weights[i][j];
				
				Viterbi vb1 = this.vb[j];

				
				final double[] ws = probs[j];
				for(int w=0; w<ws.length; w++) {
					
					final int a = i==0 ? 0 : vb1.path[w][i-1],
							b = vb1.path[w][i];
				
					if(this.exp_b.contains(i)) {
						double[] prior = this.compoundTransProbs[i].prior[a][b];
						s_a = this.statespace[a].state;
						st_b = this.statespace[b].stateperm;
						_t_ = st_b.length;

						for(int l=0; l<_t_; l++) {
							tmp = prior[l]*ws[w];
							s_b = st_b[l].state;
							for(int m=0; m<_p_; m++)
								trans_count[s_a[m]][s_b[m]] += tmp;
						}
					}
					
					Arrays.fill(emissG, 0);
					for(int k=0; k<_d_; k++) {
						Integer[] idx = dgMap.get(dosaS[k]);
						for(int c : idx)
							emissG[c] = dp1.likelihood[k]*
							ep1.probsMat[b][c]/ep1.probsDosaMat[b][k];
					}
					
					s_a = this.statespace[b].state;

					for(int c=0; c<_g_; c++) {
						A = coeff*emissG[c]*ws[w];
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
	
	private static long[] scale_times = new long[10]; 
	
	private void makeVBT() {
		// TODO Auto-generated method stub
		this.transProbs = new TP[this.M];
		transProbs[0] = new TP(this.hs, -1, true, false);
		for(int i=1; i<this.M; i++) {
			System.out.print(i);
			transProbs[i] = exp_b.contains(i) ? 
				new TP(this.hs, this.distance[i-1], 
					false, true) : 
				new TP(this.hs, this.distance[i-1], 
					false, false);
		}
		this.emissProbs = new EP[this.M];
		for(int i=0; i<this.M; i++)
			emissProbs[i] = new EP(
				this.hs, 
				this.obspace[i].allele, 
				this.bfrac[i]);

		this.compoundTransProbs = new CompoundTP[this.M];
		int _n_ = this.statespace.length;
		for(int i=0; i<this.M; i++) {
			this.compoundTransProbs[i] = new CompoundTP(
					this.transProbs[i],
					this.statespace);
		}
		
		this.compoundEmissProbs = new CompoundEP[this.M];
		for(int i=0; i<this.M; i++) {
			this.compoundEmissProbs[i] = new CompoundEP(
				this.emissProbs[i],
				this.obspace[i],
				this.statespace);
		}
		this.updateDP();
		
		int _m_ = this.dp.length;
		int _k_ = this.dp[0].length;
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
			for(Viterbi v : this.vb) probability += v.probability();
			return probability;
		}
	}
	
	public int[][][] resample(final int n){
		// TODO Auto-generated method stub
		return null;
	}
	
	
	public void print(boolean detail) {
		// TODO Auto-generated method stub
		double probability = 0;
		//for(Viterbi vb: this.vb) probability += vb.probability;
		//logger.info(" log-likelihood "+probability);
		if(iteration==0)
			myLogger.info(" hmm initialised.");
		else {
			for(Viterbi v : this.vb) probability += v.probability();
			myLogger.info(" log-likelihood "+probability);
		}

		if(detail) {
			myLogger.info("Viterbi path:");
			for(int i=0; i<this.sample.length; i++) {
				String[] path = this.vb[i].path_str[0];
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
				String[] path = this.vb[i].path_str[0];
				List<String[]> path_s = new ArrayList<String[]>();
				for(int k=0; k<path.length; k++)
					path_s.add(path[k].split("_"));
				for(int k=0; k<Constants._ploidy_H; k++) {
					out.write(("# id "+this.sample[i]+":"+(k+1)+"\t\t\t").getBytes());
					for(int s=0; s<path_s.size(); s++)
						out.write(path_s.get(s)[k].getBytes());
					out.write("\n".getBytes());
				}
			}
			
			out.putNextEntry(new ZipEntry("phasedStates/"+experiment+"_ld.txt"));
			out.write((""+this.loglik()+"\n").getBytes());
			out.write((""+this.dp.length+"\n").getBytes());
			for(int i=0; i<this.sample.length; i++) {
				Viterbi vb1 = this.vb[i];
				double ll = vb1.probability();
				int n = vb1.probability.length;
				for(int j=0; j<n; j++) {
					if(ll-vb1.probability(j)>loglik_diff) break;
					String[] path = vb1.path_str[j];
					List<String[]> path_s = new ArrayList<String[]>();
					for(int k=0; k<path.length; k++)
						path_s.add(path[k].split("_"));
					for(int k=0; k<Constants._ploidy_H; k++) {
						out.write(("# id "+this.sample[i]+":"+(k+1)+"\t\t\t").getBytes());
						for(int s=0; s<path_s.size(); s++)
							out.write(path_s.get(s)[k].getBytes());
						out.write("\n".getBytes());
					}
				}
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
