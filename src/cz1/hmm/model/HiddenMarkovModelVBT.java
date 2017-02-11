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
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import cz1.hmm.data.DataEntry;
import cz1.util.Algebra;
import cz1.util.Combination;
import cz1.util.Constants;
import cz1.util.Constants.Field;
import cz1.util.Utils;
import cz1.util.Permutation;
import cz1.util.Dirichlet;

public class HiddenMarkovModelVBT extends HiddenMarkovModel {
	private final static Logger logger = LogManager.getLogger(HiddenMarkovModelVBT.class.getName());

	public HiddenMarkovModelVBT(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse,
			boolean trainExp,
			Field field) {
		super(de, seperation, reverse, trainExp, field);
		
		this.trainAllExp();
		
		this.makeVBT();
		this.makeViterbi();
	}

	public void train() {
		iteration++;
		logger.info("###################");
		logger.info("train: "+iteration);
		long[] tic = new long[10];
		int k=0;
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
		
		this.updateDP();
		tic[k++] = System.nanoTime();
		logger.info("update dp "+(tic[k-1]-tic[k-2])+"ns");
		
		this.makeViterbi();
		tic[k++] = System.nanoTime();
		logger.info("make Viterbi "+(tic[k-1]-tic[k-2])+"ns");
		
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
			
			if(i==0 || this.exp_b.contains(i)) {
				double[][] trans_count = this.transProbs[i].pseudo();
				//clear(trans_count);
				for(int j=0;j<this.N; j++) {
					Viterbi vb1 = this.vb[j];
					final int a = i==0 ? 0 : vb1.path[i-1],
							b = vb1.path[i];
					double[] prior = this.transProbs[i].prior[a][b];
					s_a = this.statespace[a].state;
					st_b = this.statespace[b].stateperm;
					_t_ = st_b.length;
					for(int l=0; l<_t_; l++) {
						tmp = prior[l];
						s_b = st_b[l].state;
						for(int m=0; m<_p_; m++)
							trans_count[s_a[m]][s_b[m]] += tmp;
					}
				}
				this.transProbs[i].posterior();
			}
			
			EP ep1 = this.compoundEmissProbs[i];
			double[][] emiss_count = this.emissProbs[i].pseudo();
			//clear(emiss_count);
			double[] emissG = new double[_g_];
			for(int j=0;j<this.N; j++) {
				DP dp1 = this.dp[i][j];
				Viterbi vb1 = this.vb[j];
				final int a = vb1.path[i];
				Arrays.fill(emissG, 0);
				for(int k=0; k<_d_; k++) {
					Integer[] idx = dgMap.get(dosaS[k]);
					for(int c : idx)
						emissG[c] = dp1.likelihood[k]*
							ep1.probsMat[a][c]/dp1.emiss[a];
				}
				
				s_a = this.statespace[a].state;
				
				for(int c=0; c<_g_; c++) {
					A = emissG[c];
					s_l = genotype[c];
					for(int l=0; l<_p_; l++)
						emiss_count[s_a[l]][alMap.get(s_l[l])] += A;
				}
			}
			this.emissProbs[i].posterior();
		}
	}
	
	private static long[] scale_times = new long[10]; 
	
	private void makeVBT() {
		// TODO Auto-generated method stub
		this.transProbs = new TP[this.M];
		transProbs[0] = new TP(this.statespace,
				this.hs, -1, true, false);
		for(int i=1; i<this.M; i++) {
			System.out.print(i);
			transProbs[i] = exp_b.contains(i) ? 
				new TP(this.statespace, 
					this.hs, this.distance[i-1], 
					false, true) : 
				new TP(this.statespace, 
					this.hs, this.distance[i-1], 
					false, false);
		}
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
			for(Viterbi v : this.vb) probability += v.probability;
			return probability;
		}
	}
	
	public void print(boolean detail) {
		// TODO Auto-generated method stub
		double probability = 0;
		//for(Viterbi vb: this.vb) probability += vb.probability;
		//logger.info(" log-likelihood "+probability);
		if(iteration==0)
			logger.info(" hmm initialised.");
		else {
			for(Viterbi v : this.vb) probability += v.probability;
			logger.info(" log-likelihood "+probability);
		}

		if(detail) {
			logger.info("Viterbi path:");
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
	
	public void write(String output, 
			String experiment, 
			String scaff) {
		// TODO Auto-generated method stub
		String root = experiment+"."+
				scaff+"."+
				System.nanoTime()+"_1_1_all_0_200mb";
		
		try {
			ZipOutputStream out = new ZipOutputStream(new BufferedOutputStream(new 
					FileOutputStream(output+"/"+root+".zip"), 65536));
		
			out.putNextEntry(new ZipEntry(root+"/"));
			out.putNextEntry(new ZipEntry(root+"/"+"phasedStates/"));
			out.putNextEntry(new ZipEntry(root+"/"+"phasedStates/"+experiment+".txt"));
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
			
			out.putNextEntry(new ZipEntry(root+"/"+"results_hmm/"));
			out.putNextEntry(new ZipEntry(root+"/"+"results_hmm/emissionModel.txt"));
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

			out.putNextEntry(new ZipEntry(root+"/"+"results_hmm/transitionModel.txt"));
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

			out.putNextEntry(new ZipEntry(root+"/"+"stderr_true"));
			
			out.write("cz1.model.HiddenMarkovModel:\n".getBytes());
			out.write("cz1.model.HiidenMarkovModel$EM:\n".getBytes());
			out.write(("log prob is "+this.loglik()+" at "+iteration+"\n").getBytes());
			out.write(("random seed is "+Constants.seed).getBytes());
			
			out.putNextEntry(new ZipEntry(root+"/"+"snp_"+experiment+".txt"));
			
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
