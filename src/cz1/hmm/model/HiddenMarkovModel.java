package cz1.hmm.model;

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
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.log4j.Logger;

import cz1.hmm.data.DataEntry;
import cz1.hmm.tools.VCFtools;
import cz1.util.Algebra;
import cz1.util.Combination;
import cz1.util.Constants;
import cz1.util.Constants.Field;
import cz1.util.Utils;
import cz1.util.Permutation;
import cz1.util.Dirichlet;

public abstract class HiddenMarkovModel {
	protected final static Logger myLogger = Logger.getLogger(HiddenMarkovModel.class);
	
	protected final static Runtime runtime = Runtime.getRuntime();
	
	protected int iteration = 0;
	
	protected final Field field;
	protected final DataEntry de;
	protected final DP[][] dp;
	protected final ST[] statespace;
	protected final OB[] obspace;
	protected final String[] hs; // all hidden states 1,2,3,4,5,6,7,8
	protected final String[] str_statespace;
	protected final double[] bfrac;
	protected final boolean train_exp;
	protected final Set<Integer> exp_b = new HashSet<Integer>();

	protected double[] distance;
	protected TP[] transProbs;
	protected EP[] emissProbs;
	
	protected String[] sample;
	protected String[] parent;
	protected int[] parent_i;
	protected int N;
	protected int M;
	protected Viterbi vb[];
	protected boolean[][] miss;

	protected CompoundTP[] compoundTransProbs;
	protected CompoundEP[] compoundEmissProbs;
	
	protected final int[] pedigree; //0 and 1 represent two parents, and -1s represent F1 offspring
	protected final List<Integer[]> validStateSpace; //valid statespace for each sample

	protected final Map<String, int[][]> spMap;
	
	public HiddenMarkovModel(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse,
			boolean trainExp,
			Field field) {
		this.field = field;
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
	}
	
	public abstract void train();
	public abstract double loglik();
	public abstract void write(String output, 
			String experiment, 
			String contig);
	public abstract void write(String output, 
			String experiment, 
			String contig,
			int resampling,
			double loglik_diff);
	public abstract void print(boolean details);
	protected abstract void EM();
	public abstract int[][][] resample(final int n);
	public void resample() {
		this.resample(1);
	}
	public void print() {this.print(false);}
	
	protected void printViterbiPath() {
		for(int i=0; i<this.sample.length; i++) {
			String[] path = this.vb[i].path_str[0];
			List<String[]> path_s = new ArrayList<String[]>();
			for(int j=0; j<path.length; j++)
				path_s.add(path[j].split("_"));
			for(int j=0; j<Constants._ploidy_H; j++) {
				System.out.print(("# id "+this.sample[i]+":"+(j+1)+"\t\t\t"));
				for(int k=0; k<path_s.size(); k++)
					System.out.print(path_s.get(k)[j]);
				System.out.print("\n");
			}
		}
	}
	
	protected DataEntry catDE(DataEntry[] de2, 
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

	protected double[] bfrac() {
		// TODO Auto-generated method stub
		String[] sample = this.de.getSample();
		Set<Integer> p = new HashSet<Integer>();
		for(int i=0; i<sample.length; i++)
			if(!isParent(this.sample[i]))
				p.add(i);
		double[] baf = new double[this.M];
		
		switch(this.field) {
		case GT:
			List<List<String[]>> gt = this.de.getGenotype();
			if(gt==null) throw new RuntimeException("GT field not available!!! Try PL/GL (-L/--genotype-likelihood) or AD (-D/--allele-depth) option.");
			for(int i=0; i<this.M; i++) {
				double acnt = 0, bcnt = 0;
				String[] allele = this.de.getAllele().get(i);
				for(int j : p) {
					String[] g = gt.get(i).get(j);
					for(int k=0; k<g.length; k++) {
						acnt += (g[k].equals(allele[0]) ? 1 : 0);
						bcnt += (g[k].equals(allele[1]) ? 1 : 0);
					}
				}
				if(acnt+bcnt==0) baf[i] = .5;
				baf[i] = bcnt/(acnt+bcnt);
				if(baf[i]==0 || baf[i]==1) baf[i] = .5;
			}
			break;
		case PL:
		case GL:
			List<List<double[]>> gl = this.de.getGenotypeLikelihood();
			if(gl==null) throw new RuntimeException("PL/GL feild not available!!! Try GT (-G/--genotype) or AD (-D/--allele-depth) option.");
			int cp = gl.get(0).get(0).length-1;
			for(int i=0; i<this.M; i++) {
				double acnt = 0;
				double bcnt = 0;
				for(int j : p) {
					double[] ll = gl.get(i).get(j);
					if(ll[0]<0) continue;
					for(int k=0; k<ll.length; k++) {
						acnt += ll[k]*(cp-k);
						bcnt += ll[k]*k;
					}
				}
				if(acnt+bcnt==0) baf[i] = .5;
				baf[i] = bcnt/(acnt+bcnt);
				if(baf[i]==0 || baf[i]==1) baf[i] = .5;
			}
			break;

		case AD:
			List<List<int[]>> ad = this.de.getAlleleDepth();
			if(ad==null) throw new RuntimeException("AD feild not available!!! Try PL/GL (-L/--genotype-likelihood) or GT (-G/--genotype) options.");
			for(int i=0; i<this.M; i++) {
				double acnt = 0, bcnt = 0;
				for(int j : p) {
					int[] aa = ad.get(i).get(j);
					acnt += aa[0];
					bcnt += aa[1];
				}
				if(acnt+bcnt==0) baf[i] = .5;
				baf[i] = bcnt/(acnt+bcnt);
				if(baf[i]==0 || baf[i]==1) baf[i] = .5;
			}
			break;
		default:
			throw new RuntimeException("Undefined Field!!!");
		}
		return baf;
	}

	protected void trainAllExp() {
		// TODO Auto-generated method stub
		for(int i=0; i<this.M; i++)
			this.exp_b.add(i);
	}
	
	protected String[] makeSS() {
		// TODO Auto-generated method stub
		String[] str_states = new String[this.statespace.length];
		for(int i=0; i<str_states.length; i++) 
			str_states[i] = StringUtils.join(
					this.statespace[i].state_str,"_");
		return str_states;
	}

	protected Map<String, int[][]> makeSpMap() {
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

	protected List<Integer[]> validStateSpace() {
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

	protected int[] pedigree() {
		// TODO Auto-generated method stub
		int pedigree[] = new int[this.sample.length];
		this.parent_i = new int[this.parent.length];
		Arrays.fill(pedigree, -1);
		Arrays.fill(this.parent_i, -1);
		for(int i=0; i<pedigree.length; i++) {
			for(int j=0; j<this.parent.length; j++) {
				if(this.sample[i].equals(this.parent[j])) {
					pedigree[i] = j;
					parent_i[j] = i;
				}
			}
		}
		return pedigree;
	}

	protected void initialise() {
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

	protected void printExpTrain() {
		// TODO Auto-generated method stub
		for(final Integer i : this.exp_b) {
			double cumR = 0.0;
			double[][] tp_i = this.transProbs[i].probsMat;
			for(int j=0; j<tp_i.length; j++)
				for(int k=0; k<tp_i[j].length; k++)
					if(j!=k) cumR+=tp_i[j][k];
			myLogger.info("********** sep "+i+" "+cumR/(tp_i.length-1));
		}
	}

	protected void updateDP() {
		// TODO Auto-generated method stub
		for(int i=0; i<M; i++) {
			String[] dosaS = this.obspace[i].dosaS;
			int _d_ = dosaS.length;
			// Map<String, Integer[]> dgMap = this.obspace[i].dgMap;
			CompoundEP ep1 = this.compoundEmissProbs[i];
			for(int j=0; j<N; j++) {
				Integer[] vst = this.validStateSpace.get(j);
				
				DP dp1 = this.dp[i][j];
				for(int k : vst) {
					double ss = 0;
					for(int a=0; a<_d_; a++) {
						/***
						Integer[] ix = dgMap.get(dosaS[a]);
						for(int c : ix)
							ss += dp1.likelihood[a]*ep1.probsMat[k][c];
						 **/
						ss += dp1.likelihood[a]*ep1.probsDosaMat[k][a];
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

	protected double[][] normalise(double[][] mat, boolean byrow,
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

	protected double[][] logspace(double[][] mat) {
		// TODO Auto-generated method stub
		for(int i=0; i<mat.length; i++)
			mat[i] = logspace(mat[i]);
		return mat;
	}

	protected double[] logspace(double[] array) {
		// TODO Auto-generated method stub
		for(int i=0; i<array.length; i++)
			array[i] = Math.log(array[i]);
		return array;
	}

	protected double[][] normalspace(double[][] mat) {
		// TODO Auto-generated method stub
		for(int i=0; i<mat.length; i++)
			mat[i] = normalspace(mat[i]);
		return mat;
	}

	protected double[] normalspace(double[] array) {
		// TODO Auto-generated method stub
		for(int i=0; i<array.length; i++)
			array[i] = Math.exp(array[i]);
		return array;
	}

	protected double[][] transpose(double[][] mat) {
		// TODO Auto-generated method stub
		double[][] tMat = new double[mat[0].length][mat.length];
		for (int i = 0; i < mat.length; i++)
			for (int j = 0; j < mat[0].length; j++)
				tMat[j][i] = mat[i][j];
		return tMat;
	}

	protected void clear(double[][][] matrix) {
		// TODO Auto-generated method stub
		int r = matrix.length, 
				c = matrix[0].length;
		for(int i=0; i<r; i++)
			for(int j=0; j<c; j++)
				Arrays.fill(matrix[i][j], 0);

	}

	protected void clear(double[][] matrix) {
		// TODO Auto-generated method stub
		int r = matrix.length;
		for(int i=0; i<r; i++)
			Arrays.fill(matrix[i], 0);
	}

	protected void clear(int[][] matrix) {
		// TODO Auto-generated method stub
		int r = matrix.length;
		for(int i=0; i<r; i++)
			Arrays.fill(matrix[i], 0);
	}

	protected int maxIndex(int[][] trace, 
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

	protected String[] makeHS() {
		// TODO Auto-generated method stub
		String[] hs = new String[Constants._haplotype_z+1];
		hs[0] = ".";
		for(int k=1; k<hs.length; k++)
			hs[k] = k<10 ? ""+k : ""+(char)('a'+k-10);
		return hs;
	}

	protected void makeCompoundEP() {
		// TODO Auto-generated method stub
		for(final CompoundEP ep : this.compoundEmissProbs) ep.update();
	}

	protected void makeCompoundTP() {
		// TODO Auto-generated method stub
		for(int i=0; i<this.M; i++)
			if(this.exp_b.contains(i))  
				this.compoundTransProbs[i].update();
	}
	
	protected double makeViterbi() {
		// TODO Auto-generated method stub
		int _m_ = this.dp.length+1;
		int _n_ = this.dp[0].length;
		
		double probability = 0.0;
		for(int i=0; i<_n_; i++) {
			
			Integer[] vst = this.validStateSpace.get(i);
			Set<Integer> vst_copy;
			double[][] v = this.vb[i].v;
			int[][] trace = this.vb[i].trace;
			clear(v);
			clear(trace);
			v[0][0] = 1.0;
			
			double lower_bound = Double.NEGATIVE_INFINITY;
			for(int k : vst) {
				double x = this.dp[0][i].emiss[k]*
						this.compoundTransProbs[0].probsMat[0][k];
				double logscale = 0.0;
				for(int j=2; j<_m_; j++) {
					x *= this.dp[j-1][i].emiss[k]*
							this.compoundTransProbs[j-1].probsMat[k][k];
					if(x<Constants.threshMin) {
						logscale += Constants.logThreshMax;
						x /= Constants.threshMax;
					}
				}
				x = Math.log(x)+logscale;
				if(lower_bound < x) lower_bound = x;
			}
			
			for(int j=1; j<_m_; j++) {
				
				double logscale_j = this.vb[i].logscale[j];
				
				vst_copy = new HashSet<Integer>(Arrays.asList(vst));
				for(int k : vst)
					if(Math.log(v[j-1][k])+logscale_j < lower_bound)
						vst_copy.remove(k);
				
				DP dp1 = this.dp[j-1][i];
				CompoundTP tp1 = this.compoundTransProbs[j-1];
				
				for(int k : vst) {
					double ss = 0.0, ss0 = 0.0, ss_tmp;
					int st = -1;
					
					for(int a : vst_copy) {
						ss_tmp = v[j-1][a]*tp1.probsMat[a][k];
						if(ss_tmp > ss) {
							ss = ss_tmp;
							st = a;
						}
						if(a==k) ss0 = ss_tmp;
					}
					
					if(ss0==ss) st = k;
					trace[j-1][k] = st;
					v[j][k] = dp1.emiss[k]*ss;
				}
				
				this.vb[i].scale(j);
			}
			
			this.vb[i].finalise();
			this.vb[i].traceAll();
			probability += this.vb[i].probability();
		}
		
		return probability;
	}

	protected OB[] makeOB() {
		// TODO Auto-generated method stub
		List<String[]> allele = this.de.getAllele();
		OB[] obs = new OB[allele.size()];
		for(int i=0; i<obs.length; i++)
			obs[i] = new OB(allele.get(i));
		return obs;
	}

	protected ST[] makeST() {
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

	protected DP[][] makeDP() {
		// TODO Auto-generated method stub
		DP[][] dp = new DP[this.M][this.N];
		miss = new boolean[this.M][this.N];
		switch(this.field) {
		
		case GT:
			//TO-DO: calculate genotype likelihoods from genotype
			List<List<String[]>> gt = this.de.getGenotype();
			List<String[]> allele = this.de.getAllele();
			if(gt==null) throw new RuntimeException("GT field not available!!! Try PL/GL (-L/--genotype-likelihood) or AD (-D/--allele-depth) option.");
			for(int i=0; i<this.M; i++) {
				Arrays.fill(miss[i], true);
				for(int j=0; j<this.N; j++) {
					double[] ll =  VCFtools.fitGL(gt.get(i).get(j), allele.get(i)[0]);
					double l = 0;
					for(int k=1; k<ll.length; k++) {
						if(ll[k]!=l) {
							miss[i][j] = false;
							break;
						}
					}
					dp[i][j] = new DP(ll,
							this.statespace.length,
							this.obspace[i].genotype.length,
							isParent(this.sample[j]),
							false);
				}
			}
			break;
			
		case PL:
		case GL:
			List<List<double[]>> gl = this.de.getGenotypeLikelihood();
			if(gl==null) throw new RuntimeException("PL/GL feild not available!!! Try GT (-G/--genotype) or AD (-D/--allele-depth) option.");
			for(int i=0; i<this.M; i++) {
				for(int j=0; j<this.N; j++) {
					double[] ll = gl.get(i).get(j);
					double l = 0;
					for(int k=1; k<ll.length; k++) {
						if(ll[k]!=l) {
							miss[i][j] = false;
							break;
						}
					}
					dp[i][j] = new DP(ll,
							this.statespace.length,
							this.obspace[i].genotype.length,
							isParent(this.sample[j]),
							false);
				}
			}
			break;
			
		case AD:
			//TO-DO: calculate genotype likelihoods from allele depth
			List<List<int[]>> ad = this.de.getAlleleDepth();
			if(ad==null) throw new RuntimeException("AD feild not available!!! Try PL/GL (-L/--genotype-likelihood) or GT (-G/--genotype) options.");
			for(int i=0; i<this.M; i++) {
				for(int j=0; j<this.N; j++) {
					double[] ll = VCFtools.fitGL(ad.get(i).get(j), Constants._ploidy_H);
					double l = 0;
					for(int k=1; k<ll.length; k++) {
						if(ll[k]!=l) {
							miss[i][j] = false;
							break;
						}
					}
					dp[i][j] = new DP(ll, 
							this.statespace.length,
							this.obspace[i].genotype.length,
							isParent(this.sample[j]),
							false);
				}
			}
			break;
		
		default:
			throw new RuntimeException("Invalid field provided. "
					+ "Should be \"PL\", \"GL\", \"AD\", or \"GT\". ");
		}
		return dp;
	}

	protected boolean isParent(String sample) {
		// TODO Auto-generated method stub
		return new HashSet<String>(Arrays.asList(this.parent)).contains(sample);
	}

	protected long memory(String type) {
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
	
	public class OB {
		protected final String[] allele;
		protected final String[][] genotype;
		protected final String[][] dosage;
		protected final String[] dosaS;
		protected final String[] genoS;

		protected Map<String, Integer[]> dgMap; // dosage and genotype map

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

		public String[] allele() {
			// TODO Auto-generated method stub
			return this.allele;
		}
	}
	
	protected class ST {
		protected final boolean isparent; // is valid for parents only
		protected final int[] state;
		protected final String[] state_str;
		protected final boolean isdupstate;
		protected final ST[] stateperm;

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

		protected ST[] stateperm() {
			// TODO Auto-generated method stub
			if(!isdupstate||isparent) return new ST[]{
					new ST(this.state,
							this.state_str,
							this.isparent)};
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
	
	protected class DP {

		protected final boolean isparent; // is parent
		protected final double[] likelihood; // likelihood of different dosage
		protected final double[] hweDist1;
		protected boolean allowtrans; // if transition is allowed
		protected boolean logspace; // if is in log space
		protected double[] emiss; // product of likelihood and emiss for dosa
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
				Arrays.fill(likelihood, 1.0);
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
	
	public class CompoundTP {
		protected final TP tp;
		protected final ST[] statespace;
		
		protected double[][][] prior;
		protected double[][] probsMat;
		
		public CompoundTP(final TP tp,
				final ST[] statespace) {
			this.tp = tp;
			this.statespace = statespace;
			final int _m_ = this.statespace.length;
			this.probsMat = new double[_m_][_m_];
			this.prior = new double[_m_][_m_][];
			this.update();
		}

		private void update() {
			// TODO Auto-generated method stub
			int _a_ = this.statespace.length;
			for(int a=0; a<_a_; a++) {
				int[] s_a = this.statespace[a].state;
				for(int b=0; b<_a_; b++) { 
					ST[] st_b = this.statespace[b].stateperm();
					int _i_ = st_b.length;
					int _j_ = st_b[0].state.length;
					
					if(this.prior[a][b]==null) 
						this.prior[a][b] = new double[_i_];
					final double[] prior = this.prior[a][b];
				
					Arrays.fill(prior, 1.0);
					for(int i=0; i<_i_; i++) {
						int[] s_i = st_b[i].state;
						for(int j=0; j<_j_; j++)
							prior[i] *= tp.probsMat[s_a[j]][s_i[j]];
					}
					
					this.probsMat[a][b] = StatUtils.sum(prior);
					Algebra.normalize(prior);
				}
			}
		}
	}
	
	public class TP {
		protected final boolean isDotState;
		protected final double distance;
		protected final String[] str_statespace;
		protected final boolean trainExp;
		protected final double[][] probsMat;
		
		protected boolean logspace;
		protected double[][] count;
		protected double exp;
		protected double[][] alpha;

		public TP(String[] str_statespace, 
				double distance, 
				boolean isDotState,
				boolean trainExp) {
			this.str_statespace = str_statespace;
			this.distance = distance;
			this.isDotState = isDotState;
			this.prior();
			this.trainExp = trainExp;
			int _n_=this.str_statespace.length;
			this.probsMat = new double[_n_][_n_];
			this.count = new double[_n_][_n_];
			this.logspace = false;
			this.update();
		}

		public TP(String[] str_statespace, 
				double distance, 
				boolean isDotState, 
				boolean trainExp, 
				String tranProbs) {
			// TODO Auto-generated constructor stub
			this.str_statespace = str_statespace;
			this.distance = distance;
			this.isDotState = isDotState;
			this.prior();
			this.trainExp = trainExp;
			int _n_=this.str_statespace.length;
			this.probsMat = new double[_n_][_n_];
			this.count = new double[_n_][_n_];
			this.logspace = false;
			this.update(tranProbs);
		}
		
		protected void update() {
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
		}
		
		private void update(String tranProbs) {
			// TODO Auto-generated method stub
			Pattern p = Pattern.compile("\\{(.*?)\\}");
			Matcher m = p.matcher(tranProbs);
			
			int _n_ = this.str_statespace.length;
			String[] prob_str;
			double[] probs;
			for(int a=0; a<_n_; a++) {
				probs = this.probsMat[a];
				m.find();
				prob_str = m.group(1).split(",|;");
				for(int b=0; b<_n_; b++) 
					probs[b] = Double.parseDouble(prob_str[b*2+1]);
			}
		}
		
		private void prior() {
			// TODO Auto-generated method stub	
			double e = Math.exp(-this.distance*Constants._con_base_r);
			double alpha = Math.max(Constants._mu_J_e*(1-e), Constants.eps);
			double beta = Math.max(Constants._mu_J_e*e, Constants.eps);
			this.exp = new BetaDistribution(Constants.rg, alpha, beta).sample();
			// TODO should be fixed? should be a smaller number?
			// this.exp = Math.max(this.exp, 1e-3);
			int z = Constants._haplotype_z/2;
			this.alpha = new double[2][z];
			//Dirichlet diri = new Dirichlet(distz(z), Constants._mu_alpha_e);
			for(int i=0; i<2; i++)
				//this.alpha[i] = ArrayUtils.toPrimitive(diri.sample());
				Arrays.fill(this.alpha[i], 1.0/z);
			return;
		}

		protected void posterior() {
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
			
			/***
			int _n_=this.str_statespace.length;
			for(int i=0; i<_n_; i++) {
				Arrays.fill(this.count[i], .0);
			}
			***/
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

		public double[][] probs() {
			// TODO Auto-generated method stub
			return this.probsMat;
		}
	}

	public class CompoundEP {
		protected final EP ep;
		protected final OB ob;
		protected final ST[] statespace;
		protected final double[][] probsMat;
		protected final double[][] probsDosaMat;

		public CompoundEP(EP ep,
				OB obspace,
				ST[] statespace) {
			// TODO Auto-generated constructor stub
			this.ep = ep;
			this.ob = obspace;
			this.statespace = statespace;
			int _n_ = this.statespace.length;
			int _g_ = obspace.genotype.length;
			int _k_ = obspace.dosage.length;
			this.probsMat = new double[_n_][_g_];
			this.probsDosaMat = new double[_n_][_k_];
			this.update();
		}

		private void update() {
			// TODO Auto-generated method stub
			final int _s_ = this.statespace.length;
			
			for(int i=0; i<_s_; i++) {
				int[] state = this.statespace[i].state;
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
					this.probsMat[i][a] = ss;
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
						ss += this.probsMat[i][idx[l]];
					this.probsDosaMat[i][k] = ss;
				}
			}
		}
	}
	
	public class EP {

		protected final String[] statespace;
		protected final String[] allele;
		protected final double bfrac;
		protected final double[][] probsMat;
		protected boolean logspace;
		protected double[][] count;

		protected final Map<String, Integer> ssMap;
		protected final Map<String, Integer> esMap;

		public EP(String[] statespace, 
				String[] allele, 
				double bfrac) {
			// TODO Auto-generated constructor stub
			this.statespace = statespace;
			this.allele = allele;
			this.bfrac = bfrac;
			int _a_ = this.statespace.length;
			int _b_ = this.allele.length;
			this.probsMat = new double[_a_][_b_];
			this.count = new double[_a_][_b_];
			this.prior();
			this.ssMap = new HashMap<String, Integer>();
			for(int i=0; i<_a_; i++)
				this.ssMap.put(statespace[i], i);
			this.esMap = new HashMap<String, Integer>();
			for(int i=0; i<_b_; i++)
				this.esMap.put(allele[i], i);
			this.logspace = false;
		}

		public EP(String[] statespace,
				String[] allele, 
				double bfrac, 
				String emissProbs) {
			// TODO Auto-generated constructor stub
			this.statespace = statespace;
			this.allele = allele;
			this.bfrac = bfrac;
			int _a_ = this.statespace.length;
			int _b_ = this.allele.length;
			this.probsMat = new double[_a_][_b_];
			this.count = new double[_a_][_b_];
			this.prior(emissProbs);
			this.ssMap = new HashMap<String, Integer>();
			for(int i=0; i<statespace.length; i++)
				this.ssMap.put(statespace[i], i);
			this.esMap = new HashMap<String, Integer>();
			for(int i=0; i<allele.length; i++)
				this.esMap.put(allele[i], i);
			this.logspace = false;
			this.count = new double[this.statespace.length][this.allele.length];
		}

		private void prior(String emissProbs) {
			// TODO Auto-generated method stub
			Pattern p = Pattern.compile("\\{(.*?)\\}");
			Matcher m = p.matcher(emissProbs);
			
			int _a_ = this.statespace.length;
			int _b_ = this.allele.length;
			String[] prob_str;
			double[] probs;
			
			for(int a=0; a<_a_; a++) {
				probs = this.probsMat[a];
				m.find();
				prob_str = m.group(1).split(",|;");
				for(int b=0; b<_b_; b++) 
					probs[b] = Double.parseDouble(prob_str[b*2+1]);
			}
		}

		protected void prior() {
			// TODO Auto-generated method stub
			Arrays.fill(this.probsMat[0], 0);
			for(int i=1; i<this.probsMat.length; i++) {
				//Dirichlet diri = new Dirichlet(distz(this.allele.length), 
				//		Constants._mu_theta_e);
				Dirichlet diri = new Dirichlet(new double[]{1-bfrac, bfrac}, 
						Constants._mu_theta_e);
				this.probsMat[i] = ArrayUtils.toPrimitive(diri.sample());
			}
		}

		public void getCounts(double[] hittingProb) {
			int _i_ = count.length;
			double sum = 0;
			for(int i=1; i<_i_; i++) {
				hittingProb[i] = StatUtils.sum(count[i]);
				sum += hittingProb[i];
			}
			if(sum==0) 
				for(int i=1; i<_i_; i++) 
					hittingProb[i] = 1.0/(_i_-1);
			else
				for(int i=1; i<_i_; i++) 
					hittingProb[i] /= sum;
		}
		
		protected void posterior() {
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

		public double[][] probs() {
			// TODO Auto-generated method stub
			return this.probsMat;
		}
	}
	
	protected class Viterbi {
		protected double[][] v;
		protected int m;
		protected int[][] trace;
		protected String[] statespace;
		protected double[] logscale;
		protected String[][] path_str;
		protected int[][] path;
		protected double[] probability;
		protected int[] sortv;
		
		public Viterbi(int _m_, int _n_, String[] statespace) {
			// TODO Auto-generated constructor stub
			this.v = new double[_m_][_n_];
			this.sortv = new int[_n_];
			this.m = _m_;
			this.trace = new int[_m_-1][_n_];
			this.statespace = statespace;
			this.logscale = new double[_m_];
			this.path_str = new String[_n_][_m_-1];
			this.path = new int[_n_][_m_-1];
			this.probability = new double[_n_];
		}

		public double probability(final int i) {
			// TODO Auto-generated method stub
			return probability[i];
		}
		
		public double probability() {
			// TODO Auto-generated method stub
			return probability(0);
		}

		public void finalise() {
			this.sortv = IntStream.range(0, v[m-1].length)
					.boxed().sorted((i, j) -> Double.compare(v[m-1][j], v[m-1][i]))
					.mapToInt(ele -> ele).toArray();
			for(int i=0; i<probability.length; i++)
				this.probability[i] = Math.log(this.v[m-1][sortv[i]])+this.logscale[m-1];
		}
		
		protected void trace(int v) {
			int tr = sortv[v];
			this.path[v][m-2] = tr;
			this.path_str[v][m-2] = this.statespace[tr];
			for(int i=m-3; i>=0; i--) {
				tr = trace[i+1][tr];
				this.path[v][i] = tr;
				this.path_str[v][i] = maxMatch(this.statespace[tr], 
						this.path_str[v][i+1]);
			}
			return;
		}
		
		private String maxMatch(String str1, String str2) {
			// TODO Auto-generated method stub
			final String[] s1 = str1.split("_");
			final Set<String> s2 = Arrays.stream(str2.split("_")).collect(Collectors.toSet());;
			final String[] s = new String[s1.length];
			int i = 0;
			for(final String c1 : s1) {
				if(s2.contains(c1)) {
					s[i++] = c1;
					s2.remove(c1);
				}
			}
			for(final String c2 : s2) s[i++] = c2;
			return StringUtils.join(s, '_');
		}

		protected void trace() {
			this.trace(0);
		}
		
		protected void traceAll() {
			for(int i=0; i<this.path.length; i++)
				this.trace(i);
		}
		
		protected void scale(final int i) {
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
	
	protected class FB { /** forward/backward algorithm object */
		protected double probability;
		protected double[][] probsMat;
		protected boolean logspace;
		protected double[] logscale;
		protected final boolean backward;

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

	public int hs() {
		// TODO Auto-generated method stub
		return this.hs.length;
	}

	public EP ep(int i) {
		// TODO Auto-generated method stub
		return this.emissProbs[i];
	}
	
	public TP tp(int i) {
		return this.transProbs[i];
	}

	public int noSnps() {
		// TODO Auto-generated method stub
		return this.M;
	}

	public OB ob(int i) {
		// TODO Auto-generated method stub
		return this.obspace[i];
	}

	public DataEntry de() {
		// TODO Auto-generated method stub
		return this.de;
	}
	
	public String[] state_str() {
		// TODO Auto-generated method stub
		return this.str_statespace;
	}
}
