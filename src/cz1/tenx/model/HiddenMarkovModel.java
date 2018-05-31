package cz1.tenx.model;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.inference.TestUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import cz1.util.Algebra;
import cz1.util.Constants;
import cz1.util.Utils;
import cz1.util.Dirichlet;

public class HiddenMarkovModel {
	private final static Logger myLogger = LogManager.getLogger(HiddenMarkovModel.class.getName());
	
	protected final static Runtime runtime = Runtime.getRuntime();
	
	protected int iteration = 0;
	
	protected String[] hs = null;
	protected EP[] emissProbs = null;
	protected Viterbi vb[] = null;
	
	protected int N = -1;
	protected int M = -1;
	
	private FB[] forward, backward;
	
	protected String rangeChr = null;
	protected int rangeLowerBound = Integer.MIN_VALUE;
	protected int rangeUpperBound = Integer.MAX_VALUE;
	
	protected boolean simulated_annealing = true;
	
	public HiddenMarkovModel(String vcf_file,
			String dat_file,
			String rangeChr) {
		this(vcf_file, dat_file, rangeChr, true);
	}
	
	public HiddenMarkovModel(String vcf_file,
			String dat_file,
			String rangeChr,
			boolean use_sa) {
		this(vcf_file, dat_file, rangeChr, Integer.MIN_VALUE, Integer.MAX_VALUE, use_sa);
	}
	
	public HiddenMarkovModel(String vcf_file,
			String dat_file,
			String rangeChr,
			int rangeLowerBound,
			int rangeUpperBound) {
		this(vcf_file, dat_file, rangeChr, rangeLowerBound, rangeUpperBound, true);
	}
	
	public HiddenMarkovModel(String vcf_file,
			String dat_file,
			String rangeChr,
			int rangeLowerBound,
			int rangeUpperBound,
			boolean use_sa) {
		this.rangeChr = rangeChr;
		this.rangeLowerBound = rangeLowerBound;
		this.rangeUpperBound = rangeUpperBound;
		this.simulated_annealing = use_sa;
		this.setVariantDataFile(vcf_file);
		this.setDataEntryFile(dat_file);
		this.bfrac();
		this.makeBWT();
		if(this.simulated_annealing) this.makeSA();
	}
	
	public void train() {
		// TODO Auto-generated method stub
		if(this.isNullModel()) throw new RuntimeException("cannot train a null HMM model!!!");
		iteration++;
		myLogger.info("###################");
		myLogger.info("train: "+iteration);
		long[] tic = new long[10];
		int k=0;
		tic[k++] = System.nanoTime();
		this.updateDP();
		tic[k++] = System.nanoTime();
		myLogger.info("update dp "+(tic[k-1]-tic[k-2])+"ns");
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
		if(this.simulated_annealing) {
			this.SA();
			tic[k++] = System.nanoTime();
			myLogger.info("simulated annealing "+(tic[k-1]-tic[k-2])+"ns");
		}
		return;
	}
	
	// this implements a MCL to detect molecule clusters
	
	
	// <K,V> = <Index,Position>
	protected BidiMap<Integer, Integer> varidx = new DualHashBidiMap<Integer, Integer>();
	protected Variant[] variants = null;
	private void setVariantDataFile(String vcf_file) {
		// TODO Auto-generated method stub
		try {
			BufferedReader br = Utils.getBufferedReader(vcf_file);
			String line;
			String[] s;
			int position;
			final List<Variant> variants = new ArrayList<Variant>();
			while( (line=br.readLine())!=null ){
				if(line.startsWith("#")) continue;
				s = line.split("\\s+");
				if(!s[0].equals(rangeChr)) continue;
				position = Integer.parseInt(s[1]);
				if(position<rangeLowerBound) continue;
				if(position>rangeUpperBound) break;
				variants.add(new Variant(position, s[3], s[4]));
			}
			br.close();
			this.variants = new Variant[variants.size()];
			variants.toArray(this.variants);
			this.M = this.variants.length;
			for(int i=0; i<M; i++) varidx.put(i, this.variants[i].position);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	protected DataEntry[] dp;
	// <K,V>=<MarkerIndex,DP_SET>
	protected Map<Integer, Set<Integer>> dpCrossRef = new HashMap<Integer, Set<Integer>>();
	private void setDataEntryFile(String dat_file) {
		// TODO Auto-generated method stub
		try {
			BufferedReader br = Utils.getBufferedReader(dat_file);
			final List<DataEntry> dp = new ArrayList<DataEntry>();
			String line;
			String[] s, s2;
			int position;
			final List<Integer> index = new ArrayList<Integer>();
			final List<Integer> allele = new ArrayList<Integer>();
			int[] index_arr, allele_arr;
			int entryCount = 0;
			while( (line=br.readLine())!=null ) {
				s = line.split("\\s+");
				if(!s[0].equals(rangeChr)) continue;
				if(Integer.parseInt(s[2])<rangeLowerBound) continue;
				if(Integer.parseInt(s[1])>rangeUpperBound) break;
				s = s[3].split(";");
				index.clear();
				allele.clear();
				for(String entry : s) {
					s2 = entry.split(",");
					position = Integer.parseInt(s2[0]);
					if(position<rangeLowerBound||position>rangeUpperBound) continue;
					index.add(varidx.getKey(position));
					allele.add(Integer.parseInt(s2[1]));
				}
				if(index.size()<2) continue;
				index_arr  = ArrayUtils.toPrimitive( index.toArray(new Integer[ index.size()]));
				allele_arr = ArrayUtils.toPrimitive(allele.toArray(new Integer[allele.size()]));
				dp.add(new DataEntry(index_arr, allele_arr));
				entryCount++;
			}
			br.close();
			myLogger.info(entryCount+" data entry loaded.");
			this.dp = dp.toArray(new DataEntry[dp.size()]);
			this.N = this.dp.length;
			for(int i=0; i<this.M; i++) { 
				dpCrossRef.put(i, new HashSet<Integer>());
			}
			for(int i=0; i<this.N; i++) {
				final Set<Integer> indexSet = this.dp[i].index.values();
				for(int j : indexSet) dpCrossRef.get(j).add(i);
			}
			myLogger.info("Cross reference constuction finished.");
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	final Set<Integer> trainLoci = new HashSet<Integer>();
	protected double[] bfrac = null;
	protected final static int minD = 10;
	protected final static double minP = 0.05;
	
	protected void bfrac() {
		// TODO Auto-generated method stub
		long[][] depth = new long[M][2];
		double[] probs;
		for(int i=0; i<M; i++) {
			Set<Integer> crossRef = dpCrossRef.get(i);
			for(int j : crossRef) {
				probs = dp[j].probs[dp[j].index.getKey(i)];
				depth[i][probs[0]>probs[1]?0:1]++;
			}
		}
		int ploidy = Constants._ploidy_H;
		final double[][] config = new double[ploidy-1][2];
		for(int i=1; i<ploidy; i++) {
			config[i-1][0] = (double)i/ploidy;
			config[i-1][1] = 1.0-config[i-1][0];
		}
		final double[] chisq_p = new double[ploidy-1];
		int maxPval;
		bfrac = new double[M];
		for(int i=0; i<M; i++) {
			if(depth[i][0]+depth[i][1]<minD) continue;
			for(int z=0; z<ploidy-1; z++) 
				chisq_p[z] = TestUtils.chiSquareTest(config[z], depth[i]);
			maxPval = Algebra.maxIndex(chisq_p);
			if(chisq_p[maxPval]>=minP) {
				// a locus to train
				trainLoci.add(i);
				bfrac[i] = config[maxPval][1];
			}
		}
		myLogger.info("#Loci: "+trainLoci.size()+"/"+M);
		final StringBuilder os = new StringBuilder();
		for(int i=0; i<M; i++) {
			if(!trainLoci.contains(i)) continue;
			os.setLength(0);
			os.append(Utils.fixedLengthPaddingString(""+i, 8));
			os.append(": ");
			os.append(Utils.fixedLengthPaddingString(""+depth[i][0], 3));
			os.append(",");
			os.append(Utils.fixedLengthPaddingString(""+depth[i][1], 3));
			os.append("\t");
			os.append(String.format("%.3f", bfrac[i]));
			myLogger.info(os.toString());
		}
	}
	
	private void makeBWT() {
		// TODO Auto-generated method stub
		this.emissProbs = new EP[M];
		for(int i=0; i<M; i++)
			emissProbs[i] = new EP(bfrac[i]);

		this.forward = new FB[N];
		for(int i=0; i<N; i++) 
			this.forward[i] = new FB(false,
					this.dp[i].probs.length,
					Constants._ploidy_H);
		this.backward = new FB[N];
		for(int i=0; i<N; i++) 
			this.backward[i] = new FB(true,
					this.dp[i].probs.length,
					Constants._ploidy_H);
		this.vb = new Viterbi[N];
		for(int i=0; i<N; i++)
			this.vb[i] = new Viterbi(
					this.dp[i].probs.length,
					Constants._ploidy_H);
		return;
	}
	
	private DataEntry[] SAdp;
	private EP[] SAemissProbs;
	private FB[] SAforward;
	
	private void makeSA() {
		// TODO Auto-generated method stub
		// need a copy of dp
		SAdp = new DataEntry[N];
		for(int i=0; i<N;i++)
			SAdp[i] = dp[i].clone();
		// need a copy of emissProbs
		SAemissProbs = new EP[M];
		for(int i=0; i<M;i++)
			SAemissProbs[i] = emissProbs[i].clone();
		// need a copy of forward
		SAforward = new FB[N];
		for(int i=0; i<N;i++)
			SAforward[i] = forward[i].clone();
		return;
	}
	
	private void updateDP() {
		// TODO Auto-generated method stub
		this.updateDP(dp, emissProbs);
	}
	
	private void updateDP(DataEntry[] dp, EP[] emissProbs) {
		// TODO Auto-generated method stub
		for(int i=0; i<N; i++) {
			double[][] probs = dp[i].probs; // Mx2
			double[][] weightedProbs = dp[i].weightedProbs; // MxP
			BidiMap<Integer, Integer> index = dp[i].index;
			int M = probs.length;
			for(int j=0; j<M; j++) {
				double[][] emiss = 
						emissProbs[index.get(j)].probsMat; // Px2
				for(int k=0; k<Constants._ploidy_H; k++) {
					weightedProbs[j][k] = 
							probs[j][0]*emiss[k][0]+
							probs[j][1]*emiss[k][1];
				}
			}
			double[] ll = new double[Constants._ploidy_H];
			for(int j=0; j<Constants._ploidy_H; j++) {
				for(int k=0; k<M; k++) 
					ll[j] += Math.log(weightedProbs[k][j]);
			}
			Utils.print(ll);
			System.out.println("aaa");
		}
		return;
	}

	private void makeForward(DataEntry[] dp, FB[] forward) {
		// TODO Auto-generated method stub
		for(int i=0; i<N; i++) {
			int M = dp[i].probs.length;
			double[][] probs = dp[i].weightedProbs; // MxP
			double[][] probsMat = forward[i].probsMat; // MxP

			System.arraycopy(probs[0], 0, probsMat[0], 0, Constants._ploidy_H);
			for(int j=1; j<M; j++) {
				for(int k=0; k<Constants._ploidy_H; k++) 
					probsMat[j][k] = probsMat[j-1][k]*probs[j][k];	
				forward[i].scale(j);	
			}
			forward[i].probability(StatUtils.sum(probsMat[M-1]));
		}
		return;
	}
	
	private void makeForward() {
		// TODO Auto-generated method stub
		this.makeForward(dp, forward);
	}

	private void makeBackward(DataEntry[] dp, FB[] backward) {
		// TODO Auto-generated method stub
		for(int i=0; i<N; i++) {
			int M = dp[i].probs.length;
			double[][] probs = dp[i].weightedProbs; // MxP
			double[][] probsMat = backward[i].probsMat; // MxP

			System.arraycopy(probs[M-1], 0, probsMat[M-1], 0, Constants._ploidy_H);
			for(int j=M-2; j>=0; j--) {
				for(int k=0; k<Constants._ploidy_H; k++) 
					probsMat[j][k] = probsMat[j+1][k]*probs[j][k];	
				backward[i].scale(j);	
			}
			backward[i].probability(StatUtils.sum(probsMat[0]));
		}
		return;
	}
	
	private void makeBackward() {
		// TODO Auto-generated method stub
		this.makeBackward(dp, backward);
	}
	
	private void checkFW() {
		// TODO Auto-generated method stub
		if(iteration==0) return;

		System.out.println(this.loglik()+"---"+this.loglik1());
		for(int i=0; i<this.forward.length; i++) {
			double r = Math.abs(this.forward[i].probability-
					this.backward[i].probability);
			if(r>1e-6) myLogger.info(i+" | "+r+" --- FORWARD-BACKWARD PRECISION NOT RIGHT!!!");
		}
	}
	
	private void EM() {
		// TODO Auto-generated method stub
		
		double coeff, exp_c, exp;
		int z;
		FB fw1, bw1;
		DataEntry dp1;
		Set<Integer> crossRef;
		for(int i=0; i<M; i++) {
			
			if(!this.trainLoci.contains(i)) continue;
			
			double[][] emiss_count = this.emissProbs[i].pseudo(); // Px2
			crossRef = this.dpCrossRef.get(i);
			for(int j : crossRef) {
				fw1 = this.forward[j];
				bw1 = this.backward[j];
				dp1 = dp[j];
				z = dp1.index.getKey(i);
				
				exp_c = fw1.logscale[z]+
						bw1.logscale[z]-
						fw1.probability;
				
				if(exp_c>Constants.MAX_EXP_DOUBLE) { 
					// cannot calculate exponential directly
					// logarithm and then exponential 
					// time consuming
					for(int k=0; k<Constants._ploidy_H; k++) {
						coeff = fw1.probsMat[z][k]*bw1.probsMat[z][k];
						for(int w=0; w<2; w++) 
							emiss_count[k][w] += 
							Math.exp(Math.log(coeff*dp1.probs[z][w])+exp_c);
					}
				} else { 
					// exponential is safe
					// ideal way but dangerous
					exp = Math.exp(exp_c);
					for(int k=0; k<Constants._ploidy_H; k++) {
						coeff = fw1.probsMat[z][k]*bw1.probsMat[z][k];
						for(int w=0; w<2; w++) 
							emiss_count[k][w] += coeff*dp1.probs[z][w]*exp;
					}
				}
			}
			
			this.emissProbs[i].posterior();
		}
	}

	private static double temperature =  100;
	private static double coolingRate = 0.01;
	
	private void SA() {
		// TODO Auto-generated method stub
		// need to copy the emissProbs
		for(int i=0; i<M; i++) {
			double[][] emiss   = emissProbs[i].probsMat; // Px2
			double[][] SAemiss = SAemissProbs[i].probsMat; // Px2
			for(int j=0; j<Constants._ploidy_H; j++)
				System.arraycopy(emiss[j], 0, SAemiss[j], 0, 2);
			SAemissProbs[i].shrink();
		}
		
		// need to update SAdp
		this.updateDP(SAdp, SAemissProbs);
		// need to make SAforward
		this.makeForward(SAdp, SAforward);
		
		double llold = this.loglik(forward), 
				llnew = this.loglik(SAforward);
		if(llnew>llold||Math.exp((llnew-llold)/temperature)>Math.random()) {
			// accept SA local
			// need to replace emissProbs
			EP[] tmpEP = this.emissProbs;
			this.emissProbs = SAemissProbs;
			this.SAemissProbs = tmpEP;
			// need to replace dp
			DataEntry[] tmpDP = this.dp;
			this.dp = SAdp;
			this.SAdp = tmpDP;
			// need to replace forward
			FB[] tmpFB = this.forward;
			this.forward = SAforward;
			this.SAforward = tmpFB;
			// we run EM from here
			this.makeBackward();
			this.checkFW();
			this.EM();
			myLogger.info("SA local ACCEPTED at temperature "+temperature);
		} else {
			myLogger.info("SA local REJECTED at temperature "+temperature);
		}
		temperature *= 1-coolingRate;
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
	
	protected double makeViterbi() {
		// TODO Auto-generated method stub
		return -1;
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

	private static double shinkage = 0.5d;
	
	public class EP {
		protected final double bfrac;
		protected double[][] probsMat;
		protected double[][] count;
		protected boolean logspace;

		public EP(double[][] probsMat,
				double bfrac,
				boolean logspace) {
			// TODO Auto-generated constructor stub
			this.probsMat = probsMat;
			this.logspace = logspace;
			this.bfrac = bfrac;
			if( logspace ) this.setNormalspace();
			this.count = new double[Constants._ploidy_H][2];
		}
		
		public void shrink() {
			// TODO Auto-generated method stub
			for(int i=0; i<Constants._ploidy_H; i++) {
				if(this.probsMat[i][0]<0.5) {
					this.probsMat[i][0] = this.probsMat[i][0]*(1-shinkage);
				} else {
					this.probsMat[i][0] = this.probsMat[i][0]+shinkage*this.probsMat[i][1];
				}
				this.probsMat[i][1] = 1-this.probsMat[i][0];
			}
			shinkage *= 1-coolingRate;
		}
		
		public void shrink2() {
			// TODO Auto-generated method stub
			double[] probs = new double[Constants._ploidy_H];
			for(int i=0; i<Constants._ploidy_H; i++)
				probs[i] = this.probsMat[i][0];
			Arrays.sort(probs);
			int pivot = (int) Math.round(Constants._ploidy_H*bfrac);
			
			if(pivot==0&&probs[0]>0.5 ||
					pivot==Constants._ploidy_H&&
					probs[Constants._ploidy_H-1]<0.5 ||
					probs[pivot-1]<0.5&&probs[pivot]>0.5) 
				return;
			
			double pthres;
			if(pivot==0) 
				pthres = Double.NEGATIVE_INFINITY;
			else if(pivot==Constants._ploidy_H)
				pthres = Double.POSITIVE_INFINITY;
			else pthres = probs[pivot-1];
			
			for(int i=0; i<Constants._ploidy_H; i++) {
				if(this.probsMat[i][0]<=pthres) {
					this.probsMat[i][0] = this.probsMat[i][0]*(1-shinkage);
				} else {
					this.probsMat[i][0] = this.probsMat[i][0]+shinkage*this.probsMat[i][1];
				}
				this.probsMat[i][1] = 1-this.probsMat[i][0];
			}
		}
		
		public EP clone() {
			return new EP(this.bfrac, false);
		}
		
		public EP(double bfrac, boolean prior) {
			// TODO Auto-generated constructor stub
			this.bfrac = bfrac;
			this.probsMat = new double[Constants._ploidy_H][2];
			if(prior) this.prior();
			this.logspace = false;
			this.count = new double[Constants._ploidy_H][2];
		}

		public EP(double bfrac) {
			// TODO Auto-generated constructor stub
			this(bfrac, true);
		}

		protected void prior() {
			// TODO Auto-generated method stub
			for(int i=0; i<Constants._ploidy_H; i++) {
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

		protected void setProbsMat(double[][] probsMat) {
			this.probsMat = probsMat;
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
		protected double[] logscale;
		protected int[] path;
		protected double probability;

		public Viterbi(int _m_, int _n_) {
			// TODO Auto-generated constructor stub
			this.v = new double[_m_][_n_];
			this.m = _m_;
			this.trace = new int[_m_-1][_n_];
			this.logscale = new double[_m_];
			this.path = new int[_m_-1];
			this.probability = 0.0;
		}

		protected void trace() {
			this.probability = Math.log(StatUtils.max(v[m-1]))+
					this.logscale[m-1];
			int tr = Algebra.maxIndex(v[m-1]);
			this.path[m-2] = tr;
			for(int i=m-3; i>=0; i--) {
				tr = trace[i+1][tr];
				this.path[i] = tr;
			}
			return;
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

		public FB clone() {
			return new FB(this.backward, 
					this.probsMat.length, 
					this.probsMat[0].length);
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
	
	protected final double soften = 5e-3;
	protected class DataEntry { // data entry
		// <K,V>=<Index,MarkerIndex>
		final BidiMap<Integer, Integer> index;
		final double[][] probs; // Mx2 matrix
		final double[][] weightedProbs; // MxP matrix
				
		public DataEntry(final int[] index, final int[] allele) {
			this.probs = softProbs(allele);
			this.weightedProbs = new double[allele.length]
					[Constants._ploidy_H];
			this.index = new DualHashBidiMap<Integer, Integer>();
			this.fillIndexMap(index);
		}
		
		public DataEntry(final BidiMap<Integer, Integer> index,
				final double[][] probs) {
			this.index = index;
			this.probs = probs;
			this.weightedProbs = new double[index.size()]
					[Constants._ploidy_H];;
		}
		
		public DataEntry clone() {
			return new DataEntry(this.index, this.probs);
		}

		private void fillIndexMap(int[] index) {
			// TODO Auto-generated method stub
			for(int i=0; i<index.length; i++) 
				this.index.put(i, index[i]);
		}

		private double[][] softProbs(int[] allele) {
			// TODO Auto-generated method stub
			final double[][] probs = new double[allele.length][2]; 
			for(int i=0; i<allele.length; i++) {
				probs[i][  allele[i]] = 1-soften;
				probs[i][1-allele[i]] =   soften;
			}
			return probs;
		}
	}
	
	protected class Variant { // variant
		private final int position;
		private final String refAllele;
		private final String altAllele;
		
		public Variant(int position,
				String refAllele,
				String altAllele) {
			this.position  = position;
			this.refAllele = refAllele;
			this.altAllele = altAllele;
		}
		
		public int position() {
			return this.position;
		}
		
		public String refAllele() {
			return this.refAllele;
		}
		
		public String altAllele() {
			return this.altAllele;
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

	public int noSnps() {
		// TODO Auto-generated method stub
		return this.M;
	}
	
	public int noLoci() {
		// TODO Auto-generated method stub
		return this.trainLoci.size();
	}

	public DataEntry[] de() {
		// TODO Auto-generated method stub
		return this.dp;
	}
	
	public double loglik() {
		return this.loglik(forward);
	}
	
	public double loglik(FB[] forward) {
		if(iteration==0)
			return Double.NEGATIVE_INFINITY;
		else {
			double probability = 0;
			for(FB fw : forward) probability += fw.probability;
			return probability;
		}
	}

	public double loglik1() {
		return this.loglik1(backward);
	}
	
	public double loglik1(FB[] backward) {
		if(iteration==0)
			return Double.NEGATIVE_INFINITY;
		else {
			double probability = 0;
			for(FB bw : backward) probability += bw.probability;
			return probability;
		}
	}

	public void write(String out) {
		// TODO Auto-generated method stub
		try {
			BufferedWriter bw = Utils.getBufferedWriter(out);
			bw.write("##"+this.loglik()+"\n");
			for(int i=0; i<M; i++) {
				if(this.trainLoci.contains(i)) {
					double[][] probs = this.emissProbs[i].probsMat; // Px2
					Variant variant = variants[i];
					bw.write(String.format("%1$12s", variant.position)+"\t"+
							variant.refAllele+"\t"+variant.altAllele);
					for(int j=0; j<Constants._ploidy_H; j++) 
						bw.write("\t"+(probs[j][0]<0.5?0:1));
					bw.write("\n");	
				}
			}
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void print() {
		// TODO Auto-generated method stub
		myLogger.info("##"+this.loglik()+"\n");
		final StringBuilder os = new StringBuilder();
		for(int i=0; i<M; i++) {
			if(this.trainLoci.contains(i)) {
				double[][] probs = this.emissProbs[i].probsMat; // Px2
				Variant variant = variants[i];
				os.setLength(0);
				os.append(String.format("%1$12s", variant.position)+"\t"+
						variant.refAllele+"\t"+variant.altAllele);
				for(int j=0; j<Constants._ploidy_H; j++) 
					os.append("\t"+String.format("%.3f", probs[j][0]));
			}
			myLogger.info(os.toString());
		}
	}

	public boolean isNullModel() {
		// TODO Auto-generated method stub
		return this.trainLoci.isEmpty();
	}
}
