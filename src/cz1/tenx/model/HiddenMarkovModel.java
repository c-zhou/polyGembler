package cz1.tenx.model;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.IntStream;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.CollectionUtils;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.inference.TestUtils;
import org.apache.log4j.Logger;

import cz1.algebra.matrix.SparseMatrix;
import cz1.graph.cluster.KMedoids;
import cz1.graph.cluster.MarkovClustering;
import cz1.graph.cluster.SpectralClustering;
import cz1.util.Algebra;
import cz1.util.Constants;
import cz1.util.Utils;
import cz1.util.Dirichlet;

public class HiddenMarkovModel {
	protected final static Logger myLogger = Logger.getLogger(HiddenMarkovModel.class);
	
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

	protected boolean clustering = true;
	protected boolean simulated_annealing = false;

	protected int[] haps = null;
	
	protected boolean  debug = false;
	protected boolean ddebug = false;
	
	public HiddenMarkovModel(String vcf_file,
			String dat_file,
			String rangeChr) {
		this(vcf_file, dat_file, rangeChr, true, false);
	}

	public HiddenMarkovModel(String vcf_file,
			String dat_file,
			String rangeChr,
			boolean use_clus,
			boolean use_sa) {
		this(vcf_file, dat_file, rangeChr, Integer.MIN_VALUE, Integer.MAX_VALUE, use_clus, use_sa);
	}

	public HiddenMarkovModel(String vcf_file,
			String dat_file,
			String rangeChr,
			int rangeLowerBound,
			int rangeUpperBound) {
		this(vcf_file, dat_file, rangeChr, rangeLowerBound, rangeUpperBound, true, false);
	}

	public HiddenMarkovModel(String vcf_file,
			String dat_file,
			String rangeChr,
			int rangeLowerBound,
			int rangeUpperBound,
			boolean use_clus,
			boolean use_sa) {
		this(vcf_file, dat_file, rangeChr, rangeLowerBound, rangeUpperBound, use_clus, use_sa, false, false);
	}
	
	public HiddenMarkovModel(String vcf_file,
			String dat_file,
			String rangeChr,
			int rangeLowerBound,
			int rangeUpperBound,
			boolean use_clus,
			boolean use_sa,
			boolean debug,
			boolean ddebug) {
		this.rangeChr = rangeChr;
		this.rangeLowerBound = rangeLowerBound;
		this.rangeUpperBound = rangeUpperBound;
		this.clustering = use_clus;
		this.simulated_annealing = use_sa;
		this.debug = debug;
		this.ddebug = ddebug;
		this.setDataEntryFile(vcf_file, dat_file);
		if(this.isNullModel()) return;
		if(this.clustering) this.makeMCL();
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

	protected Variant[] variants = null;
	protected double[] bfrac = null;
	protected final static int minD = 2;
	protected final static double minP = 0.00;
	protected DataEntry[] dp;
	// <K,V>=<MarkerIndex,DP_SET>
	protected Map<Integer, Set<Integer>> dpCrossRef = new HashMap<Integer, Set<Integer>>();
	
	private void setDataEntryFile(String vcf_file, String dat_file) {
		// TODO Auto-generated method stub
		try {
			// read VCF file
			BufferedReader br_vcf = Utils.getBufferedReader(vcf_file);
			String line;
			String[] s;
			int position;
			final List<Variant> variant_list = new ArrayList<Variant>();
			int var_start = Integer.MAX_VALUE, var_end = Integer.MIN_VALUE;
			int ind = 0;
			final BidiMap<Integer, Integer> idxmap = new DualHashBidiMap<Integer, Integer>();
			
			while( (line=br_vcf.readLine())!=null ){
				if(line.startsWith("#")) continue;
				++ind;
				s = line.split("\\s+");
				if(!s[0].equals(rangeChr)) continue;
				position = Integer.parseInt(s[1]);
				if(position<rangeLowerBound) continue;
				if(position>rangeUpperBound) break;
				if(var_start>ind) var_start = ind;
				if(var_end<ind)   var_end   = ind;
				idxmap.put(ind, variant_list.size());
				variant_list.add(new Variant(position, s[3], s[4]));
			}
			br_vcf.close();

			// read DAT file
			BufferedReader br_dat = Utils.getBufferedReader(dat_file);
			final List<DataEntry> dp_list = new ArrayList<DataEntry>();
			
			final List<Integer> index = new ArrayList<Integer>();
			final List<Integer> allele = new ArrayList<Integer>();
			int[] index_arr, allele_arr;
			int n, starti;
			while( (line=br_dat.readLine())!=null ) {
				s = line.split("\\s+");
				//if(!s[1].split(":")[0].equals(rangeChr)) continue;
				n = s.length;
				if(Integer.parseInt(s[n-3])+s[n-2].length()-1<var_start) 
					continue;
				if(Integer.parseInt(s[5])>var_end) break;
				index.clear();
				allele.clear();
				for(int i=5; i<n-2; i+=2) {
					starti = Integer.parseInt(s[i]);
					for(char c : s[i+1].toCharArray()) {
						if(starti>=var_start&&starti<=var_end) {
							index.add(idxmap.get(starti));
							allele.add(c-'0');
						}
						++starti;
					}
				}
				if(index.size()<2) continue;
				index_arr  = ArrayUtils.toPrimitive( index.toArray(new Integer[ index.size()]));
				allele_arr = ArrayUtils.toPrimitive(allele.toArray(new Integer[allele.size()]));
				dp_list.add(new DataEntry(index_arr, allele_arr));
			}
			br_dat.close();
			
			// bfrac
			int M = variant_list.size();
			int N = dp_list.size();
			long[][] depth = new long[M][2];
			double[] probs;
			DataEntry dp1;
			for(int i=0; i<N; i++) {
				dp1 = dp_list.get(i);
				for(int j : dp1.index.values()) {
					probs = dp1.probs[dp1.index.getKey(j)];
					depth[j][probs[0]>probs[1]?0:1]++;
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
			
			final List<Variant> variants = new ArrayList<Variant>();
			final List<Double> bfrac = new ArrayList<Double>();
			final Set<Integer> noTrainLoci = new HashSet<Integer>();
			final BidiMap<Integer, Integer> oldnewidx = new DualHashBidiMap<Integer, Integer>();
			for(int i=0; i<M; i++) {
				if(depth[i][0]==0||depth[i][1]==0||depth[i][0]+depth[i][1]<minD) {
					noTrainLoci.add(i);
					continue;
				}
				oldnewidx.put(i, variants.size());
				variants.add(variant_list.get(i));
				for(int z=0; z<ploidy-1; z++) 
					chisq_p[z] = TestUtils.chiSquareTest(config[z], depth[i]);
				maxPval = Algebra.maxIndex(chisq_p);
				bfrac.add(config[maxPval][1]);
			}
			
			// need to update dp
			// remove loci not train
			final List<DataEntry> dp = new ArrayList<DataEntry>();
			int z, w, k;
			Set<Integer> vals;
			for(int i=0; i<N; i++) {
				dp1 = dp_list.get(i);
				vals = dp1.index.values();
				vals.removeAll(noTrainLoci);
				z = vals.size();
				if(z>1) {
					// add to final dp
					final int[] indexz  = new int[z];
					final int[] allelez = new int[z];
					k = 0;
					for(int j : vals) {
						w = dp1.index.getKey(j);
						indexz[k]  = oldnewidx.get(j);
						allelez[k] = dp1.probs[w][0]>dp1.probs[w][1]?0:1;
						++k;
					}
					dp.add(new DataEntry(indexz, allelez));
				}
			}

			this.variants = new Variant[variants.size()];
			variants.toArray(this.variants);
			this.M = this.variants.length;
			this.dp = dp.toArray(new DataEntry[dp.size()]);
			this.N = this.dp.length;
			for(int i=0; i<this.M; i++) { 
				dpCrossRef.put(i, new HashSet<Integer>());
			}
			for(int i=0; i<this.N; i++) {
				final Set<Integer> indexSet = this.dp[i].index.values();
				for(int j : indexSet) dpCrossRef.get(j).add(i);
			}
			
			this.bfrac = ArrayUtils.toPrimitive(bfrac.toArray(new Double[this.M]));
			this.haps   = new int[N];
			
			myLogger.info("Data entry loaded. #Loci: "+this.M+"/"+M);
			if(debug) {	
				final StringBuilder os = new StringBuilder();
				for(int i=0; i<this.M; i++) {
					os.setLength(0);
					os.append(Utils.fixedLengthPaddingString(""+this.variants[i].position, 8));
					os.append(": ");
					os.append(Utils.fixedLengthPaddingString(""+depth[oldnewidx.getKey(i)][0], 3));
					os.append(",");
					os.append(Utils.fixedLengthPaddingString(""+depth[oldnewidx.getKey(i)][1], 3));
					os.append("\t");
					os.append(String.format("%.3f", this.bfrac[i]));
					myLogger.info(os.toString());
				}
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void makeBWT() {
		// TODO Auto-generated method stub
		this.emissProbs = new EP[M];
		for(int i=0; i<M; i++)
			emissProbs[i] = new EP(i);

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

	/**
	for(int i=0; i<this.M; i++) { 
		dpCrossRef.put(i, new HashSet<Integer>());
	}
	for(int i=0; i<this.N; i++) {
		final Set<Integer> indexSet = this.dp[i].index.values();
		for(int j : indexSet) dpCrossRef.get(j).add(i);
	}
	 **/

	private final Map<Integer, Integer> initMolClus = new HashMap<Integer, Integer>();
	private final static int maxIter = 100;
	
	// this implements a MCL to detect molecule clusters
	private void makeMCL() {

		List<Set<Integer>> clusters = new ArrayList<Set<Integer>>();
		for(int i=0; i<N; i++) {
			Set<Integer> clus = new HashSet<Integer>();
			clus.add(i);
			clusters.add(clus);
		}

		final List<Map<Integer, Integer>> hapcombs = new ArrayList<Map<Integer,Integer>>();
		int N;
		SparseMatrix dpAdj = null;
		boolean recalc = true;
		
		for(int iter=0; iter<=maxIter; iter++) {
			N = clusters.size();

			hapcombs.clear();
			for(int i=0; i<N; i++) 
				hapcombs.add(getHapFromRead(clusters.get(i)));
			
			if(ddebug) { 
				myLogger.info("####haplotypes ");
				for(int i=0; i<N; i++) {
					printHaps(hapcombs.get(i));
				}
			}
			
			dpAdj = getSimularityMatrix(hapcombs);
			
			/***
			try {
				BufferedWriter bw = Utils.getBufferedWriter("c:/users/chenxi.zhou/desktop/aa.txt");
				for(int i=0; i<N; i++) 
					for(int j=i; j<N; j++) 
						if(dpAdj.get(i, j)>0) 
							bw.write(i+"\t"+j+"\t"+dpAdj.get(i, j)+"\n");
				bw.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			try {
				BufferedWriter bw = Utils.getBufferedWriter("c:/users/chenxi.zhou/desktop/bb.txt");
				for(int i=0; i<N; i++) {
					bw.write(""+dpAdj.get(i, 0));
					for(int j=1; j<N; j++) 
						bw.write("\t"+dpAdj.get(i, j));
					bw.write("\n");
				}
				bw.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			**/
			
			MarkovClustering mcl = new MarkovClustering(dpAdj, 1e-6, 2, 1, 1e-6);
			mcl.run();

			myLogger.info("MCL converged. "+mcl.progress());

			// now interpret MCL clusters
			mcl.parse();

			final List<Set<Integer>> clusterz = mcl.getCluster();
			myLogger.info("MCL parsed.");

			// only keep the #ploidy largest clusters
			if(clusterz.size()<Constants._ploidy_H) {
				myLogger.warn("initial graph based assignment failed!!! Use random initialisation. "+
						this.rangeChr+":"+this.rangeLowerBound+"-"+this.rangeUpperBound);
				return;
			}
			final List<Set<Integer>> clusters2 = new ArrayList<Set<Integer>>();
			for(int i=0; i<clusterz.size(); i++) {
				final Set<Integer> clus = clusterz.get(i);
				final Set<Integer> elements = new HashSet<Integer>();
				for(final int j : clus) {
					elements.addAll(clusters.get(j));
				}
				clusters2.add(elements);
			}
			
			int n = clusters2.size();
			myLogger.info("Error correction");
			errorCorrection(clusters2);
			myLogger.info("Error correction done. "+(clusters2.size()-n)+" singletons.");
			
			if(equals(clusters, clusters2)) {
				recalc = false;
				break;
			}
			
			clusters = clusters2;
			myLogger.info("####MCL round "+iter+", #clusters="+clusters.size());
		}

		if(recalc) {
			N = clusters.size();
			hapcombs.clear();
			for(int i=0; i<N; i++) 
				hapcombs.add(getHapFromRead(clusters.get(i)));
		}
		
		// loci covered by each cluster
		int n = clusters.size();
		final boolean[][] clusCov = new boolean[n][M];
		final int[] cluSize = new int[n];
		cluStats(clusters, clusCov, cluSize, n);
		int[] sortedClus = IntStream.range(0, n)
                .boxed().sorted((i, j) -> Integer.compare(cluSize[j],cluSize[i]))
                .mapToInt(i -> i).toArray();
		
		// a greedy algorithm to categorize the 'super blocks'
		// first find seeds
		final Set<Integer> placed = new HashSet<Integer>();
		final List<Set<Integer>> haplos = new ArrayList<Set<Integer>>();
		for(int i=0; i<Constants._ploidy_H; i++) 
			haplos.add(new HashSet<Integer>());
		final double[][] dissimilarity = this.getDissimilarityMatrix(hapcombs);
		int maxz;
		double d, maxd;
		final int[] seeds = new int[Constants._ploidy_H];
		seeds[0] = sortedClus[0];
		placed.add(seeds[0]);
		
		for(int i=1; i<Constants._ploidy_H; i++) {
			maxz = -1;
			maxd = Double.NEGATIVE_INFINITY;
			for(int j=0; j<n; j++) {
				if(placed.contains(j)) continue;
				d = 0;
				for(int k=0; k<i; k++) {
					d += dissimilarity[j][seeds[k]];
				}
				if(d>maxd) {
					maxd = d;
					maxz = j;
				}
			}
			seeds[i]  = maxz;
			placed.add(maxz);
		}
		for(int i=0; i<Constants._ploidy_H; i++) haplos.get(i).addAll(clusters.get(seeds[i]));
		
		if(debug) {
			myLogger.info("####seeding haplotypes");
			for(final int i : seeds) printHaps(hapcombs.get(i));
		}
		
		// now we have seeds
		// calculate similarities between seeds and 'super blocks' 
		final double[][] similarity = this.getSimilarityMatrix(hapcombs, seeds, placed);
		int z, maxi;
		double s, maxs;
		int merged = 0;
		while(n>placed.size()) {
			maxs = Double.NEGATIVE_INFINITY;
			maxz = -1;
			maxi = -1;
			for(int i=0; i<n; i++) {
				if(placed.contains(i)) continue;
				z = Algebra.maxIndex(similarity[i]);
				s = similarity[i][z];
				if(maxs<s) {
					maxs = s;
					maxz = z;
					maxi = i;
				}
			}
			if(maxs<0) break;
			// now we place selected maxi
			mergeHaps(hapcombs.get(seeds[maxz]), hapcombs.get(maxi));
			haplos.get(maxz).addAll(clusters.get(maxi));
			placed.add(maxi);
			this.getSimilarityMatrix(similarity, hapcombs, seeds[maxz], maxz, placed);
		
			if(ddebug) {
				myLogger.info("####haplotypes merged "+(++merged));
				for(final int i : seeds) printHaps(hapcombs.get(i));
			}
		}
		
		if(debug) {
			myLogger.info("####greedy initialisation");
			for(final int i : seeds) printHaps(hapcombs.get(i));
		}
		
		for(int i=0; i<Constants._ploidy_H; i++) {
			Set<Integer> haplo = haplos.get(i);
			for(int j : haplo) initMolClus.put(j, i);
		}
	}

	private void errorCorrection(List<Set<Integer>> clusters) {
		// TODO Auto-generated method stub
		int N = clusters.size();
		for(int i=0; i<N; i++) {
			List<Set<Integer>> clus = errorCorrection(clusters.get(i));
			for(int j=0; j<clus.size(); j++)
				clusters.add(clus.get(j));
		}
	}

	private List<Set<Integer>> errorCorrection(Set<Integer> dat) {
		// TODO Auto-generated method stub
		final List<Set<Integer>> cluster = new ArrayList<Set<Integer>>();
		
		final Map<Integer, Integer> hap = new HashMap<Integer, Integer>();
		final Set<Integer> mismatch = new HashSet<Integer>();
		DataEntry dp1;
		Set<Integer> m1;
		int h;
		
		// detect mismatch
		for(final int i : dat) {
			dp1 = this.dp[i];
			m1 = dp1.index.values();
			for(final int j : m1) {
				h = dp1.probs[dp1.index.getKey(j)][0]==soften ? 1 : 0;
				if(hap.containsKey(j)) {
					if(hap.get(j)!=h) 
						mismatch.add(j);
				} else {
					hap.put(j, h);
				}
			}
		}
		
		if(mismatch.isEmpty()) return cluster;
		
		// mismatch detected
		final int[] alleleCount = new int[2];
		final int[] markerCount = new int[2];
		final int[] misCount    = new int[2];
		final Set<Integer> singleton = new HashSet<Integer>();
		final Set<Integer> bucket = new HashSet<Integer>();
		final List<Set<Integer>> a = new ArrayList<Set<Integer>>();
		for(int i=0; i<2; i++) a.add(new HashSet<Integer>());
		int n;
		for(final int i : mismatch) {
			final Set<Integer> indvi = dpCrossRef.get(i);
			Arrays.fill(alleleCount, 0);
			a.get(0).clear();
			a.get(1).clear();
			for(final int j : indvi) {
				if(!dat.contains(j)||singleton.contains(j)) continue;
				dp1 = this.dp[j];
				h = dp1.probs[dp1.index.getKey(i)][0]==soften ? 1 : 0;
				++alleleCount[h];
				a.get(h).add(j);
			}
			
			if(alleleCount[0]==0||alleleCount[1]==0) continue;
			
			if(alleleCount[0]>alleleCount[1]) {
				singleton.addAll(a.get(1));
			} else if(alleleCount[0]<alleleCount[1]) {
				singleton.addAll(a.get(0));
			} else {
				for(int k=0; k<2; k++) {
					bucket.clear();
					for(final int j : a.get(k)) 
						bucket.addAll(this.dp[j].index.values());
					n = bucket.size();
					bucket.removeAll(mismatch);
					markerCount[k] = bucket.size();
					misCount[k]    = n-markerCount[k];
				}
				
				if(markerCount[0]>markerCount[1]) {
					singleton.addAll(a.get(1));
				} else if(markerCount[0]<markerCount[1]) {
					singleton.addAll(a.get(0));
				} else {
					if(misCount[0]>misCount[1]) {
						singleton.addAll(a.get(0));
					} else {
						singleton.addAll(a.get(1));
					}
				}
			}
		}
		
		dat.removeAll(singleton);
		for(final int i : singleton) {
			final Set<Integer> clus = new HashSet<Integer>();
			clus.add(i);
			cluster.add(clus);
		}
		return cluster;
	}

	private double[][] getDissimilarityMatrix(List<Map<Integer, Integer>> hapcombs) {
		// TODO Auto-generated method stub
		Map<Integer, Integer> m1, m2;
		Collection<Integer> comm;
		int mismatch;
		int N = hapcombs.size();
		double z;
		double[][] dpAdj = new double[N][N];
		for(int i=0; i<N; i++) {
			m1 = hapcombs.get(i);
			for(int j=i+1; j<N; j++) {
				m2 = hapcombs.get(j);
				comm = CollectionUtils.intersection(m1.keySet(), m2.keySet());
				mismatch = 0;
				for(int k : comm) 
					if(m1.get(k)!=m2.get(k))
						++mismatch;
				if(mismatch==0) 
					z = Double.NEGATIVE_INFINITY;
				else z = Math.log((double)mismatch);
				dpAdj[i][j] = z;
				dpAdj[j][i] = z;
			}
		}
		myLogger.info("Dissimulatirty matrix construction done.");
		return dpAdj;
	}

	private void getSimilarityMatrix(double[][] similarity, 
			List<Map<Integer, Integer>> hapcombs, 
			final int seed, 
			final int seedx,
			final Set<Integer> redundas) {
		// TODO Auto-generated method stub
		Map<Integer, Integer> m1, m2;
		int N = hapcombs.size();
		m2 = hapcombs.get(seed);
		for(int i=0; i<N; i++) {
			if(redundas.contains(i)) continue;
			m1 = hapcombs.get(i);
			similarity[i][seedx] = this.getHaplotypeSimularity(m1, m2);
		}
		return;
	}
	
	private double[][] getSimilarityMatrix(List<Map<Integer, Integer>> hapcombs, int[] seeds, final Set<Integer> redundas) {
		// TODO Auto-generated method stub
		Map<Integer, Integer> m1, m2;
		int N = hapcombs.size();
		int S = seeds.length;
		double[][] dpAdj = new double[N][S];
		for(int i=0; i<N; i++) {
			if(redundas.contains(i)) continue;
			m1 = hapcombs.get(i);
			for(int j=0; j<S; j++) {
				m2 = hapcombs.get(seeds[j]);
				dpAdj[i][j] = this.getHaplotypeSimularity(m1, m2);
			}
		}
		return dpAdj;
	}
	
	private SparseMatrix getSimularityMatrix(List<Map<Integer, Integer>> hapcombs) {
		Map<Integer, Integer> m1, m2;
		int N = hapcombs.size();
		SparseMatrix dpAdj = new SparseMatrix(N, N);
		double z;
		for(int i=0; i<N; i++) {
			m1 = hapcombs.get(i);
			for(int j=i+1; j<N; j++) {
				m2 = hapcombs.get(j);
				z = this.getHaplotypeSimularity(m1, m2);
				if(z>0) {
					dpAdj.set(i, j, z);
					dpAdj.set(j, i, z);
				}
			}
		}
		myLogger.info("MCL sparse matrix construction done.");
		return dpAdj;
	}
	
	private final static double penalty = 10.0;
	private double getHaplotypeSimularity(final Map<Integer, Integer> m1, final Map<Integer, Integer> m2) {
		// TODO Auto-generated method stub
		Collection<Integer> comm = CollectionUtils.intersection(m1.keySet(), m2.keySet());
		int match = 0, mismatch = 0;
		for(final int i : comm) {
			if(m1.get(i)==9 || 
					m2.get(i)==9)
				continue;
			if(m1.get(i)==m2.get(i)) 
				++match;
			else ++mismatch;
		}
		return match-penalty*mismatch;
	}
	
	private void cluStats(final List<Set<Integer>> clusters,
			final boolean[][] clusCov, 
			final int[] cluSize, 
			final int n) {
		// TODO Auto-generated method stub
		DataEntry dp1;
		for(int i=0; i<n; i++) {
			for(final int j : clusters.get(i)) {
				dp1 = dp[j];
				for(final int k : dp1.index.values())
					clusCov[i][k] = true;
			}
			for(int j=0; j<M; j++) {
				if(clusCov[i][j]) ++cluSize[i];
			}

			if(ddebug) {
				StringBuilder os = new StringBuilder();
				os.append("####cluster ");
				os.append(String.format("%1$3s", i));
				os.append(":");
				for(int j=0; j<M; j++) {
					os.append(" ");
					os.append(clusCov[i][j]?1:0);
				}
				myLogger.info(os.toString());
			}
		}
	}

	private boolean equals(List<Set<Integer>> clusters, List<Set<Integer>> clusters2) {
		// TODO Auto-generated method stub
		if(clusters.size()!=clusters2.size()) return false;
		for(int i=0; i<clusters.size(); i++) {
			if(!clusters.get(i).containsAll(clusters2.get(i)))
				return false;
		}
		return true;
	}

	private void mergeHaps(Map<Integer, Integer> dat, Map<Integer, Integer> dat2) {
		// TODO Auto-generated method stub
		
		int p, h, h2;
		for(final Map.Entry<Integer, Integer> entry : dat2.entrySet()) {
			p = entry.getKey();
			h2 = entry.getValue();
			if(dat.containsKey(p)) {
				h = dat.get(p);
				if(h!=h2) 
					dat.put(p, 9);
			} else {
				dat.put(p, h2);
			}
		}
	}
	
	private Map<Integer, Integer> getHapFromRead(Set<Integer> dat) {
		// TODO Auto-generated method stub

		final Map<Integer, Integer> hap = new HashMap<Integer, Integer>();
		final Set<Integer> mismatch = new HashSet<Integer>();
		DataEntry dp1;
		Set<Integer> m1;
		int h;
		for(final int i : dat) {
			dp1 = this.dp[i];
			m1 = dp1.index.values();
			for(final int j : m1) {
				h = dp1.probs[dp1.index.getKey(j)][0]==soften ? 1 : 0;
				if(hap.containsKey(j)) {
					if(hap.get(j)!=h) 
						mismatch.add(j);
				} else {
					hap.put(j, h);
				}
			}
		}
		for(final int i : mismatch) hap.put(i, 9);	
	
		return hap;
	}
	
	private void printHaps(Map<Integer, Integer> hap) {
		StringBuilder os = new StringBuilder();
		int mismatch = 0;
		for(final int i : hap.values()) mismatch += i==9?1:0;
		for(int i=0; i<M; i++) 
			os.append(hap.containsKey(i)?hap.get(i):"-");
		os.append("    @");
		os.append(mismatch);
		os.append(" mismatches");
		myLogger.info(os.toString());
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
		}
		return;
	}

	private void makeForward() {
		// TODO Auto-generated method stub
		this.makeForward(dp, forward);
	}

	private void makeForward(DataEntry[] dp, FB[] forward) {
		// TODO Auto-generated method stub
		DataEntry dp1;
		for(int i=0; i<N; i++) {
			dp1 = dp[i];
			int M = dp1.probs.length;
			double[][] probs = dp1.weightedProbs; // MxP
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

	private void makeBackward() {
		// TODO Auto-generated method stub
		this.makeBackward(dp, backward);
	}

	private void makeBackward(DataEntry[] dp, FB[] backward) {
		// TODO Auto-generated method stub
		DataEntry dp1;
		for(int i=0; i<N; i++) {
			dp1 = dp[i];
			int M = dp1.probs.length;
			double[][] probs = dp1.weightedProbs; // MxP
			double[][] probsMat = backward[i].probsMat; // MxP

			Arrays.fill(probsMat[M-1], 1.);
			for(int j=M-2; j>=0; j--) {
				for(int k=0; k<Constants._ploidy_H; k++) 
					probsMat[j][k] = probsMat[j+1][k]*probs[j+1][k];	
				backward[i].scale(j);	
			}
			double ll = 0.;
			for(int j=0; j<Constants._ploidy_H; j++) 
				ll+=probsMat[0][j]*probs[0][j];
			backward[i].probability(ll);
		}
		return;
	}

	private void checkFW() {
		// TODO Auto-generated method stub
		if(iteration==0) return;

		myLogger.info(this.loglik()+"---"+this.loglik1());
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
						coeff = fw1.probsMat[z][k]*bw1.probsMat[z][k]*exp;
						for(int w=0; w<2; w++) 
							emiss_count[k][w] += coeff*dp1.probs[z][w];
					}
				}
			}

			this.emissProbs[i].posterior();
		}
	}

	private double temperature =  100;
	private double coolingRate = 0.01;

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


	private void findHap() {
		// TODO Auto-generated method stub
		DataEntry dp1;
		double[][] probs1, emiss1;
		int lc;
		double[] ll = new double[Constants._ploidy_H];
		for(int i=0; i<N; i++) {
			dp1 = this.dp[i];
			probs1 = dp1.probs;
			Arrays.fill(ll, 0);
			for(Map.Entry<Integer, Integer> ent : dp1.index.entrySet()) {
				lc = ent.getKey();
				emiss1 = this.emissProbs[ent.getValue()].probsMat;
				for(int j=0; j<Constants._ploidy_H; j++) {
					ll[j] += Math.log(probs1[lc][0]*emiss1[j][0]+
							probs1[lc][1]*emiss1[j][1]);
				}
			}
			this.haps[i] = Algebra.maxIndex(ll);
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

	private double shinkage = 0.5d;

	public class EP {
		protected final double baf; // B-allele frequency
		protected double[][] probsMat;
		protected double[][] count;
		protected boolean logspace;

		public EP(double[][] probsMat,
				double bfrac,
				boolean logspace) {
			// TODO Auto-generated constructor stub
			this.probsMat = probsMat;
			this.logspace = logspace;
			this.baf = bfrac;
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
			int pivot = (int) Math.round(Constants._ploidy_H*baf);

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
			return new EP(this.baf);
		}

		public EP(double bfrac) {
			// TODO Auto-generated constructor stub
			this.baf = bfrac;
			this.probsMat = new double[Constants._ploidy_H][2];
			this.logspace = false;
			this.count = new double[Constants._ploidy_H][2];
		}

		public EP(final int i, boolean prior) {
			// TODO Auto-generated constructor stub
			this.baf = bfrac[i];
			this.probsMat = new double[Constants._ploidy_H][2];
			if(prior) this.prior(i);
			this.logspace = false;
			this.count = new double[Constants._ploidy_H][2];
		}

		public EP(final int i) {
			// TODO Auto-generated constructor stub
			this(i, true);
		}

		protected void prior(final int mInd) {
			// TODO Auto-generated method stub
			Dirichlet diri = new Dirichlet(new double[]{1-baf, baf}, 
					Constants._mu_theta_e);
			if(initMolClus.isEmpty()) {
				// random initialisation using baf
				for(int i=0; i<Constants._ploidy_H; i++) {
					this.probsMat[i] = ArrayUtils.toPrimitive(diri.sample());
				}
			} else {
				// initialisation using initial assignment
				Set<Integer> indivs = dpCrossRef.get(mInd);
				for(int i=0; i<Constants._ploidy_H; i++) {
					final double[] a = new double[2];
					for(int j : indivs) {
						if(initMolClus.containsKey(j)&&initMolClus.get(j)==i) 
							++a[dp[j].probs[dp[j].index.getKey(mInd)][1]==soften?0:1];
					}
					double s = a[0]+a[1];
					if(s==0) {
						this.probsMat[i] = ArrayUtils.toPrimitive(diri.sample());
						continue;
					}
					
					/**
					double b = a[1]/s;
					if(b==0) b = 0.01;
					if(b==1) b = 0.99;
					this.probsMat[i] = ArrayUtils.toPrimitive(new Dirichlet(new double[]{1-b, b}, 
							Constants._mu_theta_e).sample());
					 **/
					if(a[0]>a[1]) this.probsMat[i] = new double[]{0.99,0.01};
					if(a[1]>a[0]) this.probsMat[i] = new double[]{0.01,0.99};
					if(a[0]==a[1]) this.probsMat[i] = ArrayUtils.toPrimitive(diri.sample());
				}
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
				count[i][0] = (1-baf)*Constants._pseudo_[1];
				count[i][1] = baf*Constants._pseudo_[1];
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

	protected final double soften = 1e-16;
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

		public int n() {
			// TODO Auto-generated method stub
			return this.probs.length;
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

		private int softenAlleleIndex(int markerIndex) {
			// TODO Auto-generated method stub
			return this.probs[this.index.getKey(markerIndex)][0]==soften?0:1;
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
		return this.M;
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
		// stats for haplotype coverage
		this.findHap();
		final boolean[][] hapCov = new boolean[Constants._ploidy_H][M];
		boolean[] cov;
		DataEntry dp1;
		for(int i=0; i<N; i++) {
			cov = hapCov[this.haps[i]];
			dp1 = this.dp[i];
			for(final int j : dp1.index.values())
				cov[j] = true;
		}
		
		try {
			BufferedWriter bw = Utils.getBufferedWriter(out);
			bw.write("##"+this.M+" "+this.loglik()+"\n");
			for(int i=0; i<M; i++) {
				double[][] probs = this.emissProbs[i].probsMat; // Px2
				Variant variant = variants[i];
				bw.write(this.rangeChr+"\t"+String.format("%1$12s", variant.position)+"\t"+
						variant.refAllele+"\t"+variant.altAllele);
				for(int j=0; j<Constants._ploidy_H; j++)
					bw.write("\t"+(hapCov[j][i]?(probs[j][0]<0.5?0:1):"-"));
				bw.write("\n");	
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
			double[][] probs = this.emissProbs[i].probsMat; // Px2
			Variant variant = variants[i];
			os.setLength(0);
			os.append(String.format("%1$12s", variant.position)+"\t"+
					variant.refAllele+"\t"+variant.altAllele);
			for(int j=0; j<Constants._ploidy_H; j++) 
				os.append("\t"+String.format("%.3f", probs[j][0]));
			myLogger.info(os.toString());
		}
	}

	public boolean isNullModel() {
		// TODO Auto-generated method stub
		return this.M == 0;
	}
}
