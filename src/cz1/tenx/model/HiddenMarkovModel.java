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
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.stream.IntStream;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.CollectionUtils;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.inference.TestUtils;
import org.apache.log4j.Logger;

import com.google.common.collect.Range;

import cz1.algebra.matrix.SparseMatrix;
import cz1.graph.cluster.KMedoids;
import cz1.graph.cluster.MarkovClustering;
import cz1.graph.cluster.SpectralClustering;
import cz1.tenx.model.HiddenMarkovModel2.DataEntry;
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

	protected int N = -1;
	protected int M = -1;

	protected String rangeChr = null;
	protected int rangeLowerBound = Integer.MIN_VALUE;
	protected int rangeUpperBound = Integer.MAX_VALUE;

	protected boolean clustering = true;
	protected boolean simulated_annealing = false;

	protected int[] haps = null;
	protected double[][] hap_assign = null;
	
	protected boolean  debug = false;
	protected boolean ddebug = false;
	
	protected double loglik = Double.NEGATIVE_INFINITY;
	
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
		this.makeEM();
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
		this.optimise();
		tic[k++] = System.nanoTime();
		myLogger.info("EM algorithm "+(tic[k-1]-tic[k-2])+"ns");
		myLogger.info("log likelihood "+this.loglik);	
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
			this.hap_assign = new double[N][Constants._ploidy_H];
			
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

	private void makeEM() {
		// TODO Auto-generated method stub
		this.emissProbs = new EP[M];
		for(int i=0; i<M; i++)
			emissProbs[i] = new EP(i);

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
		int N = 0;
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
		

		/**
		final List<Map<Integer, Integer>> hapcounts = new ArrayList<Map<Integer,Integer>>();
		for(int i=0; i<N; i++) hapcounts.add(getCountFromRead(clusters.get(i)));
		hapcombs.clear();
		for(int i=0; i<N; i++) hapcombs.add(getHapFromRead(clusters.get(i)));
	
		this.N = hapcombs.size();
		this.dp = new DataEntry[this.N];
		for(int i=0; i<this.N; i++) {
			final Map<Integer, Integer> hap   = hapcombs.get(i);
			final Map<Integer, Integer> count = hapcounts.get(i);
			final int n1 = count.size();
			final int[] indexz  = new int[n1];
			final int[] allelez = new int[n1];
			final int[] depth   = new int[n1];
			int j=0;
			for(Map.Entry<Integer, Integer> entry : hap.entrySet()) {
				indexz[j]  = entry.getKey();
				allelez[j] = entry.getValue();
				depth[j]   = count.get(indexz[j]);
				++j;
			}
			this.dp[i] = new DataEntry(indexz, allelez, depth);
		}
		for(int i=0; i<this.M; i++) { 
			dpCrossRef.put(i, new HashSet<Integer>());
		}
		for(int i=0; i<this.N; i++) {
			final Set<Integer> indexSet = this.dp[i].index.values();
			for(int j : indexSet) dpCrossRef.get(j).add(i);
		}
		this.haps   = new int[this.N];
		this.hap_assign = new double[this.N][Constants._ploidy_H];
		**/
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
				dpAdj[i][j] = mismatch;
				dpAdj[j][i] = mismatch;
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
	
	private Map<Integer, Integer> getCountFromRead(Set<Integer> dat) {
		// TODO Auto-generated method stub

		final Map<Integer, Integer> count = new HashMap<Integer, Integer>();
		DataEntry dp1;
		Set<Integer> m1;
		for(final int i : dat) {
			dp1 = this.dp[i];
			m1 = dp1.index.values();
			for(final int j : m1) {
				if(!count.containsKey(j))
					count.put(j, 0);
				count.put(j, count.get(j)+1);
			}
		}
		
		return count;
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

	private void optimise () {
		// TODO Auto-generated method stub
		final double[] loglik = new double[Constants._ploidy_H];
		int z;
		DataEntry dp1;
		Set<Integer> crossRef;
		double exp_c;
		
		this.findHap();
		for(int i=0; i<M; i++) {
			// exception
			double[][] emiss_count = this.emissProbs[i].pseudo(); // Px2
			crossRef = this.dpCrossRef.get(i);
			for(int j : crossRef) {
				dp1 = this.dp[j];
				z = dp1.index.getKey(i);
				
				System.arraycopy(this.hap_assign[j], 0, loglik, 0, Constants._ploidy_H);
				exp_c = StatUtils.max(loglik);
				for(int k=0; k<loglik.length; k++) 
					loglik[k] = Math.exp(loglik[k]-exp_c);
				Algebra.normalize(loglik);
				
				for(int k=0; k<loglik.length; k++) {
					emiss_count[k][0] += loglik[k]*dp1.depth[z][0];
					emiss_count[k][1] += loglik[k]*dp1.depth[z][1];
				}
			}
			
			// maximisation
			this.emissProbs[i].posterior();
		}
	}
	
	private void findHap() {
		// TODO Auto-generated method stub
		DataEntry dp1;
		double[][] probs1, emiss1;
		int lc;
		double[] ll;
		this.loglik = 0;
		for(int i=0; i<N; i++) {
			dp1 = this.dp[i];
			probs1 = dp1.probs;
			ll = hap_assign[i];
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
			this.loglik += ll[this.haps[i]];
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
			for(int i=0; i<_i_; i++) {
				hittingProb[i] = StatUtils.sum(count[i]);
				sum += hittingProb[i];
			}
			if(sum==0) 
				for(int i=0; i<_i_; i++) 
					hittingProb[i] = 1.0/_i_;
			else
				for(int i=0; i<_i_; i++) 
					hittingProb[i] /= sum;
		}

		protected void posterior() {
			// TODO Auto-generated method stub
			int _i_ = count.length,
					_j_ = count[0].length;
			for(int i=0; i<_i_; i++) {
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
			for(int i=0; i<_i_; i++) {
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

	protected final double soften = 1e-16;
	protected class DataEntry { // data entry
		// <K,V>=<Index,MarkerIndex>
		final BidiMap<Integer, Integer> index;
		final double[][] depth; // Mx2 matrix
		final double[][] probs; // Mx2 matrix
		
		public DataEntry(final int[] index, final int[] allele) {
			this.probs = softProbs(allele);
			this.depth = this.nullDepth();
			this.index = new DualHashBidiMap<Integer, Integer>();
			this.fillIndexMap(index);
		}
		
		public DataEntry(final int[] index, final int[] allele, final int[] d) {
			this.probs = softProbs(allele);
			this.depth = this.depth(d);
			this.index = new DualHashBidiMap<Integer, Integer>();
			this.fillIndexMap(index);
		}

		private double[][] nullDepth() {
			// TODO Auto-generated method stub
			final int M = probs.length;
			final int[] d = new int[M]; 
			Arrays.fill(d, 1);
			return this.depth(d);
		}

		private double[][] depth(int[] d) {
			// TODO Auto-generated method stub
			final int M = probs.length;
			final double[][] depth = new double[M][2];
			for(int i=0; i<M; i++) {
				depth[i][0] = probs[i][0]*d[i];
				depth[i][1] = probs[i][1]*d[i];
			}
			return depth;
		}

		public int n() {
			// TODO Auto-generated method stub
			return this.probs.length;
		}

		public DataEntry(final BidiMap<Integer, Integer> index,
				final double[][] probs) {
			this.index = index;
			this.probs = probs;
			this.depth = this.nullDepth();
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
		// TODO Auto-generated method stub
		return this.loglik;
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
		
		final List<Integer> breakpoints = this.breakpoints();
		
		try {
			BufferedWriter bw = Utils.getBufferedWriter(out);
			bw.write("##"+this.M+" "+this.loglik()+"\n");
			for(int b=1; b<breakpoints.size(); b++) {
				int start = breakpoints.get(b-1);
				int end   = breakpoints.get(b  );
				bw.write("BLOCK: "+this.rangeChr+"\t"+this.variants[start].position+"\t"+this.variants[end-1].position+"\n");
				for(int i=start; i<end; i++) {
					double[][] probs = this.emissProbs[i].probsMat; // Px2
					Variant variant = variants[i];
					bw.write(this.rangeChr+"\t"+String.format("%1$12s", variant.position)+"\t"+
							variant.refAllele+"\t"+variant.altAllele);
					for(int j=0; j<Constants._ploidy_H; j++)
						bw.write("\t"+(hapCov[j][i]?(probs[j][0]>0.5?0:1):"-"));
					bw.write("\n");	
				}
				bw.write("********\n");
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

	public List<Integer> breakpoints() {
		// TODO Auto-generated method stub
		// aims to detect block switch errors

		// here we will create a graph for the data entry
		// for adjacent markers
		// if no path connect them
		// we then break the haplotype at that point
		final List<Set<Integer>> adjmat = new ArrayList<Set<Integer>>();
		for(int i=0; i<M; i++) adjmat.add(new HashSet<Integer>());
		Set<Integer> ms1, ms2;
		for(int i=0; i<N; i++) {
			ms1 = this.dp[i].index.values();
			for(int j=i; j<N; j++) {
				ms2 = this.dp[j].index.values();
				if(CollectionUtils.containsAny(ms1, ms2)) {
					for(int u : ms1) {
						for(int v : ms2) {
							if(u!=v) {
								adjmat.get(u).add(v);
								adjmat.get(v).add(u);
							}
						}
					}
				}
			}
		}

		// now for each pair of adjacent markers 
		// we check if there is a path to connect them
		int source, target;
		boolean reachable;
		final Set<Integer> visited = new HashSet<Integer>();
		final LinkedList<Integer> queue = new LinkedList<Integer>();
		Set<Integer> neighbours;
		final List<Integer> breakpoints = new ArrayList<Integer>();
		for(int i=1; i<M; i++) {
			source = i-1;
			target = i  ;
			queue.clear();
			visited.clear();
			queue.offer(source);
			reachable = false;
			while(!queue.isEmpty()) {
				source = queue.pop();
				neighbours = adjmat.get(source);
				if(source==target||
						neighbours.contains(target)) {
					reachable = true;
					break;
				} else {
					visited.add(source);
					for(int z : neighbours) {
						if(!visited.contains(z))
							queue.offer(z);
					}
				}
			}
			if(!reachable) breakpoints.add(i);
		}

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
		
		int[][] hapInt = new int[Constants._ploidy_H][M];
		for(int i=0; i<M; i++) {
			double[][] probs = this.emissProbs[i].probsMat; // Px2
			for(int j=0; j<Constants._ploidy_H; j++)
				hapInt[j][i] = hapCov[j][i]?(probs[j][0]>0.5?0:1):-1;
		}
		
		int h, a1, a2, m1, m2, h1, h2;
		final int[] depth = new int[M-1];
		final TreeMap<Range<Integer>, Integer> mismatch = new TreeMap<Range<Integer>, Integer>(
				new Comparator<Range<Integer>>() {

					@Override
					public int compare(Range<Integer> r0, Range<Integer> r1) {
						// TODO Auto-generated method stub
						if(r0.lowerEndpoint().equals(r1.lowerEndpoint())) {
							return Integer.compare(r0.upperEndpoint().intValue(), 
									r1.upperEndpoint().intValue());
						}
						return Integer.compare(r0.lowerEndpoint().intValue(), 
								r1.lowerEndpoint().intValue());
					}

				});
		Range<Integer> range;
		for(int i=0; i<N; i++) {
			dp1 = this.dp[i];
			h   = this.haps[i];
			final List<Integer> ms = new ArrayList<Integer>(dp1.index.values());
			Collections.sort(ms);
			m1 = ms.get(0);
			a1 = dp1.probs[dp1.index.getKey(m1)][1]==soften?0:1;
			h1 = hapInt[h][m1];
			for(int j=1; j<ms.size(); j++) {
				m2 = ms.get(j);
				a2 = dp1.probs[dp1.index.getKey(m2)][1]==soften?0:1;
				h2 = hapInt[h][m2];
				range = Range.open(m1, m2);
				if( (a1==h1)!=(a2==h2) ) {
					// so we have a switch error here
					if(!mismatch.containsKey(range)) {
						mismatch.put(range, 1);
					} else {
						mismatch.put(range, mismatch.get(range)+1);
					}
				}
				for(int k=m1; k<m2; k++) ++depth[k];
				m1 = m2;
				a1 = a2;
				h1 = h2;
			}
		}
		
		if(this.ddebug) {
			myLogger.info("#### START mismatches");
			for(Map.Entry<Range<Integer>, Integer> entry : mismatch.entrySet()) {
				range = entry.getKey();
				myLogger.info(range.toString()+": "+entry.getValue());
			}
			myLogger.info("###### END mismatches");
		}
		
		// now for each pair of adjacent markers
		// if the mismatch rate is greater that 0.5 
		// then we break the block at the position
		int bp, mismatch_count, mc;
		final Set<Range<Integer>> keys = new HashSet<Range<Integer>>(),
				tmp_keys = new HashSet<Range<Integer>>();
		
		while(true) {
			bp = -1;
			mismatch_count = 0;
			// find a point with maximum inconsistency
			for(int i=0; i<M-1; i++) {
				mc = 0;
				tmp_keys.clear();
				for(Map.Entry<Range<Integer>, Integer> entry : mismatch.entrySet()) {
					range = entry.getKey();
					if(range.upperEndpoint()<=i) continue;
					if(range.lowerEndpoint()> i) break;
					mc += entry.getValue();
					tmp_keys.add(range);
				}
				if(this.ddebug) myLogger.info("#mismatch@"+i+": "+mc+"/"+depth[i]);
				if(mc>depth[i]*0.5&&mc>mismatch_count) {
					bp = i;
					mismatch_count = mc;
					keys.clear();
					keys.addAll(tmp_keys);
				}
			}
			if(bp==-1) break;
			breakpoints.add(bp+1);
			for(final Range<Integer> key : keys) mismatch.remove(key);
		}
		
		breakpoints.add(0);
		breakpoints.add(M);

		Collections.sort(breakpoints);
		
		return breakpoints;
	}
}





