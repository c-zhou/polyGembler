package cz1.hmm.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Serializable;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.HashSet;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.text.NumberFormat;
import java.text.DecimalFormat;
import java.lang.StringBuilder;

import cz1.hmm.tools.RFUtils.FileExtraction;
import cz1.hmm.tools.RFUtils.FileLoader;
import cz1.hmm.tools.RFUtils.FileObject;
import cz1.hmm.tools.RFUtils.InputStreamObj;
import cz1.hmm.tools.RFUtils.PhasedDataCollection;
import cz1.util.Algebra;
import cz1.util.ArgsEngine;
import cz1.util.Combination;
import cz1.util.Constants;
import cz1.util.Utils;
import cz1.util.JohnsonTrotter;
import cz1.util.Permutation;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.builder.EqualsBuilder;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.correlation.KendallsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.apache.commons.math3.stat.inference.GTest;

// recombination frequency estimator from resampling
public class RFEstimatorRS3 extends RFUtils {	

	// private static long[] scale_times = new long[10]; 

	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--hap-file               Directory with input haplotype files.\n"
						+ " -o/--prefix                 Output file prefix.\n"
						+ " -ex/--experiment-id         Common prefix of haplotype files for this experiment.\n"
						+ " -f/--parent                 Parent samples (seperated by a \":\").\n"
						+ " -p/--ploidy                 Ploidy of genome (default 2).\n"
						+ " -nb/--best                  The most likely nb haplotypes will be used (default 10).\n"
						+ " -phi/--skew-phi             For a haplotype inference, the frequencies of parental \n"
						+ "                             haplotypes need to be in the interval [1/phi, phi], \n"
						+ "                             otherwise will be discared (default 2).\n"
						+ " -nd/--drop                  At least nd haplotype inferences are required for \n"
						+ "                             a contig/scaffold to be analysed (default 1).\n"
						+ " -t/--threads                Threads (default 1).\n"	
				);
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}

		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add( "-ex", "--experiment-id", true);
			myArgsEngine.add( "-i", "--hap-file", true);
			myArgsEngine.add( "-o", "--prefix", true);
			myArgsEngine.add( "-f", "--parent", true);
			myArgsEngine.add( "-p", "--ploidy", true);
			myArgsEngine.add( "-nb", "--best", true);
			myArgsEngine.add( "-t", "--threads", true);
			myArgsEngine.add( "-phi", "--skew-phi", true);
			myArgsEngine.add( "-nd", "--drop", true);
			myArgsEngine.parse(args);
		}
		
		if(myArgsEngine.getBoolean("-i")) {
			in_haps = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your input zip file.");
		}

		if(myArgsEngine.getBoolean("-o")) {
			out_prefix = myArgsEngine.getString("-o");
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify your output file prefix.");
		}
		
		if(myArgsEngine.getBoolean("-ex")) {
			expr_id = myArgsEngine.getString("-ex");
		}  else {
			expr_id = guessExperimentId();
			myLogger.warn("No experiment prefix provided, I guess it's "+expr_id+". Please\n"
					+ "specify it with -ex/--experiment-id option if it's incorrect.");
		}

		if(myArgsEngine.getBoolean("-f")) {
			founder_haps = myArgsEngine.getString("-f").split(":");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the parent samples (seperated by a \":\").");
		}
		
		if(myArgsEngine.getBoolean("-p")) {
			int ploidy = Integer.parseInt(myArgsEngine.getString("-p"));
			Constants.ploidy(ploidy);
			Constants._haplotype_z = ploidy*2;
			probs_uniform = new double[ploidy*2];
			Arrays.fill(probs_uniform, .5/ploidy);
			ploidy2 = ploidy/2;
			shift_bits2 = mask_length*ploidy;
		}
		
		if(myArgsEngine.getBoolean("-nb")) {
			best_n = Integer.parseInt(myArgsEngine.getString("-nb"));
		}
		
		if(myArgsEngine.getBoolean("-t")) {
			THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}
		
		if(myArgsEngine.getBoolean("-phi")) {
			skew_phi = Double.parseDouble(myArgsEngine.getString("-phi"));
		}
		
		if(myArgsEngine.getBoolean("-nd")) {
			drop_thres = Integer.parseInt(myArgsEngine.getString("-nd"));
		}
		
		setHaps();
	}

	private BufferedWriter rfMinimumWriter;
	
	@Override
	public void run() {
		// TODO Auto-generated method stub

		this.initialise();
		
		this.hashDC();
		
		this.initial_thread_pool();
		rfMinimumWriter = Utils.getBufferedWriter(out_prefix+".txt");
		
		Set<String> kept_scaffs = new HashSet<String>();
		for(int i=0; i<this.dc.length; i++) {
			if(this.dc[i][0]!=null)
				kept_scaffs.add(this.dc[i][0].markers[0].
						replaceAll("_[0-9]{1,}$", ""));
		}
		try {
			for(String scaff : kept_scaffs)
					rfMinimumWriter.write("##"+scaff+"\n");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		long start = System.nanoTime();
		for(int i=0; i<dc.length; i++) 
			for(int j=i+1; j<dc.length; j++) 
				executor.submit(new MyRfCalculator(i, j));
		this.waitFor();
		
		long end = System.nanoTime();
		myLogger.info("elapsed, "+(end-start)+"; R pool key-value pair, "+rf_pool.size());
		try {
			rfMinimumWriter.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		this.initial_thread_pool();
		for(int i=0; i<dc.length; i++) 
			executor.submit(new MyMapCalculator(i));
		this.waitFor();
		writeMap();
		
		myLogger.info("["+Utils.getSystemTime()+"] DONE.");
	}
	
	protected MyPhasedDataCollection[][] dc;
	
	protected class MyPhasedDataCollection extends PhasedDataCollection {

		final boolean[][] dup_flag; // if duplicated false, else ture
		final int[] non_dup_haps; // number of haplotypes not duplicated
		private boolean[][][][] data;

		public MyPhasedDataCollection(String file,
				String[] markers, 
				int[] start_end,
				boolean[] dup_haps) {
			// TODO Auto-generated constructor stub
			super(file, markers, start_end);	
			this.dup_flag = new boolean[2][Constants._ploidy_H];
			System.arraycopy(dup_haps, 0, 
					this.dup_flag[0], 0, Constants._ploidy_H);
			System.arraycopy(dup_haps, Constants._ploidy_H, 
					this.dup_flag[1], 0, Constants._ploidy_H);
			non_dup_haps = new int[2];
			for(int i=0; i<2; i++)
				for(int j=0; j<Constants._ploidy_H; j++)
					if(this.dup_flag[i][j])
						non_dup_haps[i]++;
		}

		public MyPhasedDataCollection(String file,
				String[] markers, 
				int[] start_end, 
				boolean[][] dup_haps,
				int[] non_dup_haps,
				boolean[][][][] data) {
			// TODO Auto-generated constructor stub
			super(file, markers, start_end);
			this.dup_flag = dup_haps;
			this.non_dup_haps = non_dup_haps;
			this.data = data;
		}

		protected void data() {
			// TODO Auto-generated method stub
			this.data = this.data(
					this.start_end_position[0],
					this.start_end_position[1]);
		}

		protected MyPhasedDataCollection clone() {
			return new MyPhasedDataCollection(this.file, 
					this.markers, 
					this.start_end_position,
					this.cloneDupHaps(),
					this.non_dup_haps,
					this.cloneData());
		}
		
		private boolean[][] cloneDupHaps() {
			// TODO Auto-generated method stub
			final boolean[][] duplicated_haps = new boolean[2][Constants._ploidy_H];
			duplicated_haps[0] = this.dup_flag[0].clone();
			duplicated_haps[1] = this.dup_flag[1].clone();
			return duplicated_haps;
		}

		protected boolean[][][][] cloneData() {
			final boolean[][][][] data = 
					new boolean[2][2][Constants._ploidy_H][nF1];
			for(int i=0; i<2; i++) 
				for(int j=0; j<2; j++)
					for(int k=0; k<Constants._ploidy_H; k++)
						data[i][j][k] = this.data[i][j][k].clone();
			return data;
		}

		private boolean[][][][] data(int start, int end) {
			// TODO Auto-generated method stub
			boolean[][][][] data = new boolean[2][2][Constants._ploidy_H][nF1];

			try {
				InputStreamObj isObj = new InputStreamObj(this.file);
				isObj.getInputStream("PHASEDSTATES");
				BufferedReader br = Utils.getBufferedReader(isObj.is);
				
				String line;
				String[] s;
				String stateStr;

				br.readLine();
				br.readLine();
				outerloop:
					for(int i=0; i<nF1; i++) {
						for(int k=0; k<Constants._ploidy_H; k++) {
							line = br.readLine();
							s = line.split("\\s+|:");
							if(Arrays.asList(founder_haps).contains(s[2])) {
								for(int z=0; z<Constants._ploidy_H-1; z++) br.readLine();
								i--;
								continue outerloop;
							}
							stateStr = s[s.length-1];
							int p = k<Constants._ploidy_H/2 ? 0 : 1;
							data[p][0][hap_index[stateStr.charAt(start)]][i] = true;
							data[p][1][hap_index[stateStr.charAt(end)]][i] = true;
						}
					}
				br.close();
				isObj.close();
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(1);
			}

			return data;
		}
	}
	
	private void hashDC() {
		// TODO Auto-generated method stub
		
		this.initial_thread_pool();
		for(int i=0; i<dc.length; i++) {
			for(int j=0; j<best_n; j++) {
				
				executor.submit(new Runnable() {
					private int i;
					private int j;
					
					@Override
					public void run() {
						// TODO Auto-generated method stub
						try {
							if(dc[i][j]==null) return;
							dc[i][j].data();
						} catch (Exception e) {
							Thread t = Thread.currentThread();
							t.getUncaughtExceptionHandler().uncaughtException(t, e);
							e.printStackTrace();
							executor.shutdown();
							System.exit(1);
						}
					}

					public Runnable init(int i, int j) {
						// TODO Auto-generated method stub
						this.i = i;
						this.j = j;
						return(this);
					}
					
				}.init(i, j));
			}
		}
		this.waitFor();
	}

	private final int[] hap_index = new int[256];
	private final int mask_length = 1;
	private int ploidy2 = 1;
	private int shift_bits2 = mask_length*2;
	private final Set<String> scaff_only = new HashSet<String>();	

	public RFEstimatorRS3() {}
			
	public RFEstimatorRS3 (String in_haps, 
				String out_prefix,
				String expr_id, 
				int ploidy, 
				String[] founder_haps,
				int threads,
				double skew_phi,
				int drop_thres,
				int best_n) { 
		this.in_haps = in_haps;
		this.out_prefix = out_prefix;
		this.expr_id = expr_id;
		Constants.ploidy(ploidy);
		ploidy2 = ploidy/2;
		shift_bits2 = mask_length*ploidy;
		this.founder_haps = founder_haps;
		THREADS = threads;
		this.skew_phi = skew_phi;
		this.drop_thres = drop_thres;
		this.best_n = best_n;
		probs_uniform = new double[ploidy*2];
		Arrays.fill(probs_uniform, .5/ploidy);
		setHaps();
	}
	
	private boolean[][] gamete_ref = null;
	private char[][] haps = null;
	private int factorial_ploidy = -1;
	private int combinant_ploidy = -1;
	private int[][] johnson_trotter_permutation;
	private List<Map<Integer, Integer>> johnson_trotter_swap_index;
			
	private void setHaps() {
		factorial_ploidy = (int) Permutation.factorial(Constants._ploidy_H);
		List<List<Character>> haps_list = new ArrayList<List<Character>>();
		List<Character> hap = new ArrayList<Character>();
		for(int i=1; i<=Constants._ploidy_H; i++) {
			hap.add( i<10 ? Character.toChars('0'+i)[0] : 
				Character.toChars('a'+i-10)[0] );
		}
		haps_list.add(hap);
		hap = new ArrayList<Character>();
		for(int i=Constants._ploidy_H+1; 
				i<=Constants._haplotype_z; i++) {
			hap.add( i<10 ? Character.toChars('0'+i)[0] : 
				Character.toChars('a'+i-10)[0] );
		}
		haps_list.add(hap);
		haps = new char[haps_list.size()][haps_list.get(0).size()];
		for(int i=0; i<haps.length; i++) {
			for(int j=0; j<haps[i].length; j++) {
				haps[i][j] = haps_list.get(i).get(j);
				hap_index[haps[i][j]] = j;
			}
		}
		johnson_trotter_permutation = JohnsonTrotter.perm(Constants._ploidy_H);
		initialiseRFactory();
		
		Integer[] h = new Integer[Constants._ploidy_H];
		for(int i=0; i<Constants._ploidy_H; i++) h[i] = i;
		ArrayList<List<Integer>> combs = Combination.combination(h, Constants._ploidy_H/2);
		combinant_ploidy = combs.size();
		gamete_ref = new boolean[combinant_ploidy][Constants._ploidy_H];
		for(int i=0; i<combinant_ploidy; i++) {
			List<Integer> comb = combs.get(i);
			for(Integer j : comb) gamete_ref[i][j] = true;
		}
		
		johnson_trotter_swap_index = new ArrayList<Map<Integer, Integer>>();
		boolean[][] gametes = new boolean[combinant_ploidy][Constants._ploidy_H];
		for(int i=0; i<combinant_ploidy; i++) gametes[i] = gamete_ref[i].clone();
		
		for(int i=0; i<factorial_ploidy; i++) {
			Map<Integer, Integer> map_i = new HashMap<Integer, Integer>();
			int[] jt = johnson_trotter_permutation[i];
			swapCol(gametes, jt[0], jt[1]);
			for(int j=0; j<combinant_ploidy; j++) 
				map_i.put(intKeyFromBoolArray(gametes[j]), j);
			johnson_trotter_swap_index.add(map_i);
		}
	}
	

	private void swapCol(final boolean[][] arr, 
			final int i, 
			final int j) {
		// TODO Auto-generated method stub
		for(int k=0; k<arr.length; k++) {
			boolean tmp = arr[k][i];
			arr[k][i] = arr[k][j];
			arr[k][j] = tmp;
		}
	}

	private int intKeyFromBoolArray(boolean[] bool) {
		// TODO Auto-generated method stub
		int key = 0;
		for(boolean b : bool) {
			key <<= 1;
			if(b) key += 1;
		}
		return key;
	}
	
	private int intKeyFromIntArray(int[] ints) {
		// TODO Auto-generated method stub
		boolean[] bool = new boolean[Constants._ploidy_H];
		for(int i : ints) 
			bool[i] = true;
		return intKeyFromBoolArray(bool);
	}

	private void initialiseRFactory() {
		// TODO Auto-generated method stub
		List<List<Integer>> combs = Combination.combination(
				Constants._ploidy_H, ploidy2);
		for(int c=0; c<combs.size(); c++) {
			List<Integer> comb = combs.get(c);
			boolean[] bools = new boolean[Constants._ploidy_H];
			for(int i=0; i<ploidy2; i++)
				bools[comb.get(i)] = true;
			int key = hashcode(bools);
			for(int c2=0; c2<combs.size(); c2++) {
				List<Integer> comb2  = combs.get(c2);
				boolean[] bools2 = new boolean[Constants._ploidy_H];
				for(int i=0; i<ploidy2; i++)
					bools2[comb2.get(i)] = true;
				int key2 = hashcode(bools2);
				byte rf = recombinationFreq(bools, bools2);
				rf_factory.put((key<<shift_bits2)+key2, rf);
			}
		}
		return;
	}

	private byte recombinationFreq(boolean[] bools, 
			boolean[] bools2) {
		// TODO Auto-generated method stub
		byte rf = 0;
		for(int i=0; i<Constants._ploidy_H; i++)
			if(bools[i]&&bools2[i])
				rf++;
		return rf;
	}

	private int hashcode(boolean[] bools) {
		// TODO Auto-generated method stub
		int hash = 0;
		for(int i=0; i<Constants._ploidy_H; i++)
			hash = (hash<<mask_length)+(bools[i] ? 1 : 0);
		return hash;
	}

	protected class MyFileLoader extends FileLoader {

		public MyFileLoader(String id, FileObject[] files, int i) {
			// TODO Auto-generated constructor stub
			super(id, files, i);
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub

			try {
				double[] ll = new double[this.files.length];
				for(int k=0; k<ll.length; k++) 
					ll[k]=Double.NEGATIVE_INFINITY;
				int marker_n = this.files[0].markers.length;

				if( marker_n<2 ) {
					myLogger.warn("warning: "+
							this.id +" #marker is less than 2.");
					return;
				}

				for(int k=0; k<this.files.length; k++) {
					String file = this.files[k].file;
					try {
						InputStreamObj isObj = new InputStreamObj(file);
						if( !isObj.getInputStream("PHASEDSTATES") ) {
							myLogger.warn("warning: "+file+
									" exsits, but phased states do not.");
							continue;
						} 

						BufferedReader br = Utils.getBufferedReader(isObj.is);
						br.readLine();
						String marker_str = br.readLine();
						br.close();
						isObj.closeIs();

						if( marker_str==null ) {
							myLogger.warn("warning: "+file+
									" exists, but phased states are NULL.");
							continue;
						}
						double frac = marker_n/Double.parseDouble(marker_str);

						String line;
						if( isObj.getInputStream("STDERR") ) {
							BufferedReader br2 = Utils.getBufferedReader(isObj.is);
							String lprob=null;
							while( (line=br2.readLine()) !=null) {
								if(line.startsWith("log"))
									lprob = line;
							}
							br2.close();
							isObj.closeIs();
							if( lprob!=null ) {
								String[] s0 = lprob.split("\\s+");
								ll[k] = Double.parseDouble(s0[3])*frac;
							}
						} else {
							isObj.getInputStream("PHASEDSTATES");
							BufferedReader br2 = Utils.getBufferedReader(isObj.is);
							String lprob = br2.readLine();
							br2.close();
							isObj.closeIs();
							if( lprob!=null ) 
								ll[k] = Double.parseDouble(lprob)*frac;
							myLogger.info(ll[k]);
						}
						isObj.closeIn();
					} catch (IOException e) {
						e.printStackTrace();
						System.exit(1);
					}
				}

				int[] maxN = maxN(ll);
				StringBuilder oos = new StringBuilder();
				boolean[] drop = new boolean[maxN.length];
				final boolean[][] dup_haps = new boolean[maxN.length][Constants._haplotype_z];
				int dropped = 0;
				
				for(int k=0; k<maxN.length; k++) {

					final int[][] haps_duplicated = duplicateHaplotypes(maxN[k]);
					final int[] haps_observed = readHaplotypes(maxN[k]);
					int dz = Constants._haplotype_z;
					final boolean[] dup_hapsK = new boolean[Constants._haplotype_z];
					Arrays.fill(dup_hapsK, true);
					for(int u=0; u<2; u++) { 
						if(haps_duplicated[u]!=null) {
							for(int v=0; v<haps_duplicated[u].length; v++) {
								dup_hapsK[u*Constants._ploidy_H+haps_duplicated[u][v]] = false;
								dz--;
							}
						}
					}
					dup_haps[k] = dup_hapsK;
					
					if(dz<2) {
						drop[k] = true;
						dropped++;
						continue;
					}
					
					long[] observed;
					double p;
					int w;
					switch(goodness_of_fit) {
					case "fraction":
						double[] phases = new double[dz];
						w = 0;
						for(int z=0; z<Constants._haplotype_z; z++)
							if(dup_hapsK[z])
								phases[w++] = (double) haps_observed[z];
						double expected = StatUtils.sum(phases)/Constants._ploidy_H/2;
						double maf = StatUtils.max(phases)/expected, 
								mif = StatUtils.min(phases)/expected;
						if( maf>skew_phi || mif<1/skew_phi) {
							myLogger.info(this.files[maxN[k]].file+
									" was dropped due to large haploptype frequency variance. (" +
									cat(phases, ",") +")");
							drop[k] = true;
						}
						if(drop[k]) dropped++;
						oos.append("["+(drop[k]?"drop](maf,":"keep](maf,")+maf+";mif,"+mif+") "+
								cat(phases,",")+"\t"+
								this.files[maxN[k]].file+"\n");
						break;
					case "chisq":
						observed = new long[dz];
						w = 0;
						for(int z=0; z<Constants._haplotype_z; z++) 
							if(dup_hapsK[z])
								observed[w++] = (long) haps_observed[z];
						p = new ChiSquareTest().chiSquareTest(
								uniformProbs(Constants._haplotype_z-dz), observed);
						if(p<skew_phi) drop[k] = true;
						if(drop[k]) dropped++;
						oos.append("["+(drop[k]?"drop](p,":"keep](p,")+formatter.format(p)+") "+
								cat(observed,",")+"\t"+
								this.files[maxN[k]].file+"\n");
						break;
					case "gtest":
						observed = new long[dz];
						w = 0;
						for(int z=0; z<Constants._haplotype_z; z++) 
							if(dup_hapsK[z])
								observed[w++] = (long) haps_observed[z];
						p = new GTest().gTest(
								uniformProbs(Constants._haplotype_z-dz), observed);
						if(p<skew_phi) drop[k] = true;
						if(drop[k]) dropped++;
						oos.append("["+(drop[k]?"drop](p,":"keep](p,")+formatter.format(p)+") "+
								cat(observed,",")+"\t"+
								this.files[maxN[k]].file+"\n");
						break;
					default:
						throw new RuntimeException("Goodness-of-fit test should be fraction, "
								+ "chisq or gTest.");
					}
				}
				myLogger.info(oos.toString());
				myLogger.info(this.id+" - dropped "+dropped);
				if( drop.length-dropped<drop_thres ) {
					myLogger.info("Scaffold "+this.id+" dropped.");
				} else {
					int z=0;
					for(int k=0; k<drop.length; k++) {
						if(!drop[k]) {
							FileObject fobj = this.files[maxN[k]];
							this.collectData(i, z, dup_haps[k], fobj);
							z++;
						}
						if(z>=best_n) break;
					}
				}

			} catch (Exception e) {
				// TODO Auto-generated catch block
				Thread t = Thread.currentThread();
				t.getUncaughtExceptionHandler().uncaughtException(t, e);
				e.printStackTrace();
				executor.shutdown();
				System.exit(1);
			}
		}

		private double[] uniformProbs(final int i) {
			// TODO Auto-generated method stub
			double[] probs = new double[i];
			for(int k=0; k<i; k++)
				probs[k] = 1.0/i;
			return probs;
		}

		protected void collectData(int i, 
				int z,
				boolean[] dup_haps, 
				FileObject fobj) {
			// TODO Auto-generated method stub
			dc[i][z] = new MyPhasedDataCollection(
					fobj.file,
					fobj.markers,
					fobj.start_end_position,
					dup_haps);
		}

		@Override
		protected String getPhaseFile() {
			// TODO Auto-generated method stub
			return "PHASEDSTATES";
		}
		
		private int[][] duplicateHaplotypes(final int f) {
			// TODO Auto-generated method stub
			final int[][] dup_haps = new int[2][];
			
			final int start = this.files[f].start_end_position[0];
			final int end = this.files[f].start_end_position[1];
			final int nM = end-start+1;
			final String[][] alleles = new String[nM][Constants._haplotype_z];
			try {
				InputStreamObj isObj = new InputStreamObj(this.files[f].file);
				isObj.getInputStream("EMISS");
				BufferedReader br = Utils.getBufferedReader(isObj.is);
				for(int i=1; i<start; i++) br.readLine();
				
				String[] probs_str;
				String line;
				Pattern pat = Pattern.compile("\\{(.*?)\\}");
				for(int i=start; i<=end; i++) {
					line = br.readLine();
					Matcher m = pat.matcher(line);
					m.find(); // skip dummy hidden state
					for(int j=0; j<Constants._haplotype_z; j++) {
						m.find();
						probs_str = m.group(1).split(",|;");
						String a = probs_str[0];
						double p = Double.parseDouble(probs_str[1]);
						for(int k=2; k<probs_str.length-1; k+=2) { 
							double q = Double.parseDouble(probs_str[k+1]);
							if(q>p) {
								a = probs_str[k];
								p = q;
							}
						}
						alleles[i-start][j] = a;
					}
				}
				br.close();
				isObj.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			for(int i=0; i<2; i++) {
				int sf = i*Constants._ploidy_H;
				boolean[][] iden = new boolean[Constants._ploidy_H][Constants._ploidy_H];
				for(int j=0; j<Constants._ploidy_H; j++) {
					iden[j][j] = true;
					for(int k=0; k<Constants._ploidy_H; k++) {
						int a = 0;
						for(int l=0; l<nM; l++)
							if(alleles[l][j+sf].equals(alleles[l][k+sf])) 
								a++;
							else break;
						if(a==nM) {
							iden[j][k] = true;
							iden[j][k] = true;
						}
					}
				}
				List<Integer> hap_list = new ArrayList<Integer>();
				for(int j=0; j<Constants._ploidy_H; j++) {
					int z = 0;
					for(int k=0; k<Constants._ploidy_H; k++) 
						if(iden[j][k]) z++;
					if(z>1) hap_list.add(j);
				}
				if(hap_list.size()>0) 
					dup_haps[i] = ArrayUtils.toPrimitive( hap_list.
							toArray(new Integer[hap_list.size()]) ); 
			}
			return dup_haps;
		}

		@Override
		protected void collectData(int i, int z, FileObject fobj) {
			// TODO Auto-generated method stub
			throw new RuntimeException("!!!");
		}
	}
	
	protected class MyMapCalculator extends MapCalculator {
		
		public MyMapCalculator(int i) {
			// TODO Auto-generated constructor stub
			super(i);
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub

			try {
				if(dc[i][0]==null) return;
				double[][] rfs = new double[dc[i].length][];

				for(int k=0; k<dc[this.i].length; k++) {
					MyPhasedDataCollection dc_ik = dc[i][k];
					if(dc_ik==null) break;
					myLogger.info(dc_ik.file);
					if(dc_ik !=null ) {
						rfs[k] = calcGDs(
								dc_ik.file,
								Constants._ploidy_H,
								founder_haps,
								nF1,
								dc_ik.start_end_position);
					}
				}

				mapCalc.put(dc[i][0].markers[0].replaceAll("_[0-9]{1,}$", ""), rfs);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				Thread t = Thread.currentThread();
				t.getUncaughtExceptionHandler().uncaughtException(t, e);
				e.printStackTrace();
				executor.shutdown();
				System.exit(1);
			}

			/**
			String contig = dc[i][0].markers[0].replaceAll("_[0-9]{1,}$", "");
			try {
				mapWriter.write("*"+contig+"\t"+
						median(kd_all)+"\t"+
						StatUtils.sum(kosambi)+"\t"+
						cat(kosambi, ",")+"\n");
			} catch (MathIllegalArgumentException | IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			 **/
		}
	}
		
	private class MyRfCalculator extends RfCalculator {
	
		public MyRfCalculator(int i, int j) {
			// TODO Auto-generated constructor stub
			super(i, j);
		}
		
		@Override
		public void run() {
			// TODO Auto-generated method stub
			try {
				MyPhasedDataCollection[] dc_i = dc[this.i];
				MyPhasedDataCollection[] dc_j = dc[this.j];
				String contig=null, contig2=null;
				int ii = 0;
				while(ii<dc_i.length && contig==null) {
					contig = dc_i[ii]==null ? null : 
						dc_i[ii].markers[0].replaceAll("_[0-9]{1,}$", "");
					ii++;
				}
				int jj = 0;
				while(jj<dc_j.length && contig2==null) {
					contig2 = dc_j[jj]==null ? null : 
						dc_j[jj].markers[0].replaceAll("_[0-9]{1,}$", "");
					jj++;
				}
				if(contig==null || contig2==null) return;
				if(!scaff_only.isEmpty() && 
						!scaff_only.contains(contig) && 
						!scaff_only.contains(contig2)) return;
	
				int m = dc_i.length, n = dc_j.length;
				double[][] rf_all = new double[m*n][4];
				for(int k=0; k<rf_all.length; k++)
					Arrays.fill(rf_all[k], -1);
	
				for(int k=0; k<m; k++) 
					for(int s=0; s<n; s++) 
						if(dc_i[k]!=null && dc_j[s]!=null)
							rf_all[k*m+s] = calcRFs(k, s);
	
				rf_all = Algebra.transpose(rf_all);
				StringBuilder os = new StringBuilder();
				os.append("#");
				os.append(this.i);
				os.append(" ");
				os.append(this.j);
				os.append("\n");
				for(int k=0; k<rf_all.length; k++) {
					os.append(formatter.format(rf_all[k][0]));
					for(int u=1; u<rf_all[k].length; u++) {
						os.append(" ");
						os.append(formatter.format(rf_all[k][u]));
					}
					os.append("\n");
				}
				double[] rf; 
				String w;
	
				rf = new double[4];
				for(int k=0; k<4; k++) {
					double[] rf_tmp = removeNEG(rf_all[k]);
					if(rf_tmp==null)
						rf[k] = -1;
					else
						rf[k] = StatUtils.min(rf_tmp);
				}
				w = StatUtils.min(rf)+"\t";
				for(int k=0; k<rf.length; k++) w += rf[k]+"\t";
				w += contig+"\t"+contig2+"\n";
				rfMinimumWriter.write(w);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				Thread t = Thread.currentThread();
				t.getUncaughtExceptionHandler().uncaughtException(t, e);
				e.printStackTrace();
				executor.shutdown();
				System.exit(1);
			}
		}
	
		private double[] calcRFs(int phase, int phase2) {
			// TODO Auto-generated method stub
			boolean[][][][] dc_i = dc[i][phase].data;
			boolean[][] dh_i = dc[i][phase].dup_flag;
			double[] rf = new double[4];
			boolean[][][][] data_j = dc[j][phase2].cloneData();
			boolean[][] dh_j = dc[j][phase2].cloneDupHaps();
			
			int no_haps = 0; // number of haps
			for(int p=0; p<2; p++) {
				boolean[][][] data_j_p = data_j[p];
				boolean[] dh_j_p = dh_j[p];
				final int pp = Math.min(dc[i][phase].non_dup_haps[p], 
						dc[j][phase2].non_dup_haps[p]);
				no_haps += pp*nF1/2;
				double[][] rf_all_p = new double[factorial_ploidy][4];
				
				
				for(int f=0; f<factorial_ploidy; f++) {
					next(data_j_p, dh_j_p, f);
					rf_all_p[f] = calcRFs(dc_i[p], 
							dh_i[p], 
							data_j_p, 
							dh_j_p,
							pp);
				}
				add(rf, max(rf_all_p));
			}
			return mininus(1, divide(rf, no_haps));
		}

		private double[] calcRFs(final boolean[][][] data,
				final boolean[] dh, 
				final boolean[][][] data2, 
				final boolean[] dh2,
				final int p) {
			// TODO Auto-generated method stub
			double[] rf = new double[4];
			double z;
			for(int m=0; m<2; m++) {
				boolean[][] data_m = data[m];
				for(int k=0; k<2; k++) {
					boolean[][] data_k = data2[k];
					z = 0;
					for(int w=0; w<nF1; w++)
						for(int v=0; v<Constants._ploidy_H; v++) 
							if(dh[v]&&data_m[v][w]&&dh2[v]&&data_k[v][w]) 
								z++;
					rf[m*2+k] = z;
				}
			}
			return rf;
		}

		private double[] doubles(int[] ints) {
			// TODO Auto-generated method stub
			double[] ds = new double[ints.length];
			for(int i=0; i!=ints.length; i++)
				ds[i] = (double) ints[i];
			return ds;
		}

		private void next(boolean[][][] data, boolean[] dh, int f) {
			// TODO Auto-generated method stub
			final int i = johnson_trotter_permutation[f][0],
					i2 = johnson_trotter_permutation[f][1];
			for(int m=0; m<2; m++) {
				boolean[][] tmps = data[m];
				boolean[] tmp = tmps[i];
				tmps[i] = tmps[i2];
				tmps[i2] = tmp;
			}
			boolean tmp = dh[i];
			dh[i] = dh[i2];
			dh[i2] = tmp;
		}
	
		private double[] divide(double[] ds, double deno) {
			// TODO Auto-generated method stub
			for(int i=0; i<4; i++) ds[i] /= deno;
			return ds;
		}
		
		private double[] mininus(double d, double[] minus) {
			// TODO Auto-generated method stub
			for(int i=0; i<4; i++) minus[i] = d-minus[i];
			return minus;
		}
	
		private void add(double[] ds, double[] ds2) {
			// TODO Auto-generated method stub
			for(int i=0; i<4; i++) ds[i] += ds2[i];
		}
	
		private double[] max(double[][] ds) {
			// TODO Auto-generated method stub
			double tot=Double.NEGATIVE_INFINITY, 
					rf=Double.NEGATIVE_INFINITY;
			double f;
			int r = -1;
			int n = ds.length;
			for(int i=0; i<n; i++) 
				if( (f=StatUtils.max(ds[i]))>rf ||
						f==rf && StatUtils.sum(ds[i])>tot ) {
					rf = f;
					tot = StatUtils.sum(ds[i]);
					r = i;
				}
			return ds[r];
		}
		
		private double[] min(double[][] ds) {
			// TODO Auto-generated method stub
			double tot=Double.POSITIVE_INFINITY, 
					rf=Double.POSITIVE_INFINITY;
			double f;
			int r = -1;
			int n = ds.length;
			for(int i=0; i<n; i++) 
				if( (f=StatUtils.min(ds[i]))<rf ||
						f==rf && StatUtils.sum(ds[i])<tot ) {
					rf = f;
					tot = StatUtils.sum(ds[i]);
					r = i;
				}
			return ds[r];
		}
	}

	private final Map<Integer, Byte> rf_factory = new HashMap<Integer, Byte>();
	private final Map<Long, byte[]> rf_pool = new ConcurrentHashMap<Long, byte[]>();

	@Override
	protected void initialise() {
		// TODO Auto-generated catch block
		super.initialise();
		this.initial_thread_pool();
		String[] scaff_all = new String[fileObj.keySet().size()];
		fileObj.keySet().toArray(scaff_all);
		this.dc = new MyPhasedDataCollection[scaff_all.length][best_n];
		for(int i=0; i<scaff_all.length; i++) {
			Set<FileObject> files = fileObj.get(scaff_all[i]);
			executor.submit(new MyFileLoader(scaff_all[i],
					files.toArray(new FileObject[files.size()]),
					i));
		}
		this.waitFor();

		myLogger.info("["+Utils.getSystemTime()+"] LOADING FILES DONE.");
		myLogger.info("["+Utils.getSystemTime()+"] READING LOG LIKELIHOOD DONE.");
	}
}



























