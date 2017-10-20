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
public class RFEstimatorRS2 extends RFUtils {	

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
		
		protected MyPhasedDataCollection(String file,
				String[] markers,
				int[] start_end) {
			// TODO Auto-generated constructor stub
			super(file, markers, start_end);
			this.indentical_gametes = 
					this.makeIndenticalGameteMap(start_end[0],start_end[1]);
		}

		final boolean[][][][] indentical_gametes;
		private boolean[][][][] data;
		
		protected MyPhasedDataCollection(String file,
				String[] markers, 
				int[] start_end,
				boolean[][][][] indentical_gametes,
				boolean[][][][] data) {
			// TODO Auto-generated constructor stub
			super(file, markers, start_end);
			this.indentical_gametes = indentical_gametes;
			this.data = data;
		}
		
		private boolean[][][][] makeIndenticalGameteMap(int start, int end) {
			// TODO Auto-generated method stub
			final boolean[][][][] indentical_gametes = 
					new boolean[2][combinant_ploidy][][];
			int nM = end-start+1;
			String[][] alleles = new String[nM][Constants._haplotype_z];
			try {
				InputStreamObj isObj = new InputStreamObj(this.file);
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
				boolean[][] eq = new boolean[Constants._ploidy_H][Constants._ploidy_H];
				for(int j=0; j<Constants._ploidy_H; j++) {
					eq[j][j] = true;
					for(int k=0; k<Constants._ploidy_H; k++) {
						int a = 0;
						for(int l=0; l<nM; l++)
							if(alleles[l][j+sf].equals(alleles[l][k+sf])) 
								a++;
							else break;
						if(a==nM) {
							eq[j][k] = true;
							eq[j][k] = true;
						}
					}
				}

				for(int j=0; j<combinant_ploidy; j++) {	
					boolean[] gamete_j = gamete_ref[j];
					int[] gamete = new int[Constants._ploidy_H/2];
					int z = 0;
					for(int k=0; k<Constants._ploidy_H; k++)
						if(gamete_j[k]) gamete[z++] = k;
					int[][] na = new int[Constants._ploidy_H/2][];
					int[] n = new int[Constants._ploidy_H/2];
					int N = 1;
					for(int k=0; k<Constants._ploidy_H/2; k++) {
						na[k] = getIndeticalAllele(eq[gamete[k]]);
						n[k] = na[k].length;
						N *= n[k];
					}
					int[][] indices = new int[N][Constants._ploidy_H/2];
					
					int rep_a = 1, rep_b = N;
					for(int k=0; k<Constants._ploidy_H/2; k++) {
						rep_b /= n[k];
						int i_n = 0;
						for(int l=0; l<rep_a; l++) 
							for(int w=0; w<n[k]; w++)
								for(int v=0; v<rep_b; v++)
									indices[i_n++][k] = na[k][w];
						rep_a *= n[k];
					}
					for(int k=0; k<N; k++) Arrays.sort(indices[k]);
					Set<Integer> keys = new HashSet<Integer>();
					List<Integer> keys_retained = new ArrayList<Integer>();
					outerloop:
						for(int k=0; k<N; k++) {
							for(int l=1; l<Constants._ploidy_H/2; l++)
								if(indices[k][l-1]==indices[k][l])
									continue outerloop;
							int key = intKeyFromIntArray(indices[k]);
							if(!keys.contains(key)) {
								keys.add(key);
								keys_retained.add(k);
							}
						}
					if(!keys.contains(intKeyFromIntArray(gamete)))
						throw new RuntimeException("!!!");
					
					indentical_gametes[i][j] = new boolean[Constants._ploidy_H][keys_retained.size()];
					for(int k=0; k<keys_retained.size(); k++) {
						int w = keys_retained.get(k);
						for(int l=0; l<Constants._ploidy_H/2; l++)
							indentical_gametes[i][j][indices[w][l]][k] = true;
					}
				}
			}
			
			return indentical_gametes;
		}
		
		private int[] getIndeticalAllele(boolean[] bs) {
			// TODO Auto-generated method stub
			int n = 0;
			for(boolean b : bs) if(b) n++;
			int[] na = new int[n];
			int k = 0;
			for(int i=0; i<bs.length; i++) 
				if(bs[i]) na[k++] = i;
			return na;
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
					this.cloneIndenticalGameteMap(),
					this.cloneData());
		}
		
		private boolean[][][][] cloneIndenticalGameteMap() {
			// TODO Auto-generated method stub
			final boolean[][][][] indentical_gametes = 
					new boolean[2][combinant_ploidy][Constants._ploidy_H][];
			for(int p=0; p<2; p++)
				for(int i=0; i<combinant_ploidy; i++)
					for(int j=0; j<Constants._ploidy_H; j++) 
						indentical_gametes[p][i][j] = 
						this.indentical_gametes[p][i][j].clone();
			return indentical_gametes;
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

	public RFEstimatorRS2() {}
			
	public RFEstimatorRS2 (String in_haps, 
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
		protected void collectData(int i, int z, FileObject fobj) {
			// TODO Auto-generated method stub
			dc[i][z] = new MyPhasedDataCollection(
					fobj.file,
					fobj.markers,
					fobj.start_end_position);
		}

		@Override
		protected String getPhaseFile() {
			// TODO Auto-generated method stub
			return "PHASEDSTATES";
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
			boolean[][][][] iden_i = dc[i][phase].indentical_gametes;
			double[] rf = new double[4];
			boolean[][][][] data_j = dc[j][phase2].cloneData();
			boolean[][][][] iden_j = dc[j][phase2].cloneIndenticalGameteMap();
			
			for(int p=0; p<2; p++) {
				boolean[][][] data_j_p = data_j[p];
				boolean[][][] iden_j_p = iden_j[p];
				
				double[][] rf_all_p = new double[factorial_ploidy][4];
				
				for(int f=0; f<factorial_ploidy; f++) {
					next(data_j_p, iden_j_p, f);
					rf_all_p[f] = calcRFs(dc_i[p], 
							iden_i[p], 
							johnson_trotter_swap_index.get(0), // dc_i next changed 
							data_j_p, 
							iden_j_p, 
							johnson_trotter_swap_index.get(f)); // has been swapped
				}
				add(rf, max(rf_all_p));
			}
			return mininus(1, divide(rf, Constants._ploidy_H*nF1));
		}

		private double[] calcRFs(boolean[][][] data, 
				boolean[][][] iden, 
				Map<Integer, Integer> swap_index, 
				boolean[][][] data2,
				boolean[][][] iden2, 
				Map<Integer, Integer> swap_index2) {
			// TODO Auto-generated method stub
			double[] rf = new double[4];
			for(int m=0; m<2; m++) {
				boolean[][] data_m = Algebra.transpose(data[m]);
				for(int k=0; k<2; k++) {
					boolean[][] data_k = Algebra.transpose(data2[k]);
					double z = 0;
					for(int w=0; w<nF1; w++) { 
						boolean[] data_m_w = data_m[w];
						boolean[] data_k_w = data_k[w];
						int key_m = intKeyFromBoolArray(data_m_w);
						int key_k = intKeyFromBoolArray(data_k_w);
						boolean[][] iden_m = iden[swap_index.get(key_m)];
						boolean[][] iden_k = iden2[swap_index2.get(key_k)];
						int x = iden_m[0].length, y = iden_k[0].length;
						double[] fs = new double[x*y];
						int q = 0;
						for(int u=0; u<x; u++) {
							for(int v=0; v<y; v++) {
								int c = 0;
								for(int a=0; a<Constants._ploidy_H; a++)
									if(iden_m[a][u]&&iden_k[a][v])
										c++;
								fs[q++] = c;
							}
						}
						z += StatUtils.max(fs);
					}
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

		private void next(boolean[][][] data, boolean[][][] iden, int f) {
			// TODO Auto-generated method stub
			final int i = johnson_trotter_permutation[f][0],
					i2 = johnson_trotter_permutation[f][1];
			for(int m=0; m<2; m++) {
				boolean[][] tmps = data[m];
				boolean[] tmp = tmps[i];
				tmps[i] = tmps[i2];
				tmps[i2] = tmp;
			}
			for(int m=0; m<combinant_ploidy; m++) {
				boolean[][] tmps = iden[m];
				boolean[] tmp = tmps[i];
				tmps[i] = tmps[i2];
				tmps[i2] = tmp;
			}
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



























