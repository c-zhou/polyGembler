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
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import java.text.NumberFormat;
import java.text.DecimalFormat;
import java.lang.StringBuilder;

import cz1.hmm.tools.RFUtils.PhasedDataCollection;
import cz1.util.Algebra;
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
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.apache.commons.math3.stat.inference.GTest;

public class Copy_11_of_SingleCNVHapRF extends RFUtils {	

	private static long[] scale_times = new long[10]; 

	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		
	}
	
	public static void main(String[] args) 
			throws IOException, InterruptedException {

		CommandLineParser parser = new PosixParser();

		// create the Options
		Options options = new Options();
		options.addOption( "e", "experiment", true, "experiment name." );
		options.addOption( "w", "workspace", true, "directory contains input files." );
		options.addOption( "o", "output", true, "output file prefix." );
		options.addOption( "p", "parent-obj",true, "parents; dilimited by a colon(:)." );
		options.addOption( "P", "ploidy",true, "ploidy of the genome." );
		options.addOption( "n", "n-best", true, "the best n run is used.");
		options.addOption( "t", "thread", true, "number threads.");
		options.addOption( "s", "skew", true, "skew of haplotype frequencies.");
		options.addOption( "d", "drop", true, "threshold do drop a scaffold.");
		options.addOption( "T", "test", true, "goodness-of-fit test.");
		options.addOption( "C", "scaff-only", true, "calculate RFs only for specified scaffolds.");
		String experiment=null, workspace=null, output=null, parentObj=null,
				goodnessOfFitTest="fraction";
		String[] parents = new String[2];
		Arrays.fill(parents, null);
		int ploidy=-1, n=1, t=1, d=0;
		double s=0.1;
		try {
			// parse the command line arguments
			CommandLine line = parser.parse( options, args );
			if( line.hasOption("e") ) {
				experiment = line.getOptionValue('e');
			}
			if(line.hasOption("w")) {
				workspace = line.getOptionValue("w");
			}
			if(line.hasOption("o")) {
				output = line.getOptionValue("o");
			}
			if(line.hasOption("p")) {
				parentObj = line.getOptionValue("p");
				parents = parentObj.split(":");
			}
			if(line.hasOption("P")) {
				ploidy = Integer.parseInt(line.getOptionValue("P"));
				if(ploidy<1 || ploidy%2==1)
					throw new RuntimeException("invalid ploidy!!!");
			}
			if(line.hasOption("n")) {
				n = Integer.parseInt(line.getOptionValue("n"));
			}
			if(line.hasOption("t")) {
				t = Integer.parseInt(line.getOptionValue("t"));
			}
			if(line.hasOption("s")) {
				s = Double.parseDouble(line.getOptionValue("s"));
			}
			if(line.hasOption("d")) {
				d = Integer.parseInt(line.getOptionValue("d"));
			}
			if(line.hasOption("T")) {
				goodnessOfFitTest = line.getOptionValue("T");
			}
			if(line.hasOption("C")) {
				String[] scaff_all = line.getOptionValue("C").split(",|;|:");
				for(String scaff : scaff_all) scaff_only.add(scaff);
			}
		}
		catch( ParseException exp ) {
			System.out.println( "Unexpected exception:" + exp.getMessage() );
		}

		Copy_11_of_SingleCNVHapRF sgRF = new Copy_11_of_SingleCNVHapRF(
				workspace,experiment,ploidy,parents,t,s,d,goodnessOfFitTest,n);
		sgRF.calcRFsForAll2(output);
	}
	
	private static int ploidy;
	private static int ploidy2;
	private static int mask_length;
	private static int shift_bits2;
	private static double[] probs_uniform;
	private final static Set<String> scaff_only = new HashSet<String>();	

	public Copy_11_of_SingleCNVHapRF (String workspace_, 
			String experiment_, 
			int ploidy_, 
			String[] parents_,
			int threads,
			double s,
			int d,
			String t,
			int n) {
		workspace = workspace_;
		experiment = experiment_;
		ploidy = ploidy_;
		ploidy2 = ploidy_/2;
		//mask_length = (int) Math.ceil(
		//		Math.log(ploidy2+1)/Math.log(2));
		mask_length = 1;
		shift_bits2 = mask_length*ploidy;
		shift_bits4 = mask_length*ploidy*2;
		parents = parents_;
		THREADS = threads;
		drop_thres = d;
		skew_thres = s;
		goodness_of_fit = t.toLowerCase();
		probs_uniform = new double[ploidy*2];
		Arrays.fill(probs_uniform, .5/ploidy);
		best_n = n;
		setHaps();
	}

	private static char[][] haps = null;
	private static int factorial_ploidy = -1;
	private final int[] hap_index = new int[256];
	private static int[][] johnson_trotter_permutation;
	
	private PhasedDataCollection[][] dc;
	
	protected class PhasedDataCollection extends RFUtils.PhasedDataCollection{	
		protected final boolean[][][][] data;
		protected int[][][] hashcode;

		public PhasedDataCollection(String file,
				String[] markers, 
				int[] start_end) {
			// TODO Auto-generated constructor stub
			super(file, markers, start_end);
			this.data = this.data(start_end[0], start_end[1]);
			this.hash();
		}
		
		public PhasedDataCollection(String file,
				String[] markers,
				boolean[][][][] data,
				int[] start_end) {
			// TODO Auto-generated constructor stub
			super(file, markers, start_end);
			this.data = data;
		}

		private void hash() {
			// TODO Auto-generated method stub
			this.hashcode = new int[2][2][nF1];
			for(int i=0; i<2; i++)
				for(int j=0; j<2; j++) {
					boolean[][] d = this.data[i][j];
					int[] code = this.hashcode[i][j];
					for(int k=0; k<nF1; k++) {
						int key = 0;
						for(int p=0; p<Constants._ploidy_H; p++)
							key = (key<<mask_length)+(d[p][k] ? 1 : 0);
						code[k] = key;
					}
				}
		}

		public PhasedDataCollection clone() {
			final boolean[][][][] data = new boolean[2][2][Constants._ploidy_H][nF1];
			for(int i=0; i<2; i++) 
				for(int j=0; j<2; j++)
					for(int k=0; k<Constants._ploidy_H; k++)
						data[i][j][k] = this.data[i][j][k].clone();
			return new PhasedDataCollection(this.file, 
					this.markers, data, 
					this.start_end_position);
		}

		private boolean[][][][] data(int start, int end) {
			// TODO Auto-generated method stub
			boolean[][][][] data = new boolean[2][2][Constants._ploidy_H][nF1];

			try {
				BufferedReader br = Utils.getBufferedReader(
						in_haps+"/"+
								this.file+
								"/phasedStates/"+
								expr_id+".txt");
				String line;
				String[] s;
				String stateStr;
				int n=0;

				while( (line=br.readLine())!=null ) {

					if(!line.startsWith("#")) continue;
					//if(skip++<2) continue;
					s = line.split("\\s+|:");
					if(Arrays.asList(founder_haps).contains(s[2])) continue;

					stateStr = s[s.length-1];
					data[0][0][hap_index[stateStr.charAt(start)]][n] = true;
					data[0][1][hap_index[stateStr.charAt(end)]][n] = true;

					for(byte i=1; i<Constants._ploidy_H/2; i++) {
						line = br.readLine();
						s = line.split("\\s+");
						stateStr = s[s.length-1];
						data[0][0][hap_index[stateStr.charAt(start)]][n] = true;
						data[0][1][hap_index[stateStr.charAt(end)]][n] = true;
					}
					for(byte i=0; i<Constants._ploidy_H/2; i++) {
						line = br.readLine();
						s = line.split("\\s+");
						stateStr = s[s.length-1];
						data[1][0][hap_index[stateStr.charAt(start)]][n] = true;
						data[1][1][hap_index[stateStr.charAt(end)]][n] = true;
					}

					n++;
				}
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(1);
			}

			return data;
		}
	}
	
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
	}
	

	private void initialiseRFactory() {
		// TODO Auto-generated method stub
		List<List<Integer>> combs = Combination.combination(ploidy, ploidy2);
		for(int c=0; c<combs.size(); c++) {
			List<Integer> comb = combs.get(c);
			boolean[] bools = new boolean[ploidy];
			for(int i=0; i<ploidy2; i++)
				bools[comb.get(i)] = true;
			int key = hashcode(bools);
			for(int c2=0; c2<combs.size(); c2++) {
				List<Integer> comb2  = combs.get(c2);
				boolean[] bools2 = new boolean[ploidy];
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
		for(int i=0; i<ploidy; i++)
			if(bools[i]&&bools2[i])
				rf++;
		return rf;
	}

	private int hashcode(boolean[] bools) {
		// TODO Auto-generated method stub
		int hash = 0;
		for(int i=0; i<ploidy; i++)
			hash = (hash<<mask_length)+(bools[i] ? 1 : 0);
		return hash;
	}

	private final static Map<Integer, Byte> rf_factory = new HashMap<Integer, Byte>();
	private final static Map<Long, byte[]> rf_pool = new ConcurrentHashMap<Long, byte[]>();

	private class rfCalculator implements Runnable {

		private final int i;
		private final int j;

		public rfCalculator(int i, int j) {
			this.i = i;
			this.j = j;
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub
			PhasedDataCollection[] dc_i = dc[this.i];
			PhasedDataCollection[] dc_j = dc[this.j];
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
			try {
				//long tic = System.nanoTime();
				for(int k=0; k<m; k++) 
					for(int s=0; s<n; s++) 
						if(dc_i[k]!=null && dc_j[s]!=null)
							//System.out.println(phases_i[k]);
							//System.out.println(phases_j[s]);
							rf_all[k*m+s] = calcRFs(k, s);
				//scale_times[0] += System.nanoTime()-tic;

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
				w = "NULL\t"+StatUtils.min(rf)+"\t";
				for(int k=0; k<rf.length; k++) w += rf[k]+"\t";
				w += contig+"\t"+contig2
						+"\t"+contigs.get(contig)+"\t"+contigs.get(contig2)+"\n";
				rfMinimumWriter.write(w);
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(1);
			}
		}

		private double[] calcRFs(int phase, int phase2) {
			// TODO Auto-generated method stub
			PhasedDataCollection dc_i = dc[i][phase];
			int[] key4 = new int[4], key2 = new int[4];
			long key;
			int[][] hash_i = new int[2][nF1],
					hash_j = new int[2][nF1];

			double[] rf = new double[4];
			//long tic = System.nanoTime();
			PhasedDataCollection dc_j = dc[j][phase2].clone();
			//scale_times[6] += System.nanoTime()-tic;
			boolean[][][][] data_j = dc_j.data;

			for(int p=0; p<2; p++) {

				boolean[][][] data_j_p = data_j[p];
				hash_i = dc_i.hashcode[p];

				double[][] rf_all_p = new double[factorial_ploidy][4];
				for(int f=0; f<factorial_ploidy; f++) {

					//tic = System.nanoTime();
					next(data_j_p, f, hash_j);
					//scale_times[7] += System.nanoTime()-tic;

					for(int n=0; n<nF1; n++) {
						key4[0] = hash_i[0][n];
						key4[1] = hash_i[1][n];
						key4[2] = hash_j[0][n];
						key4[3] = hash_j[1][n];

						//tic = System.nanoTime();
						key = ((((((long) key4[0]
								<<shift_bits2)+key4[1])
								<<shift_bits2)+key4[2])
								<<shift_bits2)+key4[3];
						//scale_times[2] += System.nanoTime()-tic;

						
						if(rf_pool.containsKey(key)) {
							//tic = System.nanoTime();
							add(rf_all_p[f], rf_pool.get(key));
							//scale_times[4] += System.nanoTime()-tic;
						} else {
							//tic = System.nanoTime();
							add(rf_all_p[f], key4, key2, key);
							//scale_times[5] += System.nanoTime()-tic;	
					
						}
					}
				}
				add(rf, min(rf_all_p));
			}
			return divide(rf, nF1*ploidy);
		}

		private void add(double[] ds, byte[] bs) {
			// TODO Auto-generated method stub
			for(int i=0; i<4; i++) ds[i] += bs[i];
		}

		private void add(double[] ds, int[] key4, int[] key2, long key) {
			// TODO Auto-generated method stub
			key2[0] = (key4[0]<<shift_bits2)+key4[2];
			key2[1] = (key4[0]<<shift_bits2)+key4[3];
			key2[2] = (key4[1]<<shift_bits2)+key4[2];
			key2[3] = (key4[1]<<shift_bits2)+key4[3];
			byte[] bsall = new byte[4];
			for(int i=0; i<4; i++)
				bsall[i] = rf_factory.get(key2[i]);
			add(ds, bsall);
			rf_pool.put(key, bsall);
		}

		private double sum(byte[] bs) {
			// TODO Auto-generated method stub
			double sum = 0;
			for(byte b : bs) sum += b;
			return sum;
		}

		private double min(byte[] bs) {
			// TODO Auto-generated method stub
			byte min = Byte.MAX_VALUE;
			for(byte b : bs)
				if(b<min) min = b;
			return min;
		}

		private void next(boolean[][][] data, int f, int[][] hash) {
			// TODO Auto-generated method stub
			final int i = johnson_trotter_permutation[f][0],
					i2 = johnson_trotter_permutation[f][1];
			boolean[] tmp;
			boolean[][] tmps;
			int key;
			for(int m=0; m<2; m++) {
				tmps = data[m];
				tmp = tmps[i];
				tmps[i] = tmps[i2];
				tmps[i2] = tmp;
				for(int n=0; n<nF1; n++) {
					key = 0;
					for(int p=0; p<ploidy; p++)
						key=(key<<mask_length)+(tmps[p][n] ? 1 : 0);
					hash[m][n] = key;
				}
			}
		}

		private double[] divide(double[] ds, double deno) {
			// TODO Auto-generated method stub
			for(int i=0; i<4; i++) ds[i] /= deno;
			return ds;
		}

		private void divide(double[][][] rfs, int deno) {
			// TODO Auto-generated method stub
			for(int i=0; i<rfs.length; i++) 
				for(int j=0; j<rfs[i].length; j++)
					for(int k=0; k<rfs[i][j].length; k++)
						rfs[i][j][k] /= deno;
		}

		private void print(double[][] ds) {
			// TODO Auto-generated method stub
			//AnsiConsole.systemInstall();

			for(int i=0; i<ds.length; i++) {
				double m = StatUtils.min(ds[i]);
				for(int j=0; j<ds[i].length; j++) { 
					if(ds[i][j]==m) { 
						//System.out.print(ansi().fg(RED).
						//		a(formatter.format(ds[i][j])+"\t").reset());
						System.out.print(formatter.format(ds[i][j])+"\t");
						System.out.flush();
					} else {
						System.out.print(formatter.format(ds[i][j])+"\t");
						System.out.flush();
					}
				}
				System.out.println();
				System.out.flush();
			}
			System.out.println();
			System.out.flush();
			//AnsiConsole.systemUninstall();
		}

		private double[] search(double[][][] rfs) {
			// TODO Auto-generated method stub

			//print(rfs[0]);
			//print(rfs[1]);

			double[] rf = new double[rfs[0][0].length];
			for(int i=0; i<rf.length; i++) {
				double p1=0.0, p2=0.0;
				boolean shift1 = true, shift2 = true;
				while (true) {
					if(shift1) {
						p1 = rfs[0][0][i];
						for(int j=0; j<factorial_ploidy; j++) 
							if(rfs[0][j][i]<p1) 
								p1 = rfs[0][j][i];
						for(int x=0; x<rfs[0].length; x++)
							if(rfs[0][x][i]==p1)
								rfs[0][x][i] = Double.POSITIVE_INFINITY;
					}
					if(shift2) {
						p2 = rfs[1][0][i];
						for(int j=0; j<factorial_ploidy; j++) 
							if(rfs[1][j][i]<p2) 
								p2 = rfs[1][j][i];
						for(int x=0; x<rfs[0].length; x++)
							if(rfs[1][x][i]==p2)
								rfs[1][x][i] = Double.POSITIVE_INFINITY;
					}
					double d = Math.abs(p1-p2);
					if( Math.abs(d)>.1) {
						if(p1>p2) {
							shift2 = true;
							shift1 = false;
						} else {
							shift1 = true;
							shift2 = false;
						}
					} else if(p1==Double.POSITIVE_INFINITY ||
							p2==Double.POSITIVE_INFINITY) { 
						break;
					} else {
						break;
					}
				}
				rf[i] = p1+p2;
				if(rf[i]>.5) rf[i] = .5;
			}
			return rf;
		}

		private void add(double[] ds, double[] ds2) {
			// TODO Auto-generated method stub
			for(int i=0; i<4; i++) ds[i] += ds2[i];
		}

		private double[] min2(double[][] ds) {
			// TODO Auto-generated method stub
			double[] dss = new double[ds[0].length];
			Arrays.fill(dss, Double.POSITIVE_INFINITY);
			for(int i=0; i<ds.length; i++)
				for(int j=0; j<ds[i].length; j++)
					if(ds[i][j]<dss[j])
						dss[j] = ds[i][j];
			return dss;
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


		private double[] max(double[][] ds) {
			// TODO Auto-generated method stub
			double tot=Double.NEGATIVE_INFINITY, 
					rf=Double.NEGATIVE_INFINITY;
			double f;
			int r = -1;
			for(int i=0; i<ds.length; i++) 
				if( (f=StatUtils.max(ds[i]))>rf ||
						f==rf && StatUtils.sum(ds[i])>tot ) {
					rf = f;
					tot = StatUtils.sum(ds[i]);
					r = i;
				}
			return ds[r];
		}
	}
	
	private static BufferedWriter rfMinimumWriter;

	private void calcRFsForAll2(String outputFilePath) throws IOException, InterruptedException {

		File folder = new File(in_haps);
		File[] listFiles = folder.listFiles();
		nF1 = 0;

		for(File file:listFiles) {
			String name = file.getName();
			if( name.startsWith(expr_id) ) {
				if(nF1<1) {
					try {
						InputStream is = this.getInputStream(
								file.getAbsolutePath(),
								"phasedStates");
						if( is!=null ) {
							BufferedReader br = Utils.getBufferedReader(is);
							int n = 0;
							String l;

							while( (l=br.readLine())!=null ) 
								if(l.startsWith("#")) n++;

							nF1 = (n/Constants._ploidy_H)-2;
							myLogger.info(nF1+" F1 samples in the experiment.");
							br.close();
						}
						is.close();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			if(nF1>0) break;
		}

		Map<String, List<File>> map = new HashMap<String, List<File>>();
		List<File> list;
		String[] s;
		for(File file:listFiles) {
			String name = file.getName();
			if( name.startsWith(expr_id) ) {
				name = name.replace(expr_id,"experiment");
				s = name.split("\\.");

				if(map.get(s[1])==null) {
					list = new ArrayList<File>();
					list.add(file);
					map.put(s[1], list);
				} else{
					map.get(s[1]).add(file);
				}
			}
		}
		
		String[] keys = new String[map.keySet().size()];
		map.keySet().toArray(keys);

		for(int i=0; i<keys.length; i++) {
			List<File> files = map.get(keys[i]);
			executor.submit(new FileExtraction(
					files.toArray(new File[files.size()])));
		}
		executor.shutdown();
		executor.awaitTermination(365, TimeUnit.DAYS);

		this.initial_thread_pool();
		String[] scaff_all = new String[fileObj.keySet().size()];
		fileObj.keySet().toArray(scaff_all);
		dc = new PhasedDataCollection[scaff_all.length][best_n];
		for(int i=0; i<scaff_all.length; i++) {
			Set<FileObject> files = fileObj.get(scaff_all[i]);
			executor.submit(new FileLoader(scaff_all[i],
					files.toArray(new FileObject[files.size()]),
					i));
		}
		executor.shutdown();
		executor.awaitTermination(365, TimeUnit.DAYS);

		System.out.println(map.keySet().size());
		System.out.println("["+Utils.getSystemTime()+"] LOADING FILES DONE.");
		System.out.println("["+Utils.getSystemTime()+"] READING LOG LIKELIHOOD DONE.");

		StringBuilder os = new StringBuilder();
		for(int i=0; i<dc.length; i++) {
			os.append("##");
			os.append(dc[i][0]==null ? "null" : dc[i][0].file);
			for(int j=1; j<dc[i].length; j++) {
				os.append(" ");
				os.append(dc[i][j]==null ? "null" : dc[i][j].file);
			}
			os.append("\n");
		}

		this.initial_thread_pool();
		rfMinimumWriter = Utils.getBufferedWriter(outputFilePath+".min");
		long start = System.nanoTime();
		for(int i=0; i<dc.length; i++) 
			for(int j=i+1; j<dc.length; j++) 
				//executor.submit(new rfCalculator(i, j));
				new rfCalculator(i, j).run();
		executor.shutdown();
		executor.awaitTermination(365, TimeUnit.DAYS);
		long end = System.nanoTime();
		System.out.println("elapsed, "+(end-start)+"; R pool key-value pair, "+rf_pool.size());
		Utils.print(scale_times);
		rfMinimumWriter.close();


		this.initial_thread_pool();
		mapWriter = Utils.getBufferedWriter(outputFilePath+".map");
		mapWriter.write("##id\tkosambi\thaldane\tkSum\thSum\tkAll\thAll\n");
		for(int i=0; i<dc.length; i++) 
			executor.submit(new mapCalculator(i));
		executor.shutdown();
		executor.awaitTermination(365, TimeUnit.DAYS);
		mapWriter.close();

		System.out.println("["+Utils.getSystemTime()+"] DONE.");
	}
}



























