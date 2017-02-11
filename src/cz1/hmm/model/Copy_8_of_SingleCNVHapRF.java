package cz1.hmm.model;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
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

public class Copy_8_of_SingleCNVHapRF {	

	private static final Runtime instance = Runtime.getRuntime();
	private static final double kb = 1024;
	private static final double mb = 1024*1024;
	private static final double gb = 1024*1024*1024;

	private static long[] scale_times = new long[10]; 
	
	public static void main2(String[] args) {
		ci("C:\\Users\\chenxi.zhou\\Desktop\\genetic mapping writing-up\\"
				+ "meta_final\\rf_for_r_tetrasim192.ci.gz",
				"C:\\Users\\chenxi.zhou\\Desktop\\genetic mapping writing-up\\"
						+ "meta_final\\MAP\\rf_for_r_tetrasim192.min",
						"C:\\Users\\chenxi.zhou\\Desktop\\genetic mapping writing-up\\"
								+ "meta_final\\MAP\\rf_for_r_tetrasim192.ci");
	}

	public static void main(String[] args) 
			throws IOException, InterruptedException {

		//double[] rf = calcRFs("C:\\Users\\chenxi.zhou\\Desktop\\putty\\aaaaa\\"
		//		+ "scaffold_chrom8_3430740_6654037.txt",
		//		"C:\\Users\\chenxi.zhou\\Desktop\\putty\\aaaaa\\"
		//		+ "scaffold_chrom8_7718111_10260293.txt");
		//IO.print(rf);
		//rf = calcRFs("C:\\Users\\chenxi.zhou\\Desktop\\genetic mapping writing-up\\"
		//		+ "metafile\\phased\\0.0\\scaffold_chrom8_3430740_6654037.txt",
		//		"C:\\Users\\chenxi.zhou\\Desktop\\\\genetic mapping writing-up\\"
		//		+ "metafile\\phased\\0.0\\scaffold_chrom8_7718111_10260293.txt");
		//IO.print(rf);
		//if(args[3].equals("simulation"))
		//	calcRFsForAll(args[0],args[1],args[2]);
		//else
		//	calcRFsForAll2(args[0],args[1],args[2],args[3]);
		//calcRFsForGroup("C:\\Users\\chenxi.zhou\\Desktop\\putty\\analysis\\rf_for_r.txt",
		//		"C:\\Users\\chenxi.zhou\\Desktop\\putty\\analysis\\scaffolding-gt20.txt",
		//		"C:\\Users\\chenxi.zhou\\Desktop\\contigs.txt",
		//		"C:\\Users\\chenxi.zhou\\Desktop\\putty\\analysis\\group-distance\\");

		// create the command line parser
		CommandLineParser parser = new PosixParser();

		// create the Options
		Options options = new Options();
		options.addOption( "e", "experiment", true, "experiment name." );
		options.addOption( "w", "workspace", true, "directory contains input files." );
		options.addOption( "o", "output", true, "output file prefix." );
		options.addOption( "c", "contig-file", true, "file contaions contigs information.");
		options.addOption( "p", "parent-obj",true, "parents; dilimited by a colon(:)." );
		options.addOption( "P", "ploidy",true, "ploidy of the genome." );
		options.addOption( "n", "n-best", true, "the best n run is used.");
		options.addOption( "t", "thread", true, "number threads.");
		options.addOption( "s", "skew", true, "skew of haplotype frequencies.");
		options.addOption( "d", "drop", true, "threshold do drop a scaffold.");
		options.addOption( "T", "test", true, "goodness-of-fit test.");
		options.addOption( "C", "scaff-only", true, "calculate RFs only for specified scaffolds.");
		String experiment=null, workspace=null, output=null, contigFile=null, 
				parentObj=null, goodnessOfFitTest="fraction";
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
			if(line.hasOption("c")) {
				contigFile = line.getOptionValue("c");
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

		Copy_8_of_SingleCNVHapRF sgRF = new Copy_8_of_SingleCNVHapRF(
				workspace,experiment,ploidy,parents,t,s,d,goodnessOfFitTest,n);
		sgRF.calcRFsForAll2(output, contigFile);
	}

	private final static int NUM_CORES = 
			Runtime.getRuntime().availableProcessors();
	private static String workspace;
	private static String experiment;
	private static String[] parents;
	private static int ploidy;
	private static int ploidy2;
	private static int mask_length;
	private static int shift_bits2;
	private static int shift_bits4;
	private static int THREADS = NUM_CORES-1;
	private static int drop_thres;
	private static double skew_thres;
	private static double[] probs_uniform;
	private static String goodness_of_fit;
	private static int best_n;
	private final static Set<String> scaff_only = new HashSet<String>();	

	public Copy_8_of_SingleCNVHapRF (String workspace_, 
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
		mask_length = (int) Math.ceil(
				Math.log(ploidy2+1)/Math.log(2));
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

	public static void calcRFsForGroup(String rfFilePath, 
			String groupFilePath,
			String contigFilePath,
			String outputFilePath) throws IOException {
		String[] s;
		Map<String, Integer> contigs = new HashMap<String, Integer>();
		BufferedReader br = getBufferedReader(contigFilePath);
		int c=0;
		String line;
		while( (line=br.readLine())!=null ) {
			s = line.split("\\s+");
			contigs.put(s[0], c++);
		}
		br.close();

		double[][] rf = new double[363][363];
		br = getBufferedReader(rfFilePath);
		while( (line=br.readLine())!=null ) {
			s = line.split("\\s+");
			rf[Integer.parseInt(s[6])-1][Integer.parseInt(s[7])-1] = 
					Double.parseDouble(s[1]);
			rf[Integer.parseInt(s[7])-1][Integer.parseInt(s[6])-1] = 
					Double.parseDouble(s[1]);
		}
		br.close();

		String experiment = new File(groupFilePath).getName();
		br = getBufferedReader(groupFilePath);
		while( (line=br.readLine())!=null ){
			s = line.split("-");
			if(s.length>1) {
				for(int i=0; i<s.length; i++) s[i] = s[i].trim();
				s = sort(s);
				String chr = s[0].split("_")[1];
				BufferedWriter bw = getBufferedWriter(outputFilePath
						+chr+"-"+experiment);
				bw.write("#");
				for(int i=0; i<s.length-1;i++) bw.write(s[i]+",");
				bw.write(s[s.length-1]+"\n");
				double[][] r = new double[s.length][s.length];
				for(int i=0; i<s.length; i++) 
					for(int j=i+1; j<s.length; j++) { 
						r[j][i] = rf[contigs.get(s[j])]
								[contigs.get(s[i])];
						r[i][j] = r[j][i];
					}
				for(int i=0; i<s.length; i++) {
					for(int j=0; j<s.length; j++) {
						bw.write(r[i][j]+" ");
					}
					bw.write("\n");
				}

				if(s.length>1 && s.length<11) {
					Integer[] ele = new Integer[s.length];
					for(int i=0; i<s.length; i++) ele[i] = i;
					ArrayList<List<Integer>> permutations  = 
							Permutation.permutation(ele);
					double MIN = Double.MAX_VALUE;
					for(int i=0; i<permutations.size(); i++) {
						double min = 0;
						List<Integer> permutation = permutations.get(i);
						for(int j=0; j<permutation.size()-1; j++)
							min += r[permutation.get(j)][permutation.get(j+1)];
						if(min <= MIN) {
							System.out.println("["+min+"]");
							for(int j=0; j<permutation.size()-1;j++)
								System.out.print(s[permutation.get(j)]+"-");
							System.out.println(s[permutation.get(permutation.size()-1)]+"\n");
							MIN = min;
						}
					}
				}

				bw.close();
			}
		}
		br.close();
	}

	private static void ci (final String ci_in, final String rf_in, final String out) {
		try {
			final BufferedReader br_ci = Utils.getBufferedReader(ci_in);
			String[] s;
			String line = br_ci.readLine();
			Map<String, double[]> cis = new HashMap<String, double[]>();
			Map<Integer, String> scaff = new HashMap<Integer, String>();
			int k = 0;
			while( line!=null) {
				if(line.startsWith("##")) {
					if(line.startsWith("##null")) 
						k++;
					else
						scaff.put(k++, line.split("\\.")[1]);
					line = br_ci.readLine();
					continue;
				}
				if(line.startsWith("#")) {
					int scf = Integer.parseInt(line.replace("#",""));
					double[][] rfs = new double[30*29/2][4];
					for(int i=0; i<4; i++) {
						line = br_ci.readLine();
						s = line.split("\\s+");
						for(int j=0; j<rfs.length; j++)
							rfs[j][i] = Double.parseDouble(s[j]);
					}

					double[] mins = new double[rfs.length];
					for(int i=0; i<mins.length; i++)
						mins[i] = StatUtils.min(rfs[i]);
					double std = Math.sqrt(StatUtils.variance(mins));
					double min = StatUtils.min(mins);
					double mea = StatUtils.mean(mins);
					double med = StatUtils.percentile(mins, 0.1);
					cis.put(scaff.get(scf), new double[]{std, min, mea, med});
					line = br_ci.readLine();
				}
			}
			br_ci.close();

			final BufferedWriter bw_out = Utils.getBufferedWriter(out);
			final BufferedReader br_rf = Utils.getBufferedReader(rf_in);
			while( (line=br_rf.readLine())!=null ) {
				s = line.split("\\s+");
				double[] ci_1 = cis.get(s[6]),
						ci_2 = cis.get(s[7]);
				double[] ci = new double[4];
				for(int i=0; i<ci.length; i++)
					ci[i] = ci_1[i]>ci_2[i] ? ci_1[i] : ci_2[i];
					bw_out.write(line+"\t"+cat(ci,"\t")+"\n");
			}
			br_rf.close();
			bw_out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static double[] min(double[][] ds) {
		// TODO Auto-generated method stub
		double tot=Double.POSITIVE_INFINITY, 
				rf=Double.POSITIVE_INFINITY;
		double f;
		int r = -1;
		for(int i=0; i<ds.length; i++) 
			if( (f=StatUtils.min(ds[i]))<rf ||
					f==rf && StatUtils.sum(ds[i])<tot ) {
				rf = f;
				tot = StatUtils.sum(ds[i]);
				r = i;
			}
		return ds[r];
	}

	private static String[] sort(String[] strArray) {
		// TODO Auto-generated method stub
		List<ContigString> list = new ArrayList<ContigString>();
		for(String s : strArray)
			list.add(new ContigString(s));
		Collections.sort(list);
		String[] order = new String[list.size()];
		for(int j=0; j<list.size(); j++)
			order[j] = list.get(j).contig;

		return order;
	}

	private static class ContigString implements Comparable<ContigString> {
		String contig;

		public ContigString(String contig) {
			this.contig = contig;
		}

		@Override
		public int compareTo(ContigString o) {
			// TODO Auto-generated method stub
			return Integer.parseInt(this.contig.split("_")[2])-
					Integer.parseInt(o.contig.split("_")[2]);
		}
	}


	private static char[][] haps = null;
	private static int factorial_ploidy = -1;;
	private final int[] hap_index = new int[256];
	private static int[][] johnson_trotter_permutation;

	private void setHaps() {
		factorial_ploidy = (int) Permutation.factorial(ploidy);
		List<List<Character>> haps_list = new ArrayList<List<Character>>();
		List<Character> hap = new ArrayList<Character>();
		for(int i=1; i<=ploidy; i++) {
			hap.add( i<10 ? Character.toChars('0'+i)[0] : 
				Character.toChars('a'+i-10)[0] );
		}
		haps_list.add(hap);
		hap = new ArrayList<Character>();
		for(int i=ploidy+1; 
				i<=2*ploidy; i++) {
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
		List<List<Integer>> inner_order_list = 
				Permutation.permutation(ploidy2);
		inner_order_length = inner_order_list.size();
		inner_order = new int[inner_order_length][ploidy2];
		for(int i=0; i<inner_order_length; i++)
			for(int j=0; j<ploidy2; j++)
				inner_order[i][j] = inner_order_list.get(i).get(j);
		johnson_trotter_permutation = JohnsonTrotter.perm(ploidy);
		initialiseRFactory();
	}

	private void initialiseRFactory() {
		// TODO Auto-generated method stub
		List<List<Integer>> combs = Combination.combination(ploidy, ploidy2);
		List<List<Integer>> perms = Permutation.permutation(ploidy2);
		for(int c=0; c<combs.size(); c++) {
			List<Integer> comb = combs.get(c);
			for(List<Integer> perm : perms) {
				byte[] bytes = new byte[ploidy];
				for(int i=0; i<ploidy2; i++)
					bytes[comb.get(i)] = (byte)(perm.get(i)+1);
				int key = hashcode(bytes);
				for(int c2=0; c2<combs.size(); c2++) {
					List<Integer> comb2  = combs.get(c2);
					for(List<Integer> perm2 : perms) {
						byte[] bytes2 = new byte[ploidy];
						for(int i=0; i<ploidy2; i++)
							bytes2[comb2.get(i)] = (byte)(perm2.get(i)+1);
						int key2 = hashcode(bytes2);
						byte[] rf = recombinationFreq(bytes, bytes2);
						rf_factory.put((key<<shift_bits2)+key2, rf);
					}
				}
			}
		}
	}

	private byte[] recombinationFreq(byte[] bytes, 
			byte[] bytes2) {
		// TODO Auto-generated method stub
		byte[] bytes_hap = new byte[ploidy2],
				bytes_hap2 = new byte[ploidy2];
		for(byte i=0; i<ploidy; i++)
			if(bytes[i]!=0)
				bytes_hap[bytes[i]-1] = i;
		for(byte i=0; i<ploidy; i++)
			if(bytes2[i]!=0)
				bytes_hap2[bytes2[i]-1] = i;

		byte[] rf = new byte[inner_order.length];
		for(int i=0; i<rf.length; i++) {
			byte r = 0;
			for(int j=0; j<ploidy2; j++) {
				if(bytes_hap[j]!=bytes_hap2[inner_order[i][j]])
					r++;
			}
			rf[i] = r;
		}
		return rf;
	}

	private int hashcode(byte[] bytes) {
		// TODO Auto-generated method stub
		int hash = 0;
		for(int i=0; i<ploidy; i++)
			hash = (hash<<mask_length)+bytes[i];
		return hash;
	}

	private int[][] inner_order = null; 
	private int inner_order_length = 0;
	private final int[][] order = new int[][] {
			new int[]{0,0},
			new int[]{0,1},
			new int[]{1,0},
			new int[]{1,1}};

	private class FileObject {
		private final File file;
		private final String[] markers; 
		private final int[] start_end_position;

		public FileObject(File file, 
				String[] markers,
				int[] start_end_position) {
			this.file = file;
			this.markers = markers;
			this.start_end_position = start_end_position;
		}
	}

	private final static Map<String, Set<FileObject>> fileObj = 
			new HashMap<String, Set<FileObject>>();
	private final static Object lock = new Object();

	private class FileExtraction implements Runnable {
		private final File[] files;

		public FileExtraction(File[] files) {
			this.files = files;
		}
		@Override
		public void run() {
			// TODO Auto-generated method stub
			List<String> scaff_all = new ArrayList<String>();
			String[][] markers = null;
			int[][] start_end_position = null;
			int scaff_n = 0;

			try {
				BufferedReader br = Utils.getBufferedReader(this.files[0].getAbsolutePath()+
						"/snp_"+experiment+".txt");
				List<List<String>> markers_all = new ArrayList<List<String>>();

				String marker = br.readLine().split("\\s+")[3];
				String scaff_prev = marker.replaceAll("_[0-9]{1,}$", ""),
						scaff;
				scaff_all.add(scaff_prev);
				String line;
				markers_all.add(new ArrayList<String>());
				int n_=0;
				markers_all.get(n_).add(marker);
				while( (line=br.readLine())!=null ) {
					marker = line.split("\\s+")[3];
					scaff = marker.replaceAll("_[0-9]{1,}$", "");
					if(scaff.equals(scaff_prev))
						markers_all.get(n_).add(marker);
					else {
						markers_all.add(new ArrayList<String>());
						n_++;
						markers_all.get(n_).add(marker);
						scaff_prev = scaff;
						scaff_all.add(scaff_prev);
					}
				}
				br.close();
				int cuv = 0;
				scaff_n = scaff_all.size();
				markers = new String[scaff_n][];
				start_end_position = new int[scaff_n][2];
				for(int i=0; i<scaff_n; i++) {
					markers[i] = new String[markers_all.get(i).size()];
					markers_all.get(i).toArray(markers[i]);
					int s = Integer.parseInt(markers[i][0].
							replaceAll(".*[^\\d](\\d+).*", "$1")),
							e = Integer.parseInt(markers[i][1].
									replaceAll(".*[^\\d](\\d+).*", "$1"));
					if(s<=e) {
						start_end_position[i][0] = cuv;
						cuv += markers[i].length;
						start_end_position[i][1] = cuv-1;
					} else {
						start_end_position[i][1] = cuv;
						cuv += markers[i].length;
						start_end_position[i][0] = cuv-1;
					}
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			for(int i=0; i<this.files.length; i++) 
				for(int j=0; j<scaff_n; j++) {
					String scaff = scaff_all.get(j); 
					synchronized(lock) {
						if(!fileObj.containsKey(scaff))
							fileObj.put(scaff, new HashSet<FileObject>());
						fileObj.get(scaff).add(new FileObject(
								this.files[i],
								markers[j],
								start_end_position[j]));
					}
				}
		}			
	}

	private class FileLoader implements Runnable {
		private final String id;
		private final FileObject[] files;
		private final int i;

		public FileLoader(String id, FileObject[] files, int i) {
			this.id = id;
			this.files = files;
			this.i = i;
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub

			double[] ll = new double[this.files.length];
			for(int k=0; k<ll.length; k++) 
				ll[k]=Double.NEGATIVE_INFINITY;
			int marker_n = this.files[0].markers.length;

			if( marker_n<2 ) {
				System.err.println("warning: "+
						this.id +" #marker is less than 2.");
				return;
			}

			for(int k=0; k<this.files.length; k++) {
				File file = this.files[k].file;
				if(!new File(file.getAbsolutePath()+
						"/phasedStates/"+experiment+".txt").exists()) {
					System.err.println("warning: "+
							file.getName()+
							" exsits, but phased states do not.");
					continue;
				}
				try {
					BufferedReader br = 
							Utils.getBufferedReader(file.getAbsolutePath()+
									"/phasedStates/"+experiment+".txt");
					br.readLine();
					String marker_str = br.readLine();
					br.close();
					if( marker_str==null ) {
						System.err.println("warning: "+
								file.getName()+
								" exists, but phased states are NULL.");
						continue;
					}
					double frac = marker_n/Double.parseDouble(marker_str);

					String res = file.getAbsolutePath()+"/stderr_true";
					//IO.println(res);
					File _res_hmm = new File(res);
					String line;
					if( _res_hmm.exists() ) {
						BufferedReader br2 = getBufferedReader(res);
						String lprob=null;
						while( (line=br2.readLine()) !=null) {
							if(line.startsWith("log"))
								lprob = line;
						}
						br2.close();
						if( lprob!=null ) {
							String[] s0 = lprob.split("\\s+");
							ll[k] = Double.parseDouble(s0[3])*frac;
						}
					} else {
						BufferedReader br2 = 
								getBufferedReader(file.getAbsolutePath()+
										"/phasedStates/"+
										experiment+".txt");
						String lprob=br2.readLine();
						br2.close();
						if( lprob!=null ) 
							ll[k] = Double.parseDouble(lprob)*frac;
						System.out.println(ll[k]);
					}
				} catch (IOException e) {
					e.printStackTrace();
					System.exit(1);
				}
			}

			int[] maxN = maxN(ll);
			StringBuilder oos = new StringBuilder();
			boolean[] drop = new boolean[maxN.length];
			int dropped = 0;
			for(int k=0; k<maxN.length; k++) {

				int[] haps_observed = readHaplotypes(maxN[k]);
				/**
					br_states = IO.getBufferedReader(this.files[maxN[k]]+
							"/phasedStates/"+experiment+".txt");
					String line_states;
					String[] s_states;
					double[] phases = new double[ploidy*2];
					while( (line_states=br_states.readLine())!=null ) {
						if(!line_states.startsWith("#"))
							continue;
						s_states = line_states.split("\\s+");
						phases[s_states[s_states.length-1].charAt(0)-'1'] += 1.0;
					}
				 ***/

				long[] observed;
				double p;
				switch(goodness_of_fit) {
				case "fraction":
					double[] phases = new double[ploidy*2];
					for(int z=0; z<phases.length; z++) 
						phases[z] = (double) haps_observed[z];
					double expected = StatUtils.sum(phases)/ploidy/2;
					double maf = StatUtils.max(phases)/expected, 
							mif = StatUtils.min(phases)/expected;
					if( maf>1/skew_thres || mif<skew_thres) {
						System.out.println(this.files[maxN[k]].file.getName()+
								" was dropped due to large haploptype frequency variance. (" +
								cat(phases, ",") +")");
						drop[k] = true;
					}
					if(drop[k]) dropped++;
					oos.append("["+(drop[k]?"drop](maf,":"keep](maf,")+maf+";mif,"+mif+") "+
							cat(haps_observed,",")+"\t"+
							this.files[maxN[k]].file.getName()+"\n");
					break;
				case "chisq":
					observed = new long[ploidy*2];
					for(int z=0; z<observed.length; z++) 
						observed[z] = (long) haps_observed[z];
					p = new ChiSquareTest().chiSquareTest(probs_uniform, observed);
					if(p<skew_thres) drop[k] = true;
					if(drop[k]) dropped++;
					oos.append("["+(drop[k]?"drop](p,":"keep](p,")+formatter.format(p)+") "+
							cat(haps_observed,",")+"\t"+
							this.files[maxN[k]].file.getName()+"\n");
					break;
				case "gtest":
					observed = new long[ploidy*2];
					for(int z=0; z<observed.length; z++) 
						observed[z] = (long) haps_observed[z];
					p = new GTest().gTest(probs_uniform, observed);
					if(p<skew_thres) drop[k] = true;
					if(drop[k]) dropped++;
					oos.append("["+(drop[k]?"drop](p,":"keep](p,")+formatter.format(p)+") "+
							cat(haps_observed,",")+"\t"+
							this.files[maxN[k]].file.getName()+"\n");
					break;
				default:
					System.err.println("Goodness-of-fit test should be fraction, chisq or gTest.");
					System.exit(1);
				}
			}
			System.out.print(oos.toString());
			System.err.println(this.id+" - dropped "+dropped);
			if( drop.length-dropped<drop_thres ) {
				System.err.println("Scaffold "+this.id+" dropped.");
			} else {
				int kk=0;
				for(int k=0; k<drop.length; k++) {
					if(!drop[k]) {
						FileObject fobj = this.files[maxN[k]];
						dc[i][kk] = new PhasedDataCollection(
								fobj.file.getName(),
								fobj.markers,
								fobj.start_end_position);
						kk++;
					}
					if(kk>=best_n) break;
				}
			}
		}

		private int[] readHaplotypes(final int i) {
			// TODO Auto-generated method stub
			try {
				BufferedReader br_states = Utils.getBufferedReader(
						this.files[i].file+
						"/phasedStates/"+experiment+".txt");;
						String line, stateStr;
						String[] s;
						int[] haps_observed = new int[ploidy*2];
						while( (line=br_states.readLine())!=null ) {
							if(!line.startsWith("#")) continue;
							s = line.split("\\s+|:");
							stateStr = s[s.length-1];
							for(char h : stateStr.toCharArray())
								haps_observed[h>'9'?(h-'a'+9):(h-'1')]++;
						}
						br_states.close();
						return haps_observed;
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				System.exit(1);
			}
			return null;
		}

		private int[] maxN(double[] ll) {
			double[] ll0 = Arrays.copyOf(ll, ll.length);
			int n = ll.length;//best_n_phases[0].length;
			int[] maxN = new int[n];
			Arrays.fill(maxN, -1);
			for(int k=0; k<n; k++) {
				if(k>=ll0.length) return maxN;
				int p = 0;
				double e = Double.NEGATIVE_INFINITY;
				for(int s=0; s<ll0.length; s++)
					if(ll0[s]>e) {
						e = ll0[s];
						p = s;
					}
				maxN[k] = p;
				ll0[p] = Double.NEGATIVE_INFINITY;
			}
			return maxN;
		}

		private int[] maxN(double[] ll, int N) {
			double[] ll0 = Arrays.copyOf(ll, ll.length);
			int[] maxN = new int[N];
			Arrays.fill(maxN, -1);
			for(int k=0; k<N; k++) {
				if(k>=ll0.length) return maxN;
				int p = 0;
				double e = Double.NEGATIVE_INFINITY;
				for(int s=0; s<ll0.length; s++)
					if(ll0[s]>e) {
						e = ll0[s];
						p = s;
					}
				maxN[k] = p;
				ll0[p] = Double.NEGATIVE_INFINITY;
			}
			return maxN;
		}
	}

	private final static Map<Integer, byte[]> rf_factory = new HashMap<Integer, byte[]>();
	private final static Map<Long, byte[]> rf_pool = new ConcurrentHashMap<Long, byte[]>();

	private static int counter = 0;
	
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
				
				for(int k=0; k<m; k++) 
					for(int s=0; s<n; s++) 
						if(dc_i[k]!=null && dc_j[s]!=null)
							//System.out.println(phases_i[k]);
							//System.out.println(phases_j[s]);
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
				rfAllWriter.write(os.toString());
				
				double[] rf; 
				String w;
				rf = new double[4];
				for(int k=0; k<4; k++) {
					double[] rf_tmp = removeNEG(rf_all[k]);
					if(rf_tmp==null)
						rf[k] = -1;
					else
						rf[k] = median(rf_all[k]);
				}
				w = "NULL\t"+StatUtils.min(rf)+"\t";
				for(int k=0; k<rf.length; k++) w += rf[k]+"\t";
				w += contig+"\t"+contig2
						+"\t"+contigs.get(contig)+"\t"+contigs.get(contig2)+"\n";
				rfMedianWriter.write(w);

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

				rf = new double[4];
				for(int k=0; k<4; k++) {
					double[] rf_tmp = removeNEG(rf_all[k]);
					if(rf_tmp==null)
						rf[k] = -1;
					else
						rf[k] = StatUtils.mean(rf_tmp);
				}
				w = "NULL\t"+StatUtils.min(rf)+"\t";
				for(int k=0; k<rf.length; k++) w += rf[k]+"\t";
				w += contig+"\t"+contig2
						+"\t"+contigs.get(contig)+"\t"+contigs.get(contig2)+"\n";
				rfMeanWriter.write(w);
				
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(1);
			}
		}

		private double[] calcRFs(int phase, int phase2) {
			// TODO Auto-generated method stub
			long tic;
			
			PhasedDataCollection dc_i = dc[i][phase];
			int[] key4 = new int[4], key2 = new int[4];
			long key;
			int[][] hash_i = new int[2][nF1],
					hash_j = new int[2][nF1];

			double[] rf = new double[order.length];
			PhasedDataCollection dc_j = dc[j][phase2].clone();
			int[][][][] data_j = dc_j.data;
			
			for(int p=0; p<2; p++) {

				int[][][] data_j_p = data_j[p];
				hash_i = dc_i.hashcode[p];

				double[][] rf_all_p = new double[factorial_ploidy][4];
				for(int f=0; f<factorial_ploidy; f++) {
					
					tic = System.nanoTime();
					next(data_j_p, f, hash_j);
					scale_times[0] += System.nanoTime()-tic;
					
					for(int n=0; n<nF1; n++) {
						
						tic = System.nanoTime();
						key4[0] = hash_i[0][n];
						key4[1] = hash_i[1][n];
						key4[2] = hash_j[0][n];
						key4[3] = hash_j[1][n];
						key = ((((((long) key4[0]
								<<shift_bits2)+key4[1])
								<<shift_bits2)+key4[2])
								<<shift_bits2)+key4[3];
						scale_times[1] += System.nanoTime()-tic;
						
						if(rf_pool.containsKey(key)) {
							tic = System.nanoTime();
							add(rf_all_p[f], rf_pool.get(key));
							scale_times[2] += System.nanoTime()-tic;
						}
						else {
							tic = System.nanoTime();
							add(rf_all_p[f], key4, key2, key);
							scale_times[3] += System.nanoTime()-tic;
						}
					}
				}
				tic = System.nanoTime();
				add(rf, min(rf_all_p));
				scale_times[4] += System.nanoTime()-tic;
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
			byte[][] bsall = new byte[inner_order_length][4];
			for(int i=0; i<4; i++) {
				byte[] bs = rf_factory.get(key2[i]);
				for(int j=0; j<inner_order_length; j++)
					bsall[j][i] = bs[j];
			}
			double t=Double.POSITIVE_INFINITY, 
					rf=Double.POSITIVE_INFINITY;
			double f;
			int r = -1;
			for(int i=0; i<inner_order_length; i++) 
				if( (f=min(bsall[i]))<rf ||
						f==rf && sum(bsall[i])<t ) {
					rf = f;
					t = sum(bsall[i]);
					r = i;
				}		
			add(ds, bsall[r]);
			rf_pool.put(key, bsall[r]);
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

		private void next(int[][][] data, int f, int[][] hash) {
			// TODO Auto-generated method stub
			final int i = johnson_trotter_permutation[f][0],
					i2 = johnson_trotter_permutation[f][1];
			int[] tmp;
			int[][] tmps;
			int key;
			for(int m=0; m<2; m++) {
				tmps = data[m];
				tmp = tmps[i];
				tmps[i] = tmps[i2];
				tmps[i2] = tmp;
				for(int n=0; n<nF1; n++) {
					key = 0;
					for(int p=0; p<ploidy; p++)
						key=(key<<mask_length)+tmps[p][n];
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

	private double median(double[] ds) {
		// TODO Auto-generated method stub
		double[] ds0 = removeNEG(ds);
		if(ds0==null) return -1;
		Arrays.sort(ds0);
		int n = ds0.length;
		if (n % 2 == 0)
			return (ds0[n/2] + ds0[n/2-1])/2;
		else
			return ds0[n/2];
	}

	private double[] removeNEG(double[] ds) {
		List<Double> ds0 = new ArrayList<Double>();
		for(double d : ds)
			if(d>=0)
				ds0.add(d);
		if(ds0.isEmpty()) return null;
		return ArrayUtils.toPrimitive(
				ds0.toArray(new Double[ds0.size()]));

	}

	private class mapCalculator implements Runnable {

		private final int i;

		public mapCalculator(int i) {
			this.i = i;
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub
			double[][] kosambi_all = null;
			double[] kd_all = new double[dc[i].length];
			Arrays.fill(kd_all, -1);

			for(int k=0; k<dc[this.i].length; k++) {
				PhasedDataCollection dc_ik= dc[i][k];
				if(dc_ik !=null ) {
					double[] kosambi_ = calcGDs(workspace+"/"+dc_ik.file+"/phasedStates/"+experiment+".txt",
							ploidy,
							parents,
							nF1,
							dc_ik.start_end_position);
					if(kosambi_all==null) {
						kosambi_all = new double[dc[i].length][kosambi_.length];
						fill(kosambi_all, -1);
					}

					kosambi_all[k] = kosambi_;

					kd_all[k] = calcGD(workspace+"/"+dc[i][k].file+"/phasedStates/"+experiment+".txt",
							ploidy,
							parents,
							nF1,
							dc_ik.start_end_position,
							"kosambi");
				}
			}

			if(kosambi_all==null) return;

			double[] kosambi = new double[kosambi_all[0].length];
			kosambi_all = Algebra.transpose(kosambi_all);
			for(int k=0; k<kosambi.length; k++) 
				kosambi[k] = median(kosambi_all[k]);

			if(removeNEG(kd_all)==null || removeNEG(kosambi)==null)
				return;

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
		}

	}

	private class gdCalculator implements Runnable {
		private final int i;

		public gdCalculator(int i) {
			this.i = i;
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub
			try {
				PhasedDataCollection _0 = null;
				for(PhasedDataCollection f : dc[i]) {
					_0 = f;
					if(f!=null)
						break;
				}
				if(_0==null) return;


				String[] snp_id = _0.markers;
				int[] snp_pos = new int[snp_id.length];
				for(int ii=0; ii<snp_pos.length; ii++) 
					snp_pos[ii] = Integer.parseInt(snp_id[ii].
							replaceAll(".*[^\\d](\\d+).*", "$1"));

				int nSNP = snp_id.length;
				double[][] rfAll = new double[nSNP*(nSNP-1)/2]
						[dc[i].length];
				fill(rfAll, -1);
				for(int k=0; k<dc[i].length; k++)
					if(dc[i][k]!=null)
						calcGDsAll(workspace+"/"+dc[i][k].file+
								"/phasedStates/"+experiment+".txt",
								ploidy,
								parents,
								nF1,
								dc[i][k].start_end_position,
								rfAll,
								k);
				int c=0;
				StringBuilder os = new StringBuilder();
				for(int m=0; m<nSNP; m++)
					for(int n=m+1; n<nSNP; n++) {
						double rf = median(rfAll[c++]);
						if(rf!=-1) {
							os.append(snp_id[m]+"\t"+snp_id[n]+
									"\t"+Math.abs(snp_pos[m]-snp_pos[n])+
									"\t"+rf+"\n");
						}
					}
				gdWriter.write(os.toString());
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	private void fill(double[][] dss, 
			double d) {
		// TODO Auto-generated method stub
		for(double[] ds : dss) 
			Arrays.fill(ds, d);
	}

	private static PhasedDataCollection[][] dc;
	private static int nF1;
	private static Map<String, Integer> contigs;
	private static BufferedWriter rfMedianWriter;
	private static BufferedWriter rfMeanWriter;
	private static BufferedWriter rfMinimumWriter;
	private static BufferedWriter rfAllWriter;
	private static BufferedWriter mapWriter;
	private static BufferedWriter gdWriter;
	private static BufferedWriter ciWriter;
	private static NumberFormat formatter = new DecimalFormat("#0.000");

	private class PhasedDataCollection {
		private final String file;
		private final String[] markers;
		private final int[][][][] data;
		private final int[] start_end_position;
		private int[][][] hashcode;

		public PhasedDataCollection(String file,
				String[] markers, int[] start_end) {
			// TODO Auto-generated constructor stub
			this.file = file;
			this.markers = markers;
			if(markers.length!=
					Math.abs(start_end[0]-start_end[1])+1)
				throw new RuntimeException("!!!");
			this.data = this.data(start_end[0], start_end[1]);
			this.start_end_position = start_end;
			this.hash();
		}

		private void hash() {
			// TODO Auto-generated method stub
			this.hashcode = new int[2][2][nF1];
			for(int i=0; i<2; i++)
				for(int j=0; j<2; j++) {
					int[][] d = this.data[i][j];
					int[] code = this.hashcode[i][j];
					for(int k=0; k<nF1; k++) {
						int key = 0;
						for(int p=0; p<ploidy; p++)
							key = (key<<mask_length)+d[p][k];
						code[k] = key;
					}
				}
		}

		public PhasedDataCollection(String file,
				String[] markers,
				int[][][][] data,
				int[] start_end_position) {
			// TODO Auto-generated constructor stub
			this.file = file;
			this.markers = markers;
			this.data = data;
			this.start_end_position = start_end_position;
		}

		public PhasedDataCollection clone() {
			final int[][][][] data = new int[2][2][ploidy][nF1];
			for(int i=0; i<2; i++) 
				for(int j=0; j<2; j++)
					for(int k=0; k<ploidy; k++)
						data[i][j][k] = this.data[i][j][k].clone();
			return new PhasedDataCollection(this.file, 
					this.markers, data, 
					this.start_end_position);
		}

		private int[][][][] data(int start, int end) {
			// TODO Auto-generated method stub
			int[][][][] data = new int[2][2][ploidy][nF1];

			try {
				BufferedReader br = getBufferedReader(
						workspace+"/"+
								this.file+
								"/phasedStates/"+
								experiment+".txt");
				String line;
				String[] s;
				String stateStr;
				int n=0;

				while( (line=br.readLine())!=null ) {

					if(!line.startsWith("#")) continue;
					//if(skip++<2) continue;
					s = line.split("\\s+|:");
					if(Arrays.asList(parents).contains(s[2])) continue;

					stateStr = s[s.length-1];
					data[0][0][hap_index[stateStr.charAt(start)]][n] = 1;
					data[0][1][hap_index[stateStr.charAt(end)]][n] = 1;

					for(byte i=1; i<ploidy/2; i++) {
						line = br.readLine();
						s = line.split("\\s+");
						stateStr = s[s.length-1];
						data[0][0][hap_index[stateStr.charAt(start)]][n] = i+1;
						data[0][1][hap_index[stateStr.charAt(end)]][n] = i+1;
					}
					for(int i=0; i<ploidy/2; i++) {
						line = br.readLine();
						s = line.split("\\s+");
						stateStr = s[s.length-1];
						data[1][0][hap_index[stateStr.charAt(start)]][n] = i+1;
						data[1][1][hap_index[stateStr.charAt(end)]][n] = i+1;
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

		private Map[][] map() {
			// TODO Auto-generated method stub
			Map[][] map = new HashMap[2][2];
			for(int i=0; i<2; i++)
				for(int j=0; j<2; j++)
					map[i][j] = new HashMap();
			return map;
		}
	}

	private static ExecutorService executor = 
			Executors.newFixedThreadPool(THREADS);
	private void reset() {
		executor = Executors.newFixedThreadPool(THREADS);
	}

	private void calcRFsForAll2(String outputFilePath,
			String contigFilePath) throws IOException, InterruptedException {

		String[] s;
		contigs = new HashMap<String, Integer>();
		BufferedReader br = getBufferedReader(contigFilePath);
		int c=0;
		String line;
		while( (line=br.readLine())!=null ) {
			s = line.split("\\s+");
			contigs.put(s[0], ++c);
			if(!s[0].startsWith("chr"))
				contigs.put("chr"+s[0], c);
		}
		br.close();

		File folder = new File(workspace);
		File[] listFiles = folder.listFiles();
		nF1 = 0;

		for(File file:listFiles) {
			String name = file.getName();
			if(file.isDirectory() &&
					name.startsWith(experiment)) {
				if(nF1<1) {
					String phasedStates = file.getAbsolutePath()+
							"/phasedStates/"+experiment+".txt";
					if(new File(phasedStates).exists()) {
						br = getBufferedReader(phasedStates);
						int n = 0;
						String l;
						while( (l=br.readLine())!=null ) 
							if(l.startsWith("#")) n++;
						nF1 = (n/ploidy)-2;
						System.out.println(nF1+" F1 samples in the experiment.");
						br.close();
					}
				}
			}
			if(nF1>0) break;
		}

		Map<String, List<File>> map = new HashMap<String, 
				List<File>>();
		List<File> list;
		for(File file:listFiles) {
			String name = file.getName();
			if(file.isDirectory() &&
					name.startsWith(experiment)) {
				name = name.replace(experiment,"experiment");
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

		reset();
		String[] scaff_all = new String[fileObj.keySet().size()];
		fileObj.keySet().toArray(scaff_all);
		dc = new PhasedDataCollection[scaff_all.length][best_n];
		for(int i=0; i<scaff_all.length; i++) {
			Set<FileObject> files = fileObj.get(scaff_all[i]);
			executor.submit(new FileLoader(scaff_all[i],
					files.toArray(new FileObject[files.size()]),
					i));
			//new FileLoader(keys[i],
					//		files.toArray(new File[files.size()]),
			//		i).run();
		}
		executor.shutdown();
		executor.awaitTermination(365, TimeUnit.DAYS);

		//System.exit(1);
		System.out.println(map.keySet().size());
		System.out.println("["+getSystemTime()+"] LOADING FILES DONE.");

		//if( (size=map.get(key).size()) != 10) {
		//	System.out.println(key+"\t"+size);
		//}
		System.out.println("["+getSystemTime()+"] READING LOG LIKELIHOOD DONE.");

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

		/**
		reset();
		ciWriter = getGZIPBufferedWriter(outputFilePath+".ci.gz");
		ciWriter.write(os.toString());
		for(int i=0; i<dc.length; i++) 
			executor.submit(new ciCalculator(i));
			//new ciCalculator(i).run();
		executor.shutdown();
		executor.awaitTermination(365, TimeUnit.DAYS);
		ciWriter.close();

		System.exit(0);
		 **/

		reset();
		rfMedianWriter = getBufferedWriter(outputFilePath+".med");
		rfMeanWriter = getBufferedWriter(outputFilePath+".mea");
		rfMinimumWriter = getBufferedWriter(outputFilePath+".min");
		rfAllWriter = getGZIPBufferedWriter(outputFilePath+".all.gz");
		//rfAllWriter = getBufferedWriter(outputFilePath+".all");


		rfAllWriter.write(os.toString());
		long start = System.nanoTime();
		for(int i=0; i<dc.length; i++) 
			for(int j=i+1; j<dc.length; j++) 
				executor.submit(new rfCalculator(i, j));
		//new rfCalculator(i, j).run();
		executor.shutdown();
		executor.awaitTermination(365, TimeUnit.DAYS);
		long end = System.nanoTime();
		System.out.println("elapsed, "+(end-start)+
				"; \n R pool key-value pair, "+rf_pool.size() +
				"; \n R factory key-value pair, "+rf_factory.size());
		//for(long key : rf_pool.keySet()) System.out.println(key);
		Utils.print(scale_times);
		System.out.println(counter);
		
		rfMedianWriter.close();
		rfMeanWriter.close();
		rfMinimumWriter.close();
		rfAllWriter.close();

		//System.exit(1);

		reset();
		mapWriter = getBufferedWriter(outputFilePath+".map");
		mapWriter.write("##id\tkosambi\thaldane\tkSum\thSum\tkAll\thAll\n");
		for(int i=0; i<dc.length; i++) 
			executor.submit(new mapCalculator(i));
		//new mapCalculator(i).run();
		executor.shutdown();
		executor.awaitTermination(365, TimeUnit.DAYS);
		mapWriter.close();

		/**
		reset();
		gdWriter = getBufferedWriter(outputFilePath+".gds");
		for(int i=0; i<dc.length; i++) 
			executor.submit(new gdCalculator(i));
			//new gdCalculator(i).run();
		executor.shutdown();
		executor.awaitTermination(365, TimeUnit.DAYS);
		gdWriter.close();
		 **/

		System.out.println("["+getSystemTime()+"] DONE.");
	}

	private static String cat(double[] array, String sep) {
		String s = ""+array[0];
		for(int i=1; i<array.length; i++)
			s += sep+array[i];
		return s;
	}

	private static String cat(int[] array, String sep) {
		String s = ""+array[0];
		for(int i=1; i<array.length; i++)
			s += sep+array[i];
		return s;
	}

	private static class State implements Serializable {

		private static final long serialVersionUID = 7181643075632419109L;
		/**
		 * 
		 */
		private List<Character> state;
		private String str_state;

		public State(List<Character> state) {
			this.state = state;
			this.str_state = state.toString();
		}

		public int hashCode() {
			return new HashCodeBuilder(17, 31). // two randomly chosen prime numbers
					// if deriving: appendSuper(super.hashCode()).
					append(str_state).
					toHashCode();
		}

		public boolean equals(Object obj) {
			if (!(obj instanceof State))
				return false;
			if (obj == this)
				return true;
			return new EqualsBuilder().append(this.str_state,
					((State) obj).str_state).isEquals();
		}

		public String getState2String() {
			// TODO Auto-generated method stub
			return this.str_state;
		}

		public List<Character> getState() {
			return this.state;
		}

		public int size() {
			return this.state.size();
		}
	}

	private static BufferedReader getBufferedReader(String path) throws IOException {
		BufferedReader br = null;
		if(path.endsWith(".gz")){
			br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new
					FileInputStream(path))), 65536);
		}else{
			br = new BufferedReader(new FileReader(path), 65536);
		}
		return br;
	}

	private static BufferedWriter getBufferedWriter(String path) throws IOException {
		return new BufferedWriter(new FileWriter(new File(path)));
	}

	private static BufferedWriter getGZIPBufferedWriter(String path) throws IOException {
		return new BufferedWriter(new OutputStreamWriter(new
				GZIPOutputStream(new FileOutputStream(path))));
	}

	private static String getSystemTime(){
		return new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").
				format(Calendar.getInstance().getTime());
	}

	public void calcGDsAll(String phasedStates, int ploidy, 
			String[] parents, int nF1, int[] start_end,  
			double[][] rfAll, int s) {
		// TODO Auto-generated method stub
		char[][] h = readHaplotypes(phasedStates, 
				ploidy, parents, nF1);
		int c = 0;
		if(start_end[0]<=start_end[1]) {
			for(int i=start_end[0]; i<=start_end[1]; i++) { 
				for(int j=i+1; j<=start_end[1]; j++) {
					double r = 0;
					for(int k=0; k<h.length; k++) 
						r += h[k][i]==h[k][j] ? 0 : 1;
					rfAll[c++][s] = r/h.length;
				}
			}
		} else {
			for(int i=start_end[0]; i>=start_end[1]; i--) { 
				for(int j=i-1; j>=start_end[1]; j--) {
					double r = 0;
					for(int k=0; k<h.length; k++) 
						r += h[k][i]==h[k][j] ? 0 : 1;
					rfAll[c++][s] = r/h.length;
				}
			}
		}
	}

	public static double[] calcGDs(String phasedStates, int ploidy, 
			String[] parents, int nF1, int[] start_end) {
		// TODO Auto-generated method stub
		char[][] h = readHaplotypes(phasedStates, ploidy, parents, nF1);
		if(start_end[0]<=start_end[1]) {
			double[] d = new double[start_end[1]-start_end[0]+1];
			for(int i=start_end[0]; i<start_end[1]; i++) {
				double c = 0;
				for(int j=0; j<h.length; j++) 
					c += h[j][i]==h[j][i+1] ? 0 : 1;
				//d[i] = geneticDistance( c/h.length, mapFunc);
				d[i-start_end[0]] = c/h.length;
			}
			return d;
		} else {
			double[] d = new double[start_end[0]-start_end[1]+1];
			for(int i=start_end[0]; i>start_end[1]; i--) {
				double c = 0;
				for(int j=0; j<h.length; j++) 
					c += h[j][i]==h[j][i-1] ? 0 : 1;
				//d[i] = geneticDistance( c/h.length, mapFunc);
				d[start_end[0]-i] = c/h.length;
			}
			return d;
		}
	}

	public static double calcGD(String phasedStates, int ploidy,
			String[] parents, int nF1, int[] start_end, String mapFunc) {
		char[][] h = readHaplotypes(phasedStates, ploidy, parents, nF1);
		double c = 0;
		for(int i=0; i<h.length; i++)
			c += h[i][start_end[0]]==h[i][start_end[1]] ? 0 : 1;
		//return geneticDistance( c/h.length, mapFunc);
		return c/h.length;	
	}

	private static double geneticDistance(double r, String mapFunc) {
		// TODO Auto-generated method stub
		switch(mapFunc.toUpperCase()) {
		case "KOSAMBI":
			return .25*Math.log((1+2*r)/(1-2*r));
		case "HALDANE":
			return -.5*Math.log(1-2*r);	
		default:
			System.err.println("Error - Undefined genetic mapping function.");
			System.exit(1);
		}
		return -1;
	}

	private static char[][] readHaplotypes(String phasedStates, int ploidy,
			String[] parents, int nF1) {
		// TODO Auto-generated method stub
		try {
			BufferedReader br = Utils.getBufferedReader(phasedStates);
			br.readLine();
			int m = Integer.parseInt(br.readLine());
			char[][] h = new char[nF1*ploidy][m];
			String line, stateStr;
			String[] s;
			int c = 0;
			while( (line=br.readLine())!=null ) {
				if(!line.startsWith("#")) continue;
				//if(skip++<2) continue;
				s = line.split("\\s+|:");
				if(Arrays.asList(parents).contains(s[2])) continue;
				stateStr = s[s.length-1];
				h[c++] = stateStr.toCharArray();
			}
			br.close();
			return h;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.exit(1);
		}

		return null;
	}

}

