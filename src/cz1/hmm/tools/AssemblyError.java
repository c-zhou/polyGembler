package cz1.hmm.tools;

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
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.apache.commons.math3.stat.inference.GTest;

import cz1.util.Algebra;

public class AssemblyError {
	public static void main(String[] args) 
			throws IOException, InterruptedException {

		// create the command line parser
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
		String experiment=null, workspace=null, output=null,
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
		}
		catch( ParseException exp ) {
			System.out.println( "Unexpected exception:" + exp.getMessage() );
		}

		AssemblyError sgRF = new AssemblyError(
				workspace,experiment,ploidy,parents,t,s,d,goodnessOfFitTest,n);
		sgRF.calcMisScaffold(output);
	}

	private final static int NUM_CORES = 
			Runtime.getRuntime().availableProcessors();
	private static String workspace;
	private static String experiment;
	private static String[] parents;
	private static int ploidy;
	private static int THREADS = NUM_CORES-1;
	private static int drop_thres;
	private static double skew_thres;
	private static double[] probs_uniform;
	private static String goodness_of_fit;
	private static int best_n;
	
	public AssemblyError (String workspace_, 
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
		parents = parents_;
		THREADS = threads;
		drop_thres = d;
		skew_thres = s;
		goodness_of_fit = t.toLowerCase();
		probs_uniform = new double[ploidy*2];
		Arrays.fill(probs_uniform, .5/ploidy);
		best_n = n;
	}
	
	private class FileLoader implements Runnable {
		private final String id;
		private final File[] files;
		private final int i;

		public FileLoader(String id, File[] files, int i) {
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

			for(int k=0; k<this.files.length; k++) {

				File file = this.files[k];
				if(!new File(file.getAbsolutePath()+
						"/phasedStates/"+experiment+".txt").exists()) {
					System.err.println("warning: "+
							file.getName()+
							" exsits, but phased states do not.");
					continue;
				}
				try {
					BufferedReader br = 
							getBufferedReader(file.getAbsolutePath()+
									"/phasedStates/"+experiment+".txt");
					br.readLine();
					String mak = br.readLine();
					br.close();
					if( mak==null ) {
						System.err.println("warning: "+
								file.getName()+
								" exists, but phased states are NULL.");
						continue;
					}
					if( Integer.parseInt(mak)<2 ) {
						System.err.println("warning: "+
								file.getName()+
								" exists, but #marker is less than 2.");
						continue;
					}

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
							ll[k] = Double.parseDouble(s0[3]);
						}
					} else {
						BufferedReader br2 = 
								getBufferedReader(file.getAbsolutePath()+
										"/phasedStates/"+
										experiment+".txt");
						String lprob=br2.readLine();
						br2.close();
						if( lprob!=null ) 
							ll[k] = Double.parseDouble(lprob);
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
						System.out.println(this.files[maxN[k]].getName()+
								" was dropped due to large haploptype frequency variance. (" +
								cat(phases, ",") +")");
						drop[k] = true;
					}
					if(drop[k]) dropped++;
					oos.append("["+(drop[k]?"drop](maf,":"keep](maf,")+maf+";mif,"+mif+") "+
							cat(haps_observed,",")+"\t"+this.files[maxN[k]].getName()+"\n");
					break;
				case "chisq":
					observed = new long[ploidy*2];
					for(int z=0; z<observed.length; z++) 
						observed[z] = (long) haps_observed[z];
					p = new ChiSquareTest().chiSquareTest(probs_uniform, observed);
					if(p<skew_thres) drop[k] = true;
					if(drop[k]) dropped++;
					oos.append("["+(drop[k]?"drop](p,":"keep](p,")+formatter.format(p)+") "+
							cat(haps_observed,",")+"\t"+this.files[maxN[k]].getName()+"\n");
					break;
				case "gtest":
					observed = new long[ploidy*2];
					for(int z=0; z<observed.length; z++) 
						observed[z] = (long) haps_observed[z];
					p = new GTest().gTest(probs_uniform, observed);
					if(p<skew_thres) drop[k] = true;
					if(drop[k]) dropped++;
					oos.append("["+(drop[k]?"drop](p,":"keep](p,")+formatter.format(p)+") "+
							cat(haps_observed,",")+"\t"+this.files[maxN[k]].getName()+"\n");
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
						dc[i][kk] = this.files[maxN[k]]
								.getName();
						
						kk++;
					}
					if(kk>=best_n) break;
				}
			}
		}
		
		private int[] readHaplotypes(final int i) {
			// TODO Auto-generated method stub
			try {
				BufferedReader br_states = getBufferedReader(this.files[i]+
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
	}

	private class mapCalculator implements Runnable {
		private final int i;
		public mapCalculator(int i) {
			this.i = i;
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub
			double[][] rf_all = new double[dc[this.i].length][];

			for(int k=0; k<dc[this.i].length; k++) 
				if( dc[i][k]!=null ) 
					rf_all[k] = calcGDs(workspace+"/"+dc[i][k]+"/phasedStates/"+experiment+".txt",
							ploidy,
							parents,
							nF1);

			String contig = dc[i][0].replace(experiment,"experiment")
					.split("\\.")[1];
			rf_all = Algebra.transpose(rf_all);
			
			try {
				StringBuilder os = new StringBuilder();
				os.append("*"+contig+"\n");
				for(int k=0; k<rf_all.length; k++)
					if(rf_all[k]!=null)
						os.append(cat(rf_all[k], ",")+"\n");
				mapWriter.write(os.toString());
			} catch (MathIllegalArgumentException | IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	private static int nF1;
	private static BufferedWriter mapWriter;
	private static NumberFormat formatter = new DecimalFormat("#0.000");
	private static String[][] dc = null;
	private static ExecutorService executor = 
			Executors.newFixedThreadPool(THREADS);
	private void reset() {
		executor = Executors.newFixedThreadPool(THREADS);
	}

	private void calcMisScaffold(String outputFilePath) throws IOException, InterruptedException {
		
		String[] s;
		BufferedReader br;
		File folder = new File(workspace);
		//File folder = new File("C:\\Users\\chenxi.zhou\\Desktop\\console out");
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

		dc = new String[map.keySet().size()][best_n];
		int ii = 0;
		for(String contig : map.keySet()) {
			List<File> files = map.get(contig);
			executor.submit(new FileLoader(contig,
					files.toArray(new File[files.size()]),
					ii++));
		}
		executor.shutdown();
		executor.awaitTermination(365, TimeUnit.DAYS);

		System.out.println(map.keySet().size());
		System.out.println("["+getSystemTime()+"] LOADING FILES DONE.");
		System.out.println("["+getSystemTime()+"] READING LOG LIKELIHOOD DONE.");
		
		reset();
		mapWriter = getBufferedWriter(outputFilePath+".map");
		for(int i=0; i<dc.length; i++) 
			executor.submit(new mapCalculator(i));
		executor.shutdown();
		executor.awaitTermination(365, TimeUnit.DAYS);
		mapWriter.close();

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
			String[] parents, int nF1, int nSNP, 
			double[][] rfAll, int s) {
		// TODO Auto-generated method stub
		char[][] h = readHaplotypes(phasedStates, 
				ploidy, parents, nF1);
		int c = 0;
		for(int i=0; i<nSNP; i++) { 
			for(int j=i+1; j<nSNP; j++) {
				double r = 0;
				for(int k=0; k<h.length; k++) 
					r += h[k][i]==h[k][j] ? 0 : 1;
				rfAll[c++][s] = r/h.length;
			}
		}
	}

	public static double[] calcGDs(String phasedStates, int ploidy, 
			String[] parents, int nF1) {
		// TODO Auto-generated method stub
		char[][] h = readHaplotypes(phasedStates, ploidy, parents, nF1);
		double[] d = new double[h[0].length-1];
		for(int i=0; i<d.length; i++) {
			double c = 0;
			for(int j=0; j<h.length; j++) 
				c += h[j][i]==h[j][i+1] ? 0 : 1;
			d[i] = c/h.length;
		}
		return d;
	}

	public static double calcGD(String phasedStates, int ploidy,
			String[] parents, int nF1) {
		char[][] h = readHaplotypes(phasedStates, ploidy, parents, nF1);
		double c = 0;
		for(int i=0; i<h.length; i++)
			c += h[i][0]==h[i][h[i].length-1] ? 0 : 1;
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
			BufferedReader br = getBufferedReader(phasedStates);
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
