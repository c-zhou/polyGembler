package cz1.hmm.tools;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.HashSet;

import cz1.hmm.model.ModelReader;
import cz1.math.JohnsonTrotter;
import cz1.math.Permutation;
import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Utils;

import org.apache.commons.math.stat.StatUtils;
import org.apache.log4j.Logger;

public class TwoPointAnalysis extends RFUtils {
	private final static Logger myLogger = Logger.getLogger(TwoPointAnalysis.class);
	
	protected String out_prefix;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--hap-file               Directory with input haplotype files.\n"
						+ " -o/--prefix                 Output file prefix.\n"
						+ " -ex/--experiment-id         Common prefix of haplotype files for this experiment.\n"
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
		
		if(myArgsEngine.getBoolean("-nb")) {
			best_n = Integer.parseInt(myArgsEngine.getString("-nb"));
		}
		
		if(myArgsEngine.getBoolean("-t")) {
			THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}
		
		if(myArgsEngine.getBoolean("-phi")) {
			skew_phi = Integer.parseInt(myArgsEngine.getString("-phi"));
		}
		
		if(myArgsEngine.getBoolean("-nd")) {
			drop_thres = Integer.parseInt(myArgsEngine.getString("-nd"));
		}
	}

	private BufferedWriter rfWriter;
	private int factorial_ploidy;
	private int[][] jt_permutation;
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		this.initialise();
		factorial_ploidy = (int) Permutation.factorial(ploidy);
		jt_permutation = JohnsonTrotter.perm(ploidy);
		if(jt_permutation[0][0]!=0||jt_permutation[0][1]!=0)
			throw new RuntimeException("!!!");
		rfWriter = Utils.getBufferedWriter(out_prefix+".txt");
		
		Set<String> scaffs = new HashSet<String>();
		for(int i=0; i<this.dc.length; i++) {
			if(this.dc[i][0]!=null)
				scaffs.add(this.dc[i][0].scaff);
		}
		try {
			for(String scaff : scaffs)
				rfWriter.write("##"+scaff+"\n");

			this.initial_thread_pool();
			for(int i=0; i<dc.length; i++) 
				for(int j=i+1; j<dc.length; j++) 
					executor.submit(new RfCalculator(i, j));
			this.waitFor();
			rfWriter.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		myLogger.info("["+Utils.getSystemTime()+"] DONE.");
	}

	public TwoPointAnalysis() {}
			
	public TwoPointAnalysis (String in_haps, 
				String out_prefix,
				String expr_id, 
				int threads,
				double skew_phi,
				int drop_thres,
				int best_n) { 
		this.in_haps = in_haps;
		this.out_prefix = out_prefix;
		this.expr_id = expr_id;
		THREADS = threads;
		this.skew_phi = skew_phi;
		this.drop_thres = drop_thres;
		this.best_n = best_n;
	}
	
	private class RfCalculator implements Runnable {
		private PhasedDataCollection[] data_i, data_j;
		private String scaf_i, scaf_j;
		private final boolean conj;
		
		public RfCalculator(int i, int j) {
			// TODO Auto-generated constructor stub
			data_i = dc[i];
			data_j = dc[j];
			scaf_i = data_i[0].scaff;
			scaf_j = data_j[0].scaff;
			conj = conjPair.contains(scaf_i+Constants.collapsed_str+scaf_j);
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub
			try {
				int maxc = 0, sumc = 0, max, sum;
				int[] comn = new int[4], stat;
				int[] a0 = new int[4], a1 = new int[4], a2 = new int[4], 
						a3 = new int[4], a4 = new int[4], a5 = new int[4];

				PhasedDataCollection di, dj;
				for(int i = 0; i<best_n; i++) {
					di = data_i[i];
					if(di==null) continue;
					for(int j=0; j<best_n; j++) {
						dj = data_j[j];
						if(dj==null) continue;
						if(conj) {
							comns(di.p0, di.p1, dj.p0, dj.p1, a0, false);
							comns(di.m0, di.m1, dj.m0, dj.m1, a1, false);
							sum(a0, a1, a4);
							stat = a4;
						} else {
							comns(di.p0, di.p1, dj.p0, dj.p1, a0, true);
							comns(di.m0, di.m1, dj.m0, dj.m1, a1, true);
							comns(di.p0, di.p1, dj.m0, dj.m1, a2, true);
							comns(di.m0, di.m1, dj.p0, dj.p1, a3, true);
							sum(a0, a1, a4);
							sum(a2, a3, a5);
							stat = max(a4, a5);
						}
						
						if((max=max(stat))>maxc) {
							System.arraycopy(stat, 0, comn, 0, 4);
							maxc = max;
							sumc = sum(stat);
						} else if(max==maxc&&(sum=sum(stat))>sumc) {
							System.arraycopy(stat, 0, comn, 0, 4);
							maxc = max;
							sumc = sum;
						}
					}
				}

				double[] rfs = new double[4];
				double deno = nF1*ploidy;
				for(int i=0; i<4; i++) rfs[i] = 1-comn[i]/deno;

				synchronized(lock) {
					try {
						rfWriter.write(StatUtils.min(rfs)+"\t");
						for(double rf : rfs) rfWriter.write(rf+"\t");
						rfWriter.write(scaf_i+"\t");
						rfWriter.write(scaf_j+"\n");
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			} catch (Exception e) {
				Thread t = Thread.currentThread();
				t.getUncaughtExceptionHandler().uncaughtException(t, e);
				e.printStackTrace();
				executor.shutdown();
				System.exit(1);
			}
		}
		
		private void comns(Map<Integer, Set<Integer>> a0, Map<Integer, Set<Integer>> a1,
				Map<Integer, Set<Integer>> a2, Map<Integer, Set<Integer>> a3, 
				final int[] a, boolean permutate) {
			// TODO Auto-generated method stub
			int[][][] comns = new int[4][][];
			comns[0] = comns(a0, a2);
			comns[1] = comns(a0, a3);
			comns[2] = comns(a1, a2);
			comns[3] = comns(a1, a3);
			
			int maxc = 0, sumc = 0, max, sum;
				int[] stat = new int[4];
			int[][] comn;
			for(int i=0; i<4; i++) {
				comn = comns[i];
				for(int j=0; j<ploidy; j++)
					stat[i] += comn[j][j];
			}
			
			System.arraycopy(stat, 0, a, 0, 4);
			maxc = max(stat);
			sumc = sum(stat);
			
			if(!permutate) return;
			
			int[] swap1;
			int jt0, jt1, s;
			for(int i=1; i<factorial_ploidy; i++) {
				jt0 = jt_permutation[i][0];
				jt1 = jt_permutation[i][1];
				for(int j=0; j<4; j++) {
					comn = comns[j];
					s = stat[j];
					s -= comn[jt0][jt0];
					s -= comn[jt1][jt1];
					swap1 = comn[jt0];
					comn[jt0] = comn[jt1];
					comn[jt1] = swap1;
					s += comn[jt0][jt0];
					s += comn[jt1][jt1];
					stat[j] = s;
				}
				
				if((max=max(stat))>maxc) {
					System.arraycopy(stat, 0, a, 0, 4);
					maxc = max;
					sumc = sum(stat);
				} else if(max==maxc&&(sum=sum(stat))>sumc) {
					System.arraycopy(stat, 0, a, 0, 4);
					maxc = max;
					sumc = sum;
				}
			}
			
			return;
		}

		private int[][] comns(Map<Integer, Set<Integer>> a1, Map<Integer, Set<Integer>> a2) {
			// TODO Auto-generated method stub
			int[][] comns = new int[ploidy][ploidy];
			Set<Integer> s1, s2;
			for(int i=0; i<ploidy; i++) {
				s1 = a1.get(i);
				for(int j=0; j<ploidy; j++) {
					s2 = a2.get(j);
					int n = 0;
					for(int k : s1) 
						if(s2.contains(k)) 
							++n;
					comns[j][i] = n;
				}
			}
			return comns;
		}
	}
	
	@Override
	protected void initialise() {
		// TODO Auto-generated catch block
		super.initialise();
		this.collectData();
	}

	private class PhasedDataCollection {
		private final String scaff;
		private final Map<Integer, Set<Integer>> p0, p1, m0, m1;
		
		public PhasedDataCollection(String scaff, FileObject file) {
			// TODO Auto-generated constructor stub
			this.scaff = scaff;
			p0 = new HashMap<>();
			p1 = new HashMap<>();
			m0 = new HashMap<>();
			m1 = new HashMap<>();
			for(int i=0; i<ploidy; i++) {
				p0.put(i, new HashSet<>());
				p1.put(i, new HashSet<>());
				m0.put(i, new HashSet<>());
				m1.put(i, new HashSet<>());
			}
			this.collectData(file);
		}

		private void collectData(FileObject file) {
			// TODO Auto-generated method stub
			ModelReader modelReader = new ModelReader(file.file);
			Map<String, char[][]> haps = modelReader.getHaplotypeByPosition(file.position, ploidy);
			modelReader.close();
			for(String f : parents) haps.remove(f);
			if(haps.size()!=nF1) throw new RuntimeException("!!!");
			int H = ploidy/2;
			int h;
			
			for(int f=0; f<nF1; f++) {
				char[][] hap = haps.get(progeny[f]);
				for(int i=0; i<H; i++) {
					h = getIntInd(hap[i][0]);
					p0.get(h).add(f);
					h = getIntInd(hap[i][1]);
					p1.get(h).add(f);
				}
				for(int i=H; i<ploidy; i++) {
					h = getIntInd(hap[i][0])-ploidy;
					m0.get(h).add(f);
					h = getIntInd(hap[i][1])-ploidy;
					m1.get(h).add(f);
				}	
			}
		}

		private int getIntInd(char hs) {
			// TODO Auto-generated method stub
			return hs>57?(hs-'a'+9):hs-'1';
		}
	}
	
	PhasedDataCollection[][] dc;
	
	private void collectData() {
		// TODO Auto-generated method stub
		this.initial_thread_pool();
		final String[] scaffs = new String[fileObj.keySet().size()];
		fileObj.keySet().toArray(scaffs);
		dc = new PhasedDataCollection[scaffs.length][best_n];
		for(int i=0; i<scaffs.length; i++) {
			List<FileObject> objs = fileObj.get(scaffs[i]);
			for(int j=0; j<objs.size(); j++) {
				executor.submit(new Runnable() {
					private int i;
					private int j;
					
					@Override
					public void run() {
						// TODO Auto-generated method stub
						try {
							String scaff = scaffs[i];
							FileObject obj = objs.get(j);
							dc[i][j] = new PhasedDataCollection(scaff, obj);
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
						return this;
					}
					
				}.init(i, j));	
			}
		}
		this.waitFor();

		myLogger.info("["+Utils.getSystemTime()+"] LOADING FILES DONE.");
		myLogger.info("["+Utils.getSystemTime()+"] READING LOG LIKELIHOOD DONE.");
	}
}



























