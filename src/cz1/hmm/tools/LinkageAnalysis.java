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

public class LinkageAnalysis extends RFUtils {
	private final static Logger myLogger = Logger.getLogger(LinkageAnalysis.class);
	
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

	private BufferedWriter twoPointRfWriter;
	private BufferedWriter singlePointRfWriter;
	
	private int factorial_ploidy;
	private int[][] jt_permutation;
	private int[][] permutations;
	private PhasedDataCollection[][][] dc;
	private String[] scaffolds;
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		super.initialise();
		
		factorial_ploidy = (int) Permutation.factorial(ploidy);
		jt_permutation = JohnsonTrotter.perm(ploidy);
		if(jt_permutation[0][0]!=0||jt_permutation[0][1]!=0)
			throw new RuntimeException("!!!");
		permutations = new int[factorial_ploidy][ploidy];
		for(int i=0; i<ploidy; i++) permutations[0][i] = i;
		int swap;
		for(int i=1; i<factorial_ploidy; i++) {
			System.arraycopy(permutations[i-1], 0, permutations[i], 0, ploidy);
			swap = permutations[i][jt_permutation[i][0]];
			permutations[i][jt_permutation[i][0]] = permutations[i][jt_permutation[i][1]];
			permutations[i][jt_permutation[i][1]] = swap;
		}
		
		scaffolds = new String[fileObj.keySet().size()];
		fileObj.keySet().toArray(scaffolds);
		dc = new PhasedDataCollection[scaffolds.length][][];
		
		singlePointRfWriter = Utils.getBufferedWriter(out_prefix+".map");
		try {
			this.initial_thread_pool();
			for(int i=0; i<scaffolds.length; i++) {
				executor.submit(new SinglePointAnalysis(i));
			}
			this.waitFor();
			singlePointRfWriter.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		twoPointRfWriter = Utils.getBufferedWriter(out_prefix+".txt");
		try {
			for(String scaff : scaffolds)
				twoPointRfWriter.write("##"+scaff+"\n");

			this.initial_thread_pool();
			for(int i=0; i<scaffolds.length; i++) 
				for(int j=i+1; j<scaffolds.length; j++) 
					executor.submit(new TwoPointAnalysis(i, j));
			this.waitFor();
			twoPointRfWriter.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		myLogger.info("["+Utils.getSystemTime()+"] DONE.");
	}

	public LinkageAnalysis() {}
			
	public LinkageAnalysis (String in_haps, 
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
	
	private class SinglePointAnalysis implements Runnable {
		private final int scfi;
		
		public SinglePointAnalysis(int scfi) {
			// TODO Auto-generated constructor stub
			this.scfi = scfi;
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub
			try {
				String scaff = scaffolds[scfi];
				List<FileObject> objs = fileObj.get(scaff);
				final int n = objs.size();
				dc[scfi] = new PhasedDataCollection[n][2];
				ModelReader modelReader = new ModelReader(objs.get(0).file);
				int[] distance = modelReader.getDistance();
				modelReader.close();
				int lb = objs.get(0).position[0], 
						ub = objs.get(0).position[1];
				final int m = ub-lb+1;
				int[] position = new int[m];
				for(int i=0; i<m; i++) position[i] = lb+i;
				int[] dists = new int[m-1];
				System.arraycopy(distance, lb, dists, 0, m-1);
				double[][] jumps = new double[n][];
				PhasedDataCollection d0 = null, d1 = null;
				
				int[] a0 = new int[2], a1 = new int[2];
				for(int i=0; i<n; i++) {
					FileObject obj = objs.get(i);
					modelReader = new ModelReader(obj.file);
					Map<String, char[][]> haps = modelReader.getHaplotypeByPosition(position, ploidy);
					modelReader.close();
					for(String f : parents) haps.remove(f);

					double[] jump = new double[m-1];
					d1 = new PhasedDataCollection(haps, 0);
					dc[scfi][i][0] = d1;
					for(int j=1; j<m; j++) {						
						d0 = d1;
						d1 = new PhasedDataCollection(haps, j);
						comns(d0.f0, d1.f0, a0);
						comns(d0.f1, d1.f1, a1);
						reallocate(d1.f0, a0[0]);
						reallocate(d1.f1, a1[0]);
						jump[j-1] = a0[1]+a1[1];
					}
					jumps[i] = jump;
					dc[scfi][i][1] = d1;
				}
				int deno = nF1*ploidy;
				for(int i=0; i<n; i++)
					for(int j=0; j<m-1; j++)
						jumps[i][j] = 1-jumps[i][j]/deno;
				synchronized(lock) {
					singlePointRfWriter.write("C "+scaff+"\n");
					singlePointRfWriter.write("D "+Utils.cat(dists, ",")+"\n");
					for(int i=0; i<n; i++)
						singlePointRfWriter.write(Utils.cat(jumps[i], ",")+"\n");
				}
			} catch (Exception e) {
				Thread t = Thread.currentThread();
				t.getUncaughtExceptionHandler().uncaughtException(t, e);
				e.printStackTrace();
				executor.shutdown();
				System.exit(1);
			}
		}
	}
	
	private class TwoPointAnalysis implements Runnable {
		private final int scfi, scfj;
		private final boolean conj;
		
		public TwoPointAnalysis(int scfi, int scfj) {
			// TODO Auto-generated constructor stub
			this.scfi = scfi;
			this.scfj = scfj;
			this.conj = conjPair.contains(scaffolds[scfi]+Constants.collapsed_str+scaffolds[scfj]);
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub
			try {
				int maxc = 0, sumc = 0, max, sum;
				int[] comn = new int[4], stat;
				int[] a0 = new int[4], a1 = new int[4], a2 = new int[4], 
						a3 = new int[4], a4 = new int[4], a5 = new int[4];
				
				PhasedDataCollection[] di, dj;
				for(int i = 0; i<dc[scfi].length; i++) {
					di = dc[scfi][i];
					for(int j=0; j<dc[scfj].length; j++) {
						dj = dc[scfj][j];
						if(conj) {
							comns(di[0].f0, di[1].f0, dj[0].f0, dj[1].f0, a0, true);
							comns(di[0].f1, di[1].f1, dj[0].f1, dj[1].f1, a1, true);
							sum(a0, a1, a4);
							stat = a4;
						} else {
							comns(di[0].f0, di[1].f0, dj[0].f0, dj[1].f0, a0, true);
							comns(di[0].f1, di[1].f1, dj[0].f1, dj[1].f1, a1, true);
							comns(di[0].f0, di[1].f0, dj[0].f1, dj[1].f1, a2, true);
							comns(di[0].f1, di[1].f1, dj[0].f0, dj[1].f0, a3, true);
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
						twoPointRfWriter.write(StatUtils.min(rfs)+"\t");
						for(double rf : rfs) twoPointRfWriter.write(rf+"\t");
						twoPointRfWriter.write(scaffolds[scfi]+"\t");
						twoPointRfWriter.write(scaffolds[scfj]+"\n");
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
	}
	
	private int comns(Map<Integer, Set<Integer>> a0, Map<Integer, Set<Integer>> a1,
			Map<Integer, Set<Integer>> a2, Map<Integer, Set<Integer>> a3, final int[] a, boolean permutate) {
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
		int perm = 0;
		
		if(!permutate) return perm;
		
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
				perm = i;
			} else if(max==maxc&&(sum=sum(stat))>sumc) {
				System.arraycopy(stat, 0, a, 0, 4);
				maxc = max;
				sumc = sum;
				perm = i;
			}
		}
		
		return perm;
	}

	private void comns(Map<Integer, Set<Integer>> f0, Map<Integer, Set<Integer>> f1, final int[] f) {
		// TODO Auto-generated method stub
		int[][] comns = comns(f0, f1);
		
		int maxc = 0;
		for(int i=0; i<ploidy; i++)
			maxc += comns[i][i];
		int perm = 0;
		
		int[] swap1;
		int jt0, jt1, max = maxc, s;
		for(int i=1; i<factorial_ploidy; i++) {
			jt0 = jt_permutation[i][0];
			jt1 = jt_permutation[i][1];
			
			s = max;
			s -= comns[jt0][jt0];
			s -= comns[jt1][jt1];
			swap1 = comns[jt0];
			comns[jt0] = comns[jt1];
			comns[jt1] = swap1;
			s += comns[jt0][jt0];
			s += comns[jt1][jt1];
			max = s;
			
			if(max>maxc) {
				maxc = max;
				perm = i;
			}
		}
		f[0] = perm;
		f[1] = maxc;
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
	
	private int sum(int[] arr) {
		// TODO Auto-generated method stub
		int s = 0;
		for(int a : arr) s+=a;
		return s;
	}

	private int max(int[] arr) {
		// TODO Auto-generated method stub
		int m = Integer.MIN_VALUE;
		for(int a : arr) 
			if(a>m) m = a;
		return m;
	}
	
	private int[] max(int[] arr1, int[] arr2) {
		// TODO Auto-generated method stub
		int m1 = max(arr1), m2 = max(arr2);
		if(m1>m2) return arr1;
		else if(m1<m2) return arr2;
		else {
			int s1 = sum(arr1), s2 = sum(arr2);
			if(s1>=s2) return arr1;
			else return arr2;
		}
	}
	
	private void sum(int[] arr1, int[] arr2, int[] arr) {
		// TODO Auto-generated method stub
		for(int i=0; i<arr.length; i++)
			arr[i] = arr1[i]+arr2[i];
	}


	protected void reallocate(Map<Integer, Set<Integer>> f, final int p) {
		// TODO Auto-generated method stub
		if(p==0) return;
		Map<Integer, Set<Integer>> tmp = new HashMap<>();
		for(Map.Entry<Integer, Set<Integer>> entry : f.entrySet())
			tmp.put(entry.getKey(), entry.getValue());
		for(int k=0; k<ploidy; k++) f.put(k, tmp.get(permutations[p][k])); 
	}
	
	private class PhasedDataCollection {
		private final Map<Integer, Set<Integer>> f0, f1;
		
		public PhasedDataCollection(final Map<String, char[][]> haps, final int i) {
			// TODO Auto-generated constructor stub
			f0 = new HashMap<>();
			f1 = new HashMap<>();
			for(int p=0; p<ploidy; p++) {
				f0.put(p, new HashSet<>());
				f1.put(p, new HashSet<>());
			}
			
			int H = ploidy/2;
			int h;
			
			for(int f=0; f<nF1; f++) {
				char[][] hap = haps.get(progeny[f]);
				for(int p=0; p<H; p++) {
					h = getIntInd(hap[p][i]);
					f0.get(h).add(f);
				}
				for(int p=H; p<ploidy; p++) {
					h = getIntInd(hap[p][i])-ploidy;
					f1.get(h).add(f);
				}	
			}
		}

		private int getIntInd(char hs) {
			// TODO Auto-generated method stub
			return hs>57?(hs-'a'+9):hs-'1';
		}
	}
}



























