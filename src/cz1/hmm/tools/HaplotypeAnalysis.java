package cz1.hmm.tools;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import cz1.hmm.model.ModelReader;
import cz1.math.JohnsonTrotter;
import cz1.math.Permutation;
import cz1.util.ArgsEngine;
import cz1.util.Utils;

public class HaplotypeAnalysis extends RFUtils {
	private final static Logger myLogger = Logger.getLogger(HaplotypeAnalysis.class);
	
	private String phase_file;
	private final Map<String, Map<String, int[][]>> phased_haps = new HashMap<>();
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--hap-file               Directory with input haplotype files.\n"
						+ " -p/--phase-file             Phased haplotype file.\n"
						+ " -ex/--experiment-id         Common prefix of haplotype files for this experiment.\n"
						+ " -phi/--skew-phi             For a haplotype inference, the frequencies of parental \n"
						+ "                             haplotypes need to be in the interval [1/phi, phi], \n"
						+ "                             otherwise will be discared (default 2).\n"
						+ " -t/--threads                #threads (default 1).\n"	
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
			myArgsEngine.add( "-i", "--hap-file", true);
			myArgsEngine.add( "-p", "--phased-file", true);
			myArgsEngine.add( "-ex", "--experiment-id", true);
			myArgsEngine.add( "-phi", "--skew-phi", true);
			myArgsEngine.add( "-t", "--threads", true);
			myArgsEngine.parse(args);
		}
		
		if(myArgsEngine.getBoolean("-i")) {
			in_haps = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your input zip file.");
		}

		if(myArgsEngine.getBoolean("-p")) {
			phase_file = myArgsEngine.getString("-p");
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify your phased haplotype file.");
		}
		
		if(myArgsEngine.getBoolean("-ex")) {
			expr_id = myArgsEngine.getString("-ex");
		}  else {
			expr_id = guessExperimentId();
			myLogger.warn("No experiment prefix provided, I guess it's "+expr_id+". Please\n"
					+ "specify it with -ex/--experiment-id option if it's incorrect.");
		}
		
		if(myArgsEngine.getBoolean("-phi")) {
			skew_phi = Integer.parseInt(myArgsEngine.getString("-phi"));
		}
		
		this.best_n = 1;
		this.skew_phi = Double.POSITIVE_INFINITY;
		
		if(myArgsEngine.getBoolean("-t")) {
			THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}
	}
	
	public HaplotypeAnalysis (String in_haps, 
			String phase_file,
			String expr_id, 
			int threads,
			double skew_phi) { 
		this.in_haps = in_haps;
		this.phase_file = phase_file;
		this.expr_id = expr_id;
		THREADS = threads;
		this.skew_phi = skew_phi;
	}
	
	public HaplotypeAnalysis() {
		// TODO Auto-generated constructor stub
		super();
	}

	private int factorial_ploidy;
	private int[][] jt_permutation;
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		super.initialise();
		this.collectData();
		this.parsePhasedHaplotypes();
		factorial_ploidy = (int) Permutation.factorial(ploidy);
		jt_permutation = JohnsonTrotter.perm(ploidy);
		if(jt_permutation[0][0]!=0||jt_permutation[0][1]!=0)
			throw new RuntimeException("!!!");
		this.initial_thread_pool();
		for(String scaff : dc.keySet())
			executor.submit(new Evaluator(scaff));
		this.waitFor();
		
		myLogger.info("     Total sites: "+totals);
		myLogger.info("   Matched sites: "+matches);
		myLogger.info("Phasing accuracy: "+(double)matches/totals);
	}
	
	final Map<String, PhasedDataCollection> dc = new HashMap<>();
	
	private class PhasedDataCollection {
		private final String[] markers;
		private final Map<String, int[][]> data;
		
		public PhasedDataCollection(FileObject file) {
			// TODO Auto-generated constructor stub
			ModelReader modelReader = new ModelReader(file.file);
			this.markers = modelReader.getSnpIdByPositionRange(file.position);
			this.data = modelReader.getGenotypeByPositionRange(file.position, ploidy);
			modelReader.close();
		}
	}
	
	private void collectData() {
		// TODO Auto-generated method stub
		this.initial_thread_pool();
		final String[] scaffs = fileObj.keySet().toArray(new String[fileObj.size()]);
		for(int i=0; i<scaffs.length; i++) {
			executor.submit(new Runnable() {
				String scaff;
				
				@Override
				public void run() {
					// TODO Auto-generated method stub
					try {
						dc.put(scaff, new PhasedDataCollection(fileObj.get(scaff).get(0)));
					} catch (Exception e) {
						Thread t = Thread.currentThread();
						t.getUncaughtExceptionHandler().uncaughtException(t, e);
						e.printStackTrace();
						executor.shutdown();
						System.exit(1);
					}
				}

				public Runnable init(String scaff) {
					// TODO Auto-generated method stub
					this.scaff = scaff;
					return this;
				}

			}.init(scaffs[i]));	
		}
		this.waitFor();
	}

	private void parsePhasedHaplotypes() {
		// TODO Auto-generated method stub
		try {
			BufferedReader br = Utils.getBufferedReader(phase_file);
			String line;
			String[] s, s1;
			String[] header = null;
			String scaff = null;
			List<String[]> lineBuff = new ArrayList<>();
			
			line = br.readLine();
			while(line!=null) {
				if(line.startsWith("##")) {
					line = br.readLine();
					continue;
				}
				if(line.startsWith("#")) {
					header = line.split("\\s+");
					line = br.readLine();
					continue;
				}
				
				s = line.split("\\s+");
				scaff = s[0];
				lineBuff.add(s);
				
				while((line=br.readLine())!=null) {
					s = line.split("\\s+");
					if(scaff.equals(s[0]))
						lineBuff.add(s);
					else
						break;
				}
				
				if(dc.containsKey(scaff)) {
					String[] snpIds = dc.get(scaff).markers;
					List<Integer> selected = new ArrayList<>();
					int d = 0;
					int n = lineBuff.size();
					String snpId;
					for(int i=0; i<snpIds.length; i++) {
						while(d<n) {
							snpId = lineBuff.get(d)[0]+"_"+lineBuff.get(d)[1];
							if(snpIds[i].equals(snpId)) {
								selected.add(d);
								++d;
								break;
							} else {
								++d;
							}
						}
					}
					int m = selected.size();
					if(m!=snpIds.length) throw new RuntimeException("!!!");

					n = lineBuff.get(0).length;
					Map<String, int[][]> phases = new HashMap<>();
					for(int i=9; i<n; i++) {
						int[][] phase = new int[ploidy][m];
						for(int j=0; j<m; j++) {
							s1 = lineBuff.get(selected.get(j))[i].split("\\|");
							for(int k=0; k<ploidy; k++)
								phase[k][j] = s1[k].equals(".") ? -1 : Integer.parseInt(s1[k]);
						}
						phases.put(header[i], phase);
					}
					phased_haps.put(scaff, phases);
				}
				
				lineBuff.clear();
			}
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private int totals = 0;
	private int matches = 0;
	
	private class Evaluator implements Runnable {
		private final String scaff;
		
		public Evaluator(String scaff) {
			// TODO Auto-generated constructor stub
			this.scaff = scaff;
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub
			try {
				Map<String, int[][]> phased = phased_haps.get(scaff);
				Map<String, int[][]> genuine = dc.get(scaff).data;
				for(String sample : phased.keySet())
					evaluate(phased.get(sample), genuine.get(sample));
			} catch (Exception e) {
				Thread t = Thread.currentThread();
				t.getUncaughtExceptionHandler().uncaughtException(t, e);
				e.printStackTrace();
				executor.shutdown();
				System.exit(1);
			}
		}

		private void evaluate(int[][] hap1, int[][] hap2) {
			// TODO Auto-generated method stub
			
			int[][] total_table = new int[ploidy][ploidy];
			int[][] match_table = new int[ploidy][ploidy];
			for(int i=0; i<ploidy; i++) {
				int[] h1 = hap1[i];
				for(int j=0; j<ploidy; j++) {
					int h2[] = hap2[j];
					int s = 0;
					int m = 0;
					for(int k=0; k<h1.length; k++) {
						if(h1[k]!=-1&&h2[k]!=-1) {
							++s;
							if(h1[k]==h2[k]) ++m;
						}
					}
					total_table[i][j] = s;
					match_table[i][j] = m;
				}
			}
			
			int m_total = 0, m_match = 0, total, match; 
			double m_acc = 0, acc;
			int jt0, jt1;
			int[] swap1;
			
			for(int i=0; i<factorial_ploidy; i++) {
				jt0 = jt_permutation[i][0];
				jt1 = jt_permutation[i][1];
				
				swap1 = total_table[jt0];
				total_table[jt0] = total_table[jt1];
				total_table[jt1] = swap1;
				swap1 = match_table[jt0];
				match_table[jt0] = match_table[jt1];
				match_table[jt1] = swap1;
				
				total = 0;
				match = 0;
				for(int j=0; j<ploidy; j++) {
					total += total_table[j][j];
					match += match_table[j][j];
				}
				acc = (double)match/total;
				
				if(acc>m_acc) {
					m_acc = acc;
					m_total = total;
					m_match = match;
				}
			}
			
			synchronized(lock) {
				totals += m_total;
				matches += m_match;
			}
		}
	}
}
