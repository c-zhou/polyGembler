package cz1.hmm.tools;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import cz1.hmm.model.ModelReader;
import cz1.util.ArgsEngine;
import cz1.util.Utils;

import org.apache.log4j.Logger;

public class SinglePointAnalysis extends RFUtils {
	private final static Logger myLogger = Logger.getLogger(SinglePointAnalysis.class);
	
	private String out_prefix;
	
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
		
		best_n = 30; // use as much as possible up to 30
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
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		this.initialise();
		rfWriter = Utils.getBufferedWriter(out_prefix+".map");

		final String[] scaffs = new String[fileObj.keySet().size()];
		fileObj.keySet().toArray(scaffs);
		
		this.initial_thread_pool();
		for(String scaff : scaffs) {
			executor.submit(new RfCalculator(scaff));
		}
		this.waitFor();
		
		try {
			rfWriter.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public SinglePointAnalysis() {}
			
	public SinglePointAnalysis (String in_haps, 
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
		private final String scaff;
		
		public RfCalculator(String scaff) {
			// TODO Auto-generated constructor stub
			this.scaff = scaff;
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub
			try {
				List<FileObject> objs = fileObj.get(scaff);
				final int n = objs.size();
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
				for(int i=0; i<n; i++) {
					FileObject obj = objs.get(i);
					modelReader = new ModelReader(obj.file);
					Map<String, char[][]> haps = modelReader.getHaplotypeByPosition(position, ploidy);
					modelReader.close();
					for(String f : parents) haps.remove(f);
					double[] jump = new double[m-1];
					char[] h;
					for(char[][] hap : haps.values()) {
						for(int j=0; j<ploidy; j++) {
							h = hap[j];
							for(int k=0; k<m-1; k++)
								if(h[k]!=h[k+1]) jump[k] += 1;
						}
					}
					jumps[i] = jump;
				}
				int deno = nF1*ploidy;
				for(int i=0; i<n; i++)
					for(int j=0; j<m-1; j++)
						jumps[i][j] /= deno;
				synchronized(lock) {
					rfWriter.write("C "+scaff+"\n");
					rfWriter.write("D "+Utils.cat(dists, ",")+"\n");
					for(int i=0; i<n; i++)
						rfWriter.write(Utils.cat(jumps[i], ",")+"\n");
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
}



























