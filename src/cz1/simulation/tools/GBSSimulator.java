package cz1.simulation.tools;

import java.util.concurrent.TimeUnit;

import cz1.simulation.model.GBS;
import cz1.util.ArgsEngine;
import cz1.util.Executor;

public class GBSSimulator extends Executor {

	private GBS gbs;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -f/--fasta-file     Directory contains genome fasta files to be sequenced. \n"
						+ " -e/--enzyme         Enzyme(default PstI). \n"
						+ " -l/--library        GBS protocol library preparation file (default null). \n"
						+ " -t/--threads        Number of threads (default 1).\n"
						+ " -b/--barcode-file   GBS protocol barcode file (default null).\n"
						+ " -m/--avg-depth      Depth of coverage (default 5).\n"
						+ " -s/--sdev           Standard deviation of depth of coverage (default 5).\n"
						+ " -S/--random-seed    Random seed (default system nano time). \n"
						+ " -q/--quality-file   Markov chain parameter file for quality scores (default null). \n"
						+ " -o/--output-prefix  Output directory (defult current directory).\n\n");
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
			myArgsEngine.add("-f", "--fasta-file", true);
			myArgsEngine.add("-e", "--enzyme", true);
			myArgsEngine.add("-l", "--library", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-b", "--barcode-file", true);
			myArgsEngine.add("-m", "--avg-depth", true);
			myArgsEngine.add("-s", "--sdev", true);
			myArgsEngine.add("-S", "--random-seed", true);
			myArgsEngine.add("-q", "--quality-file", true);
			myArgsEngine.add("-o", "--output-file", true);
			myArgsEngine.parse(args);
		}

		String fastaFileDir, 
			enzymeName = "PstI", 
			libPrepFilePath = null, 
			barcodeFilePath = null,
			params = null,
			outputDir = "./";
		double avg = 5, sd = 5;
		long RANDOM_SEED = System.nanoTime();
		
		if (myArgsEngine.getBoolean("-f")) {
			fastaFileDir = myArgsEngine.getString("-f");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the FASTA files.");
		}

		if (myArgsEngine.getBoolean("-e")) {
			enzymeName = myArgsEngine.getString("-e");
		}

		if (myArgsEngine.getBoolean("-l")) {
			libPrepFilePath = myArgsEngine.getString("-l");
		}

		if (myArgsEngine.getBoolean("-b")) {
			barcodeFilePath = myArgsEngine.getString("-b");
		}

		if (myArgsEngine.getBoolean("-t")) {
			THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}
		
		if (myArgsEngine.getBoolean("-m")) {
			avg = Double.parseDouble(myArgsEngine.getString("-m"));
		}
		
		if (myArgsEngine.getBoolean("-s")) {
			sd = Double.parseDouble(myArgsEngine.getString("-s"));
		}

		if (myArgsEngine.getBoolean("-S")) {
			RANDOM_SEED = Long.parseLong(myArgsEngine.getString("-S"));
		}
		
		if (myArgsEngine.getBoolean("-q")) {
			params = myArgsEngine.getString("-q");
		}
		
		if (myArgsEngine.getBoolean("-o")) {
			outputDir = myArgsEngine.getString("-o");
		}
		
		gbs = new GBS(fastaFileDir,
				enzymeName,
				avg,
				sd,
				params,
				libPrepFilePath,
				barcodeFilePath,
				outputDir,
				RANDOM_SEED);
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		this.initial_thread_pool();
		for(int f=0; f<gbs.size(); f++) {
			executor.submit(new Runnable() {
				private int f;

				@Override
				public void run() {
					try {
						gbs.simulate(f);
					} catch (Exception e) {
						Thread t = Thread.currentThread();
						t.getUncaughtExceptionHandler().uncaughtException(t, e);
						e.printStackTrace();
						executor.shutdown();
						System.exit(1);
					}
				}
				public Runnable init(int f) {
					this.f = f;
					return(this);
				}
			}.init(f) );
		}
		this.waitFor();
	}	
}
