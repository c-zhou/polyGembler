package cz1.simulation.tools;

import java.util.concurrent.TimeUnit;

import cz1.simulation.model.Population;
import cz1.util.ArgsEngine;
import cz1.util.Executor;

public class PopulationSimulator extends Executor {

	private Population pop;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -r/--reference		Reference (fasta file format). \n"
						+ " -n/--pop-size		Population size including parents (default 96). \n"
						+ " -p/--ploidy  		Copy number of chromosomes (default 2). \n"
						+ " -t/--threads 		Number of threads (default 1).\n"
						+ " -s/--run-id			Unique run id (default Sc1).\n" 
						+ " -o/--prefix 	 	Output directory (defult current directory).\n\n");
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
			myArgsEngine.add("-r", "--reference", true);
			myArgsEngine.add("-n", "--pop-size", true);
			myArgsEngine.add("-p", "--ploidy", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-s", "--run-id", true);
			myArgsEngine.add("-o", "--prefix", true);
			myArgsEngine.parse(args);
		}

		String reference, id="Sc1", output="./";
		int f1=94, ploidy=2;
		if (myArgsEngine.getBoolean("-r")) {
			reference = myArgsEngine.getString("-r");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the reference FASTA file.");
		}


		if (myArgsEngine.getBoolean("-n")) {
			f1 = Integer.parseInt(myArgsEngine.getString("-n"))-2;
		}

		if (myArgsEngine.getBoolean("-p")) {
			ploidy = Integer.parseInt(myArgsEngine.getString("-p"));
		}

		if (myArgsEngine.getBoolean("-t")) {
			THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}

		if (myArgsEngine.getBoolean("-s")) {
			id = myArgsEngine.getString("-s");
		}

		if (myArgsEngine.getBoolean("-o")) {
			output = myArgsEngine.getString("-o");
		}
		
		pop = new Population(reference, f1, ploidy, id, output);
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		this.initial_thread_pool();

		final int M = pop.indivCount();

		for (int i=0; i<M; i++) {

			executor.submit(new Runnable() {
				private int i;

				@Override
				public void run() {
					try {
						pop.writeGenomeFile(i);
					} catch (Exception e) {
						Thread t = Thread.currentThread();
						t.getUncaughtExceptionHandler().uncaughtException(t, e);
						e.printStackTrace();
						executor.shutdown();
						System.exit(1);
					}
				} 

				public Runnable init(int i) {
					this.i = i;
					return(this);
				}
			}.init(i) );
		}
		
		this.waitFor();
	}

}
