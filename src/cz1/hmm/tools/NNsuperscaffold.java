package cz1.hmm.tools;

import java.io.File;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

public class NNsuperscaffold extends Executor {

	private String rf_file = null;
	private String out_prefix = null;
	private String RLibPath = null;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " Common:\n"
						+ "     -i/--rf                     Recombination frequency file.\n"
						+ "     -o/--prefix                 Output file prefix.\n"
						+ "     -rlib/--R-external-libs     R external library path.\n\n"
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
			myArgsEngine.add("-i", "--rf", true);
			myArgsEngine.add("-o", "--prefix", true);
			myArgsEngine.add("-rlib", "--R-external-libs", true);
			myArgsEngine.parse(args);
		}
		
		if(myArgsEngine.getBoolean("-i")) {
			rf_file = myArgsEngine.getString("-i");
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify your recombinatio frequency file.");
		}
		
		if(myArgsEngine.getBoolean("-o")) {
			out_prefix = myArgsEngine.getString("-o");
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify your output file prefix.");
		}
		
		if (myArgsEngine.getBoolean("-rlib")) {
			RLibPath = myArgsEngine.getString("-rlib");
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		final String temfile_prefix = ".tmp/";
		final boolean tmpdirCreated = Utils.makeOutputDir(new File(temfile_prefix));
		final String concorde_path = 
				RFUtils.makeExecutable("cz1/hmm/executable/concorde", temfile_prefix);
		final String nnssR_path = 
				RFUtils.makeExecutable("cz1/hmm/scripts/make_nnsuperscaffold.R", temfile_prefix);
		RFUtils.makeExecutable("cz1/hmm/scripts/include.R", temfile_prefix);
		new File(concorde_path).setExecutable(true, true);
		RFUtils.makeRMatrix(rf_file, out_prefix+".RData");
		String command = "Rscript "+nnssR_path+" "
				+ "-i "+out_prefix+".RData "
				+ "-n 2 "
				+ "-o "+out_prefix+".nnss "
				+ "--concorde "+new File(concorde_path).getParent()
				+ (RLibPath==null ? "" : " --include "+RLibPath);
		this.consume(this.bash(command));
		
		if(tmpdirCreated) Utils.deleteDirectory(new File(temfile_prefix));
	}

}
