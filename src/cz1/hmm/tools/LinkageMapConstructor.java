package cz1.hmm.tools;

import java.io.File;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

public class LinkageMapConstructor extends Executor {

	private String rf_file = null;
	private String map_file = null;
	private String out_prefix = null;
	private String RLibPath = null;
	private boolean one = false;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " Common:\n"
						+ "     -i/--rf                     Recombination frequency file.\n"
						+ "     -m/--map                    Recombination map file.\n"
						+ "     -1/--one-group              One group. \n"
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
			myArgsEngine.add("-m", "--map", true);
			myArgsEngine.add("-1", "--one-group", false);
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
		
		if(myArgsEngine.getBoolean("-m")) {
			map_file = myArgsEngine.getString("-m");
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify your recombination map file.");
		}
		
		if(myArgsEngine.getBoolean("-1")) {
			one = true;
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
		Utils.makeOutputDir(temfile_prefix);
		final String concorde_path = 
				RFUtils.makeExecutable("cz1/hmm/executable/concorde", temfile_prefix);
		new File(concorde_path).setExecutable(true, true);
		final String mklgR_path = 
				RFUtils.makeExecutable("cz1/hmm/scripts/make_geneticmap.R", temfile_prefix);
		RFUtils.makeExecutable("cz1/hmm/scripts/include.R", temfile_prefix);
		RFUtils.makeRMatrix(rf_file, out_prefix+".RData");
		final String command = one ?
				"Rscript "+mklgR_path+" "
				+ "-i "+out_prefix+".RData "
				+ "-m "+map_file+" "
				+ "-1 "
				+ "-o "+out_prefix+" "
				+ "--concorde "+new File(concorde_path).getParent()
				+ (RLibPath==null ? "" : " --include "+RLibPath) 
				:
					"Rscript "+mklgR_path+" "
					+ "-i "+out_prefix+".RData "
					+ "-m "+map_file+" "
					+ "-o "+out_prefix+" "
					+ "--concorde "+new File(concorde_path).getParent()
					+ (RLibPath==null ? "" : " --include "+RLibPath);
		this.consume(this.bash(command));
		
		new File(temfile_prefix).delete();
	}

}
