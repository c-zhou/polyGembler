package cz1.hmm.tools;

import java.io.File;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

public class TwoPointConstructor extends Executor {

	private String rf_file = null;
	private String out_prefix = null;
	private double fixed_rf = -1;
	private String twopoint_file = null;
	private String RLibPath = null;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " Common:\n"
						+ "     -i/--rf                     Recombination frequency file.\n"
						+ "     -o/--prefix                 Output file prefix.\n"
						+ "     -r/--fixed-rf               Fix recombination frequency. \n"
						+ "     -2/--two-point              File contains two-point groups. \n"
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
			myArgsEngine.add("-r", "--fixed-rf", true);
			myArgsEngine.add("-2", "--two-point", true);
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
		
		if(myArgsEngine.getBoolean("-2")) {
			twopoint_file = myArgsEngine.getString("--two-point");
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify your two-point file prefix.");
		}
		
		if (myArgsEngine.getBoolean("-r")) {
			fixed_rf = Double.parseDouble(myArgsEngine.getString("-r"));
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
		final String tpR_path = 
				RFUtils.makeExecutable("cz1/hmm/scripts/make_twopoint.R", temfile_prefix);
		RFUtils.makeExecutable("cz1/hmm/scripts/include.R", temfile_prefix);
		new File(concorde_path).setExecutable(true, true);
		RFUtils.makeRMatrix(rf_file, out_prefix+".RData");
		String command = "Rscript "+tpR_path+" "
				+ "-i "+out_prefix+".RData "
				+ "-2 "+twopoint_file+" "
				+ (fixed_rf<0 ? "" : "-r "+fixed_rf+" ")
				+ "-o "+out_prefix+".tps "
				+ "--concorde "+new File(concorde_path).getParent()
				+ (RLibPath==null ? "" : " --include "+RLibPath);
		this.consume(this.bash(command));
		new File(temfile_prefix).delete();
	}

}
