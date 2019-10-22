package cz1.hmm.tools;

import java.io.File;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

public class MappingAnalysis extends Executor {

	private String rf_file = null;
	private String map_file = null;
	private String RLibPath = null;
	private String out_prefix = null;
	private double lod_thres = 3;
	private double rf_thres = RFUtils.inverseGeneticDistance(0.5, "kosambi");
	private boolean one = false;
	private boolean check_chimeric = false;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--input                  Recombination frequency file.\n"
						+ " -m/--map                    Recombination map file.\n"
						+ " -l/--lod                    LOD score threshold (default: 3).\n"
						+ " -r/--rf                     Recombination frequency threshold (default: 0.38).\n"
						+ " -1/--one-group              Keep all in one group. Do not do clustering (default: false). \n"
						+ " -c/--check-chimeric         Check chimeric joins during grouping (default: false). \n"
						+ " -rlib/--R-external-libs     R external library path.\n"
						+ " -t/--threads                #threads (default 1).\n"
						+ " -o/--prefix                 Output file prefix.\n\n"
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
			myArgsEngine.add("-i", "--input", true);
			myArgsEngine.add("-m", "--map", true);
			myArgsEngine.add("-l", "--lod", true);
			myArgsEngine.add("-r", "--rf", true);
			myArgsEngine.add("-1", "--one-group", false);
			myArgsEngine.add("-c", "--check-chimeric", false);
			myArgsEngine.add("-rlib", "--R-external-libs", true);
			myArgsEngine.add( "-t", "--threads", true);
			myArgsEngine.add("-o", "--prefix", true);
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

		if(myArgsEngine.getBoolean("-l")) {
			lod_thres = Double.parseDouble(myArgsEngine.getString("-l"));
		}
		
		if(myArgsEngine.getBoolean("-r")) {
			rf_thres = Double.parseDouble(myArgsEngine.getString("-r"));
		}
		
		if(myArgsEngine.getBoolean("-1")) {
			one = true;
		}
		
		if(myArgsEngine.getBoolean("-c")) {
			check_chimeric = true;
		}

		if (myArgsEngine.getBoolean("-rlib")) {
			RLibPath = myArgsEngine.getString("-rlib");
		}
		
		if(myArgsEngine.getBoolean("-t")) {
			THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}
		
		if(myArgsEngine.getBoolean("-o")) {
			out_prefix = myArgsEngine.getString("-o");
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify your output file prefix.");
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		final String temfile_prefix = Utils.makeTempDir();
		final boolean tmpdirCreated = Utils.makeOutputDir(new File(temfile_prefix));
		final String concorde_path =
				RFUtils.makeExecutable("cz1/hmm/executable/concorde", temfile_prefix);
		new File(concorde_path).setExecutable(true, true);
		final String mklgR_path =
				RFUtils.makeExecutable("cz1/hmm/scripts/make_geneticmap.R", temfile_prefix);
		RFUtils.makeExecutable("cz1/hmm/scripts/include.R", temfile_prefix);
		RFUtils.makeRMatrix(rf_file, out_prefix+".RData");
		myLogger.info("Using LOD score threshold: "+lod_thres+".");
		myLogger.info("Using recombination frequency threshold: "+rf_thres+".");
		final String command =
				"Rscript "+mklgR_path+" "
						+ "-i "+out_prefix+".RData "
						+ "-m "+map_file+" "
						+ "-r "+rf_thres+" "
						+ "-l "+lod_thres+" "
						+ (one?"-1 ":" ")
						+ (check_chimeric?"-x ":" ")
						+ "-p "+THREADS+" "
						+ "-o "+out_prefix+" "
						+ "--concorde "+new File(concorde_path).getParent()+" "
						+ (RLibPath==null ? "" : "--include "+RLibPath+" ")
						+ "--tmpdir "+new File(temfile_prefix).getAbsolutePath();
		this.consume(this.bash(command));

		if(tmpdirCreated) Utils.deleteDirectory(new File(temfile_prefix));
	}
}




