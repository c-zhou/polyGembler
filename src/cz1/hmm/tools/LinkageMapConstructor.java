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
	private double lod_thres = 3;
	private int ns;
	private boolean one = false;
	private boolean two = false;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " Common:\n"
						+ "     -i/--rf                     Recombination frequency file.\n"
						+ "     -m/--map                    Recombination map file.\n"
						+ "     -n/--hap-size               #haplotypes (popsize*ploidy). \n"
						+ "     -l/--lod                    LOD score threshold (default: 3).\n"
						+ "     -1/--one-group              One group. \n"
						+ "     -2/--check-group            Re-check linkage groups. \n"
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
			myArgsEngine.add("-n", "--pop-size", true);
			myArgsEngine.add("-l", "--lod", true);
			myArgsEngine.add("-1", "--one-group", false);
			myArgsEngine.add("-2", "--check-group", false);
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
		
		if(myArgsEngine.getBoolean("-n")) {
			ns = Integer.parseInt(myArgsEngine.getString("-n"));
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the number of haplotypes (popsize*ploidy).");
		}
		
		if(myArgsEngine.getBoolean("-l")) {
			lod_thres = Double.parseDouble(myArgsEngine.getString("-l"));
		}
		
		if(myArgsEngine.getBoolean("-1")) {
			one = true;
		}
		
		if(myArgsEngine.getBoolean("-2")) {
			two = true;
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
		final String command = 
				"Rscript "+mklgR_path+" "
				+ "-i "+out_prefix+".RData "
				+ "-m "+map_file+" "
				+ "-r "+calcRfFromLOD(lod_thres, ns)+" "
				+ (one?"-1 ":" ")
				+ (two?"-2 ":" ")
				+ "-o "+out_prefix+" "
				+ "--concorde "+new File(concorde_path).getParent()
				+ (RLibPath==null ? "" : " --include "+RLibPath);
		this.consume(this.bash(command));
		
		new File(temfile_prefix).delete();
	}
	
	private static double calcRfFromLOD(double lod_thres, int n) {
		// TODO Auto-generated method stub
		double lb = 0, ub = 0.5;
		double rf;
		double lod;
		while(true) {
			rf = (lb+ub)/2;
			lod = calcLODFromRF(rf, n);
			if(lod<lod_thres) {
				ub = rf;
			} else if(lod-lod_thres<1e-6) {
				return rf;
			} else {
				lb = rf;
			}
		}
	}
	
	private static double calcLODFromRF(double theta, int n) {
		// TODO Auto-generated method stub
		double r = n*theta;
		return (n-r)*Math.log10(1-theta)+r*Math.log10(theta)-n*Math.log10(0.5);
	}
}
