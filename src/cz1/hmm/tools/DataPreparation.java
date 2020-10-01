package cz1.hmm.tools;

import java.io.File;

import cz1.util.ArgsEngine;
import cz1.util.Executor;

public class DataPreparation extends Executor {

	private int min_depth = 0;
	private int max_depth = Integer.MAX_VALUE;
	private int min_qual = 0;
	private double min_maf = 0;
	private double max_missing = 1.0;
	private String vcf_in = null;
	String id = null;
	String out_file = null;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--vcf            Input VCF file.\n"
						+ " -s/--id             Unique id of this run (default: input VCF file name prefix).\n"
						+ " -l/--min-depth      Minimum depth to keep a SNP (0).\n"
						+ " -u/--max-depth      Maximum depth to keep a SNP ("+Integer.MAX_VALUE+").\n"
						+ " -q/--min-qual       Minimum quality to keep a SNP (0).\n"
						+ " -f/--min-maf        Minimum minor allele frequency to keep a SNP (default 0).\n"
						+ " -m/--max-missing    Maximum proportion of missing data to keep a SNP (default 1.0).\n"
						+ " -o/--prefix         Prefix for output files (default: input VCF file folder).\n\n"
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
			myArgsEngine.add("-i", "--vcf", true);
			myArgsEngine.add("-s", "--id", true);
			myArgsEngine.add("-l", "--min-depth", true);
			myArgsEngine.add("-u", "--max-depth", true);
			myArgsEngine.add("-q", "--min-qual", true);
			myArgsEngine.add("-f", "--min-maf", true);
			myArgsEngine.add("-m", "--max-missing", true);
			myArgsEngine.add("-o", "--prefix", true);
		}
		myArgsEngine.parse(args);
		
		if (myArgsEngine.getBoolean("-i")) {
			vcf_in = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your VCF file.");
		}
	
		if (myArgsEngine.getBoolean("-s")) {
			id = myArgsEngine.getString("-s");
		} else {
			id = new File(vcf_in).getName().replaceAll(".vcf.gz$", "").
					replaceAll(".vcf$", "")+".recode";
		}
		
		if (myArgsEngine.getBoolean("-l")) {
			min_depth = Integer.parseInt(myArgsEngine.getString("-l"));
		}
		
		if (myArgsEngine.getBoolean("-u")) {
			max_depth = Integer.parseInt(myArgsEngine.getString("-u"));
		}
		
		if (myArgsEngine.getBoolean("-q")) {
			min_qual = Integer.parseInt(myArgsEngine.getString("-q"));
		}
		
		if (myArgsEngine.getBoolean("-f")) {
			min_maf = Double.parseDouble(myArgsEngine.getString("-f"));
		}
		
		if (myArgsEngine.getBoolean("-m")) {
			max_missing = Double.parseDouble(myArgsEngine.getString("-m"));
		}
		
		if (myArgsEngine.getBoolean("-o")) {
			out_file = myArgsEngine.getString("-o");
		} else {
			out_file = new File(vcf_in).getParent();
			if(out_file==null) out_file = "./";
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		new VCFtools(min_depth, 
				max_depth, min_qual, 
				min_maf, max_missing, 
				vcf_in, out_file+"/"+id+".vcf")
		.run();
		
		new ZipDataCollection(
				out_file+"/"+id+".vcf",
				id,
				out_file)
		.run();
	}
}
