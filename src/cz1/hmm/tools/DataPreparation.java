package cz1.hmm.tools;

import java.io.File;

import cz1.util.ArgsEngine;
import cz1.util.Executor;

public class DataPreparation extends Executor {

	private int ploidy = 2;
	private int min_depth = 0;
	private int max_depth = Integer.MAX_VALUE;
	private int min_qual = 0;
	private double min_maf = 0.1;
	private double max_missing = 0.5;
	private String vcf_in = null;
	String id = null;
	String out_file = null;
	
	public DataPreparation() {}
	
	public DataPreparation(String vcf_in,
			int ploidy, 
			int min_depth, 
			int max_depth, 
			int min_qual, 
			double min_maf, 
			double max_missing,
			String id,
			String out_file) {
		this.vcf_in = vcf_in;
		this.ploidy = ploidy;
		this.min_depth = min_depth;
		this.max_depth = max_depth;
		this.min_qual = min_qual;
		this.min_maf = min_maf;
		this.max_missing = max_missing;
		this.id = id;
		this.out_file = out_file;
	}
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--vcf			Input VCF file.\n"
						+ " -s/--id             Unique id of this run (default: input VCF file name prefix).\n"
						+ " -p/--ploidy			Ploidy of genome (default 2). \n"
						+ "                     NOTE: If you called variant as diploid, then the program will \n"
						+ "                     fit a binomial model to call genotypes and genotype qualities \n"
						+ "                     from allele depth with the ploidy specified here.\n"
						+ " -l/--min-depth		Minimum depth to keep a SNP (DP).\n"
						+ " -u/--max-depth		Maximum depth to keep a SNP (DP).\n"
						+ " -q/--min-qual  		Minimum quality to keep a SNP (QUAL).\n"
						+ " -f/--min-maf		Minimum minor allele frequency to keep a SNP (default 0.1).\n"
						+ " -m/--max-missing	Maximum proportion of missing data to keep a SNP (default 0.5).\n"
						+ " -o/--prefix			Prefix for output files (default: input VCF file folder).\n\n");
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
			myArgsEngine.add("-p", "--ploidy", true);
			myArgsEngine.add("-u", "--max-depth", true);
			myArgsEngine.add("-q", "--min-qual", true);
			myArgsEngine.add("-f", "--min-maf", true);
			myArgsEngine.add("-m", "--max-missing", true);
			myArgsEngine.add("-o", "--prefix", true);
			myArgsEngine.parse(args);
		}
		
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
		
		if (myArgsEngine.getBoolean("-p")) {
			ploidy = Integer.parseInt(myArgsEngine.getString("-p"));
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
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		new VCFtools(ploidy, min_depth, 
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
