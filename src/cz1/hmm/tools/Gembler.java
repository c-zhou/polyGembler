package cz1.hmm.tools;

import java.io.File;
import java.util.Arrays;

import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Executor;
import cz1.util.IO;
import cz1.util.Constants.Field;

public class Gembler extends Executor {

	private String in_vcf = null;
	private String out_prefix = null;
	private boolean vbt = false;
	private int max_iter = 100;
	private int ploidy = 2;
	private Field field = Field.PL;
	
	private int min_snpc = 5;
	private int min_depth = 0;
	private int max_depth = Integer.MAX_VALUE;
	private int min_qual = 0;
	private double min_maf = 0.1;
	private double max_missing = 0.5;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " Common:\n"
						+ "		-i/--input					Input VCF file.\n"
						+ "		-o/--prefix					Output file location, create the directory if not exist.\n"
						+ "		-p/--ploidy					Ploidy of genome (default 2).\n"
						+ " 	-S/--random-seed			Random seed for this run.\n"
						+ " 	-t/--threads				Threads (default 1).\n"	
						
						+ " Data preparation:\n"
						+ " 	-l/--min-depth				Minimum depth to keep a SNP (DP).\n"
						+ " 	-u/--max-depth				Maximum depth to keep a SNP (DP).\n"
						+ " 	-q/--min-qual  				Minimum quality to keep a SNP (QUAL).\n"
						+ " 	-mf/--min-maf				Minimum minor allele frequency to keep a SNP (default 0.1).\n"
						+ " 	-mm/--max-missing			Maximum proportion of missing data to keep a SNP (default 0.5).\n"
						
						+ " Haplotype inferring:\n"
						+ "		-x/--max-iter				Maxmium rounds for EM optimization (default 100).\n"
						+ "		-f/--parent					Parent samples (seperated by a \":\").\n"
						+ "		-G/--genotype				Use genotypes to infer haplotypes. Mutually exclusive with"
						+ "									option -D/--allele-depth and -L/--genetype likelihood.\n"
						+ "		-D/--allele-depth			Use allele depth to infer haplotypes. Mutually exclusive "
						+ "									with option -G/--genotype and -L/--genetype likelihood.\n"
						+ "		-L/--genotype-likelihood	Use genotype likelihoods to infer haplotypes. Mutually "
						+ "									exclusive with option -G/--genotype and -L/--allele-depth "
						+ "									(default).\n"
						+ " 	-b/--segmental-kmeans		Use Vitering training instead of Baum-Welch algorithm.\n"
						+ " 	-c/--min-snp-count			Minimum number of SNPs on a scaffold to run.\n"
				);
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		// create the command line parser

		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}

		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-i", "--input", true);
			myArgsEngine.add("-o", "--prefix", true);
			myArgsEngine.add("-x", "--max-iter", true);
			myArgsEngine.add("-p", "--ploidy", true);
			myArgsEngine.add("-f", "--parent", true);
			myArgsEngine.add("-G", "--genotype", false);
			myArgsEngine.add("-D", "--allele-depth", false);
			myArgsEngine.add("-L", "--genotype-likelihood", false);
			myArgsEngine.add("-b", "--segmental-kmeans", false);
			myArgsEngine.add("-S", "--random-seed", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-l", "--min-depth", true);
			myArgsEngine.add("-u", "--max-depth", true);
			myArgsEngine.add("-q", "--min-qual", true);
			myArgsEngine.add("-mf", "--min-maf", true);
			myArgsEngine.add("-mm", "--max-missing", true);
			myArgsEngine.add("-c", "--min-snp-count", true);
			myArgsEngine.parse(args);
		}
		
		if(myArgsEngine.getBoolean("-i")) {
			in_vcf = myArgsEngine.getString("-i");
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
		
		if(myArgsEngine.getBoolean("-x")) {
			max_iter = Integer.parseInt(myArgsEngine.getString("-x"));
		}
		
		if(myArgsEngine.getBoolean("-p")) {
			ploidy = Integer.parseInt(myArgsEngine.getString("-p"));
			Constants._ploidy_H = ploidy;
			Constants._haplotype_z = ploidy*2;
		}
		
		if(myArgsEngine.getBoolean("-f")) {
			Constants._founder_haps = myArgsEngine.getString("-f");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the parent samples (seperated by a \":\").");
		}
		
		int i = 0;
		if(myArgsEngine.getBoolean("-G")) {
			field = Field.GT;
			i++;
		}
		
		if(myArgsEngine.getBoolean("-D")) {
			field = Field.AD;
			i++;
		}
		
		if(myArgsEngine.getBoolean("-L")) {
			field = Field.PL;
			i++;
		}
		if(i>1) throw new RuntimeException("Options -G/--genotype, "
				+ "-D/--allele-depth, and -L/--genotype-likelihood "
				+ "are exclusive!!!");
		
		if(myArgsEngine.getBoolean("-b")) {
			vbt = true;
			throw new RuntimeException("Viterbi training not supported yet!!!");
		}
		
		if(myArgsEngine.getBoolean("-S")) {
			Constants.seed = Long.parseLong(myArgsEngine.getString("-S"));
			Constants.setRandomGenerator();
		}
		
		if(myArgsEngine.getBoolean("-t")) {
			THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
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
		
		if (myArgsEngine.getBoolean("-mf")) {
			min_maf = Double.parseDouble(myArgsEngine.getString("-mf"));
		}
		
		if (myArgsEngine.getBoolean("-mm")) {
			max_missing = Double.parseDouble(myArgsEngine.getString("-mm"));
		}
		
		if (myArgsEngine.getBoolean("-c")) {
			min_snpc = Integer.parseInt(myArgsEngine.getString("-c"));
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
	
		IO.makeOutputDir(out_prefix);
		String prefix_vcf = new File(in_vcf).getName().
				replaceAll(".vcf$", "");
		
		//#### STEP 01 filter SNPs and create ZIP file
		IO.makeOutputDir(out_prefix+"/data");
		VCFtools vcftools = new VCFtools(ploidy, min_depth, 
				max_depth, min_qual, 
				min_maf, max_missing, 
				in_vcf, out_prefix+"/data/"+prefix_vcf+".recode.vcf");
		vcftools.run();
		
		DataPreparation datapreparation = new DataPreparation(
				out_prefix+"/data/"+prefix_vcf+".recode.vcf",
				prefix_vcf+".recode",
				out_prefix+"/data/");
		datapreparation.run();
		
		//#### STEP 02 single-point haplotype inferring
		String in_zip = out_prefix+"/data/"+prefix_vcf+".recode.zip";
		String[] scaffs = IO.topLevelFolder(in_zip);
		String out = out_prefix+"/single_hap_infer";
		
		IO.makeOutputDir(out);
		for(int i=0; i<2; i++) {
			new Haplotyper(in_zip,
					out,
					new String[]{scaffs[2]},
					ploidy,
					Constants.Field.GL).run();
		}
		
		
		//#### STEP 03 assembly errors
		
		//#### STEP 04 single-point haplotype inferring (assembly errors) 
	
		//#### STEP 05 recombination frequency estimation
		
		//#### STEP 06 building superscaffolds (nearest neighbour joining)
		
		//#### STEP 07 multi-point hapotype inferring
		
		//#### STEP 08 recombination frequency estimation
		
		//#### STEP 09 genetic mapping
		
		//#### STEP 10 genetic map refinement
		
		//#### STEP 11 pseudo melecules construction
		
	}

}
