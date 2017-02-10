package cz1.hmm.tools;

import java.util.Arrays;

import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Executor;
import cz1.util.IO;
import cz1.util.Constants.Field;

public class Gembler extends Executor {

	private static String in_vcf = null;
	private static String out_prefix = null;
	private static boolean vbt = false;
	private static int max_iter = 100;
	private static int ploidy = 2;
	private static Field field = Field.PL;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
							+" -i/--input					Input VCF file.\n"
							+" -o/--prefix					Output file location, create the directory if not exist.\n"
							+" -x/--max-iter				Maxmium rounds for EM optimization (default 100).\n"
							+" -p/--ploidy					Ploidy of genome (default 2).\n"
							+" -f/--parent					Parent samples (seperated by a \":\").\n"
							+" -G/--genotype				Use genotypes to infer haplotypes. Mutually exclusive with"
							+"								option -D/--allele-depth and -L/--genetype likelihood.\n"
							+" -D/--allele-depth			Use allele depth to infer haplotypes. Mutually exclusive "
							+"								with option -G/--genotype and -L/--genetype likelihood.\n"
							+" -L/--genotype-likelihood		Use genotype likelihoods to infer haplotypes. Mutually "
							+"								exclusive with option -G/--genotype and -L/--allele-depth "
							+"								(default).\n"
							+" -b/--segmental-kmeans		Use Vitering training instead of Baum-Welch algorithm.\n"
							+" -t/--threads					Threads (default 1).\n"			
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
			myArgsEngine.add("-t", "--threads", false);
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
			if(!out_prefix.endsWith(".zip")) out_prefix += ".zip";
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
		
		if(myArgsEngine.getBoolean("-t")) {
			THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
	
		IO.makeOutputDir(out_prefix);
		
		//#### STEP 01 filter SNPs and create ZIP file
		
		//#### STEP 02 single-point haplotype inferring
		
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
