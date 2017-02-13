package cz1.hmm.tools;

import java.io.File;
import java.util.Arrays;
import java.util.Map;
import java.util.concurrent.TimeUnit;

import cz1.hmm.data.DataCollection;
import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Executor;
import cz1.util.Utils;
import cz1.util.Constants.Field;

public class Gembler extends Executor {

	private String in_vcf = null;
	private String out_prefix = null;
	private boolean vbt = false;
	private int max_iter = 100;
	private Field field = Field.PL;
	
	private int min_snpc = 5;
	private int min_depth = 0;
	private int max_depth = Integer.MAX_VALUE;
	private int min_qual = 0;
	private double min_maf = 0.1;
	private double max_missing = 0.5;
	
	private int[] repeat = new int[]{30,30,10};
	private int refine_repeat = 10;
	
	private int nB = 10;
	private int phi = 2;
	
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
						+ " 	-r/--repeat					Repeat haplotype inferring for multiple times as EM algorithm \n"
						+ "									could be trapped in local optima. The program takes three values \n"
						+ "									from here, i.e., for scaffold, for superscaffold and for genetic \n"
						+ "									linkage map refinement. Three values should be seperated by ',' \n"
						+ "									(default 30,30,10 ).\n"
						+ "		-rr/--refinement-round		Number of rounds to refine pseudomelecules (default 10.)\n"
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
			myArgsEngine.add("-r", "--repeat", true);
			myArgsEngine.add("-rr", "--refinement-round", true);
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
			Constants.ploidy(Integer.parseInt(myArgsEngine.getString("-p")));
		}
		
		if(myArgsEngine.getBoolean("-f")) {
			Constants._founder_haps = myArgsEngine.getString("-f");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the parent samples (seperated by a \":\").");
		}
		
		int c = 0;
		if(myArgsEngine.getBoolean("-G")) {
			field = Field.GT;
			c++;
		}
		
		if(myArgsEngine.getBoolean("-D")) {
			field = Field.AD;
			c++;
		}
		
		if(myArgsEngine.getBoolean("-L")) {
			field = Field.PL;
			c++;
		}
		if(c>1) throw new IllegalArgumentException("Options -G/--genotype, "
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
		
		if (myArgsEngine.getBoolean("-r")) {
			String[] rs = myArgsEngine.getString("-r").split(",");
			if(rs.length!=3) 
				throw new IllegalArgumentException("Option -r/--repeat should "
						+ "be 3 integers seperated by \",\".");
			for(int i=0; i<3; i++)
				repeat[i] = Integer.parseInt(rs[i]);
		}
		
		if (myArgsEngine.getBoolean("-rr")) {
			refine_repeat = Integer.parseInt(myArgsEngine.getString("-rr"));
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
	
		Utils.makeOutputDir(out_prefix);
		String prefix_vcf = new File(in_vcf).getName().
				replaceAll(".vcf$", "");
		
		//#### STEP 01 filter SNPs and create ZIP file
		Utils.makeOutputDir(out_prefix+"/data");
		VCFtools vcftools = new VCFtools(Constants._ploidy_H, min_depth, 
				max_depth, min_qual, 
				min_maf, max_missing, 
				in_vcf, out_prefix+"/data/"+prefix_vcf+".recode.vcf");
		//vcftools.run();
		
		DataPreparation datapreparation = new DataPreparation(
				out_prefix+"/data/"+prefix_vcf+".recode.vcf",
				prefix_vcf+".recode",
				out_prefix+"/data/");
		//datapreparation.run();
		
		//#### STEP 02 single-point haplotype inferring
		final String in_zip = out_prefix+"/data/"+prefix_vcf+".recode.zip";
		final String out = out_prefix+"/single_hap_infer";
		final Map<String, Integer> scaffs = DataCollection.readScaff(in_zip);
		Utils.makeOutputDir(out);
		
		this.initial_thread_pool();
		
		//for(final String scaff : scaffs.keySet()) {
		for(final String scaff : new String[]{"scaffold7_size520567"}) {
			if(scaffs.get(scaff)<min_snpc) continue;
			
			for(int i=0; i<repeat[0]; i++) {
				executor.submit(new Runnable(){

					@Override
					public void run() {
						// TODO Auto-generated method stub
						try {
							new Haplotyper(in_zip,
									out,
									new String[]{scaff},
									Constants._ploidy_H,
									field).run();
						} catch (Exception e) {
							Thread t = Thread.currentThread();
							t.getUncaughtExceptionHandler().uncaughtException(t, e);
							e.printStackTrace();
							executor.shutdown();
							System.exit(1);
						}
					}
				});
			}
		}
		
		try {
			executor.shutdown();
			executor.awaitTermination(365, TimeUnit.DAYS);
		} catch (InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
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
