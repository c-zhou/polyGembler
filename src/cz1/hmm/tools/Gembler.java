package cz1.hmm.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
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
	private String[] founder_haps;
	
	private int[] repeat = new int[]{30,30,10};
	private int refine_round = 10;
	
	private int drop = 1;
	private double phi = 2;
	private int nB = 10;
	
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
						
						+ " Recombination frequency estimation:\n"
						+ "		-nb/--best					The most likely nb haplotypes will be used (default 10).\n"
						+ " 	-phi/--skew-phi				For a haplotype inference, the frequencies of parental \n"
						+ "									haplotypes need to be in the interval [1/phi, phi], \n"
						+ "									otherwise will be discared (default 2).\n"
						+ " 	-nd/--drop					At least nd haplotype inferences are required for \n"
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
			myArgsEngine.add("-nb", "--best", true);
			myArgsEngine.add("-phi", "--skew-phi", true);
			myArgsEngine.add("-nd", "--drop", true);
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
			founder_haps = Constants._founder_haps.split(":");
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
			refine_round = Integer.parseInt(myArgsEngine.getString("-rr"));
		}
		
		if (myArgsEngine.getBoolean("-nb")) {
			nB = Integer.parseInt(myArgsEngine.getString("-nb"));
		}
		
		if (myArgsEngine.getBoolean("-phi")) {
			phi = Double.parseDouble(myArgsEngine.getString("-phi"));
		}
		
		if (myArgsEngine.getBoolean("-nd")) {
			drop = Integer.parseInt(myArgsEngine.getString("-nd"));
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
		/**
		new VCFtools(Constants._ploidy_H, min_depth, 
				max_depth, min_qual, 
				min_maf, max_missing, 
				in_vcf, out_prefix+"/data/"+prefix_vcf+".recode.vcf")
		.run();
		
		new DataPreparation(
				out_prefix+"/data/"+prefix_vcf+".recode.vcf",
				prefix_vcf+".recode",
				out_prefix+"/data/")
		.run();
		**/
		//#### STEP 02 single-point haplotype inferring
		final String in_zip = out_prefix+"/data/"+prefix_vcf+".recode.zip";
		final String out = out_prefix+"/single_hap_infer";
		final Map<String, Integer> scaffs = DataCollection.readScaff(in_zip);
		final String expr_id = new File(in_zip).getName().
				replaceAll(".zip$", "").
				replace(".", "").
				replace("_", "");
		Utils.makeOutputDir(out);
		//this.runHaplotyper(scaffs, expr_id, in_zip, out);
		
		//#### STEP 03 assembly errors
		final String metafile_prefix = out_prefix+"/meta/";
		Utils.makeOutputDir(metafile_prefix);
		final String ass_err_map = metafile_prefix+prefix_vcf;
		AssemblyError assemblyError = new AssemblyError (out, 
				ass_err_map,
				expr_id, 
				Constants._ploidy_H,
				founder_haps,
				THREADS,
				phi,
				drop,
				nB);
		/**
		assemblyError.run();
		final String prefix_vcf_assError = prefix_vcf+".recode.assError";
		Set<String> scaff_breakage = 
				assemblyError.split(out_prefix+"/data/"+prefix_vcf+".recode.vcf",
						metafile_prefix+prefix_vcf_assError+".vcf");
		if(!scaff_breakage.isEmpty()) {
			this.move(scaff_breakage, out, expr_id);

			new DataPreparation(
					metafile_prefix+prefix_vcf_assError+".vcf",
					prefix_vcf_assError,
					metafile_prefix)
			.run();
			String in_zip_assError = metafile_prefix+prefix_vcf_assError+".zip";
			final Map<String, Integer> scaffs_assError = DataCollection.readScaff(
					in_zip_assError);
			this.runHaplotyper(scaffs_assError, expr_id, in_zip_assError, out);
		}
		**/
		final String rf_prefix = metafile_prefix+prefix_vcf;
		RecombinationFreqEstimator recombinationFreqEstimator = 
				new RecombinationFreqEstimator (out, 
				rf_prefix,
				expr_id, 
				Constants._ploidy_H,
				founder_haps,
				THREADS,
				phi,
				drop,
				nB);
		recombinationFreqEstimator.run();
	
		//#### STEP 04 recombination frequency estimation
		
		
		//#### STEP 05 building superscaffolds (nearest neighbour joining)
		
		//#### STEP 06 multi-point hapotype inferring
		
		//#### STEP 07 recombination frequency estimation
		
		//#### STEP 08 genetic mapping
		
		//#### STEP 09 genetic map refinement
		
		//#### STEP 10 pseudo melecules construction
		
	}

	private void move(final Set<String> scaff_breakage, 
			final String out, 
			final String expr_id) {
		// TODO Auto-generated method stub
		String out_err = out+"/assembly_error";
		Utils.makeOutputDir(out_err);
		for(final String scaff : scaff_breakage) {
			File[] files = new File(out).listFiles(
					new FilenameFilter() {
						@Override
						public boolean accept(File dir, String name) {
							return name.matches("^"+expr_id+"."+scaff+".*");    
						}
					});
			for(File f : files)
				try {
					Files.move(f.toPath(), 
							Paths.get(out_err+"/"+f.getName()),
							StandardCopyOption.REPLACE_EXISTING,
							StandardCopyOption.ATOMIC_MOVE);
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
		}
	}

	private void runHaplotyper(final Map<String, Integer> scaffs,
			final String expr_id,
			final String in_zip,
			final String out) {
		// TODO Auto-generated method stub
		this.initial_thread_pool();
		for(final String scaff : scaffs.keySet()) {
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
									field,
									expr_id).run();
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
		this.waitFor();
	}

	
}
