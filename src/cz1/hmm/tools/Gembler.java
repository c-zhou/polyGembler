package cz1.hmm.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;

import cz1.hmm.data.DataCollection;
import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Executor;
import cz1.util.Utils;

public class Gembler extends Executor {
	// universal
	private String in_vcf = null;
	private String out_prefix = null;
	private int ploidy = 2;
	private String RLibPath = null;
	
	// data preparation
	private int min_depth = 0;
	private int max_depth = Integer.MAX_VALUE;
	private int min_qual = 0;
	private double min_maf = 0;
	private double max_missing = 1.0;
	
	// haplotype inferring and pseudomolecule refinement
	private String parents = null;
	private int max_iter = 1000;
	private String field = "-D";
	private int min_snpc = 5;
	private int[] nr = new int[]{30,30,10};
	private int refine_round = 3;
	
	// recombination frequency estimation and assembly error detection
	private double err_rf = 0.1;
	private int wbp = 30000;
	private int wnm = 30;
	private int nb = 30;	
	private int drop = 1;
	private double phi = 2;

	// genetic linkage map construction
	private double lod_thresh = 3.0;
	private double lg_rf = RFUtils.inverseGeneticDistance(0.5, "kosambi");
	private boolean check_chimeric = false;
	
	// uperscaffold construction by nearest neighbor joining
	private boolean ss = true;
	private double ss_rf = RFUtils.inverseGeneticDistance(0.3, "kosambi");
	
	// pseudomolecule construction
	private String contig_file = null;
	private int n_gap = 1000;
	
	public Gembler() {
		//this.require("Rscript");
	}
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " Common:\n"
						+ "     -i/--input-vcf              Input VCF file.\n"
						+ "     -o/--prefix                 Output file location, create the directory if not exist.\n"
						+ "     -p/--ploidy                 Ploidy of genome (default 2).\n"
						+ "     -t/--threads                Threads (default 1).\n"
						+ "     -rlib/--R-external-libs     External library paths that you want R to search for packages.\n"
						+ "                                 This could be useful if you are not root users and install R \n"
						+ "                                 packages in directories other than default. \n"
						+ "                                 Multiple paths separated by ':' could be provided.\n"
						+ "     -seed/--random-seed         Random seed for this run.\n"
						+ "\n"
						+ " Data preparation:\n"
						+ "     -l/--min-depth              Minimum depth to keep a SNP (0).\n"
						+ "     -u/--max-depth              Maximum depth to keep a SNP ("+Integer.MAX_VALUE+").\n"
						+ "     -q/--min-qual               Minimum quality to keep a SNP (0).\n"
						+ "     -f/--min-maf                Minimum minor allele frequency to keep a SNP (default 0).\n"
						+ "     -m/--max-missing            Maximum proportion of missing data to keep a SNP (default 1.0).\n"
						+ "\n"
						+ " Haplotype inferring and pseudomolecule refinement:\n"
						+ "     -x/--max-iter               Maxmium rounds for EM optimization (default 1000).\n"
						+ "     -parent/--parent            Parent samples (separated by a \":\").\n"
						+ "     -G/--genotype               Use genotypes to infer haplotypes. Mutually exclusive with \n"
						+ "                                 option -D/--allele-depth.\n"
						+ "     -D/--allele-depth           Use allele depth to infer haplotypes. Mutually exclusive \n"
						+ "                                 with option -G/--genotype (default).\n"
						+ "     -c/--min-snp-count          Minimum number of SNPs on a scaffold to run (default 5).\n"
						+ "     -r/--repeat                 Repeat haplotype inferring for multiple times as EM algorithm \n"
						+ "                                 could be trapped in local optima. The program takes three values \n"
						+ "                                 from here, i.e., for scaffold, for superscaffold and for genetic \n"
						+ "                                 linkage map refinement. Three values should be separated by ',' \n"
						+ "                                 (default 30,30,10).\n"
						+ "     -rr/--refinement-round      Number of rounds to refine pseudomelecules (default 3).\n"
						+ "\n"
						+ " Recombination frequency estimation and assembly error detection:\n"
						+ "     -asmr/--asmr-thresh         Recombination frequency threshold for assembly error detection (default 0.1).\n"
						+ "     -wbp/--windows-bp           Window size (#basepairs) for RF estimation for marker pairs on \n"
						+ "                                 same scaffold (default 30000).\n"
						+ "     -wnm/--windows-nm           Window size (#markers) for RF estimation for marker pairs on \n" 
						+ "                                 same scaffold (default 30).\n"
						+ "     -nb/--best                  The most likely nb haplotypes will be used (default 30). \n"
						+ "     -phi/--skew-phi             For a haplotype inference, the frequencies of parental \n"
						+ "                                 haplotypes need to be in the interval [1/phi, phi], \n"
						+ "                                 otherwise will be discared (default 2).\n"
						+ "     -nd/--drop                  At least nd haplotype inferences are required for \n"
						+ "                                 a contig/scaffold to be analysed (default 1).\n"
						+ "\n"
						+ " Genetic linkage map construction:\n"
						+ "     -l/--lod                    LOD score threshold (default: 3).\n"
						+ "     -lr/--lr-thresh             Recombination frequency threshold for linkage mapping (default: 0.38).\n"
						+ "     -c/--check-chimeric         Check chimeric joins during grouping (default: false). \n"
						+ "\n"
						+ " Superscaffold construction by nearest neighbor joining:\n"
						+ "     -ns/--no-superscaffold      Do NOT construct superscaffold for RF refinement.\n"
                        + "     -sr/--sr-thresh             Recombination frequency threshold for superscaffold (default: 0.27).\n"
                        + "\n"
						+ " Pseudomolecule construction:\n"
						+ "     -a/--contig-file            Input assembly FASTA file.\n"
						+ "     -g/--gap                    Gap size between sequences (default 1000). \n"
						+ "                                 The gaps will be filled with character 'n'.\n"
						+ "\n"
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
			myArgsEngine.add("-i", "--input-vcf", true);
			myArgsEngine.add("-o", "--prefix", true);
			myArgsEngine.add("-p", "--ploidy", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-rlib", "--R-external-libs", true);
			myArgsEngine.add("-seed", "--random-seed", true);
			
			myArgsEngine.add("-l", "--min-depth", true);
			myArgsEngine.add("-u", "--max-depth", true);
			myArgsEngine.add("-q", "--min-qual", true);
			myArgsEngine.add("-f", "--min-maf", true);
			myArgsEngine.add("-m", "--max-missing", true);
			
			myArgsEngine.add("-x", "--max-iter", true);
			myArgsEngine.add("-parent", "--parent", true);
			myArgsEngine.add("-G", "--genotype", false);
			myArgsEngine.add("-D", "--allele-depth", false);
			myArgsEngine.add("-c", "--min-snp-count", true);
			myArgsEngine.add("-r", "--repeat", true);
			myArgsEngine.add("-rr", "--refinement-round", true);
			
			myArgsEngine.add("-asmr", "--asmr-thresh", true);
			myArgsEngine.add("-wbp", "--windows-bp", true);
			myArgsEngine.add("-wnm", "--windows-nm", true);
			myArgsEngine.add("-nb", "--best", true);
			myArgsEngine.add("-phi", "--skew-phi", true);
			myArgsEngine.add("-nd", "--drop", true);
			
			myArgsEngine.add("-l", "--lod", true);
			myArgsEngine.add("-lr", "--lr-thresh", true);
			myArgsEngine.add("-c", "--check-chimeric", true);
			
			myArgsEngine.add("-sr", "--sr-thresh", true);
			
			myArgsEngine.add("-a", "--contig-file", true);
			myArgsEngine.add("-g", "--gap", true);
			
			myArgsEngine.parse(args);
		}
		
		if(myArgsEngine.getBoolean("-i")) {
			in_vcf = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your input VCF file.");
		}
		
		if(myArgsEngine.getBoolean("-o")) {
			out_prefix = myArgsEngine.getString("-o");
			if(new File(out_prefix).exists()) {
				throw new RuntimeException("Output dir "+out_prefix+" existed.");
			}
		}  else {
			String timeStamp = new SimpleDateFormat("yyyy-MM-dd-HH:mm:ss").
					format(Calendar.getInstance().getTime());
			out_prefix = "results_"+timeStamp;
			myLogger.info("Using output dir: "+out_prefix);
		}
		
		if(myArgsEngine.getBoolean("-p")) {
			this.ploidy = Integer.parseInt(myArgsEngine.getString("-p"));
		}
		
		if(myArgsEngine.getBoolean("-t")) {
			THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}
		
		if (myArgsEngine.getBoolean("-rlib")) {
			RLibPath = myArgsEngine.getString("-rlib");
		}
		
		if(myArgsEngine.getBoolean("-seed")) {
			Constants.seed = Long.parseLong(myArgsEngine.getString("-seed"));
			Constants.setRandomGenerator();
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
		
		if(myArgsEngine.getBoolean("-x")) {
			max_iter = Integer.parseInt(myArgsEngine.getString("-x"));
		}
		
		if(myArgsEngine.getBoolean("-parent")) {
			parents = myArgsEngine.getString("-parent");
		}
		
		int c = 0;
		if(myArgsEngine.getBoolean("-G")) {
			field = "-G";
			c++;
		}
		
		if(myArgsEngine.getBoolean("-D")) {
			field = "-D";
			c++;
		}
		
		if(c>1) throw new IllegalArgumentException("Options -G/--genotype and "
				+ "-D/--allele-depth are mutually exclusive.");
		
		if (myArgsEngine.getBoolean("-c")) {
			min_snpc = Integer.parseInt(myArgsEngine.getString("-c"));
		}
		
		if (myArgsEngine.getBoolean("-r")) {
			String[] rs = myArgsEngine.getString("-r").split(",");
			if(rs.length!=3) 
				throw new IllegalArgumentException("Option -r/--repeat should "
						+ "be three integers seperated by \",\".");
			for(int i=0; i<3; i++)
				nr[i] = Integer.parseInt(rs[i]);
		}
		
		if (myArgsEngine.getBoolean("-rr")) {
			refine_round = Integer.parseInt(myArgsEngine.getString("-rr"));
		}
		
		if (myArgsEngine.getBoolean("-asmr")) {
			err_rf = Double.parseDouble(myArgsEngine.getString("-asmr"));
		}
		
		if (myArgsEngine.getBoolean("-wbp")) {
			wbp = Integer.parseInt(myArgsEngine.getString("-wbp"));
		}
		
		if (myArgsEngine.getBoolean("-wnm")) {
			wnm = Integer.parseInt(myArgsEngine.getString("-wnm"));
		}
		
		if (myArgsEngine.getBoolean("-nb")) {
			nb = Integer.parseInt(myArgsEngine.getString("-nb"));
		}
		
		if (myArgsEngine.getBoolean("-phi")) {
			phi = Double.parseDouble(myArgsEngine.getString("-phi"));
		}
		
		if (myArgsEngine.getBoolean("-nd")) {
			drop = Integer.parseInt(myArgsEngine.getString("-nd"));
		}
		
		if (myArgsEngine.getBoolean("-l")) {
			lod_thresh = Double.parseDouble(myArgsEngine.getString("-l"));
		}
		
		if (myArgsEngine.getBoolean("-lr")) {
			lg_rf = Double.parseDouble(myArgsEngine.getString("-lr"));
		}
		
		if (myArgsEngine.getBoolean("-c")) {
			check_chimeric = true;
		}
		
		if (myArgsEngine.getBoolean("-ns")) {
			ss = false;
		}
		
		if (myArgsEngine.getBoolean("-sr")) {
			ss_rf = Double.parseDouble(myArgsEngine.getString("-sr"));
		}
		
		if (myArgsEngine.getBoolean("-a")) {
			contig_file = myArgsEngine.getString("-a");
		}
		
		if (myArgsEngine.getBoolean("-g")) {
			n_gap = Integer.parseInt(myArgsEngine.getString("-g"));
		}
		
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		
		Utils.makeOutputDir(new File(out_prefix));
		String datPref = "out1";

		//#### STEP 01 prepare data
		myLogger.info("STEP 01 prepare data");
		
		DataPreparation dp = new DataPreparation();
		dp.setParameters(new String[] {
				"-i", in_vcf,
				"-s", datPref,
				"-l", String.valueOf(min_depth),
				"-u", String.valueOf(max_depth),
				"-q", String.valueOf(min_qual),
				"-f", String.valueOf(min_maf),
				"-m", String.valueOf(max_missing),
				"-o", out_prefix 
		});
		dp.run();
		
		//#### STEP 02 infer single-point haplotypes
		myLogger.info("STEP 02 infer single-point haplotypes");
		
		final String expr_id = "zzz";
		String in_zip = out_prefix+"/"+datPref+".zip";
		Map<String, Integer> scaffStats = DataCollection.readScaff(in_zip);
		String outs = out_prefix+"/h1";
		Utils.makeOutputDir(new File(outs));
		this.runHaplotyper(scaffStats, expr_id, in_zip, nr[0], outs);
		
		//#### STEP 03 detect assembly errors
		myLogger.info("STEP 03 detect assembly errors");
		
		datPref = "out2";
		String outPref = out_prefix+"/"+datPref;
		AssemblyError asmerr = new AssemblyError();
		asmerr.setParameters(new String[] {
				"-i", outs,
				"-o", outPref,
				"-vi", in_vcf,
				"-r", String.valueOf(err_rf),
				"-wbp", String.valueOf(wbp),
				"-wnm", String.valueOf(wnm),
				"-ex", expr_id,
				"-nb", String.valueOf(nb),
				"-phi", String.valueOf(phi),
				"-nd", String.valueOf(drop),
				"-t", String.valueOf(THREADS)
		});
		asmerr.run();
		
		// move haplotype phasing results for old scaffs if misassembled
		Set<String> oldScaffs = asmerr.errs.keySet();
		if(!oldScaffs.isEmpty()) {
			String outs_err = out_prefix+"/herr";
			Utils.makeOutputDir(new File(outs_err));
			move(oldScaffs, outs_err, expr_id);
		}
		
		// run haplotyper for new scaffs
		Set<String> newScaffs = asmerr.newScaffs;
		if(!newScaffs.isEmpty()) {
			// create a new zip file to include new scaffs
			// first add contents of the old VCF to new VCF
			String tmp_vcf = outPref+"_tmp.vcf";
			Utils.renameFile(outPref+".vcf", tmp_vcf);
			pendVCFContent(in_vcf, tmp_vcf);
			dp.setParameters(new String[] {
					"-i", tmp_vcf,
					"-s", datPref,
					"-l", String.valueOf(min_depth),
					"-u", String.valueOf(max_depth),
					"-q", String.valueOf(min_qual),
					"-f", String.valueOf(min_maf),
					"-m", String.valueOf(max_missing),
					"-o", out_prefix 
			});
			dp.run();
			Utils.deleteFile(tmp_vcf);
			
			// now run haplotyper for new scaffs
			in_zip = outPref+".zip";
			Map<String, Integer> newScaffStats = DataCollection.readScaff(
					in_zip, newScaffs);
			this.runHaplotyper(newScaffStats, expr_id, in_zip, nr[0], outs);
		}

		//#### STEP 04 build superscaffolds using nearest neighbor joining
		if(ss) {
			myLogger.info("STEP 04 build superscaffolds using nearest neighbor joining");
			
			// recombination frequency estimation
			outPref = out_prefix+"/out3";
			TwoPointAnalysis tp = new TwoPointAnalysis();
			tp.setParameters(new String[] {
					"-i", outs,
					"-o", outPref,
					"-ex", expr_id,
					"-nb", String.valueOf(nb),
					"-phi", String.valueOf(phi),
					"-nd", String.valueOf(drop),
					"-t", String.valueOf(THREADS)
			});
			tp.run();
			
			// build superscaffolds
			NNsuperscaffold nns = new NNsuperscaffold();
			nns.setParameters(new String[] {
					"-i", outPref+".txt",
					"-r", String.valueOf(ss_rf),
					"-o", outPref
			});
			nns.run();
			
			// multi-point hapotype inferring
			outs = out_prefix+"/h2";
			Utils.makeOutputDir(new File(outs));
			Set<String> nn_scaffs = new HashSet<String>();
			Map<String, String> nn_separation = new HashMap<String, String>();
			Map<String, String> nn_reverse = new HashMap<String, String>();
			this.readSS(outPref+".nns", nn_scaffs, nn_separation, nn_reverse);
			this.runHaplotyper(nn_scaffs, nn_separation, nn_reverse,
					expr_id, in_zip, nr[1], outs);
		}
		
		//#### STEP 05 estimate recombination frequencies
		myLogger.info("STEP 05 estimate recombination frequencies");
		
		outPref = out_prefix+"/out4";
		SinglePointAnalysis sp = new SinglePointAnalysis();
		sp.setParameters(new String[] {
				"-i", outs,
				"-o", outPref,
				"-wbp", String.valueOf(wbp),
				"-wnm", String.valueOf(wnm),
				"-ex", expr_id,
				"-nb", String.valueOf(nb),
				"-phi", String.valueOf(phi),
				"-nd", String.valueOf(drop),
				"-t", String.valueOf(THREADS)
		});
		sp.run();
		
		TwoPointAnalysis tp = new TwoPointAnalysis();
		tp.setParameters(new String[] {
				"-i", outs,
				"-o", outPref,
				"-ex", expr_id,
				"-nb", String.valueOf(nb),
				"-phi", String.valueOf(phi),
				"-nd", String.valueOf(drop),
				"-t", String.valueOf(THREADS)
		});
		tp.run();
		
		//#### STEP 06 genetic mapping
		myLogger.info("STEP 06 genetic mapping");
		
		MappingAnalysis ma = new MappingAnalysis();
		String[] args1 = new String[] {
				"-i", outPref+".txt",
				"-m", outPref+".map",
				"-l", String.valueOf(lod_thresh),
				"-r", String.valueOf(lg_rf),
				"-rlib", RLibPath, 
				"-t", String.valueOf(THREADS),
				"-o", outPref
		};
		if(check_chimeric) {
			String[] args2 = args1;
			args1 = new String[args2.length+1];
			System.arraycopy(args2, 0, args1, 0, args2.length);
			args1[args1.length-1] = "-c";
		}
		ma.setParameters(args1);
		ma.run();
		
		//#### STEP 07 refine genetic maps
		if(refine_round>0) {
			myLogger.info("STEP 07 refine genetic maps");
			
			String lgOutDir = out_prefix+"/hlg";
			Utils.makeOutputDir(new File(lgOutDir));
			
			Set<String> nn_scaffs = new HashSet<String>();
			Map<String, String> nn_separation = new HashMap<String, String>();
			Map<String, String> nn_reverse = new HashMap<String, String>();
			this.readSS(outPref+".par", nn_scaffs, nn_separation, nn_reverse);

			final int lgN = nn_scaffs.size();
			for(int i=0; i<lgN; i++) {
				String lgOutDir_i = lgOutDir+"/"+i;
				Utils.makeOutputDir(new File(lgOutDir_i));
				for(int j=0; j<this.refine_round; j++) {
					String lgOutDir_ij = lgOutDir_i+"/"+j;
					String lgOutDir_ijh = lgOutDir_i+"/"+j+"/h";
					Utils.makeOutputDir(new File(lgOutDir_ij));
					Utils.makeOutputDir(new File(lgOutDir_ijh));
				}
			}

			final List<String> scaff_list = new ArrayList<String>(nn_scaffs);
			Collections.sort(scaff_list, new Comparator<String>() {
				@Override
				public int compare(String scaff0, String scaff1) {
					// TODO Auto-generated method stub
					return StringUtils.countMatches(scaff1, ":")-
							StringUtils.countMatches(scaff0, ":");
				}
			});

			final String final_zip = in_zip;
			final int[][] task_table = new int[refine_round][lgN];
			final int[] task_progress = new int[lgN];

			for (int[] i : task_table)
				Arrays.fill(i, nr[2]);

			this.initial_thread_pool();

			for(int i=0; i<lgN; i++) { 
				String scaff_i = scaff_list.get(i);
				String lgOutDir_i = lgOutDir+"/"+i;
				String lgOutDir_ij = lgOutDir_i+"/0";
				String lgOutDir_ijh = lgOutDir_ij+"/h";
				for(int j=0; j<nr[2]; j++) {
					executor.submit(new Runnable(){
						private int i;
						private int j;
						@Override
						public void run() {
							// TODO Auto-generated method stub
							try {
								Haplotyper haplo = new Haplotyper();
								haplo.setParameters(new String[] {
										"-i", final_zip,
										"-o", lgOutDir_ijh,
										"-ex", expr_id,
										"-c", scaff_i,
										"-x", String.valueOf(max_iter),
										"-p", String.valueOf(ploidy),
										"-f", parents,
										"-s", nn_separation.get(scaff_i),
										"-r", nn_reverse.get(scaff_i),
										field
								});
								haplo.run();
								
								synchronized (task_table) {
									task_table[i][j]--;
									task_table.notify();
								}
							} catch (Exception e) {
								Thread t = Thread.currentThread();
								t.getUncaughtExceptionHandler().uncaughtException(t, e);
								e.printStackTrace();
								executor.shutdown();
								System.exit(1);
							}
						}

						public Runnable init(final int i, final int j) {
							// TODO Auto-generated method stub
							this.i = i;
							this.j = j;
							return (this);
						}
					}.init(0, i));
				}
			}

			while(true) {
				for(int i=0; i<lgN; i++) {

					if(task_progress[i]<refine_round &&
							task_table[task_progress[i]][i]==0) {

						String lgOutDir_i = lgOutDir+"/"+i;
						String lgOutDir_ij = lgOutDir_i+"/"+task_progress[i];
						String lgOutDir_ijh = lgOutDir_i+"/"+task_progress[i]+"/h";
						String lgOutPrex = lgOutDir_ij+"/out";

						
						SinglePointAnalysis sp1 = new SinglePointAnalysis();
						sp1.setParameters(new String[] {
								"-i", lgOutDir_ijh,
								"-o", lgOutPrex,
								"-wbp", String.valueOf(wbp),
								"-wnm", String.valueOf(wnm),
								"-ex", expr_id,
								"-nb", String.valueOf(nb),
								"-phi", String.valueOf(phi),
								"-nd", String.valueOf(drop),
								"-t", String.valueOf(THREADS)
						});
						sp1.run();
						
						TwoPointAnalysis tp1 = new TwoPointAnalysis();
						tp1.setParameters(new String[] {
								"-i", lgOutDir_ijh,
								"-o", lgOutPrex,
								"-ex", expr_id,
								"-nb", String.valueOf(nb),
								"-phi", String.valueOf(phi),
								"-nd", String.valueOf(drop),
								"-t", String.valueOf(THREADS)
						});
						tp1.run();
						
						//#### STEP 06 genetic mapping
						MappingAnalysis ma1 = new MappingAnalysis();
						ma1.setParameters(new String[] {
								"-i", lgOutPrex+".txt",
								"-m", lgOutPrex+".map",
								"-l", String.valueOf(lod_thresh),
								"-r", String.valueOf(lg_rf),
								"-1",
								"-rlib", RLibPath, 
								"-t", String.valueOf(THREADS),
								"-o", lgOutPrex
						});
						ma1.run();
						
						task_progress[i]++;

						if(task_progress[i]<refine_round) {

							final String lgOutDir_ijh1 = lgOutDir_i+"/"+task_progress[i]+"/h";

							Set<String> nn_scaffs_i = new HashSet<String>();
							Map<String, String> nn_separation_i = new HashMap<String, String>();
							Map<String, String> nn_reverse_i = new HashMap<String, String>();
							this.readSS(lgOutPrex+".par", nn_scaffs_i, nn_separation_i, nn_reverse_i);
							String scaff_i = nn_scaffs_i.iterator().next();

							for(int j=0; j<nr[2]; j++) {
								executor.submit(new Runnable(){
									private int i;
									private int j;
									@Override
									public void run() {
										// TODO Auto-generated method stub
										try {
											Haplotyper haplo = new Haplotyper();
											haplo.setParameters(new String[] {
													"-i", final_zip,
													"-o", lgOutDir_ijh1,
													"-ex", expr_id,
													"-c", scaff_i,
													"-x", String.valueOf(max_iter),
													"-p", String.valueOf(ploidy),
													"-f", parents,
													"-s", nn_separation_i.get(scaff_i),
													"-r", nn_reverse_i.get(scaff_i),
													field
											});
											haplo.run();
											
											synchronized (task_table) {
												task_table[i][j]--;
												task_table.notify();
											}
										} catch (Exception e) {
											Thread t = Thread.currentThread();
											t.getUncaughtExceptionHandler().uncaughtException(t, e);
											e.printStackTrace();
											executor.shutdown();
											System.exit(1);
										}
									}

									public Runnable init(final int i, final int j) {
										// TODO Auto-generated method stub
										this.i = i;
										this.j = j;
										return (this);
									}
								}.init(task_progress[i], i));
							}
						}
					}
				}
				boolean done = true;
				for(int i=0; i<lgN; i++)
					if(task_progress[i]<refine_round)
						done = false;
				if(done) break;
			}
			this.waitFor();
			
			// merge and clean file
			try {
				for(int j=0; j<refine_round; j++) {
					String outDir = lgOutDir+"/round_"+(j+1);
					Utils.makeOutputDir(new File(outDir));
					outPref = out_prefix+"/out"+(5+j);
					BufferedWriter bw_par = Utils.getBufferedWriter(outPref+".par");
					BufferedWriter bw_mct = Utils.getBufferedWriter(outPref+".mct");
					String line;
					for(int i=0; i<lgN; i++) {
						String lgOutDir_i = lgOutDir+"/"+i;
						String lgOutDir_ij = lgOutDir_i+"/"+j;
						String lgOutDir_ijh = lgOutDir_i+"/"+j+"/h";
						String lgOutPrex = lgOutDir_ij+"/out";

						BufferedReader br_par = Utils.getBufferedReader(lgOutPrex+".par");
						while((line=br_par.readLine())!=null) {
							bw_par.write(line);
							bw_par.write("\n");
						}
						br_par.close();

						bw_mct.write("group\tLG"+StringUtils.leftPad(""+(i+1), 2, '0')+"\n");
						BufferedReader br_mct = Utils.getBufferedReader(lgOutPrex+".mct");
						br_mct.readLine();
						while((line=br_mct.readLine())!=null) {
							bw_mct.write(line);
							bw_mct.write("\n");
						}
						br_mct.close();

						File[] fs = new File(lgOutDir_ijh).listFiles((File f) 
								-> f.getName().endsWith(".zip"));
						for(File f : fs) {
							Files.move(Paths.get(lgOutDir_ijh+"/"+f.getName()), 
									Paths.get(outDir+"/"+f.getName()),
									StandardCopyOption.REPLACE_EXISTING,
									StandardCopyOption.ATOMIC_MOVE);
						}
					}
					bw_par.close();
					bw_mct.close();
				}

				for(int i=0; i<lgN; i++) {
					Files.walk(Paths.get(lgOutDir+"/"+i))
					.map(Path::toFile)
					.sorted((o1, o2) -> -o1.compareTo(o2))
					.forEach(File::delete);
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			// calculate recombination frequencies for the last refinement round
			sp.setParameters(new String[] {
					"-i", lgOutDir+"/round_"+(refine_round+1),
					"-o", out_prefix+"/out"+(5+refine_round),
					"-wbp", String.valueOf(wbp),
					"-wnm", String.valueOf(wnm),
					"-ex", expr_id,
					"-nb", String.valueOf(nb),
					"-phi", String.valueOf(phi),
					"-nd", String.valueOf(drop),
					"-t", String.valueOf(THREADS)
			});
			sp.run();
			
			tp.setParameters(new String[] {
					"-i", lgOutDir+"/round_"+(refine_round+1),
					"-o", out_prefix+"/out"+(5+refine_round),
					"-ex", expr_id,
					"-nb", String.valueOf(nb),
					"-phi", String.valueOf(phi),
					"-nd", String.valueOf(drop),
					"-t", String.valueOf(THREADS)
			});
			tp.run();
		}
		
		//#### STEP 8 prepare results
		myLogger.info("STEP 8 prepare results");
		
		try {
			int h = 0;
			File[] fs = new File(out_prefix).listFiles((File f) 
					-> f.getName().endsWith(".mct"));
			for(File f : fs) {
				int h1 = Integer.parseInt(f.getName().
						replaceAll("^out", "").replaceAll(".mct$", ""));
				if(h1>h) h = h1;
			}
			Files.createSymbolicLink(Paths.get(out_prefix, "final.err"), Paths.get("out2.err"));
			Files.createSymbolicLink(Paths.get(out_prefix, "final.txt"), Paths.get("out"+h+".txt"));
			Files.createSymbolicLink(Paths.get(out_prefix, "final.map"), Paths.get("out"+h+".map"));
			Files.createSymbolicLink(Paths.get(out_prefix, "final.par"), Paths.get("out"+h+".par"));
			Files.createSymbolicLink(Paths.get(out_prefix, "final.mct"), Paths.get("out"+h+".mct"));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		//#### STEP 9 construct pseudomolecules
		myLogger.info("STEP 9 construct pseudomolecules");
		
		if(contig_file==null) myLogger.info("No assembly file provided, "
				+ "pseudomolecule construction module skipped.");
		
		Pseudomolecule pseudom = new Pseudomolecule();
		pseudom.setParameters(new String[] {
				"-i", out_prefix+"/final.mct",
				"-a", contig_file,
				"-e", out_prefix+"/final.err",
				"-n", String.valueOf(n_gap),
				"-o", out_prefix+"/final"
		});
		pseudom.run();
	}

	private void pendVCFContent(String sourceVcf, String targetVcf) {
		// TODO Auto-generated method stub
		try {
			BufferedWriter bw = Utils.getBufferedWriter(targetVcf, true);
			BufferedReader br = Utils.getBufferedReader(sourceVcf);
			String line;
			while((line=br.readLine())!=null) {
				if(!line.startsWith("#")) {
					bw.write(line);
					bw.write("\n");
				}
			}			
			br.close();
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void makeSelectedLG(String[] selected, String out) {
		// TODO Auto-generated method stub
		try {
			BufferedWriter bw_log = Utils.getBufferedWriter(out+"genetic_linkage_map.log");
			BufferedReader[] br_logs = new BufferedReader[selected.length];
			for(int i=0; i<br_logs.length; i++) 
				br_logs[i] = Utils.getBufferedReader(selected[i]+"1.log");
			String line = null;
			for(int i=0; i<br_logs.length; i++) 
				while( (line=br_logs[i].readLine())!=null &&
				!line.startsWith("$") ) {}
			while( line!=null ) {
				if(line.startsWith("$")) {
					bw_log.write(line+"\n");
					for(int i=0; i<br_logs.length; i++) 
						while( (line=br_logs[i].readLine())!=null &&
						!line.startsWith("$") )
							if(line.length()>0) bw_log.write(line+"\n");
					bw_log.write("\n");
				}
			}
			for(int i=0; i<br_logs.length; i++) br_logs[i].close();
			bw_log.close();
			
			BufferedWriter bw_mct = Utils.getBufferedWriter(out+"genetic_linkage_map.mct");
			for(int i=0; i<selected.length; i++) { 
				BufferedReader br_mct = Utils.getBufferedReader(selected[i]+"1.mct");
				br_mct.readLine();
				bw_mct.write("group\tLG"+StringUtils.leftPad(""+i, 2, '0')+"\n");
				while( (line=br_mct.readLine())!=null )
					if(line.length()>0) bw_mct.write(line+"\n");
				bw_mct.write("\n");
				br_mct.close();
			}
			bw_mct.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private double lgCM(String in_file) {
		// TODO Auto-generated method stub
		try {
			BufferedReader br = Utils.getBufferedReader(in_file);
			String line;
			while( (line=br.readLine())!=null &&
					!line.startsWith("$cm") ){}
			line = br.readLine();
			String[] s = line.split(",");
			br.close();
			double cm = 0.0;
			for(int i=0; i<s.length; i++)
				cm += Double.parseDouble(s[i]);
			return cm;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return -1;
	}

	private void readSS(String ss_in, 
			Set<String> nn_scaffs,
			Map<String, String> nn_separation,
			Map<String, String> nn_reverse) {
		// TODO Auto-generated method stub
		nn_scaffs.clear();
		nn_separation.clear();
		nn_reverse.clear();
		try {
			BufferedReader nn_br = Utils.getBufferedReader(ss_in);
			String line;
			String[] s;
			while( (line=nn_br.readLine())!=null ) {
				if(!line.startsWith("-c")) continue;
				s = line.split("\\s+");
				nn_scaffs.add(s[1]);
				if(s.length>2) {
					nn_separation.put(s[1], s[3]);
					nn_reverse.put(s[1], s[5]);
				} else {
					nn_separation.put(s[1], "0");
					nn_reverse.put(s[1], "false");
				}
			}
			nn_br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void move(final Set<String> scaffs, 
			final String out,
			final String expr_id) {
		// TODO Auto-generated method stub
		for(final String scaff : scaffs) {
			File[] files = new File(out).listFiles(
					new FilenameFilter() {
						@Override
						public boolean accept(File dir, String name) {
							return name.matches("^"+expr_id+"\\."+scaff+"\\..*");    
						}
					});
			for(File f : files)
				try {
					Files.move(f.toPath(),
							Paths.get(out),
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
			final int repeat,
			final String out) {
		// TODO Auto-generated method stub
		this.initial_thread_pool();
		for(final String scaff : scaffs.keySet()) {
			
			if(scaffs.get(scaff)<this.min_snpc) continue;
			
			for(int i=0; i<repeat; i++) {
				executor.submit(new Runnable(){
					@Override
					public void run() {
						// TODO Auto-generated method stub
						try {
							Haplotyper haplo = new Haplotyper();
							haplo.setParameters(new String[] {
									"-i", in_zip,
									"-o", out,
									"-ex", expr_id,
									"-c", scaff,
									"-x", String.valueOf(max_iter),
									"-p", String.valueOf(ploidy),
									"-f", parents,
									field
							});
							haplo.run();
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
	
	private void runHaplotyper(final Set<String> scaffs,
			final Map<String, String> separation,
			final Map<String, String> reverse,
			final String expr_id,
			final String in_zip,
			final int repeat,
			final String out) {
		// TODO Auto-generated method stub
		this.initial_thread_pool();
		for(final String scaff : scaffs) {
			for(int i=0; i<repeat; i++) {
				executor.submit(new Runnable(){

					@Override
					public void run() {
						// TODO Auto-generated method stub
						try {
							Haplotyper haplo = new Haplotyper();
							haplo.setParameters(new String[] {
									"-i", in_zip,
									"-o", out,
									"-ex", expr_id,
									"-c", scaff,
									"-x", String.valueOf(max_iter),
									"-p", String.valueOf(ploidy),
									"-f", parents,
									"-s", separation.get(scaff),
									"-r", reverse.get(scaff),
									field
							});
							haplo.run();
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
