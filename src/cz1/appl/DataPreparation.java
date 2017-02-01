package cz1.appl;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.BufferedReader;
import java.lang.StringBuilder;
import java.util.Random;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.RejectedExecutionHandler;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPInputStream;
import java.io.OutputStreamWriter;
import java.util.zip.GZIPOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.FileInputStream;
import java.io.Reader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.HashMap;
import java.util.Calendar;
import java.util.List;
import java.util.Map;
import java.text.SimpleDateFormat;

import PedigreeSim.PopulationData;
import PedigreeSim.PedigreeSimulate;
import PedigreeSim.Chromosome;
import PedigreeSim.Locus;
import PedigreeSim.Individual;

import org.apache.commons.io.FileUtils;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;

import cz1.util.ArgsEngine;
import cz1.util.Utils;

//TODO
// multi-threading has not been tested yet

public class DataPreparation {
	private final static Logger myLogger = 
			Logger.getLogger(DataPreparation.class);
	static {
		BasicConfigurator.configure();
	}
	private ArgsEngine myArgsEngine = null;
	public final static String NLS = System.getProperty("line.separator");
	public final static String SEP = System.getProperty("file.separator");
	public final static long HGB = 3079843747L; //human genome base pairs
	//mean and sd of base pairs across human chromosomes
	public final static long MEAN = 128326823L; 
	public final static long SD = 58070837L;
	//estimated variations in human genome
	public final static long VARIATIONS = 55757749L;
	//decide the percentage of VARIATIONS (SNPs) that the two parents carries
	private final static double[] POISSON_MIXTURE_LAMBDA = 
			new double[] {100, 200, 500, 1000, 2000, 5000, 10000};
	private final static double[] POISSON_MIXTURE_PROBS = 
			new double[] {.05, .1, .2, .3, .2, .1, .05};
	private final static double PROBS_PRECISION = 1e-4;
	private final static Map<Integer, Double> POISSON_MEAN_MAP = 
			new HashMap<Integer, Double>();
	static {
		int p = 0;
		int n = (int) (1/PROBS_PRECISION);
		for(int i=0; i<POISSON_MIXTURE_PROBS.length; i++) {
			int k = (int) (n*POISSON_MIXTURE_PROBS[i]);
			for(int j=0; j<k; j++) POISSON_MEAN_MAP.put(p+j, 
					POISSON_MIXTURE_LAMBDA[i]);
			p += k;
		}
	}
	private final static double PERCENT_SNPS_TO_ALL = 0.1;
	private int PLOIDY = 2;
	private int MAX_ALLELES = 2;
	private int CHROM_NUMBER = 1;
	private int PROGENY_NUMBER = 94;
	private String REFERENCE = null;
	private int THREADS = 1;
	private final int CHUNK_SIZE = 100; //when writing fasta file
	private final double CHROM_BASE_LENGTH = 100.0;
	private final double CHROM_BASE_CENTROMERE_POS = 40.0;
	private long RANDOM_SEED = 112019890314L;
	//private final long RANDOM_SEED = 314198901120L;
	private final Random random= new Random(RANDOM_SEED);
	private final double REF_ALLELE_FRQ = 0.75;
	private final char[] BASE = new char[] {'A','C','G','T'};
	private final int BASE_NUM = BASE.length;
	private final double[][] SNP_TRANS_MAT = new double[][]{
			/***To*/ /******      A      C       G       T   */
			/***From*/     /***A*/ { 0.0000, 1862.0, 6117.0, 1363.0 },
			/***C*/ { 1919.0, 0.0000, 2423.0, 5623.0 },
			/***G*/ { 6726.0, 2985.0, 0.0000, 2030.0 },
			/***T*/ { 1938.0, 7090.0, 7112.0, 0.0000 }
	};
	private Genome genome;
	private Pedigree pedigree;
	private Chromosome0 chromosome;
	private ArrayList<SNP> parentSNPs;
	private String scenario = "Sc1";
	private String filePath = ".";
	private PopulationData popData;
	private HashMap<String, Integer> SNPs2PosMap;
	private HashMap<String, String[]> Chromsome2SNPsMap;
	private static ExecutorService executor;
	private BlockingQueue<Runnable> tasks = null;

	private void initial_thread_pool() {
		tasks = new ArrayBlockingQueue<Runnable>(THREADS);
		executor = new ThreadPoolExecutor(THREADS, 
				THREADS, 
				1, 
				TimeUnit.SECONDS, 
				tasks, 
				new RejectedExecutionHandler(){
			@Override
			public void rejectedExecution(Runnable task,
					ThreadPoolExecutor arg1) {
				// TODO Auto-generated method stub
				try {
					tasks.put(task);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		});
	}

	public DataPreparation(String fastaFilePath, int chromNumber, 
			int offspring, int ploidy, String scenario, String filePath) 
					throws Exception{
		this.scenario = scenario;
		this.filePath = filePath;
		this.PLOIDY = ploidy;
		this.genome = generatePrimitiveCHROM(fastaFilePath,chromNumber);
		this.pedigree = generatePedigree(offspring);
		this.chromosome = generateChromosome();
		this.parentSNPs = generateParentSNPs();
		this.SNPs2PosMap = buildSNPs2PosMap();
		this.Chromsome2SNPsMap = buildChromsome2SNPsMap();
		writeFiles();
		this.popData = PedigreeSimulate.simulate(filePath+SEP+scenario+".par");
	}

	public DataPreparation(String fastaFilePath, int chromNumber,
			int offspring, int ploidy, String scenario) throws Exception{
		this.scenario = scenario;
		this.PLOIDY = ploidy;
		this.genome = generatePrimitiveCHROM(fastaFilePath,chromNumber);
		this.pedigree = generatePedigree(offspring);
		this.chromosome = generateChromosome();
		this.parentSNPs = generateParentSNPs();
		this.SNPs2PosMap = buildSNPs2PosMap();
		this.Chromsome2SNPsMap = buildChromsome2SNPsMap();
		createFiles();
		writeFiles();
		this.popData = PedigreeSimulate.simulate(filePath+SEP+scenario+".par");
	}

	public DataPreparation(String fastaFilePath, int chromNumber,
			int offspring, int ploidy) throws Exception{
		this.PLOIDY = ploidy;
		this.genome = generatePrimitiveCHROM(fastaFilePath,chromNumber);
		this.pedigree = generatePedigree(offspring);
		this.chromosome = generateChromosome();
		this.parentSNPs = generateParentSNPs();
		this.SNPs2PosMap = buildSNPs2PosMap();
		this.Chromsome2SNPsMap = buildChromsome2SNPsMap();
		createFiles();
		writeFiles();
		this.popData = PedigreeSimulate.simulate(filePath+SEP+scenario+".par");
	}

	public DataPreparation() throws Exception {
		this.genome = generatePrimitiveCHROM(this.REFERENCE,this.CHROM_NUMBER);
		this.pedigree = generatePedigree(this.PROGENY_NUMBER);
		this.chromosome = generateChromosome();
		this.parentSNPs = generateParentSNPs();
		this.SNPs2PosMap = buildSNPs2PosMap();
		this.Chromsome2SNPsMap = buildChromsome2SNPsMap();
		if(this.filePath.equals(".")) createFiles();
		this.writeFiles();
		this.popData = PedigreeSimulate.simulate(filePath+SEP+scenario+".par");
	}

	private void printUsage() {
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -r  Reference (fasta file format). \n"
						+ " -c  Number of chromosomes (default 1). \n"
						+ " -n  Number of offspring (F1, default 94). \n"
						+ " -p  Copy of chromosomes (ploidy, default 2). \n"
						+ " -t  Threads (default 1).\n"
						+ " -s	Scenario (default Sc1).\n" 
						+ " -o  Output directory (defult current directory).\n\n");
	}

	public void setParameters(String[] args) {
		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}

		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-r", "--reference", true);
			myArgsEngine.add("-c", "--chromosomes", true);
			myArgsEngine.add("-n", "--offspring", true);
			myArgsEngine.add("-p", "--ploidy", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-s", "--scenario", true);
			myArgsEngine.add("-o", "--output-file", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-r")) {
			this.REFERENCE = myArgsEngine.getString("-r");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the reference FASTA file.");
		}

		if (myArgsEngine.getBoolean("-c")) {
			this.CHROM_NUMBER = Integer.parseInt(myArgsEngine.getString("-c"));
		}

		if (myArgsEngine.getBoolean("-n")) {
			this.PROGENY_NUMBER = Integer.parseInt(myArgsEngine.getString("-n"));
		}

		if (myArgsEngine.getBoolean("-p")) {
			this.PLOIDY = Integer.parseInt(myArgsEngine.getString("-p"));
		}

		if (myArgsEngine.getBoolean("-t")) {
			this.THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}

		if (myArgsEngine.getBoolean("-s")) {
			this.scenario = myArgsEngine.getString("-s");
		}

		if (myArgsEngine.getBoolean("-o")) {
			this.filePath = myArgsEngine.getString("-o");
		}

	}

	public static void main(String[] args) throws Exception {
		DataPreparation dp = new DataPreparation();
		dp.writeGenomeFile();
		dp.writeGenomeFileAll();
	}

	private void writeFiles() {
		writePedigreeFile();
		writeChromosomeFile();
		writeGeneFile();
		writeMapFile();
		writeParameterFile();
		return;
	}

	private void createFiles() throws SecurityException, IOException {
		File theDir = new File(filePath);
		if (!theDir.exists() || theDir.exists() && !theDir.isDirectory())
			theDir.mkdir();
		theDir = new File(filePath+SEP+scenario);
		if (theDir.exists() && theDir.isDirectory()) 
			FileUtils.deleteDirectory(theDir);
		theDir.mkdir();
		return;
	}

	private HashMap<String, Integer> buildSNPs2PosMap() {
		HashMap<String, Integer> map = new HashMap<String, Integer>();
		SNP snp;
		for(int i=0; i<parentSNPs.size(); i++) {
			snp = parentSNPs.get(i);
			map.put(snp.getName(),snp.getPosBasePairs());
			//debug
			//println(">>>"+snp.getName()+"..."+snp.getPosBasePairs()+"..."+map.get(snp.getName()));
		}
		return map;
	}

	private HashMap<String, String[]> buildChromsome2SNPsMap() {
		HashMap<String, String[]> map = new HashMap<String, String[]>();
		String names[] = chromosome.getName();
		SNP snp;
		ArrayList<String> snpslist;
		String[] snps;
		for(int i=0; i<names.length; i++) {
			snpslist = new ArrayList<String>();
			for(int j=0; j<parentSNPs.size(); j++) {
				snp = parentSNPs.get(j);
				if(snp.getChromosome().equals(names[i]))
					snpslist.add(snp.getName());
			}
			snps = new String[snpslist.size()];
			snpslist.toArray(snps);
			map.put(names[i], snps);
		}
		return map;
	}

	/***
	 * generate the primitive chromsomes from a reference genome;
	 * the reference genome is first merged together, filtering out
	 * the fasta file identifiers and missing bases 'N';
	 * then the concatenated genome (of length L) is splitted into N 
	 * (number of chromosomes) segments, following a norm distribution 
	 * similar to human genome;
	 * bases from the same contig are not seperated.
	 */
	public Genome generatePrimitiveCHROM(String fastaFilePath, int N) {
		int[] cutsite = null;
		String ref = null;
		try{
			BufferedReader br = getBufferedReader(fastaFilePath);
			StringBuilder sb = new StringBuilder();
			String line;
			ArrayList<Integer> count = new ArrayList<Integer>();
			String contig;
			while( (line=br.readLine())!=null) {
				if(!line.startsWith(">")) {
					contig = line.replaceAll("N","");
					sb.append(contig);
					count.add(contig.length());
				}
			}
			br.close();
			ref = sb.toString();
			int L = ref.length();
			double mean = (double)L/N;
			double sd = (double)L/HGB*SD;
			long[] size = new long[N];
			int sum=0;
			for(int i=0; i<(N-1); i++) {
				size[i] = Math.round(random.nextGaussian()*sd+mean);
				sum += size[i];
			}
			size[N-1] = ref.length()-sum;

			//debug
			//println("###### "+ref.length());
			//println(""+mean);
			//println(""+sd);
			//for(int i=0; i<N; i++) System.out.println(size[i]);

			cutsite = findCutSite(count,size,ref.length());
			//debug
			//print(cutsite);
		} catch (IOException e) {
			e.printStackTrace();
		}
		return new Genome(ref,cutsite);
	}

	public void writeGenomeFile() {
		println("Writing genome file...");
		int[] cutsite = genome.getCutsite();
		int N = cutsite.length-1;
		String ref = genome.getReference();
		try(BufferedWriter bw = getBufferedWriter(filePath+SEP+scenario+
				"reference.fasta")){
			String chrom;
			for(int i=0; i<N; i++) {
				chrom = ref.substring(cutsite[i],cutsite[i+1]);
				bw.write(">CHROM"+(i+1)+NLS);
				int j=0;
				while(j+CHUNK_SIZE <= chrom.length()) {
					bw.write(chrom.substring(j,j+CHUNK_SIZE)+NLS);
					j += CHUNK_SIZE;
				}
				if(j<chrom.length()) bw.write(chrom.substring(j)+NLS);
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void writeGenomeFileAll() {
		this.initial_thread_pool();

		final int M = popData.indivCount();
		final int N = popData.chromCount();
		final String[] chromNames = chromosome.getName();
		final int[] cutsite = genome.getCutsite();
		final String ref = genome.getReference();

		for (int i=0; i<M; i++) {

			executor.submit(new Runnable() {
				private int i;

				@Override
				public void run() {
					String[][] indall;
					//StringBuilder chr;
					ArrayList<Locus> loci;
					int pos;
					String str2write;

					Individual indiv =  popData.getIndiv(i);
					try(BufferedWriter bw = getBufferedWriter(filePath+SEP+scenario+
							indiv.getIndivName()+".fasta")) {
						for (int c=0; c<N; c++) {
							Chromosome chrom = popData.getChrom(c);
							loci = chrom.getLocus();
							indall = new String[PLOIDY][loci.size()];
							for (int l=0; l<loci.size(); l++) {
								String[] inds = indiv.getLocusAllele(c,l);
								for(int h=0; h<PLOIDY; h++) {
									indall[h][l] = inds[h];
								}
							}

							//debug
							/*
                    		println("building snp map done...");
                    		println(""+indall.length+"..."+indall[0].length);
                    		println("loci size - "+loci.size());
                    		println("locus name - "+loci.get(10).getLocusName());
                    		println("position - "+SNPs2PosMap.get(loci.get(10).getLocusName()));
							 */

							for(int h=0; h<PLOIDY; h++) {
								//debug
								//println(getSystemTime()+" >>> StringBulder starting...");
								/*** String builder is slow
                        chr = new StringBuilder(ref.substring(cutsite[c],cutsite[c+1]));
                        for(int l=0; l<loci.size(); l++) {
                            pos = SNPs2PosMap.get(loci.get(l).getLocusName());
                            chr.replace(pos,pos+indall[h][l].length(),indall[h][l]);            
                        }
								 **/
								//debug
								//println(getSystemTime()+" >>> Done. BufferedWriter Start writing...");

								bw.write(">"+chromNames[c]+"|ploidy "+PLOIDY+"|genome "+(h+1)+" of "+PLOIDY+NLS);
								int j0 = cutsite[c], j1 = cutsite[c+1];
								for(int l=0; l<loci.size(); l++) {
									pos = SNPs2PosMap.get(loci.get(l).getLocusName())+cutsite[c];
									str2write = indall[h][l];
									while(j0+CHUNK_SIZE <= pos) {
										bw.write(ref.substring(j0,j0+CHUNK_SIZE)+NLS);
										j0 += CHUNK_SIZE;
									}
									//debug
									//println(">>>"+ref.length()+"..."+j0+"..."+pos+"..."+str2write);
									bw.write(ref.substring(j0,pos-str2write.length()+1));
									bw.write(str2write+NLS);
									j0 = pos+str2write.length();
								}
								if(j0 < j1) {
									while(j0+CHUNK_SIZE <= j1) {
										bw.write(ref.substring(j0,j0+CHUNK_SIZE)+NLS);
										j0 += CHUNK_SIZE;
									}
									if(j0 < j1) bw.write(ref.substring(j0,j1)+NLS);
								}

								//debug
								myLogger.info(getSystemTime()+" >>> done.");
							}
						}
						bw.close();
					} catch (Exception e) {
						e.printStackTrace();
						System.exit(1);
					}
				} 

				public Runnable init(int i) {
					this.i = i;
					return(this);
				}
			}.init(i) );
		}

		executor.shutdown();
		try {
			executor.awaitTermination(365, TimeUnit.DAYS);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public void writeParameterFile() {
		println("Writing parameter file...");
		try(BufferedWriter bw = getBufferedWriter(filePath+SEP+scenario+".par")) {
			bw.write("; This is a parameter file for PedigreeSim" + NLS +
					"; to run this example:" + NLS +
					"; - open a console window and go to the folder containing these example files" + NLS +
					"; - give the following command (with for <path> the path where PedigreeSim.jar" + NLS +
					";   and folder lib containing jsci-core.jar are located):" + NLS +
					";   java -jar <path>PedigreeSim.jar " + scenario + ".par" + NLS + NLS);
			bw.write("PLOIDY = " + PLOIDY + NLS +
					"MAPFUNCTION = HALDANE" +NLS +
					"MISSING = NA" + NLS +
					"CHROMFILE = " + filePath + SEP + scenario + ".chrom" + NLS +
					"PEDFILE = " + filePath + SEP + scenario + ".ped" + NLS +
					"MAPFILE = " + filePath + SEP + scenario + ".map" + NLS +
					"FOUNDERFILE = " + filePath + SEP + scenario + ".gen" + NLS +
					"OUTPUT = " + filePath + SEP + scenario + "_out" + NLS + NLS);
			bw.write("; the following parameters are here set to their default values," + NLS +
					"; these lines may therefore be omitted:"  + NLS + NLS +
					"ALLOWNOCHIASMATA = 1" + NLS +
					"NATURALPAIRING = 1 ; Note that this overrules the \"quadrivalent\" column" + NLS +
					"                   ; in the CHROMFILE" + NLS +
					"PARALLELQUADRIVALENTS = 0.0" + NLS +
					"PAIREDCENTROMERES = 0.0");
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void writeGeneFile() {
		println("Writing gene file...");
		try(BufferedWriter bw = getBufferedWriter(filePath+SEP+scenario+".gen")) {
			bw.write("marker\t");
			for(int i=1; i<=PLOIDY; i++) bw.write("P1_"+i+"\t");
			for(int i=1; i<=PLOIDY; i++) bw.write("P2_"+i+"\t");
			bw.write(NLS);
			SNP snp;
			char[] phase;
			for(int i=0; i<parentSNPs.size(); ) {
				bw.write(parentSNPs.get(i).getName()+"\t");
				for(int p=0; p<2; p++) {
					snp = parentSNPs.get(i++);
					phase = snp.getPhase();
					for(int j=0; j<PLOIDY; j++) bw.write(phase[j]+"\t");
				}
				bw.write(NLS);
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public void writeMapFile() {
		println("Writing map file...");
		try(BufferedWriter bw = getBufferedWriter(filePath+SEP+scenario+".map")) {
			bw.write("marker\tchromosome\tposition"+NLS);
			for(int i=0; i<parentSNPs.size(); ) {
				bw.write(parentSNPs.get(i).getName()+"\t"+
						parentSNPs.get(i).getChromosome()+"\t"+
						parentSNPs.get(i).getPosCentiMorgan()+NLS);
				i += 2;
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public ArrayList<SNP> generateParentSNPs() {
		int[] cutsite = genome.getCutsite();
		String ref = genome.getReference();
		int L = ref.length();
		int N = chromosome.getNumber();
		String[] chrom = chromosome.getName();
		int[] LBasePairs = chromosome.getLBasePairs();
		double[] LCentiMorgan = chromosome.getLCentiMorgan();

		int amount = (int)Math.round((double)VARIATIONS/HGB*L*PERCENT_SNPS_TO_ALL);

		//debug
		//println(""+VARIATIONS+"..."+HGB+"..."+L+"..."+amount);

		char[][] alleles = getRandomAlleles(amount);
		char[][][] phase = getRandomPhase(alleles); //for 2 parents

		ArrayList<SNP> snps = new ArrayList<SNP>(); //even for parent 1 and odd for parent 2
		int next, t, n=0;
		//double lambda = L/amount;
		double lambda;
		for(int i=0; i<N; i++) {
			double base = LCentiMorgan[i]/LBasePairs[i];
			t = 0;
			lambda = POISSON_MEAN_MAP.get(random.nextInt((int) (1/PROBS_PRECISION)));
			next = poisson(lambda);
			// System.out.println(lambda+"\t\t"+next);
			t += next;
			while(t<LBasePairs[i] && n<phase.length) {
				snps.add(new SNP(chrom[i]+"."+t, chrom[i], alleles[n], phase[n][0], t, t*base));
				snps.add(new SNP(chrom[i]+"."+t, chrom[i], alleles[n], phase[n][1], t, t*base));
				//debug
				//println(">>>"+t+"..."+snps.get(snps.size()-1).getPosBasePairs()+"..."+snps.get(snps.size()-1).getPosCentiMorgan());
				n++;
				lambda = POISSON_MEAN_MAP.get(random.nextInt((int) (1/PROBS_PRECISION)));
				next = poisson(lambda);
				// System.out.println(lambda+"\t\t"+next);
				t += next;
			}
		}
		return snps;
	}

	public int poisson(double lambda) {
		double L = -lambda;
		double p = 0;
		int k = 0;
		do {
			k++;
			p += Math.log(random.nextDouble());
		} while (p > L);
		return k - 1;
	}

	/*** only for biallelic SNPs currently,
	 * return a N*M*P char matrix,
	 * where N is the number of SNPs,
	 * M is the number of parents (i.e. 2 here),
	 * and P is the ploidy.
	 */
	public char[][][] getRandomPhase(char[][] alleles) {
		char[][] phaseMerge = new char[alleles.length][PLOIDY*2];
		int c;
		double r;
		for(int i=0; i<phaseMerge.length; i++) {
			c = 0;
			for(int j=0; j<PLOIDY*2; j++) {
				r = random.nextDouble();
				if(r <= REF_ALLELE_FRQ) {
					phaseMerge[i][j] = alleles[i][0];
					c++;
				} else {
					phaseMerge[i][j] = alleles[i][1];
				}
			}

			// if both the parents are homozygous,
			//  change the allele of a randomly chosen position
			if(c==0 | c==PLOIDY*2) {
				char alt = c==0 ? alleles[i][0] : alleles[i][1];
				r = random.nextDouble();
				int k = 0;
				while(r>=(double)(k+1)/PLOIDY/2) {
					k++;
				}
				//debug
				//println(""+r+"..."+alleles[i][0]+alleles[i][1]+"..."+c+"..."+alt+"..."+k);
				phaseMerge[i][k] = alt;
			}
		}
		char[][][] phase = new char[phaseMerge.length][2][PLOIDY];
		for(int i=0; i<phase.length; i++) {
			for(int j=0; j<PLOIDY; j++) {
				phase[i][0][j] = phaseMerge[i][j];
				phase[i][1][j] = phaseMerge[i][j+PLOIDY];
			}
		}
		return phase;
	}

	public char[][] getRandomAlleles(int amount) {
		char[][] alleles = new char[amount][MAX_ALLELES];
		double sum = 0;
		double[] cusumfrom = new double[BASE_NUM],
				from = new double[BASE_NUM];
		double[][] cusumto = new double[BASE_NUM][BASE_NUM];
		for(int i=0; i<BASE_NUM; i++){
			for(int j=0; j<BASE_NUM; j++) from[i]+=SNP_TRANS_MAT[i][j];
			sum += from[i];
		}
		cusumfrom[0] = from[0]/sum;
		for(int i=1; i<BASE_NUM; i++){
			cusumfrom[i] = cusumfrom[i-1]+from[i]/sum;
		}
		for(int i=0; i<BASE_NUM; i++){
			cusumto[i][0] = SNP_TRANS_MAT[i][0]/from[i];
			for(int j=1; j<BASE_NUM; j++){
				cusumto[i][j] = SNP_TRANS_MAT[i][j]/from[i]+cusumto[i][j-1];
			}
		}

		//debug
		//print(from);
		//print(cusumfrom);
		//print(cusumto);

		int ref, alt;
		for(int i=0; i<amount; i++) {
			ref = getBaseFromCusumProbs(cusumfrom,random.nextDouble());
			alt = getBaseFromCusumProbs(cusumto[ref],random.nextDouble());
			alleles[i] = new char[] {BASE[ref],BASE[alt]};
		}
		return alleles;
	}

	public int getBaseFromCusumProbs(double[] cusum, double r) {
		int base;
		for(base=0; base<BASE_NUM; base++)
			if(cusum[base]>=r) break;
		return base;
	}

	public Chromosome0 generateChromosome() {
		int[] cutsite = genome.getCutsite();
		int N = cutsite.length-1;
		double mean = (double)cutsite[N]/N;
		String[] name = new String[N];
		int[] LBasePairs = new int [N];
		double[] LCentiMorgan = new double[N], centromere = new double[N], 
				prefPairing = new double[N], quadrivalent = new double[N];

		for(int i=0; i<N; i++) {
			name[i] = "CHROM"+(i+1);
			LBasePairs[i] = cutsite[i+1]-cutsite[i];
			LCentiMorgan[i] = LBasePairs[i]/mean*CHROM_BASE_LENGTH;
			centromere[i] = CHROM_BASE_CENTROMERE_POS/CHROM_BASE_LENGTH*LCentiMorgan[i];
			prefPairing[i] = (random.nextFloat()/2+0.5);
			quadrivalent[i] = 0.0;
		}

		return new Chromosome0(N,name,LBasePairs,LCentiMorgan,centromere,prefPairing,quadrivalent);
	}

	public void writeChromosomeFile() {
		println("Writing chromosome file...");
		int[] cutsite = genome.getCutsite();
		int N = chromosome.getNumber();
		String[] name = chromosome.getName();
		double[] LCentiMorgan = chromosome.getLCentiMorgan(), 
				centromere = chromosome.getCentromere(),
				prefPairing = chromosome.getPrefPairing(), 
				quadrivalent = chromosome.getQuadrivalent();

		try(BufferedWriter bw = getBufferedWriter(filePath+SEP+scenario+".chrom")) {
			bw.write("chromosome\tlength\tcentromere\tprefPairing\tquadrivalents"+NLS);
			for(int i=0; i<N; i++) bw.write(name[i]+"\t"+LCentiMorgan[i]+"\t"+centromere[i]+
					"\t"+prefPairing[i]+"\t"+quadrivalent[i]+NLS);
			bw.close();
		} catch (IOException e){
			e.printStackTrace();
		}
	}

	public Pedigree generatePedigree(int offspring) {
		return new Pedigree(2,offspring);
	}

	public void writePedigreeFile() {
		println("Writing pedigree file...");
		int offspring = pedigree.getOffspring();
		try(BufferedWriter bw = getBufferedWriter(filePath+SEP+scenario+".ped")){
			bw.write("Name\tParent1\tParent2"+NLS);
			bw.write("P1\tNA\tNA"+NLS);
			bw.write("P2\tNA\tNA"+NLS);
			for(int i=0; i<offspring; i++) bw.write("F"+(i+1)+"\tP1\tP2"+NLS);
			bw.close();
		} catch (IOException e){
			e.printStackTrace();
		}
	}

	public int[] findCutSite(ArrayList<Integer> count, long[] size, int length) {
		int N = size.length, L=0, l;
		int[] cutsite = new int[N+1];
		Iterator<Integer> ii = count.iterator();
		for(int i=0; i<N; i++) {
			l = 0;
			while(l<size[i] && ii.hasNext()) {
				l += ii.next();
			}
			cutsite[i+1] = L+l;
			L += l;
		}
		return cutsite;

	}

	public BufferedReader getBufferedReader(String path) throws IOException {
		BufferedReader br = null;

		if(path.endsWith(".gz")){
			InputStream inputFileStream = new FileInputStream(path);
			InputStream gzipStream = new GZIPInputStream(inputFileStream);
			Reader decoder  = new InputStreamReader(gzipStream);
			br = new BufferedReader(decoder, 65536);
		}else{
			br = new BufferedReader(new FileReader(path), 65536);
		}
		return br;
	}

	public BufferedWriter getBufferedWriter(String path) throws IOException {
		return new BufferedWriter(new FileWriter(new File(path)));
	}

	public BufferedWriter getGZIPBufferedWriter(String path) throws IOException {
		return new BufferedWriter(new OutputStreamWriter(new 
				GZIPOutputStream(new FileOutputStream(path))));
	}

	public static String getSystemTime(){
		return new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").
				format(Calendar.getInstance().getTime());
	}

	public void println() {
		System.out.println();
	}

	public void print(char character) {
		System.out.print(character);
	}

	public void println(char character) {
		print(character+"\n");
	}

	public void print(String message) {
		System.out.print(message);
	}


	public void println(String message) {
		print(message+"\n");
	}

	public void print(ArrayList<Integer> array) {
		for(int i=0; i<array.size(); i++)
			print(array.get(i)+"\t");
		println();
	}

	public void print(int[][] matrix) {
		for(int i=0; i<matrix.length; i++)
			print(matrix[i]);
	}


	public void print(int[] array) {
		for(int i=0; i<array.length; i++)
			print(array[i]+"\t");
		println();
	}

	public void print(double[] array) {
		for(int i=0; i<array.length; i++)
			print(array[i]+"\t");
		println();
	}

	public void print(double[][] matrix) {
		for(int i=0; i<matrix.length; i++)
			print(matrix[i]);
	}

	public void print(long[] array) {
		for(int i=0; i<array.length; i++)
			print(array[i]+"\t");
		println();
	}
}

class Genome {
	private String reference;
	private int[] cutsite;

	public Genome(String reference, int[] cutsite) {
		this.reference = reference;
		this.cutsite = cutsite;
	}

	public String getReference() {
		return reference;
	}

	public int[] getCutsite() {
		return cutsite;
	}
}

class Pedigree {
	private int parents;
	private int offspring;

	public Pedigree(int parents, int offspring) {
		this.parents = parents;
		this.offspring = offspring;
	}

	public int getParents() {
		return parents;
	}
	public int getOffspring() {
		return offspring;
	}

}

class Chromosome0 {
	private int N;
	private String[] name;
	private int[] LBasePairs;
	private double[] LCentiMorgan;
	private double[] centromere;
	private double[] prefPairing;
	private double[] quadrivalent;

	public Chromosome0(int N, String[] name, int[] LBasePairs,
			double[] LCentiMorgan, double[] centromere,
			double[] prefPairing, double[] quadrivalent) {
		this.N = N;
		this.name = name;
		this.LBasePairs = LBasePairs;
		this.LCentiMorgan = LCentiMorgan;
		this.centromere = centromere;
		this.prefPairing = prefPairing;
		this.quadrivalent = quadrivalent;
	}

	public int getNumber() {
		return N;
	}

	public String[] getName() {
		return name;
	}

	public int[] getLBasePairs() {
		return LBasePairs;
	}

	public double[] getLCentiMorgan() {
		return LCentiMorgan;
	}

	public double[] getCentromere() {
		return centromere;
	}

	public double[] getPrefPairing() {
		return prefPairing;
	}

	public double[] getQuadrivalent() {
		return quadrivalent;
	}
}

class SNP {
	private String name;
	private String chromosome;
	private char[] alleles;
	private char[] phase;
	private int posBasePairs;
	private double posCentiMorgan;

	public SNP(String name, String chromosome, char[] alleles, char[] phase, 
			int posBasePairs, double posCentiMorgan) {
		this.name = name;
		this.chromosome = chromosome;
		this.alleles = alleles;
		this.phase = phase;
		this.posCentiMorgan = posCentiMorgan;
		this.posBasePairs = posBasePairs;
	}

	public String getName() {
		return name;
	}

	public String getChromosome() {
		return chromosome;
	}

	public char[] getAlleles() {
		return alleles;
	}

	public char[] getPhase() {
		return phase;
	}

	public int getPosBasePairs() {
		return posBasePairs;
	}

	public double getPosCentiMorgan() {
		return posCentiMorgan;
	}
}












