package cz1.simulation.model;

import cz1.util.ArgsEngine;

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
import java.util.regex.Matcher;
import java.util.regex.Pattern;
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

//TODO
// multi-threading has not been tested yet

public class Population {
	private final static Logger myLogger = 
			Logger.getLogger(Population.class);
	
	public final static String NLS = System.getProperty("line.separator");
	public final static String SEP = System.getProperty("file.separator");
	private final static double[] POISSON_MIXTURE_LAMBDA = 
			new double[] {50, 100, 200, 500, 1000, 2000, 5000};
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
	private final long RANDOM_SEED = System.nanoTime(); //112019890314L;
	private final Random random= new Random(RANDOM_SEED);
	private final double REF_ALLELE_FRQ = 0.75;
	private final char[] BASE = new char[] {'A','C','G','T'};
	private final int BASE_NUM = BASE.length;
	private double GENETIC_LENGTH = -1.0;
	private double CENTROMERE_POS = 0.4;
	
	private final double[][] SNP_TRANS_MAT = new double[][]{
			/***To*/ /******      A      C       G       T   */
			/***From*/     /***A*/ { 0.0000, 1862.0, 6117.0, 1363.0 },
			/***C*/ { 1919.0, 0.0000, 2423.0, 5623.0 },
			/***G*/ { 6726.0, 2985.0, 0.0000, 2030.0 },
			/***T*/ { 1938.0, 7090.0, 7112.0, 0.0000 }
	};
	//private Genome genome;
	private Pedigree pedigree;
	private Chromosome0[] chromosome;
	private ArrayList<SNP> parentSNPs;
	private String scenario = "Sc1";
	private String filePath = ".";
	private PopulationData popData;
	private HashMap<String, Integer> SNPs2PosMap;
	private HashMap<String, String[]> Chromsome2SNPsMap;
	private static ExecutorService executor;
	private BlockingQueue<Runnable> tasks = null;


	public Population(String fastaFilePath, int offspring, 
			int ploidy, double cM, String scenario, String filePath) {
		this.scenario = scenario;
		this.filePath = filePath;
		this.PLOIDY = ploidy;
		this.GENETIC_LENGTH = cM;
		this.pedigree = generatePedigree(offspring);
		this.chromosome = generateChromosome(fastaFilePath);
		this.parentSNPs = generateParentSNPs();
		this.SNPs2PosMap = buildSNPs2PosMap();
		this.Chromsome2SNPsMap = buildChromsome2SNPsMap();
		createFiles();
		writeFiles();
		try {
			this.popData = PedigreeSimulate.simulate(
					filePath+SEP+scenario+".par");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void writeFiles() {
		writePedigreeFile();
		writeChromosomeFile();
		writeGeneFile();
		writeMapFile();
		writeParameterFile();
		return;
	}

	private void createFiles() {
		File theDir = new File(filePath);
		if (!theDir.exists() || theDir.exists() && !theDir.isDirectory())
			theDir.mkdir();
		theDir = new File(filePath+SEP+scenario);
		if (!theDir.exists() || theDir.exists() && !theDir.isDirectory()) 
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
		SNP snp;
		ArrayList<String> snpslist;
		String[] snps;
		for(int i=0; i<chromosome.length; i++) {
			String name = chromosome[i].getName();
			snpslist = new ArrayList<String>();
			for(int j=0; j<parentSNPs.size(); j++) {
				snp = parentSNPs.get(j);
				if(snp.getChromosome().equals(name))
					snpslist.add(snp.getName());
			}
			snps = new String[snpslist.size()];
			snpslist.toArray(snps);
			map.put(name, snps);
		}
		return map;
	}

	public void writeGenomeFile(int sim_sample_index) {

		final int N = popData.chromCount();

		String[][] indall;
		//StringBuilder chr;
		ArrayList<Locus> loci;
		int pos;
		String str2write;

		Individual indiv =  popData.getIndiv(sim_sample_index);
		try(BufferedWriter bw = getGZIPBufferedWriter(filePath+SEP+
				scenario+SEP+
				indiv.getIndivName()+".fasta.gz")) {
			for (int c=0; c<N; c++) {
				Chromosome chrom = popData.getChrom(c);
				Chromosome0 chr = chromosome[c];
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
					String header = ">"+chr.getName()+"|"+(h+1)+":"+PLOIDY;
					bw.write(header+NLS);
					final StringBuilder oos = new StringBuilder();
					int j0 = 0, j1 = chr.getLBasePairs();
					for(int l=0; l<loci.size(); l++) {
						pos = SNPs2PosMap.get(loci.get(l).getLocusName());
						str2write = indall[h][l];
						oos.append( chr.getDnaSEQ().substring(j0, pos-str2write.length()+1) );
						oos.append(str2write);
						j0 = pos+str2write.length();
					}
					if(j0 < j1) 
						oos.append(chr.getDnaSEQ().substring(j0,j1));
					bw.write(this.formatOutput(oos.toString()));
					//debug
					myLogger.info(getSystemTime()+" "+header+" >>> done.");
				}
			}
			bw.close();
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	private String formatOutput(String dnaSEQ) {
		// TODO Auto-generated method stub
		Pattern p = Pattern.compile("(.{" + CHUNK_SIZE + "})", Pattern.DOTALL);
	    Matcher m = p.matcher(dnaSEQ);
	    StringBuilder os = new StringBuilder();
	    os.append(m.replaceAll("$1" + "\n"));
	    if(os.charAt(os.length()-1)!='\n') os.append("\n");
	    return os.toString();
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
		int N = chromosome.length;

		ArrayList<SNP> snps = new ArrayList<SNP>(); //even for parent 1 and odd for parent 2
		int next, t;
		double lambda;
		for(int i=0; i<N; i++) {
			double base = chromosome[i].getLCentiMorgan()/
					chromosome[i].getLBasePairs();
			t = 0;
			lambda = POISSON_MEAN_MAP.get(random.nextInt((int) (1/PROBS_PRECISION)));
			next = poisson(lambda);
			// System.out.println(lambda+"\t\t"+next);
			t += next;

			int lb = chromosome[i].getLBasePairs();
			String name = chromosome[i].getName();
			while(t<lb) {
				
				char[] alleles = getRandomAlleles();
				char[][] phase = getRandomPhase(alleles); //for 2 parents
				
				snps.add(new SNP(name+"."+t, name, alleles, phase[0], t, t*base));
				snps.add(new SNP(name+"."+t, name, alleles, phase[1], t, t*base));
				//debug
				//println(">>>"+t+"..."+snps.get(snps.size()-1).getPosBasePairs()+"..."+snps.get(snps.size()-1).getPosCentiMorgan());
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
	 * return a M*P char matrix,
	 * M is the number of parents (i.e. 2 here),
	 * and P is the ploidy.
	 */
	public char[][] getRandomPhase(char[] alleles) {
		char[] phaseMerge = new char[PLOIDY*2];
		int c;
		double r;
		c = 0;
		for(int i=0; i<PLOIDY*2; i++) {
			r = random.nextDouble();
			if(r <= REF_ALLELE_FRQ) {
				phaseMerge[i] = alleles[0];
				c++;
			} else {
				phaseMerge[i] = alleles[1];
			}
		}
		// if both the parents are homozygous,
		//  change the allele of a randomly chosen position
		if(c==0 | c==PLOIDY*2) {
			char alt = c==0 ? alleles[0] : alleles[1];
			r = random.nextDouble();
			int k = 0;
			while(r>=(double)(k+1)/PLOIDY/2) {
				k++;
			}
			//debug
			//println(""+r+"..."+alleles[i][0]+alleles[i][1]+"..."+c+"..."+alt+"..."+k);
			phaseMerge[k] = alt;
		}

		char[][] phase = new char[2][PLOIDY];
		for(int i=0; i<PLOIDY; i++) {
			phase[0][i] = phaseMerge[i];
			phase[1][i] = phaseMerge[i+PLOIDY];
		}
		return phase;
	}

	public char[] getRandomAlleles() {
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

		int ref = getBaseFromCusumProbs(cusumfrom,random.nextDouble()),
				alt = getBaseFromCusumProbs(cusumto[ref],random.nextDouble());
		return new char[] {BASE[ref],BASE[alt]};
	}

	public int getBaseFromCusumProbs(double[] cusum, double r) {
		int base;
		for(base=0; base<BASE_NUM; base++)
			if(cusum[base]>=r) break;
		return base;
	}

	public Chromosome0[] generateChromosome(String fastaFilePath) {

		List<String> chroms = new ArrayList<String>();
		double physicaL = 0.0;
		try {
			BufferedReader br = getBufferedReader(fastaFilePath);

			StringBuilder sb = new StringBuilder();
			String line = br.readLine();
			
			while( line!=null ) {
				if(line.startsWith(">")) {
					sb.setLength(0);
					while( (line=br.readLine())!=null && 
							!line.startsWith(">"))
						sb.append(line);
					chroms.add(sb.toString());
					physicaL += sb.length();
				}
			}
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		double CHROM_BASE_LENGTH = GENETIC_LENGTH>0 ? GENETIC_LENGTH/physicaL : 1e-6;
		
		int N = chroms.size();
		
		Chromosome0[] chromosome = new Chromosome0[N];
		for(int i=0; i<N; i++) {
			double cm = chroms.get(i).length()*CHROM_BASE_LENGTH;
			chromosome[i] = new Chromosome0("CHROM"+(i+1), chroms.get(i),
					chroms.get(i).length(), 
					cm, 
					CENTROMERE_POS*cm,
					random.nextFloat()/2+0.5, 0.0);
		}
		return chromosome;
	}

	public void writeChromosomeFile() {
		println("Writing chromosome file...");

		try(BufferedWriter bw = getBufferedWriter(filePath+SEP+scenario+".chrom")) {
			bw.write("chromosome\tlength\tcentromere\tprefPairing\tquadrivalents"+NLS);
			for(int i=0; i<chromosome.length; i++) {
				Chromosome0 chr = chromosome[i];
				bw.write(chr.getName()+"\t"+chr.getLCentiMorgan()+"\t"+chr.getCentromere()+
						"\t"+chr.getPrefPairing()+"\t"+chr.getQuadrivalent()+NLS);
			}
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

	public BufferedReader getBufferedReader(String path) throws IOException {
		BufferedReader br = null;

		if(path.endsWith(".gz")){
			br = new BufferedReader(new InputStreamReader(
					new GZIPInputStream(new FileInputStream(path))), 65536);
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

	private class Pedigree {
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

	private class Chromosome0 {
		private String name;
		private String dnaSEQ;
		private int LBasePairs;
		private double LCentiMorgan;
		private double centromere;
		private double prefPairing;
		private double quadrivalent;

		public Chromosome0(String name, String dnaSEQ, int LBasePairs,
				double LCentiMorgan, double centromere,
				double prefPairing, double quadrivalent) {
			this.name = name;
			this.dnaSEQ = dnaSEQ;
			this.LBasePairs = LBasePairs;
			this.LCentiMorgan = LCentiMorgan;
			this.centromere = centromere;
			this.prefPairing = prefPairing;
			this.quadrivalent = quadrivalent;
		}

		public String getName() {
			return name;
		}

		public String getDnaSEQ() {
			return dnaSEQ;
		}

		public int getLBasePairs() {
			return LBasePairs;
		}

		public double getLCentiMorgan() {
			return LCentiMorgan;
		}

		public double getCentromere() {
			return centromere;
		}

		public double getPrefPairing() {
			return prefPairing;
		}

		public double getQuadrivalent() {
			return quadrivalent;
		}
	}

	private class SNP {
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

	public int indivCount() {
		// TODO Auto-generated method stub
		return this.popData.indivCount();
	}
}











