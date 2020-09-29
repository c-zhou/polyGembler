package cz1.simulation.model;

import java.util.HashMap;
import java.util.ArrayList;
import java.util.List;
import java.io.File;
import java.io.FileReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.BufferedReader;
import java.lang.StringBuilder;
import java.util.Random;
import java.util.zip.GZIPInputStream;
import java.io.OutputStreamWriter;
import java.util.zip.GZIPOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.FileInputStream;
import java.io.Reader;
import java.io.InputStreamReader;
import java.util.Calendar;
import java.text.SimpleDateFormat;
import java.util.Collections;

import org.apache.commons.lang3.RandomStringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class GBS {

	private final static Logger myLogger = LogManager.getLogger(GBS.class);
	public final static String NLS = System.getProperty("line.separator");
	public final static String SEP = System.getProperty("file.separator");
	private final int readLength = 101;
	private final int MIN_FRAGMENT_SIZE = 100;
	private double meanDepth;
	private double sdDepth;
	private final double[] probabilityBaseMissingInitialization;
	private final double[] probabilityBaseMissingSucceeding;
	private final String QUAL_SCORE = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ";
	private final HashMap<Integer, HashMap<Character, HashMap<Integer, Character>>> 
	qualityTransitionProbabilityMap;
	private final double baseSubIndelErrorRate = 0.003;
	/*** determine the base substitution probability for substitution errors*/
	private final HashMap<Character, HashMap<Integer, Character>> 
	baseSubstitutionErrorProbabilityMap;
	private final double[][] baseSubstitutionErrorProbabilityMatrix = new double[][]{
			/***         A    T    C    G  */
			/*** A */ { 0.0, 0.6, 0.2, 0.2 },
			/*** T */ { 0.6, 0.0, 0.2, 0.2 },
			/*** C */ { 0.2, 0.2, 0.0, 0.6 },
			/*** G */ { 0.2, 0.2, 0.6, 0.0 },
	};
	/*** determine the length of an indel
	 * assuming an exponential distribution with rate (inverse scale) lambda
	 * rounding up to the nearest integer
	 * the insertion sequence is randomly generated using insertionSeqFactory
	 * */
	private final double baseIndelErrorProbabilityLambda = 0.5;
	private final double PROBS_PRECISION = 1e-4;
	private final int INV_PROBS_PRECISION = (int)Math.round(1.0/PROBS_PRECISION);
	private final char[] BASE = new char[]{'A','T','C','G'};
	private final Character MISSING_BASE_SYMBOL = 'N';
	private final Character MISSING_QUAL_SCORE = '#';
	private final ArrayList<String> fastaFileList = 
			new ArrayList<String>();
	private final HashMap<String, String> fastaFileBarcodeMap = 
			new HashMap<String, String>();
	private final String COMMON_ADAPTER = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG";
	private final Enzyme enzyme;
	private final String[][] enzymeOverhang;
	private final HashMap<Character, Character> baseComplementaryMap;
	private final Random random;
	private String GBSOutputDir;
	private String GBSFastqFilePath;
	//private BufferedWriter GBSFastqFileBufferedWriter;
	private String flowcell;
	private String lane;
	private String libraryPlate;
	private String libraryPrepID;
	private String libraryPlateID;
	private String dnaPlate;
	private String genus;
	private String species;

	
	public GBS(String fastaFileDir,
			String enzymeName,
			double avg,
			double sd,
			String parameterFilePath,
			String libPrepFilePath,
			String barcodeFilePath,
			String outputDir,
			long RANDOM_SEED) {
		/***
        meanDepth = 6.0;
        sdDepth = 6.0;
        probabilityBaseMissing = new double[readLength];
        for(int i=0; i<readLength; i++) {
            probabilityBaseMissing[i] = 0.0001;
        }

        int[] firstBaseQualityCumsum = new int[QUAL_SCORE.length()];
        int[][] baseTransitionCumsum = 
            new int[QUAL_SCORE.length()][QUAL_SCORE.length()];
        for(int i=0; i<QUAL_SCORE.length(); i++) {
            firstBaseQualityCumsum[i] = (int)((i+1.0)/
                    QUAL_SCORE.length()/PROBS_PRECISION);
            for(int j=0; j<QUAL_SCORE.length(); j++) {
                baseTransitionCumsum[i][j] = 
                    (int)((j+1.0)/QUAL_SCORE.length()/PROBS_PRECISION);
            }
        }
		 **/
		this.random = new Random(RANDOM_SEED);
		this.setLibarayPreparation(libPrepFilePath);
		this.GBSOutputDir = outputDir;
		this.GBSFastqFilePath = outputDir+SEP+flowcell+"_"+lane+"_fastq";
		this.probabilityBaseMissingInitialization = new double[readLength];
		this.probabilityBaseMissingSucceeding = new double[readLength];
		int[][] baseQualityCumsum = new int[readLength][QUAL_SCORE.length()];
		int[][][] baseTransitionCumsum =
				new int[readLength][QUAL_SCORE.length()][QUAL_SCORE.length()];
		this.meanDepth = avg;
		this.sdDepth = sd;
		
		try {
			InputStream inputSteam;
			if(parameterFilePath==null) {
				inputSteam = GBS.class.getResourceAsStream("confs/QualS_Markov_chain.cfg");
			} else {
				inputSteam = new FileInputStream(parameterFilePath);
			}
			BufferedReader br = getBufferedReader(inputSteam);
			String line;
			String[] stringSplit;
			while((line=br.readLine())!=null && line.startsWith("#")){}
			stringSplit = line.split(",");
			for(int i=0; i<readLength; i++) 
				this.probabilityBaseMissingInitialization[i] = Double.parseDouble(stringSplit[i]);
			while((line=br.readLine())!=null && line.startsWith("#")){}
			stringSplit = line.split(",");
			for(int i=0; i<readLength; i++)
				this.probabilityBaseMissingSucceeding[i] = Double.parseDouble(stringSplit[i]);
			double[] w;
			for(int l=0; l<readLength; l++) {
				while((line=br.readLine())!=null && line.startsWith("#")){}
				stringSplit = line.split(",");
				w = new double[QUAL_SCORE.length()];
				for(int i=0; i<QUAL_SCORE.length(); i++) 
					w[i] = Double.parseDouble(stringSplit[i]);
				cumsum(w);
				for(int i=0; i<QUAL_SCORE.length(); i++) 
					baseQualityCumsum[l][i] = (int)Math.round(w[i]/PROBS_PRECISION);
				for(int i=0; i<QUAL_SCORE.length(); i++) {
					while((line=br.readLine())!=null && line.startsWith("#")){}
					stringSplit = line.split(",");
					w = new double[QUAL_SCORE.length()];
					for(int j=0; j<QUAL_SCORE.length(); j++) 
						w[j] = Double.parseDouble(stringSplit[j]);
					cumsum(w);
					for(int j=0; j<QUAL_SCORE.length(); j++)
						baseTransitionCumsum[l][i][j] = 
								(int)Math.round(w[j]/PROBS_PRECISION);
					if(baseTransitionCumsum[l][i][QUAL_SCORE.length()-1]==0)
						baseTransitionCumsum[l][i] = baseQualityCumsum[l];
				}
			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}

		this.qualityTransitionProbabilityMap =
				new HashMap<Integer, HashMap<Character, HashMap<Integer, Character>>>();
		// this is to record the probability of the quality score of the fisrt base
		// the key is null
		HashMap<Integer, Character> baseTransitionBin = new HashMap<Integer, Character>();
		HashMap<Character, HashMap<Integer, Character>> baseTransitionMap =
				new HashMap<Character, HashMap<Integer, Character>>();
		int slider;

		for(int l=0; l<readLength; l++) {
			baseTransitionMap = new HashMap<Character, HashMap<Integer, Character>>();

			baseTransitionBin = new HashMap<Integer, Character>();
			slider = 0;
			for(int i=0; i<QUAL_SCORE.length(); i++) {
				while(slider<baseQualityCumsum[l][i])
					baseTransitionBin.put(slider++,QUAL_SCORE.charAt(i));
			}
			baseTransitionMap.put(null, baseTransitionBin);

			for(int i=0; i<QUAL_SCORE.length(); i++) {
				baseTransitionBin = new HashMap<Integer, Character>();    
				slider = 0;
				for(int j=0; j<QUAL_SCORE.length(); j++) {
					while(slider<baseTransitionCumsum[l][i][j])
						baseTransitionBin.put(slider++, QUAL_SCORE.charAt(j));
				}
				baseTransitionMap.put(QUAL_SCORE.charAt(i), baseTransitionBin);
			}
			this.qualityTransitionProbabilityMap.put(l, baseTransitionMap);
		}

		myLogger.info("Restriction Enzyme - "+enzymeName);
		this.enzyme = new Enzyme(enzymeName);
		this.baseComplementaryMap = getBaseComplementaryMap();
		this.baseSubstitutionErrorProbabilityMap = getBaseSubstitutionErrorProbabilityMap();
		this.enzymeOverhang = enzyme.getOverhang();
		this.getFastaFileList(fastaFileDir);
		this.getFastaFileBarcodeMap(fastaFileList, enzymeName, this.setBarcode(barcodeFilePath));
	}

	private InputStream setBarcode(String barcodeFilePath) {
		// TODO Auto-generated method stub
		if(barcodeFilePath!=null)
			try {
				return new FileInputStream(barcodeFilePath);
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		switch(this.enzyme.getName().toLowerCase()) {
		case "apeki":
			myLogger.info("Using barcode file 768-barcodes-ApeKI-size-4-9-v2.0.csv "
					+ "generated by http://www.deenabio.com/services/gbs-adapters.");
			return GBS.class.getResourceAsStream(
					"confs/768-barcodes-ApeKI-size-4-9-v2.0.csv");
		case "psti":
			myLogger.info("Using barcode file 768-barcodes-PstI-size-4-9-v2.0.csv "
					+ "generated by http://www.deenabio.com/services/gbs-adapters.");
			return GBS.class.getResourceAsStream(
					"confs/768-barcodes-PstI-size-4-9-v2.0.csv");
		case "tsei":
			myLogger.info("Using barcode file 768-barcodes-TseI-size-4-9-v2.0.csv "
					+ "generated by http://www.deenabio.com/services/gbs-adapters.");
			return GBS.class.getResourceAsStream(
					"confs/768-barcodes-TseI-size-4-9-v2.0.csv");	
		default:
			throw new RuntimeException("No barcode file available for enzyme "
		+this.enzyme.getName()+".\n"
				+ "You need to provide a customer barcode file with -b/--barcode-file."
				+ "You may generate barcodes from http://www.deenabio.com/services/gbs-adapters.\n");
		}
	}

	private void cumsum(double[] w) {
		// TODO Auto-generated method stub
		double s = 0;
		for(int i=0; i<w.length; i++) {
			s += w[i];
			w[i] = s;
		}
		if(s==0) return;
		for(int i=0; i<w.length; i++) w[i] /= s;
	}

	private void setLibarayPreparation(String libPrepFilePath) {
		if(libPrepFilePath==null) {
			setLibarayPreparation();
			return;
		}
		String line;
		String[] s;
		try(BufferedReader br = getBufferedReader(libPrepFilePath)) {
			while( (line=br.readLine())!=null && !line.startsWith("#")) {
				s = line.split("\\s:\\s");
				switch(s[0].trim().toUpperCase()) {
				case "FLOWCELL":
					this.flowcell = s[1].trim();
					break;
				case "LANE":
					this.lane = s[1].trim();
					break;
				case "LIBRARYPLATE":
					this.libraryPlate = s[1].trim();
					break;
				case "LIBRARYPREPID":
					this.libraryPrepID = s[1].trim();
					break;
				case "LIBRARYPLATEID":
					this.libraryPlateID = s[1].trim();
					break;
				case "DNA_PLATE":
					this.dnaPlate = s[1].trim();
					break;
				case "GENUS":
					this.genus = s[1].trim();
					break;
				case "SPECIES":
					this.species = s[1].trim();
					break;
				default:
					System.out.println("Invalid library preparation id. Program exit...");
					System.exit(1);
				}
			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void setLibarayPreparation() {
		// TODO Auto-generated method stub
		this.flowcell = RandomStringUtils.randomAlphanumeric(8).toUpperCase();
		this.lane = ""+1;
		this.libraryPlate = "Simulation";
		this.libraryPrepID = "R";
		this.libraryPlateID = "Random";
		this.dnaPlate = "Simulation";
		this.genus = "G";
		this.species = "S";
	}

	private HashMap<Character, Character> getBaseComplementaryMap() {
		HashMap<Character, Character> map = new HashMap<Character, Character>();
		map.put('A','T'); 
		map.put('T','A'); 
		map.put('C','G'); 
		map.put('G','C');
		return map;
	}

	private HashMap<Character, HashMap<Integer, Character>> getBaseSubstitutionErrorProbabilityMap(){
		HashMap<Character, HashMap<Integer, Character>> map = 
				new HashMap<Character, HashMap<Integer, Character>>();
		HashMap<Integer, Character> bin;
		int cumsum, slider;
		for(int i=0; i<4; i++) {
			cumsum=0;
			slider=0;
			bin = new HashMap<Integer, Character>();
			for(int j=0; j<4; j++) {
				cumsum+=Math.round(baseSubstitutionErrorProbabilityMatrix[i][j]/PROBS_PRECISION);
				while(slider<cumsum) bin.put(slider++,BASE[j]);
			}
			map.put(BASE[i],bin);
		}
		return map;
	}

	private void getFastaFileBarcodeMap(ArrayList<String> fastaFileList, 
			String enzymeName, InputStream barcodeFileStream) {
		fastaFileBarcodeMap.clear();
		try{
			//BarcodeGenerator barcodeGenerator = new BarcodeGenerator(new String[]
			//        {"-b", ""+fastaFileList.size(), "-e", enzymeName});
			//barcodeGenerator.runBarcodeGenerator();
			//BufferedReader br = getBufferedReader("barcode_list.txt");
			String line;
			BufferedReader br = getBufferedReader(barcodeFileStream);
			for(int i=0; i<16; i++) br.readLine();
			List<String> barcodes = new ArrayList<String>();
			while( (line=br.readLine())!=null && line.length()>0)
				barcodes.add(line.split(",")[1]);
			br.close(); 
			Collections.shuffle(barcodes);
			
			BufferedWriter bw = getBufferedWriter(GBSOutputDir+SEP+flowcell+"_"+lane+"_key"+".txt");
			bw.write("Flowcell\tLane\tBarcode\tDNASample\tLibraryPlate\tRow\tCol\t"+
					"LibraryPrepID\tLibraryPlateID\tEnzyme\tBarcodeWell\tDNA_Plate\t"+
					"SampleDNA_Well\tGenus\tSpecies\tPedigree\tPopulation\tSeedLot\t"+
					"FullSampleName\n");
			String dnaSample, paddedcol, pedigree;
			char row; int col;
			for(int i=0; i<this.size(); i++) {
				fastaFileBarcodeMap.put(fastaFileList.get(i),barcodes.get(i));
				dnaSample = parseSampleName(new File(fastaFileList.get(i)).getName());
				pedigree = dnaSample.contains("P1") || dnaSample.contains("P2") ? "NA" : "P1xP2";
				row = (char) (i/12+'A');
				col = i%12+1;
				paddedcol = col<10 ? "0"+col : ""+col;
				bw.write(flowcell+"\t");
				bw.write(lane+"\t");
				bw.write(barcodes.get(i)+"\t");
				bw.write(dnaSample+"\t");
				bw.write(libraryPlate+"\t");
				bw.write(row+"\t");
				bw.write(col+"\t");
				bw.write(libraryPrepID+System.nanoTime()+paddedcol+"\t");
				bw.write(libraryPlateID+"\t");
				bw.write(enzymeName+"\t");
				bw.write(row+paddedcol+"\t");
				bw.write(dnaPlate+"\t");
				bw.write(row+paddedcol+"\t");
				bw.write(genus+"\t");
				bw.write(species+"\t");
				bw.write(pedigree+"\t");
				bw.write(pedigree+"\t");
				bw.write("\t");
				bw.write(dnaSample+":"+flowcell+":"+lane+":"+libraryPrepID+paddedcol+"\n");
			}
			bw.close();
			//} catch (IOException | StopExcecutionException e) {
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	private void getFastaFileList(String fastaFileDir) {
		fastaFileList.clear();
		File folder = new File(fastaFileDir);
		File[] listOfFiles = folder.listFiles();
		for (File file : listOfFiles) {
			if (file.isFile() && isFastaFile(file.getName())) {
				fastaFileList.add(file.getPath());
			}
		}
	}

	public boolean isFastaFile(String fileName) {
		return fileName.endsWith("fasta") | fileName.endsWith("fasta.gz")
				| fileName.endsWith("fa") | fileName.endsWith("fa.gz") ? true : false;
	}


	public String parseSampleName(String fileName) {
		return new File(fileName).getName().replace(".fasta.gz","").
				replace(".fasta","").replace(".fa.gz","").replace(".fa","");
	}

	public void simulate(final int sim_sample_index) {

		try {
			//GBSFastqFileBufferedWriter = getGZIPBufferedWriter(GBSFastqFilePath);
			String name, line, fastaFilePath;
			ArrayList<Integer> cut, recognization;
			StringBuilder chromosome;
			FastqRead fastq;
			int l;
			Digest digestion;
			fastaFilePath = fastaFileList.get(sim_sample_index);
			StringBuilder oos = new StringBuilder();
			//try{
			BufferedWriter GBSFastqFileBufferedWriter = getGZIPBufferedWriter(
					GBSFastqFilePath+"_"+sim_sample_index+".gz", 65536);
			BufferedReader br = getBufferedReader(fastaFilePath);
			line = br.readLine();
			while( line != null ) {
				name = parseSampleName(fastaFilePath)+":"+line.replaceFirst(">","").replaceAll(" ", "_");
				chromosome = new StringBuilder();
				while( (line=br.readLine()) !=null && !line.startsWith(">") ) {
					chromosome.append(line);
				}
				l = chromosome.length();

				//System.out.println(getSystemTime()+">>> digest starting...");
				digestion = new Digest(enzyme, chromosome.toString());
				//System.out.println(getSystemTime()+">>> done.");
				cut = digestion.getCutsite();
				recognization = digestion.getRecognization();

				//System.out.println(getSystemTime()+">>> simulate starting...");
				// simulate 5'-3' end genome
				int coverage;
				for(int i=2; i<cut.size()-1; i++) {
					if(cut.get(i)-cut.get(i-1)<MIN_FRAGMENT_SIZE) continue;
					coverage = (int) Math.round(meanDepth+random.nextGaussian()*sdDepth);
					for(int j=0; j<coverage; j++) {
						fastq = generateFastqRead(chromosome.substring(cut.get(i-1),cut.get(i)), 
								recognization.get(i-1),
								recognization.get( i ),
								fastaFileBarcodeMap.get(fastaFilePath), 
								"@"+name+":"+cut.get(i-1)+":0:"+j);
						//writeFastqRead(fastq);
						oos.setLength(0);
						oos.append(fastq.identifier);
						oos.append(NLS);
						oos.append(fastq.sequence);
						oos.append(NLS);
						oos.append(fastq.plus);
						oos.append(NLS);
						oos.append(fastq.quality);
						oos.append(NLS);
						GBSFastqFileBufferedWriter.write(oos.toString());
					}
				}

				// simulate reverse complementary genome
				chromosome.reverse();
				for(int i=0; i<l; i++) 
					chromosome.setCharAt(i,baseComplementaryMap.get(chromosome.charAt(i)));
				
				Collections.reverse(recognization);
				int r;
				if(enzyme.getRecognizationSequence().length>1) {
					if(enzyme.getRecognizationSequence().length%2==1)
						throw new RuntimeException("!!!");
					for(int i=1; i<recognization.size()-1; i++) {
						r = recognization.get(i);
						recognization.set(i,r+(r%2==0?1:-1));
					}
				}
				
				Collections.reverse(cut);
				for(int i=1; i<cut.size()-1; i++) {
					r = recognization.get(i);
					cut.set(i,l-cut.get(i)+enzymeOverhang[r][0].length()-enzymeOverhang[r][1].length());
				}
				
				for(int i=2; i<cut.size()-1; i++) {
					if(cut.get(i)-cut.get(i-1)<MIN_FRAGMENT_SIZE) continue;
					coverage = (int) Math.round(meanDepth+random.nextGaussian()*sdDepth);
					for(int j=0; j<coverage; j++) {
						fastq = generateFastqRead(chromosome.substring(cut.get(i-1),cut.get(i)), 
								recognization.get(i-1),
								recognization.get( i ),
								fastaFileBarcodeMap.get(fastaFilePath), 
								"@"+name+":"+cut.get(i-1)+":1:"+j);
						//writeFastqRead(fastq);
						oos.setLength(0);
						oos.append(fastq.identifier);
						oos.append(NLS);
						oos.append(fastq.sequence);
						oos.append(NLS);
						oos.append(fastq.plus);
						oos.append(NLS);
						oos.append(fastq.quality);
						oos.append(NLS);
						GBSFastqFileBufferedWriter.write(oos.toString());
					}
				}
				//System.out.println(getSystemTime()+">>> done.");
				System.out.println(getSystemTime()+">>> "+name+" done.");
			}
			br.close();
			GBSFastqFileBufferedWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}

		//GBSFastqFileBufferedWriter.close();
		//} catch (IOException | InterruptedException e) {
		//} catch (InterruptedException e) {
		//	e.printStackTrace();
		//	System.exit(1);
		//}

	}

	public String insertionSeqFactory(int length) {
		StringBuilder insertion = new StringBuilder();
		for(int i=0; i<length; i++) {
			insertion.append(BASE[random.nextInt(4)]);
		}
		return insertion.toString();
	}

	public FastqRead generateFastqRead(String dnaInsert, int recog5, int recog3, 
			String barcode, String identifier) {
		String toSequence = barcode+
				enzymeOverhang[recog5][0]+
				dnaInsert+
				enzymeOverhang[recog3][1]+
				COMMON_ADAPTER;
		StringBuilder sequence = new StringBuilder();
		StringBuilder quality = new StringBuilder();
		Character prev_base = null, prev_qual = null;
		int b = 0, p = 0, indel, l = toSequence.length(), u;
		int sub = 0, ins = 0, del = 0;
		while(b<readLength && p<l) {
			if(random.nextFloat() < baseSubIndelErrorRate) {
				indel = random.nextInt(3);
				if(indel==0) { /*** substitution */
					prev_base = baseSubstitutionErrorProbabilityMap
							.get(toSequence.charAt(p)).get(random.nextInt(INV_PROBS_PRECISION));
					prev_qual = qualityTransitionProbabilityMap.get(b).get(prev_qual)
							.get(random.nextInt(INV_PROBS_PRECISION));
					sequence.append(prev_base);
					quality.append(prev_qual);
					p++; 
					b++;
					sub++;
				} else if (indel==1) { /*** insertion */
					u = (int) Math.ceil(-Math.log(1.0-random.nextFloat())/baseIndelErrorProbabilityLambda);
					u = Math.max(1, Math.min(u, readLength-b));
					sequence.append(insertionSeqFactory(u));
					for(int t=0; t<u; t++) {
						prev_qual = qualityTransitionProbabilityMap.get(b).get(prev_qual)
								.get(random.nextInt(INV_PROBS_PRECISION));
						quality.append(prev_qual);
						b++;
					}
					prev_base = sequence.charAt(b-1);
					ins++;
				} else { /*** deletion */
					p += Math.ceil(-Math.log(1.0-random.nextFloat())/baseIndelErrorProbabilityLambda);
					del++;
				}
			} else {
				if(prev_base==MISSING_BASE_SYMBOL && 
						random.nextFloat()<probabilityBaseMissingSucceeding[b]) {
					sequence.append(MISSING_BASE_SYMBOL);
					quality.append(MISSING_QUAL_SCORE);
				} else if (random.nextFloat()<probabilityBaseMissingInitialization[b]) {
					prev_base = MISSING_BASE_SYMBOL;
					prev_qual = null;
					sequence.append(MISSING_BASE_SYMBOL);
					quality.append(MISSING_QUAL_SCORE);
				} else {
					prev_base = toSequence.charAt(p);
					prev_qual = qualityTransitionProbabilityMap.get(b).get(prev_qual)
							.get(random.nextInt(INV_PROBS_PRECISION));
					sequence.append(prev_base);
					quality.append(prev_qual);
				}
				p++; 
				b++;
			}
		}
		while(b++<readLength) {
			sequence.append(MISSING_BASE_SYMBOL);
			quality.append(MISSING_QUAL_SCORE);
		}
		return new FastqRead(identifier+":"+sub+":"+ins+":"+del, sequence.toString(),"+", quality.toString());
	}

	public void writeFastqRead(FastqRead fastq, BufferedWriter writer) throws IOException {
		writer.write(fastq.identifier+NLS);
		writer.write(fastq.sequence+NLS);
		writer.write(fastq.plus+NLS);
		writer.write(fastq.quality+NLS);
	}

	//public void writeFastqRead(FastqRead fastq) throws IOException {
	//	GBSFastqFileBufferedWriter.write(fastq.identifier+NLS);
	//	GBSFastqFileBufferedWriter.write(fastq.sequence+NLS);
	//	GBSFastqFileBufferedWriter.write(fastq.plus+NLS);
	//	GBSFastqFileBufferedWriter.write(fastq.quality+NLS);
	//}

	public BufferedReader getBufferedReader(String path) throws IOException {
		BufferedReader br = null;
		if(path.endsWith(".gz")){
			br = new BufferedReader(new InputStreamReader(new GZIPInputStream(
					new FileInputStream(path))), 65536);
		}else{
			br = new BufferedReader(new FileReader(path), 65536);
		}
		return br;
	}
	
	public BufferedReader getBufferedReader(InputStream inputStream) throws IOException {
		return new BufferedReader(new InputStreamReader(inputStream), 65536);
	}

	public BufferedWriter getBufferedWriter(String path) throws IOException {
		return new BufferedWriter(new FileWriter(new File(path)));
	}

	public BufferedWriter getGZIPBufferedWriter(String path) throws IOException {
		return new BufferedWriter(new OutputStreamWriter(new
				GZIPOutputStream(new FileOutputStream(path))));
	}
	
	public BufferedWriter getGZIPBufferedWriter(String path, int size) throws IOException {
		return new BufferedWriter(new OutputStreamWriter(new
				GZIPOutputStream(new FileOutputStream(path))), size);
	}


	public static String getSystemTime(){
		return new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").
				format(Calendar.getInstance().getTime());
	}

	private class FastqRead {
		String identifier;
		String sequence;
		String plus;
		String quality;

		public FastqRead(String identifier, String sequence,
				String plus, String quality) {
			this.identifier = identifier;
			this.sequence = sequence;
			this.plus = plus;
			this.quality = quality;
		}
	}

	public int size() {
		// TODO Auto-generated method stub
		return this.fastaFileList.size();
	}
}





