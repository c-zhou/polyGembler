package cz1.gbs.tools;

import cz1.gbs.core.BaseEncoder;
import cz1.gbs.core.ReadBarcodeResult;
import cz1.gbs.model.ParseBarcodeRead;
import cz1.util.ArgsEngine;
import cz1.util.DirectoryCrawler;
import cz1.util.Executor;
import cz1.util.ObjectSizeFetcher;
import cz1.util.Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.RejectedExecutionHandler;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPOutputStream;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;

public class FastqToTagSequence extends Executor {
	
	private String myInputDirName = null;
	private String myKeyfile = null;
	private String[] myEnzyme = null;
	private String myOutputDir = "./";
	private int myMinQualS = 10;
	private int[] myLeadingTrim = new int[]{0};; 

	public FastqToTagSequence(String myInputDirName,
			String myKeyfile, String[] myEnzyme, String myOutputDir,
			int myMinQualS, int[] myLeadingTrim, int threads) {
		this.myInputDirName = myInputDirName;
		this.myKeyfile = myKeyfile;
		this.myEnzyme = myEnzyme;
		this.myOutputDir = myOutputDir;
		this.myMinQualS = myMinQualS;
		this.myLeadingTrim = myLeadingTrim;
		this.THREADS = threads;
		this.makeOutputDir();
	}
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--input-fastq		Input directory containing FASTQ files in text or gzipped text.\n"
						+ "     					NOTE: Directory will be searched recursively and should\n"
						+ "     					be written WITHOUT a slash after its name.\n\n"
						+ " -k/--key-file			Key file listing barcodes distinguishing the samples\n"
						+ " -e/--enzyme  			Enzyme used to create the GBS library, if it differs from the one listed in the key file.\n"
						+ " -q/--min-qualS			Minimum quality score (default is 10).\n"
						+ " -t/--threads			Threads (default is 1).\n"
						+ " -T/--trim-leading		The length of leading fragments to trim off.\n"
						+ " -o/--prefix				Output directory to contain .cnt files (one per FASTQ file, defaults to input directory).\n\n");
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
			myArgsEngine.add("-i", "--input-fastq", true);
			myArgsEngine.add("-k", "--key-file", true);
			myArgsEngine.add("-e", "--enzyme", true);
			myArgsEngine.add("-q", "--min-qualS", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-T", "--trim-leading", true);
			myArgsEngine.add("-b", "--unassgined-reads", true);
			myArgsEngine.add("-o", "--prefix", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-i")) {
			myInputDirName = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the location of your FASTQ files.");
		}

		if (myArgsEngine.getBoolean("-k")) {
			myKeyfile = myArgsEngine.getString("-k");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify a barcode key file.");
		}

		if (myArgsEngine.getBoolean("-e")) {
			myEnzyme = myArgsEngine.getString("-e").split("-");
		} else {
			myLogger.warn("No enzyme specified.  Using enzyme listed in key file.");
			try {
				BufferedReader br = Utils.getBufferedReader(myKeyfile);
				String[] s = br.readLine().split("\\s+");
				int k = -1;
				for(int i=0; i<s.length; i++) 
					if(s[i].toLowerCase().equals("enzyme")) 
						k=i;
				if(k<0) throw new IllegalArgumentException("No enzyme found in the key file. "
						+ "Please specify the enzyme with -e option.\n\n");
				s = br.readLine().split("\\s+");
				myEnzyme = s[k].split("-");
				br.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		if (myArgsEngine.getBoolean("-q")) {
			myMinQualS = Integer.parseInt(myArgsEngine.getString("-q"));
		}
		
		if (myArgsEngine.getBoolean("-t")) {
			THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}
		
		if (myArgsEngine.getBoolean("-T")) {
			int leading = Integer.parseInt(myArgsEngine.getString("-T"));
			
			if(leading>0) {
				List<Integer> leadings = new ArrayList<Integer>();
				leadings.add(leading);
				for(int i=1; i<4; i++) {
					if(leading-i>=0) leadings.add(leading-i);
					leadings.add(leading+i);
				}
				myLeadingTrim = new int[leadings.size()];
				for(int i=0; i<myLeadingTrim.length; i++)
					myLeadingTrim[i] = leadings.get(i);
			}
		}
		
		if (myArgsEngine.getBoolean("-o")) {
			myOutputDir = myArgsEngine.getString("-o");
		}
		
		this.makeOutputDir();
	}
	
	private void makeOutputDir() {
		// TODO Auto-generated method stub
		File out = new File(myOutputDir);
		if(!out.exists() || out.exists()&&!out.isDirectory()) {
			out.mkdir();
		} else {
			String[] tags = out.list(new FilenameFilter() {
				public boolean accept(File dir, String name) {
					String lowercaseName = name.toLowerCase();
					if (lowercaseName.endsWith(".cnt.gz")) {
						return true;
					} else {
						return false;
					}
				}
			});
			if(tags.length>0) {
				StringBuilder msg = new StringBuilder("Tag file(s) ");
				for(int i=0; i<tags.length-1; i++) {
					msg.append(tags[i]);
					msg.append(", ");
				}
				msg.append(tags[tags.length-1]);
				msg.append(" already in the output directory. "
						+ "Please make sure "
						+ (tags.length>1 ? "they're all" : "it's")
						+ " part of this project and not stale.");
				myLogger.warn(msg.toString());
			}
		}
	}

	/**
     * Derives a tagCount list for each fastq file in the fastqDirectory.
     *
     * @param keyFileS A key file (a sample key by barcode, with a plate map
     * included).
     * @param enzyme The enzyme used to create the library (currently ApeKI or
     * PstI).
     * @param fastqDirectory Directory containing the fastq files.
     * @param outputDir Directory to which the tagCounts files.
     */

	private final static Map<BitSet, short[]> tagCounts = 
			new HashMap<BitSet, short[]>();
	private Object lock = new Object();
	private static long allReads = 0;
    private static long goodBarcodedReads = 0;
    private static String os = null;
    private static int volume = 0; 
    private static long tags = 0;
    private static ParseBarcodeRead[] thePBR;  // this reads the key file and store the expected barcodes for this lane
    private static String[] taxa; // taxa names
    private static int n; // number of taxa
    private static final double load = 0.9;
    
    @Override
    public void run() {
    	
    	String[] countFileNames = null;

        File inputDirectory = new File(this.myInputDirName);
        File[] fastqFiles = inputDirectory.listFiles(
        		new FilenameFilter() {
        			@Override
        		    public boolean accept(File dir, String name) {
        		        return name.matches("(?i).*\\.fq$|.*\\.fq\\.gz$|.*\\.fastq$|.*_fastq\\.txt$|.*_fastq\\.gz$|.*_fastq\\.txt\\.gz$|.*_sequence\\.txt$|.*_sequence\\.txt\\.gz$");
        		        //                   (?i) denotes case insensitive;                 \\. denotes escape . so it doesn't mean 'any char' & escape the backslash
        			}
        		});
        if (fastqFiles.length == 0 || fastqFiles == null) {
            myLogger.warn("Couldn't find any files that end with \".fq\", \".fq.gz\", \".fastq\", \"_fastq.txt\", \"_fastq.gz\", \"_fastq.txt.gz\", \"_sequence.txt\", or \"_sequence.txt.gz\" in the supplied directory.");
            return;
        } else {
            myLogger.info("Using the following FASTQ files:");
            countFileNames = new String[fastqFiles.length];
            for (int i = 0; i < fastqFiles.length; i++) {
                countFileNames[i] = fastqFiles[i].getName().replaceAll("(?i)\\.fq$|\\.fq\\.gz$|\\.fastq$|_fastq\\.txt$|_fastq\\.gz$|_fastq\\.txt\\.gz$|_sequence\\.txt$|_sequence\\.txt\\.gz$", "");
                //                                                                  \\. escape . so it doesn't mean 'any char' & escape the backslash
                myLogger.info(fastqFiles[i].getAbsolutePath());
            }
        }
        
        for (int laneNum = 0; laneNum < fastqFiles.length; laneNum++) {
        //for (int laneNum = 2; laneNum < 3; laneNum++) {
        	if(new File(myOutputDir+
					System.getProperty("file.separator")+
					countFileNames[laneNum]+".cnt.gz").exists()) {
        		myLogger.info("Fastq file "+fastqFiles[laneNum]+" skipped.");
        		continue;
        	}
            
        	File outputFile = new File(this.myOutputDir + File.separator + countFileNames[laneNum]);
            if (outputFile.isFile()) {
                myLogger.warn("An output file " + countFileNames[laneNum] + "\n"
                        + " already exists in the output directory for file " + fastqFiles[laneNum] + ".  Skipping.");
                continue;
            }
            myLogger.info("Reading FASTQ file: " + fastqFiles[laneNum]);
            String[] filenameField = fastqFiles[laneNum].getName().split("_");
            thePBR = new ParseBarcodeRead[this.myEnzyme.length];
            for(int i=0; i<this.myEnzyme.length; i++) {
            	if (filenameField.length == 3) {
            		thePBR[i] = new ParseBarcodeRead(this.myKeyfile, this.myEnzyme[i], filenameField[0], filenameField[1]);
            	} else if (filenameField.length == 4) {
            		thePBR[i] = new ParseBarcodeRead(this.myKeyfile, this.myEnzyme[i], filenameField[0], filenameField[2]);
            	} // B08AAABXX_s_1_sequence.txt.gz
            	else if (filenameField.length == 5) {
            		thePBR[i] = new ParseBarcodeRead(this.myKeyfile, this.myEnzyme[i], filenameField[1], filenameField[3]);
            	} else {
            		myLogger.error("Error in parsing file name: " + fastqFiles[laneNum]);
            		myLogger.error("   The filename does not contain either 3, 4, or 5 underscore-delimited values.");
            		myLogger.error("   Expect: flowcell_lane_fastq.txt.gz OR flowcell_s_lane_fastq.txt.gz OR code_flowcell_s_lane_fastq.txt.gz");
            		continue;
            	}
            }
            
            taxa = thePBR[0].getSortedTaxaNames();
            n = taxa.length;
        	os = countFileNames[laneNum];
        	volume = 0;
        	
            myLogger.info("Total barcodes found in lane:" + thePBR[0].getBarCodeCount());
            if (thePBR[0].getBarCodeCount() == 0) {
                myLogger.warn("No barcodes found.  Skipping this flowcell lane.");
                continue;
            }
            String[] taxaNames = new String[thePBR[0].getBarCodeCount()];
            for (int i = 0; i < taxaNames.length; i++) {
                taxaNames[i] = thePBR[0].getTheBarcodes(i).getTaxaName();
            }
            long start = System.currentTimeMillis();
            //executor = Executors.newFixedThreadPool(THREADS);
            initial_thread_pool();
            try {
                BufferedReader br = Utils.getBufferedReader(fastqFiles[laneNum], 65536);
              
                int block = 10000;
                String[][] Qs = new String[block][2];
                int k = 0;
                allReads = 0;
                goodBarcodedReads = 0;
                tags = 0;
                String temp = br.readLine();
                while ( temp != null ) {
                	try {
                    	Qs[k][0] = br.readLine();
                    	br.readLine();
                    	Qs[k][1] = br.readLine();
                    } catch (NullPointerException e) {
                        myLogger.error("Unable to correctly parse the sequence from fastq file.  "
                        		+ "Your fastq file may have been corrupted.");
                        System.exit(1);
                    }
                    k++;
                    temp = br.readLine();
                    if(k==block || temp==null) {
                    	
                    	if( usedMemory()/maxMemory()>load ) {
                    		writeHardDisk();
                    	}
                    	
                    	executor.submit(new Runnable() {
                    		private String[][] fastq;
							@Override
							public void run() {
								// TODO Auto-generated method stub
								try {
									final Map<BitSet, short[]> block_tagCounts = 
											new HashMap<BitSet, short[]>();
											int block_allReads = 0, block_goodBarcodedReads = 0;
											ReadBarcodeResult rr = null;
											BitSet key;
											for(int i=0; i<fastq.length; i++) {

												if(fastq[i][0]==null)
													break;
												//synchronized(lock) {
												//	allReads++;
												//}
												block_allReads++;

												outerloop:
													for(int j=0; j<myLeadingTrim.length; j++) {
														for(int k=0; k<thePBR.length; k++) {
															rr = thePBR[k].parseReadIntoTagAndTaxa (
																	fastq[i][0].substring(myLeadingTrim[j]), 
																	fastq[i][1].substring(myLeadingTrim[j]), 
																	myMinQualS);
															if(rr!=null) break outerloop;
														}
													}

												if(rr!=null) {
													key = rr.read;
													/**
											synchronized(lock) {
												goodBarcodedReads++;
												if (allReads % 1000000 == 0) {
													myLogger.info("Total Reads:" + allReads + 
															" Reads with barcode and cut site overhang:" + 
															goodBarcodedReads);
												}
												if(!tagCounts.containsKey(key)) {
													tagCounts.put(key, new short[n]);
													tags++;
												}
												tagCounts.get(key)[rr.taxonId]++;
											}
													 **/
													block_goodBarcodedReads++;
													if(!block_tagCounts.containsKey(key)) {
														block_tagCounts.put(key, new short[n]);
													}
													block_tagCounts.get(key)[rr.taxonId]++;
												}

											}

											synchronized(lock) {
												allReads += block_allReads;
												goodBarcodedReads += block_goodBarcodedReads;
												//myLogger.info("Total Reads:" + allReads + 
												//		" Reads with barcode and cut site overhang:" + 
												//		goodBarcodedReads);
												for(BitSet bs : block_tagCounts.keySet()) {
													if(!tagCounts.containsKey(bs)) {
														tags++;
														tagCounts.put(bs, block_tagCounts.get(bs));
													} else {
														short[] copy = tagCounts.get(bs);
														short[] block_copy = block_tagCounts.get(bs);
														for(int i=0; i<n; i++) copy[i] += block_copy[i];
													}
												}
											}
								} catch (Exception e) {
									Thread t = Thread.currentThread();
									t.getUncaughtExceptionHandler().uncaughtException(t, e);
									e.printStackTrace();
									executor.shutdown();
									System.exit(1);
								}
							}
							
							public Runnable init(String[][] fastq) {
						        this.fastq = fastq;
						        return(this);
						    }
                    	}.init(Qs));
                    	k=0;
                    	Qs = new String[block][2];
                    }
                }
                if(!tagCounts.isEmpty())
                	writeHardDisk();
                br.close();
                myLogger.info("Total number of reads in lane=" + allReads);
        		myLogger.info("Total number of good barcoded reads=" + goodBarcodedReads);
        		myLogger.info("Total number of tags=" + tags);
                myLogger.info("Process took " + (System.currentTimeMillis() - start)/1000 + " seconds.");
            } catch (Exception e) {
                e.printStackTrace();
            }
            myLogger.info("Finished reading " + (laneNum + 1) + " of " + fastqFiles.length + " sequence files.");

            tagCounts.clear();
        }
    }

	private void writeHardDisk() throws InterruptedException, IOException {
		// TODO Auto-generated method stub
		executor.shutdown();
		executor.awaitTermination(365, TimeUnit.DAYS);
		double load_this = 0;
		if( (load_this=usedMemory()/maxMemory())>load ) {
			myLogger.info("Memory usage exceeds"+load*100+"% of the GC limit. Load "+load_this);
		}
		myLogger.info("Volume "+volume);
		myLogger.info("Total number of tags=" + tags);
		usage();
		List<BitSet> keys = new ArrayList<BitSet>();
		keys.addAll(tagCounts.keySet());
		Collections.sort(keys, new Comparator<BitSet>() {

			@Override
			public int compare(BitSet bs, BitSet bs2) {
				// TODO Auto-generated method stub
				if(false) {
					long[] l_bs = bs.toLongArray(),
							l_bs2 = bs2.toLongArray();
					if(l_bs.length!=l_bs2.length)
						return l_bs.length>l_bs2.length ? 1 : -1;
						for(int i=0; i<l_bs.length; i++)
							if(l_bs[i]!=l_bs2[i])
								return l_bs[i]>l_bs2[i] ? 1 : -1;
								return 0;
				}
				if(bs.length()!=bs2.length())
					return bs.length()<bs2.length() ? 1 : -1;
				BitSet c_bs = (BitSet) bs.clone();
				c_bs.xor(bs2);
				if(c_bs.length()==0) return 0;
				return bs.get(c_bs.nextSetBit(0)) ? -1 : 1;
			}
		});

		String output = os+(volume==0?"":"_"+volume+"")+".cnt.gz";
		volume++;
		myLogger.info("Writing TagCount to file to "+output);
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new
				GZIPOutputStream(new FileOutputStream(
						myOutputDir+
						System.getProperty("file.separator")+
						output))), 65536);
		//BufferedWriter bw = new BufferedWriter(new FileWriter(countFileNames[laneNum]));
		bw.write("#Tag");
		for(String taxon : taxa) bw.write("\t"+taxon);
		bw.write("\n");
		//int s = 0;
		for(BitSet key : keys) {
			bw.write(os(key));
			short[] count = tagCounts.get(key);
			for(int i=0; i<count.length; i++) {
				bw.write("\t"+count[i]);
				//s += count[i];
			}
			bw.write("\n");
		}
		//System.out.println(s);
		bw.close();
		//myLogger.info("tagCounts size "+ObjectSizeFetcher.getObjectSize(tagCounts));
		tagCounts.clear();
		System.gc();
		initial_thread_pool();
	}

	private String os(BitSet bs) {
		// TODO Auto-generated method stub
		StringBuilder sb = new StringBuilder();
		for(int i=0; i<bs.length(); i++)
			sb.append(bs.get(i) ? 1 : 0);
		return sb.toString();
	}
	
}
