package cz1.gbs.tools;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.RejectedExecutionHandler;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader.ValidationStringency;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;

public class SamToTaxa extends Executor {
	
	private String mySamFileName = null;
	private String myIndexFileName = null;
	private boolean mySamIsSorted = false;
	private String myOutputDir = "./";

	public SamToTaxa(String mySamFileName, String myIndexFileName,
			boolean mySamIsSorted, String myOutputDir, int threads) {
		this.mySamFileName = mySamFileName;
		this.myIndexFileName = myIndexFileName;
		this.mySamIsSorted = mySamIsSorted;
		this.myOutputDir = myOutputDir;
		this.THREADS = threads;
		this.makeOutputDir();
	}
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -s/--sam-file		Input sam/bam file.\n"
						+ " -t/--threads  		Threads (default is 1).\n"
						+ " -S/--is-sorted		The sam/bam file is sorted by name.\n" 
						+ " -i/--index-file		Tag index file.\n"
						+ " -o/--prefix			Output directory. \n\n");
	}

	private void makeOutputDir() {
		// TODO Auto-generated method stub
		File out = new File(myOutputDir);
		if(!out.exists() || out.exists()&&!out.isDirectory()) {
			out.mkdir();
		}
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
			myArgsEngine.add("-i", "--index-file", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-o", "--prefix", true);
			myArgsEngine.add("-s", "--sam-file", true);
			myArgsEngine.add("-S", "--sorted", false);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-s")) {
			mySamFileName = myArgsEngine.getString("-s");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the name of your SAM/BAM file.");
		}

		if (myArgsEngine.getBoolean("-i")) {
			myIndexFileName = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the index your FASTQ file.");
		}
		
		if (myArgsEngine.getBoolean("-S")) {
			mySamIsSorted = true;
		}

		if (myArgsEngine.getBoolean("-t")) {
			THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}

		if (myArgsEngine.getBoolean("-o")) {
			myOutputDir = myArgsEngine.getString("-o");
		}
		
		this.makeOutputDir();
	}

	/**
	 * Derives a tagCount list for each fastq file in the fastqDirectory.
	 *
	 * @param keyFileS A key file (a sample key by barcode, with a plate map
	 * included).
	 * @param enzyme The enzyme used to create the library (currently ApeKI or
	 * PstI).
	 * @param fastqDirectory Directory containing the fastq files (will be
	 * recursively searched).
	 * @param outputDir Directory to which the tagCounts files (one per fastq
	 * file) will be written.
	 * @param maxGoodReads The maximum number of barcoded reads expected in a
	 * fastq file
	 * @param minCount The minimum number of occurrences of a tag in a fastq
	 * file for it to be included in the output tagCounts file
	 */

	private class Tag {
		private final String[] taxa;
		private final int[] count;

		public Tag(String[] taxa, int[] count) {
			this.taxa = taxa;
			this.count = count;
		}
	}

	private final static Map<String, BufferedWriter> sam_writers = 
			new HashMap<String, BufferedWriter>();
	private final static Map<String, SAMFileWriter> bam_writers = 
			new HashMap<String, SAMFileWriter>();
	private final static Map<Long, Tag> indexMap = 
			new ConcurrentHashMap<Long, Tag>();
	private final Map<String, Object> locks = new HashMap<String, Object>();
	private final Object lock = new Object();
	private static BufferedReader indexReader = null;
	private static long cursor = 0;
	private static long allReads = 0;
	private static String header_str = null;
	private static SAMFileHeader header_sam = null;
	private static int is = 0;

	private void cache() {
		// TODO Auto-generated method stub
		try {
			indexMap.clear();
			String temp = "";
			int block = 10000;
			String[] Qs = new String[block];
			int k = 0;
			int buffer = 10000000, cached = 0;
			this.initial_thread_pool();
			while( cached<buffer && temp!=null ) {
				temp = indexReader.readLine();
				if(temp!=null) {
					Qs[k] = temp;
					k++;
					cached++;
					cursor++;
				}
				if(k==block || temp==null) {
					executor.submit(new Runnable() {
						private String[] index;
						@Override
						public void run() {
							// TODO Auto-generated method stub
							try{
								long tag;
								String[] s, taxa;
								int[] count;
								for(String x : index) {
									if(x==null) break;
									s = x.split("\\s+");
									tag = Long.parseLong(s[0]);
									s = s[1].split("#|:");
									taxa = new String[s.length/2];
									count = new int[s.length/2];
									for(int i=0; i<s.length/2; i++) {
										taxa[i] = s[i*2];
										count[i] = Integer.parseInt(s[i*2+1]);
									}
									indexMap.put(tag, new Tag(taxa, count));
								}
							} catch (Exception e) {
								Thread t = Thread.currentThread();
								t.getUncaughtExceptionHandler().uncaughtException(t, e);
								e.printStackTrace();
								executor.shutdown();
								System.exit(1);
							}
						}
						public Runnable init(String[] index) {
							this.index = index;
							return(this);
						}
					}.init(Qs));
					k=0;
					Qs = new String[block];
				}
			}
			this.waitFor();
			myLogger.info("Cached indices. "+cached);
			myLogger.info("Tags processed. "+allReads);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private void distributeSortedBam() {
		// TODO Auto-generated method stub
		long start = System.currentTimeMillis();
		try {
			File indexFile = new File(myIndexFileName);
			myLogger.info("Using the following index file:");
			myLogger.info(indexFile.getAbsolutePath());
			File samFile = new File(mySamFileName);
			myLogger.info("Using the following BAM file:");
			myLogger.info(samFile.getAbsolutePath());
			
			final SAMFileReader inputSam = new SAMFileReader(samFile);
			header_sam = inputSam.getFileHeader();
			header_sam.setSortOrder(SortOrder.unsorted);
			
			indexReader = Utils.getBufferedReader(indexFile, 65536);
			String line = indexReader.readLine();
			String[] samples = line.split("\\s+");
			for(int i=1; i<samples.length; i++) { 
				String taxa = samples[i];
				SAMFileHeader header = header_sam.clone();
				SAMReadGroupRecord rg = new SAMReadGroupRecord(taxa);
				rg.setSample(taxa);
				header.addReadGroup(rg);
				bam_writers.put(taxa, 
						new SAMFileWriterFactory().
						makeSAMOrBAMWriter(header,
								true, new File(myOutputDir+
										System.getProperty("file.separator")+
										taxa+".bam")));
				locks.put(taxa, new Object());
			}
			
			inputSam.setValidationStringency(ValidationStringency.SILENT);
			SAMRecordIterator iter=inputSam.iterator();
			
			cache();
			this.initial_thread_pool();
			
			final int block = 10000;
			int bS = 0;
			SAMRecord[] Qs = new SAMRecord[block];
			SAMRecord temp = iter.next();
			int k = 0;
			allReads = 0;
			long tag;
			while ( temp!=null ) {
				allReads++;
				Qs[k] = temp;
				temp = iter.hasNext() ? iter.next() : null;
				k++;
				tag = temp==null ? 0 : 
					Long.parseLong(temp.getReadName());
				
				if(k==block || temp==null || tag>=cursor) {
					executor.submit(new Runnable() {
						private SAMRecord[] sam;
						private int block_i;
						
						@Override
						public void run() {
							// TODO Auto-generated method stub
							long tag;
							Tag tagObj;
							String taxa;
							for(int i=0; i<sam.length; i++) {
								try {
									if(sam[i]==null) break;
									tag = Long.parseLong(sam[i].getReadName());
									tagObj = indexMap.get(tag);
									int l = tagObj.taxa.length;
									for(int t=0; t<l; t++) {
										taxa = tagObj.taxa[t];
										sam[i].setAttribute("RG", taxa);
										int p = tagObj.count[t];
										synchronized(locks.get(taxa)) {
											for(int r=0; r<p;r++)
												bam_writers.get(taxa).addAlignment(sam[i]);
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
							
							if (block*(block_i+1) % 1000000 == 0) 
								myLogger.info("Tag processed: " + block*(block_i+1));
						}

						public Runnable init(SAMRecord[] sam, int bS) {
							this.sam = sam;
							this.block_i = bS;
							return(this);
						}
					}.init(Qs, bS));
					
					k=0;
					Qs = new SAMRecord[block];
					bS++;
					
					if(temp==null || tag>=cursor) {
						this.waitFor();
						if(temp!=null) {
							cache();
							this.initial_thread_pool();
						}
					}
				}
			}
			iter.close();
			inputSam.close();
			indexReader.close();
			//executor.shutdown();
			//executor.awaitTermination(365, TimeUnit.DAYS);
			for(String key : bam_writers.keySet())
				bam_writers.get(key).close();
			myLogger.info("Total number of reads in lane=" + allReads);
			myLogger.info("Process took " + (System.currentTimeMillis() - start)/1000 + " seconds.");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	private void distributeUnsortedBam() {
		// TODO Auto-generated method stub
		throw new RuntimeException("Unimplemented function!!!");
	}
	
	private void distributeSortedSam() {
		// TODO Auto-generated method stub
		long start = System.currentTimeMillis();
		try {
			File indexFile = new File(myIndexFileName);
			myLogger.info("Using the following index file:");
			myLogger.info(indexFile.getAbsolutePath());
			File samFile = new File(mySamFileName);
			myLogger.info("Using the following SAM file:");
			myLogger.info(samFile.getAbsolutePath());

			indexReader = Utils.getBufferedReader(indexFile, 65536);
			BufferedReader br = Utils.getBufferedReader(samFile, 65536);
			
			cache();
			this.initial_thread_pool();
			
			String temp;
			int block = 10000;
			String[] Qs = new String[block];
			int k = 0;
			allReads = 0;
			StringBuilder sb = new StringBuilder();
			while( (temp= br.readLine())!=null && 
					temp.startsWith("@")) 
				sb.append(temp+"\n");
			header_str = sb.toString();
			long tag;
			while ( temp != null ) {
				Qs[k] = temp;
				temp = br.readLine();
				k++;
				tag = temp==null ? 0 : 
					Long.parseLong(temp.substring(0,10).
							split("\\s+",2)[0]);
				if(k==block || temp==null || tag>=cursor) {

					executor.submit(new Runnable() {
						private String[] sam;
						@Override
						public void run() {
							// TODO Auto-generated method stub

							long tag;
							Tag tagObj;
							String taxa;
							for(int i=0; i<sam.length; i++) {
								try {
									if(sam[i]==null)
										break;
									synchronized(lock) {
										allReads++;
									}

									tag = Long.parseLong(sam[i].substring(0,16).
											split("\\s+",2)[0]);
									if (allReads % 1000000 == 0) 
										myLogger.info("Total Reads:" + allReads);
									synchronized(lock) {
										tagObj = indexMap.get(tag);
										//if(tagObj==null) continue;
										for(int t=0; t<tagObj.taxa.length; t++) {
											taxa = tagObj.taxa[t];
											if(!sam_writers.containsKey(taxa)) { 
												sam_writers.put(taxa, Utils.getBufferedWriter(
														myOutputDir+
														System.getProperty("file.separator")+
														taxa+".sam"));
												sam_writers.get(taxa).write(header_str);
											}

											for(int r=0; r++<tagObj.count[t];)
												sam_writers.get(taxa).write(sam[i]+"\n");
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
						}

						public Runnable init(String[] sam) {
							this.sam = sam;
							return(this);
						}
					}.init(Qs));
					
					k=0;
					Qs = new String[block];
					if(temp==null || tag>=cursor) {
						this.waitFor();
						if(temp!=null) {
							cache();
							this.initial_thread_pool();
						}
					}
				}
			}
			br.close();
			indexReader.close();
			//executor.shutdown();
			//executor.awaitTermination(365, TimeUnit.DAYS);
			for(String key : sam_writers.keySet())
				sam_writers.get(key).close();
			myLogger.info("Total number of reads in lane=" + allReads);
			myLogger.info("Process took " + (System.currentTimeMillis() - start)/1000 + " seconds.");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void distributeUnsortedSam() {

		long start = System.currentTimeMillis();
		//executor = Executors.newFixedThreadPool(THREADS);
		try {
			this.initial_thread_pool();
			File indexFile = new File(myIndexFileName);
			myLogger.info("Using the following index file:");
			myLogger.info(indexFile.getAbsolutePath());
			BufferedReader br = Utils.getBufferedReader(indexFile, 65536);
			String temp;
			int block = 10000;
			String[] Qs = new String[block];
			int k = 0;
			temp = br.readLine();
			while( temp!=null ) {
				Qs[k] = temp;
				temp = br.readLine();
				k++;
				if(k==block || temp==null) {
					executor.submit(new Runnable() {
						private String[] index;
						@Override
						public void run() {
							// TODO Auto-generated method stub
							long tag;
							String[] s, taxa;
							int[] count;
							for(String x : index) {
								s = x.split("\\s+");
								tag = Long.parseLong(s[0]);
								s = s[1].split("#|:");
								taxa = new String[s.length/2];
								count = new int[s.length/2];
								for(int i=0; i<s.length/2; i++) {
									taxa[i] = s[i*2];
									count[i] = Integer.parseInt(s[i*2+1]);
								}
								indexMap.put(tag, new Tag(taxa, count));
								synchronized(lock) {
									is++;
								}
							}
						}
						public Runnable init(String[] index) {
							this.index = index;
							return(this);
						}
					}.init(Qs));
					k=0;
					Qs = new String[block];
				}
			}
			br.close();
			this.waitFor();
			myLogger.info("Loading indices done. "+is+" tags.");

			this.initial_thread_pool();
			File samFile = new File(mySamFileName);
			myLogger.info("Using the following SAM file:");
			myLogger.info(samFile.getAbsolutePath());

			br = Utils.getBufferedReader(samFile, 65536);

			Qs = new String[block];
			k = 0;
			allReads = 0;
			StringBuilder sb = new StringBuilder();
			while( (temp= br.readLine())!=null && 
					temp.startsWith("@")) 
				sb.append(temp+"\n");
			header_str = sb.toString();

			while ( temp != null ) {
				Qs[k] = temp;
				temp = br.readLine();
				k++;
				if(k==block || temp==null) {

					executor.submit(new Runnable() {
						private String[] sam;
						@Override
						public void run() {
							// TODO Auto-generated method stub

							long tag;
							Tag tagObj;
							String taxa;
							for(int i=0; i<sam.length; i++) {
								try {
									if(sam[i]==null)
										break;
									synchronized(lock) {
										allReads++;
									}

									tag = Long.parseLong(sam[i].substring(0,16).
											split("\\s+",2)[0]);
									if (allReads % 1000000 == 0) 
										myLogger.info("Total Reads:" + allReads);
									synchronized(lock) {
										tagObj = indexMap.get(tag);
										//if(tagObj==null) continue;
										for(int t=0; t<tagObj.taxa.length; t++) {
											taxa = tagObj.taxa[t];
											if(!sam_writers.containsKey(taxa)) { 
												sam_writers.put(taxa, Utils.getBufferedWriter(
														myOutputDir+
														System.getProperty("file.separator")+
														taxa+".sam"));
												sam_writers.get(taxa).write(header_str);
											}

											for(int r=0; r++<tagObj.count[t];)
												sam_writers.get(taxa).write(sam[i]+"\n");
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
						}

						public Runnable init(String[] sam) {
							this.sam = sam;
							return(this);
						}
					}.init(Qs));
					k=0;
					Qs = new String[block];
				}
			}
			br.close();
			indexReader.close();
			this.waitFor();
			for(String key : sam_writers.keySet())
				sam_writers.get(key).close();
			myLogger.info("Total number of reads in lane=" + allReads);
			myLogger.info("Process took " + (System.currentTimeMillis() - start)/1000 + " seconds.");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		if(this.mySamIsSorted)
			if(this.mySamFileName.endsWith(".sam"))
				distributeSortedSam();
			else if(this.mySamFileName.endsWith(".bam"))
				distributeSortedBam();
			else
				throw new RuntimeException(
						"input file should be sorted sam or bam file.");
		else
			if(this.mySamFileName.endsWith(".sam"))
				distributeUnsortedSam();
			else if(this.mySamFileName.endsWith(".bam"))
				distributeUnsortedBam();
			else
				throw new RuntimeException(
						"input file should be sam or bam file.");
	}
}
