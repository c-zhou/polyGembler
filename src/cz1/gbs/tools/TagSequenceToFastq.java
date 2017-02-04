package cz1.gbs.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import cz1.gbs.core.BaseEncoder;
import cz1.util.ArgsEngine;
import cz1.util.DirectoryCrawler;
import cz1.util.Executor;
import cz1.util.IO;
import cz1.util.Utils;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;

public class TagSequenceToFastq extends Executor {

	private String myInputDir = null;
	private String myOutputDir = "./";
	private final static String poly5 = "55555555555555555555555555555"
			+ "5555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555";

	public TagSequenceToFastq(String myInputDir,
			String myOutputDir, int threads) {
		this.myInputDir = myInputDir;
		this.myOutputDir = myOutputDir;
		this.THREADS = threads;
		IO.makeOutputDir(this.myOutputDir);
	}
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--input-tag  	Input directory containing tag sequence files.\n"
						+ " -t/--threads		Threads (default is 1).\n"
						+ " -o/--prefix  		Output directory.\n\n");
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
			myArgsEngine.add("-i", "--input-tag", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-o", "--prefix", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-i")) {
			myInputDir = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the location of your FASTQ files.");
		}

		if (myArgsEngine.getBoolean("-t")) {
			THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}
		
		if (myArgsEngine.getBoolean("-o")) {
			myOutputDir = myArgsEngine.getString("-o");
		}
		
		IO.makeOutputDir(this.myOutputDir);
	}

	BufferedWriter bw, bw2;
	private final List<Long> synchedBlock = 
			Collections.synchronizedList(new LinkedList<Long>());
	
	public void write(final long i, final String o) 
			throws InterruptedException, IOException {
		synchronized (synchedBlock) {
			while (synchedBlock.get(0)!=i)
				synchedBlock.wait();
			bw2.write(o);
			synchedBlock.remove(0);
			synchedBlock.notifyAll();
		}
	}
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		File inputDirectory = new File(this.myInputDir);
		File[] countFileNames = inputDirectory.listFiles(
        		new FilenameFilter() {
        			@Override
        		    public boolean accept(File dir, String name) {
        		        return name.matches("(?i).*\\.cnt$|.*\\.cnt\\.gz$|.*_cnt\\.txt$|.*_cnt\\.txt\\.gz$");
        		        //                   (?i) denotes case insensitive;                 \\. denotes escape . so it doesn't mean 'any char' & escape the backslash
        			}
        		});
		if (countFileNames.length == 0 || countFileNames == null) {
			myLogger.warn("Couldn't find any files that end with \".cnt\", \".cnt.gz\", \"_cnt.txt\", or \"_cnt.txt.gz\" in the supplied directory.");
			return;
		} else {
			myLogger.info("Using the following TagCount files:");
			for (int i = 0; i < countFileNames.length; i++) 
				myLogger.info(countFileNames[i].getAbsolutePath());
		}
		
		bw = Utils.getBufferedWriter(
				myOutputDir+
				System.getProperty("file.separator")+
				"master.fastq.gz");
		bw2 = Utils.getBufferedWriter(
				myOutputDir+
				System.getProperty("file.separator")+
				"master.index.gz");
		
		long start = System.currentTimeMillis();
		this.initial_thread_pool();
		try{
			for(int i=0; i<countFileNames.length; i++) {
				BufferedReader br = Utils.getBufferedReader(
						countFileNames[i], 65536);
				String line = br.readLine();
				String[] s = line.split("\\s+");
				final String[] taxa = Arrays.copyOfRange(s, 1, s.length);
				final int block = 10000;
				String[] Qs = new String[block];
				int k = 0;
				String temp = br.readLine();
				long bS = 0;
				while( temp!=null ) {
					Qs[k] = temp;
					k++;
					temp = br.readLine();
					
					if(k==block || temp==null) {
						
						synchronized (synchedBlock) {synchedBlock.add(bS);}
						
						executor.submit(new Runnable() {
							private String[] tag;
							private long block_i;
							@Override
							public void run() {
								// TODO Auto-generated method stub
								
								final long start_i = block_i*block;
								final StringBuilder fastq = new StringBuilder();
								final StringBuilder index = new StringBuilder();
								int c;
								String[] s;
								String sequence;
								for(int i=0; i<tag.length; i++) {
									if(tag[i]==null) break;
									
									s = tag[i].split("\\s+");
									sequence = BaseEncoder.getSequenceFromBitSet(str2BitSet(s[0]));
									
									fastq.append("@");
									fastq.append(start_i+i);
									fastq.append("\n");
									fastq.append(sequence);
									fastq.append("\n+\n");
									fastq.append(poly5.substring(0, sequence.length()));
									fastq.append("\n");
									
									index.append(start_i+i);
									index.append("\t");
									for(int k=1; k<s.length; k++) {
										c = Integer.parseInt(s[k]);
										if(c>0) {
											index.append(taxa[k-1]);
											index.append('#');
											index.append(c);
											index.append(':');
										}
									}
									index.setLength(index.length()-1);
									index.append("\n");
								}
								try {
									bw.write(fastq.toString());
									write(block_i, index.toString());
								} catch (InterruptedException | IOException e) {
									// TODO Auto-generated catch block
									e.printStackTrace();
								}
								
								if(start_i%1000000==0)
									myLogger.info(start_i+" tags processed");
							}
							
							public Runnable init(final String[] Qs, final long bS) {
						        this.tag = Qs;
						        this.block_i = bS;
						        return(this);
						    }
						}.init(Qs, bS));
						bS++;
                    	k = 0;
                    	Qs = new String[block];
					}
				}
				br.close();
			}
			executor.shutdown();
			executor.awaitTermination(365, TimeUnit.DAYS);
			bw.close();
			bw2.close();
		} catch(IOException | InterruptedException e) {
			e.printStackTrace();
			System.exit(1);
		}
		myLogger.info("Process took " + (System.currentTimeMillis() - start)/1000 + " seconds.");
	}

	
	public void single_thread_run() {
		// TODO Auto-generated method stub
		File inputDirectory = new File(this.myInputDir);
		File[] countFileNames = DirectoryCrawler.listFiles("(?i).*\\.cnt$|.*\\.cnt\\.gz$|.*_cnt\\.txt$|.*_cnt\\.txt\\.gz$", inputDirectory.getAbsolutePath());
		//                                              (?i) denotes case insensitive;                 \\. denotes escape . so it doesn't mean 'any char' & escape the backslash
		if (countFileNames.length == 0 || countFileNames == null) {
			myLogger.warn("Couldn't find any files that end with \".cnt\", \".cnt.gz\", \"_cnt.txt\", or \"_cnt.txt.gz\" in the supplied directory.");
			return;
		} else {
			myLogger.info("Using the following TagCount files:");
			for (int i = 0; i < countFileNames.length; i++) 
				myLogger.info(countFileNames[i].getAbsolutePath());
		}

		long start = System.currentTimeMillis();
		try{
			for(int i=0; i<countFileNames.length; i++) {
				BufferedReader br = Utils.getBufferedReader(
						countFileNames[i], 65536);
				String line = br.readLine();
				String[] s = line.split("\\s+");
				String[] taxa = Arrays.copyOfRange(s, 1, s.length);
				StringBuilder read = new StringBuilder();
				String sequence;
				int c;
				long ss = 0;
				while( (line=br.readLine())!=null ) {
					s = line.split("\\s+");
					sequence = BaseEncoder.getSequenceFromBitSet(str2BitSet(s[0]));
					read.setLength(0);
					read.append("@");
					read.append(ss);
					read.append("\n");
					read.append(sequence);
					read.append("\n+\n");
					read.append(poly5.substring(0, sequence.length()));
					read.append("\n");
					bw.write(read.toString());
					
					read.setLength(0);
					read.append(ss);
					read.append("\t");
					for(int k=1; k<s.length; k++) {
						c = Integer.parseInt(s[k]);
						if(c>0) {
							read.append(taxa[k-1]);
							read.append('#');
							read.append(c);
							read.append(':');
						}
					}
					read.setLength(read.length()-1);
					read.append("\n");
					bw2.write(read.toString());
					ss++;
					if(ss%1000000==0)
						myLogger.info(ss+" tags processed");
				}
				br.close();
			}
			bw.close();
			bw2.close();
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
		myLogger.info("Process took " + (System.currentTimeMillis() - start)/1000 + " seconds.");
	}

	private BitSet str2BitSet(String seq) {
		// TODO Auto-generated method stub
	    BitSet bs = new BitSet(seq.length());

	    for (int i=seq.length()-1; i >= 0; i--) {
	        if (seq.charAt(i)=='1'){
	            bs.set(i);                            
	        }               
	    }
		return bs;
	}
}
