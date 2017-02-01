package cz1.appl;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.BitSet;

import cz1.model.BaseEncoder;
import cz1.util.ArgsEngine;
import cz1.util.DirectoryCrawler;
import cz1.util.Utils;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;

public class TagToFastqPlugin {

	private final static Logger myLogger = 
			Logger.getLogger(TagToFastqPlugin.class);
	private ArgsEngine myArgsEngine = null;
	private String myInputDirName = null;
	private String myOutputDir = null;
	private int myMinCount = 1;
	private final static int mb = 1024*1024;
	private final static Runtime instance = Runtime.getRuntime();
	private final static String poly5 = "55555555555555555555555555555"
			+ "5555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555";
	
	static {
		BasicConfigurator.configure();
	}

	public static void main(String[] args) {
		TagToFastqPlugin ttf = new TagToFastqPlugin();
		ttf.setParameters(args);
		ttf.callFastq();
	}

	private void printUsage() {
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i  Input directory containing FASTQ files in text or gzipped text.\n"
						+ "     NOTE: Directory will be searched recursively and should\n"
						+ "     be written WITHOUT a slash after its name.\n\n"
						+ " -c  Minimum count of a tag to output (default is 1).\n"
						+ " -o  Output directory to contain .cnt files (one per FASTQ file, defaults to input directory).\n\n");

	}

	public void setParameters(String[] args) {
		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}

		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-i", "--input-directory", true);
			myArgsEngine.add("-c", "--min-count", true);
			myArgsEngine.add("-o", "--output-file", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-i")) {
			myInputDirName = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the location of your FASTQ files.");
		}

		if (myArgsEngine.getBoolean("-c")) {
			myMinCount = Integer.parseInt(myArgsEngine.getString("-c"));
		}

		if (myArgsEngine.getBoolean("-o")) {
			myOutputDir = myArgsEngine.getString("-o");
		} else {
			myOutputDir = myInputDirName;
		}

	}

	private void callFastq() {
		// TODO Auto-generated method stub
		File inputDirectory = new File(this.myInputDirName);
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

		BufferedWriter bw = Utils.getBufferedWriter(
				myOutputDir+
				System.getProperty("file.separator")+
				"master.fastq.gz"),
				bw2 = Utils.getBufferedWriter(
						myOutputDir+
						System.getProperty("file.separator")+
						"master.index.gz");
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
				int ss = 0;
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

	private static void usage() {
		myLogger.info("Max Memory: "+instance.maxMemory() / mb);
		myLogger.info("Total Memory: "+instance.totalMemory() / mb);
		myLogger.info("Free Memory: "+instance.freeMemory() / mb);
		myLogger.info("Used Memory: "+(instance.totalMemory()-instance.freeMemory()) / mb);
	}
}
