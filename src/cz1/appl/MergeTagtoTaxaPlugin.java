package cz1.appl;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import cz1.util.ArgsEngine;
import cz1.util.DirectoryCrawler;
import cz1.util.Utils;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;

public class MergeTagtoTaxaPlugin {

	private final static Logger myLogger = 
			Logger.getLogger(MergeTagtoTaxaPlugin.class);
	private ArgsEngine myArgsEngine = null;
	private String myInputDirName = null;
	private String myOutputDir = null;
	private int myMinCount = 1;
	private final static int mb = 1024*1024;
	private final static Runtime instance = Runtime.getRuntime();

	static {
		BasicConfigurator.configure();
	}

	public static void main(String[] args) {
		MergeTagtoTaxaPlugin ftt = new MergeTagtoTaxaPlugin();
		ftt.setParameters(args);
		ftt.mergeTags();
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

	private void mergeTags() {
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

		BufferedReader[] brs = new BufferedReader[countFileNames.length];
		BufferedWriter bw = Utils.getBufferedWriter(
				myOutputDir+
				System.getProperty("file.separator")+
				"tag.cnt.gz");
		Map<String, int[]> taxa = new HashMap<String, int[]>();
		String[] s, s0;
		int n = brs.length;
		long start = System.currentTimeMillis();
		try{
			for(int i=0; i<n; i++) {
				brs[i] = Utils.getBufferedReader(countFileNames[i], 65536);
				
				s = brs[i].readLine().split("\\s+");
				
				for(int j=1; j<s.length; j++) {
					s0 = s[j].split(":");
					if(taxa.containsKey(s0[0])) 
						taxa.get(s0[0])[i] = (taxa.get(s0[0])[i]<<8)+j;
					else {
						int[] ii = new int[n];
						ii[i] = j;
						taxa.put(s0[0], ii);
					}
				}
			}
			
			String[] taxaName = taxa.keySet().
					toArray(new String[taxa.keySet().size()]);
			int[][] c = new int[n][97];
			for(int i=0; i<taxaName.length; i++) {
				int[] ii = taxa.get(taxaName[i]);
				for(int j=0; j<n; j++) {
					int jj = ii[j];
					while( jj>0 ) {
						int kk = jj & 127;
						c[j][kk] = i;
						jj = jj>>8;
					}
				}
			}
			bw.write("#Tag");
			for(String t : taxaName) bw.write("\t"+t);
			bw.write("\n");
			
			String[][] tags = new String[n][]; 
			for(int i=0; i<n; i++) 
				tags[i] = brs[i].readLine().split("\\s+");
			String tag, line;
			int nTaxa = taxa.keySet().size();
			int all = 0;
			int C=0;
			while(true) {
				tag = max(tags, 0);
				if(tag==null) break;
				int[] count = new int[nTaxa];
				for(int i=0; i<n; i++) {
					if(tags[i][0]!=null &&
							tag.compareTo(tags[i][0])==0) {
						for(int j=1; j<tags[i].length; j++) {
							count[c[i][j]] += 
									Integer.parseInt(tags[i][j]);
						}
						
						line = brs[i].readLine();
						if(line==null)
							tags[i][0] = null;
						else
							tags[i] = line.split("\\s+");
					}
				}
				bw.write(tag);
				for(int cc : count) {
					bw.write("\t"+cc);
					C += cc;
				}
				bw.write("\n");
				all++; 
				if(all%1000000==0)
					myLogger.info(all+" tags");
			}

			myLogger.info(C+" tags processed");
			bw.close();
			for(int i=0; i<brs.length; i++) brs[i].close();
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
		myLogger.info("Process took " + (System.currentTimeMillis() - start)/1000 + " seconds.");
	}

	private String max(String[][] tags, int ii) {
		// TODO Auto-generated method stub
		String tag = tags[0][ii];
		for(int i=1; i<tags.length; i++) {
			if(tags[i][ii]==null)
				continue;
			if(tag==null) {
				tag = tags[i][ii];
				continue;
			}
			if(tags[i][ii].length()>tag.length()|| 
					tags[i][ii].length()==tag.length() &&
					tags[i][ii].compareTo(tag)>0)
				tag = tags[i][ii];
		}
		return tag;
	}

	private static void usage() {
		myLogger.info("Max Memory: "+instance.maxMemory() / mb);
		myLogger.info("Total Memory: "+instance.totalMemory() / mb);
		myLogger.info("Free Memory: "+instance.freeMemory() / mb);
		myLogger.info("Used Memory: "+(instance.totalMemory()-instance.freeMemory()) / mb);
	}
}
