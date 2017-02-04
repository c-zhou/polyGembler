package cz1.gbs.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import cz1.util.ArgsEngine;
import cz1.util.DirectoryCrawler;
import cz1.util.Executor;
import cz1.util.Utils;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;

public class MergeTagSequence extends Executor{

	private String myInputDirName = null;
	private String myOutputDir = "./";

	public MergeTagSequence(String myInputDirName,
			String myOutputDir) {
		this.myInputDirName = myInputDirName;
		this.myOutputDir = myOutputDir;
		this.makeOutputDir();
	}
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--input-tags		Input directory containing tag sequence files.\n"
						+ " -o/--prefix			Output directory.\n\n");
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
			myArgsEngine.add("-i", "--input-tags", true);
			myArgsEngine.add("-o", "--prefix", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-i")) {
			myInputDirName = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the location of your FASTQ files.");
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
						+ "This is might be messed up in the next step.");
				myLogger.warn(msg.toString());
			}
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		File inputDirectory = new File(this.myInputDirName);
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

		if(countFileNames.length == 1) {
			try {
				myLogger.info("Only one tag file found. Copying file...");
				Files.copy(countFileNames[0].toPath(), 
						new File(myOutputDir+"/"+"tag.cnt.gz").toPath());
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			return;
		}
		
		BufferedReader[] brs = new BufferedReader[countFileNames.length];
			
		BufferedWriter bw = Utils.getBufferedWriter(
				myOutputDir+
				System.getProperty("file.separator")+
				"tag.cnt.gz");
		Map<String, List<List<Integer>>> taxa = 
				new HashMap<String, List<List<Integer>>>();
		
		String[] s;
		String s0;
		int n = brs.length;
		long start = System.currentTimeMillis();
		int maxS = 0;
		try{
			for(int i=0; i<n; i++) {
				brs[i] = Utils.getBufferedReader(countFileNames[i], 65536);
				
				s = brs[i].readLine().split("\\s+");
				
				for(int j=1; j<s.length; j++) {
					s0 = s[j].split(":")[0];
					if(taxa.containsKey(s0)) 
						taxa.get(s0).get(i).add(j);
					else {
						List<List<Integer>> ii = new ArrayList<List<Integer>>();
						for(int k=0; k<n; k++)
							ii.add(new ArrayList<Integer>());
						ii.get(i).add(j);
						taxa.put(s0, ii);
					}
				}
				if(s.length>maxS) maxS = s.length;
			}
			
			String[] taxaName = taxa.keySet().
					toArray(new String[taxa.keySet().size()]);
			final int[][] c = new int[n][maxS];
			for(int i=0; i<taxaName.length; i++) {
				List<List<Integer>> ii = taxa.get(taxaName[i]);
				for(int j=0; j<n; j++) {
					List<Integer> jj = ii.get(j);
					for(Integer k : jj)
						c[j][k] = i;
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
			long all = 0;
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
				for(int cc : count) 
					bw.write("\t"+cc);
				bw.write("\n");
				all++; 
				if(all%1000000==0) myLogger.info(all+" tags");
			}

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
}
