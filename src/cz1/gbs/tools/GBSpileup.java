package cz1.gbs.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

public class GBSpileup extends Executor {

	private String myInputDirName = null;
	private String myKeyfile = null;
	private String[] myEnzyme = null;
	private String myOutputDir = "./";
	private int myMinQualS = 10;
	private int[] myLeadingTrim = new int[]{0};

	public GBSpileup() {
		require("freebayes");
		require("samtools");
		require("bwa");
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
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		File out = new File(this.myOutputDir);
		if(!out.exists() || out.exists()&&!out.isDirectory())
			out.mkdir();
		
		FastqToTagSequence fastq2TagSequence = 
				new FastqToTagSequence(myInputDirName, myKeyfile, 
						myEnzyme, myOutputDir+"/tags", myMinQualS, myLeadingTrim, THREADS);
		fastq2TagSequence.run();
		
		MergeTagSequence mergeTagSequence = new MergeTagSequence(
				myOutputDir+"/tags",
				myOutputDir+"/mergedTags");
		mergeTagSequence.run();
		
		TagSequenceToFastq tagSequence2Fastq = new TagSequenceToFastq(
				myOutputDir+"/mergedTags",
				myOutputDir+"/tagFastq", THREADS);
		tagSequence2Fastq.run();
		
		SamToTaxa sam2Taxa = new SamToTaxa(myOutputDir+"/sortedBam/master.sorted.bam", 
				myOutputDir+"/tagFastq/master.index.gz",
				true, myOutputDir+"/bam", THREADS);
		sam2Taxa.run();
	}
	
}
