package cz1.gbs.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

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
	private String myReference = null;
	private int myPloidy = 2;
	
	public GBSpileup() {
		this.require("freebayes");
		this.require("samtools");
		this.require("bwa");
	}

	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--input-fastq        Input directory containing FASTQ files in text or gzipped text.\n"
						+ "                         NOTE: Directory will be searched recursively and should\n"
						+ "                         be written WITHOUT a slash after its name.\n"
						+ " -k/--key-file           Key file listing barcodes distinguishing the samples\n"
						+ " -e/--enzyme             Enzyme used to create the GBS library, if it differs from the one listed in the key file.\n"
						+ " -q/--min-qualS          Minimum quality score (default is 10).\n"
						+ " -p/--ploidy             Ploidy of the genome (default is 2).\n"
						+ " -t/--threads            Threads (default is 1).\n"
						+ " -T/--trim-leading       The length of leading fragments to trim off.\n"
						+ " -f/--reference          The reference genome (in fasta format).\n"
						+ " -o/--prefix             Output directory to contain .cnt files (one per FASTQ file, defaults to input directory).\n\n");
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
			myArgsEngine.add("-p", "--ploidy", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-T", "--trim-leading", true);
			myArgsEngine.add("-b", "--unassgined-reads", true);
			myArgsEngine.add("-f", "--reference", true);
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
		
		if (myArgsEngine.getBoolean("-f")) {
			myReference = myArgsEngine.getString("-f");
			if(!new File(myReference+".amb").exists() || 
					!new File(myReference+".ann").exists() || 
					!new File(myReference+".bwt").exists() ||
					!new File(myReference+".pac").exists() ||
					!new File(myReference+".sa").exists()) {
				String index = "bwa index -p "+myReference+" -a bwtsw "+myReference;
				try {
					this.bash(index).waitFor();
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the reference.");
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
		
		if (myArgsEngine.getBoolean("-p")) {
			myPloidy = Integer.parseInt(myArgsEngine.getString("-p"));
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
		
		File out = new File(myOutputDir);
		if (out.exists()&&out.isDirectory()) {
			myLogger.warn("Output directory "+myOutputDir+" exsits. "
					+ "We strongly recommend a new location.");
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		
		Utils.makeOutputDir(this.myOutputDir);
		
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
		
		Utils.makeOutputDir(myOutputDir+"/tagBam");
		String bwaMem = "bwa mem -t "+THREADS+" "+myReference+
				" "+myOutputDir+"/tagFastq/master.fastq.gz > "+myOutputDir+"/tagBam/master.sam";
		String sam2bam = "samtools view -S -b "+myOutputDir+"/tagBam/master.sam > "+
				myOutputDir+"/tagBam/master.bam";
		String clearSam = "rm "+myOutputDir+"/tagBam/master.sam";
		String sortBam = "samtools sort -n -o "+myOutputDir+"/tagBam/master.sorted.bam "+ 
				myOutputDir+"/tagBam/master.bam";
		try {
			this.bash(bwaMem).waitFor();
			this.bash(sam2bam).waitFor();
			this.bash(clearSam).waitFor();
			this.bash(sortBam).waitFor();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		SamToTaxa sam2Taxa = new SamToTaxa(myOutputDir+"/tagBam/master.sorted.bam", 
				myOutputDir+"/tagFastq/master.index.gz",
				true, myOutputDir+"/bam", THREADS);
		sam2Taxa.run();
		
		File[] bams = new File(myOutputDir+"/bam").listFiles();
		String[] commands = new String[bams.length];
		for(int i=0; i<commands.length; i++) {
			String bam = bams[i].getName();
			String newBam = bam.replaceAll(".bam$", ".sorted.bam");
			commands[i] = "samtools sort -o "+myOutputDir+"/bam/"+newBam+" "+ 
					myOutputDir+"/bam/"+bam+" && rm "+myOutputDir+"/bam/"+bam;
		}
		this.bash(commands);
		
		this.generateSplitReference();
		
		SamFileSplit samFileSplit = 
				new SamFileSplit(myOutputDir+"/bam", myOutputDir+"/bed", 
						myOutputDir+"/splitBam", THREADS);
		samFileSplit.run();
		
		List<String> comm = new ArrayList<String>();
		File[] splitBam = new File(myOutputDir+"/splitBam").listFiles();
		for(File f : splitBam) {
			File[] sB = f.listFiles();
			for(File f2 : sB)
				comm.add("samtools index "+f2.getAbsolutePath());
		}
		commands = new String[comm.size()];
		comm.toArray(commands);
		this.bash(commands);
		
		Utils.makeOutputDir(myOutputDir+"/freebayes");
		Utils.makeOutputDir(myOutputDir+"/freebayes/bam_list");
		Utils.makeOutputDir(myOutputDir+"/freebayes/vcf_list");
		commands = new String[THREADS];
		for(int i=0; i<THREADS; i++)
			commands[i] = "ls "+myOutputDir+"/splitBam/"+i+"/*.bam > "
					+ myOutputDir+"/freebayes/bam_list/"+i+".list && "
					+ "freebayes -f "+myOutputDir+"/ref/"+i
					+".fa -L "+myOutputDir+"/freebayes/bam_list/"+i
					+".list --genotype-qualities -p "+myPloidy+" -v "
					+myOutputDir+"/freebayes/vcf_list/____"+i+".vcf";
		this.bash(commands);
		
		String command = "grep -v '^##contig=' "+myOutputDir+"/freebayes/vcf_list/____0.vcf > "
				+myOutputDir+"/freebayes/out.vcf";
		for(int i=1; i<THREADS; i++) 
			command += " && grep -v '^#' "+myOutputDir+"/freebayes/vcf_list/____"+i+".vcf >> "
					+myOutputDir+"/freebayes/out.vcf";
		try {
			this.bash(command).waitFor();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void generateSplitReference() {
		// TODO Auto-generated method stub
		Utils.makeOutputDir(myOutputDir+"/bed");
		Utils.makeOutputDir(myOutputDir+"/ref");
		final Map<Long, Set<String>> sizes = new HashMap<Long, Set<String>>();
		final Map<String, String> scaffs = new HashMap<String, String>();
		try {
			final BufferedReader br = Utils.getBufferedReader(this.myReference);
			String line = br.readLine();
			while( line!=null ) {
				if(line.startsWith(">")) {
					String id = line;
					long l = 0; 
					StringBuilder scaff = new StringBuilder();
					while( (line=br.readLine())!=null && 
							!line.startsWith(">") ) {
						l += line.length();
						scaff.append(line);
						scaff.append("\n");
					}
					//if(l<300) continue;
					if(!sizes.containsKey(l)) sizes.put(l, new HashSet<String>());
					sizes.get(l).add(id);
					scaffs.put(id, scaff.toString());
				}
			}
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		List<Long> longs = new ArrayList<Long>();
		for(Entry<Long, Set<String>> scaff : sizes.entrySet()) {
			int r = scaff.getValue().size();
			long key = scaff.getKey();
			for(int i=0; i<r; i++) longs.add(key);
		}
		Collections.sort(longs);
		Collections.reverse(longs);
		final List<Map<Long, Integer>> Bs = new ArrayList<Map<Long, Integer>>(); 
		for(int i=0; i<THREADS; i++) Bs.add(new HashMap<Long, Integer>());
		long[] sums = new long[THREADS];
		for(Long l : longs) {
			int min_i = 0;
			long min = sums[0];
			for(int i=1; i<THREADS; i++) {
				if(sums[i]<min) {
					min_i = i;
					min = sums[i];
				}
			}
			sums[min_i] += l;
			Map<Long, Integer> B = Bs.get(min_i);
			if(!B.containsKey(l)) B.put(l, 1);
			else B.put(l, B.get(l)+1);
		}
		try {
			final BufferedWriter[] bw_ref = new BufferedWriter[THREADS];
			final BufferedWriter[] bw_bed = new BufferedWriter[THREADS];
			for(int i=0; i<THREADS; i++) {
				bw_ref[i] = new BufferedWriter(new FileWriter(myOutputDir+"/ref/"+i+".fa"));
				bw_bed[i] = new BufferedWriter(new FileWriter(myOutputDir+"/bed/"+i+".bed"));
			}
			for(Long l : sizes.keySet()) {
				String[] scaff = new String[sizes.get(l).size()];
				sizes.get(l).toArray(scaff);
				int p = 0;
				for(int i=0; i<THREADS; i++) {
					if(!Bs.get(i).containsKey(l)) continue;
					int r = Bs.get(i).get(l);
					for(int j=0; j<r; j++) {
						bw_ref[i].write(scaff[p+j]);
						bw_ref[i].write("\n");
						bw_ref[i].write(scaffs.get(scaff[p+j]));
						
						bw_bed[i].write(scaff[p+j].substring(1));
						bw_bed[i].write("\t0\t");
						bw_bed[i].write(""+l);
						bw_bed[i].write("\n");
					}
					p += r;
				}
			}
			for(int i=0; i<THREADS; i++) {
				bw_ref[i].close();
				bw_bed[i].close();
			}	
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
