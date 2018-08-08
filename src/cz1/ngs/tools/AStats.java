package cz1.ngs.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Executor;
import cz1.util.Utils;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class AStats extends Executor {

	private int numContigsForInitialEstimate = 20;
	private int numIterations = 3;
	private int singleCopyThreshold = 30;
	private boolean bKeepDuplicates = true;
	private int minLength = 0;
	private long genomeSize = 0;
	private double arrivalRate = 0;
	private String[] bamList;
	private int num_threads = Runtime.getRuntime().availableProcessors();
	private boolean debug  = false;
	private boolean ddebug = false;
	private String out_prefix = null;

	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -a/--align              Alignment file(s). Multiple file are separated by ':'. \n"
						+ " -aL/--align-list        Alignment file list.\n"
						+ " -m/--min-len            Only compute a-stat for contigs at least the specified size.\n"
						+ " -b/--best               Number of longest contigs to perform the initial estimate \n"
						+ "                         of the arrival rate (default: 20).\n"
						+ " -n/--bootstrap          Number of bootstrap iterations to perform for the estimate.\n"
						+ " -g/--genome-size        The genome size. User this instead of estimating it.\n"
						+ " -no-dup/--no-duplicates Do not use duplicate reads to calculate statistics.\n"
						+ " -t/--threads            Number of threads to use (default 16).\n"
						+ " -d/--debug              Debugging mode will have extra information printed out.\n"
						+ " -dd/--debug-debug       Debugging mode will have more information printed out than -d mode.\n"
						+ " -o/--out-prefix         Prefix of the output files.\n"
						+ "\n");	
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-a", "--align", true);
			myArgsEngine.add("-aL", "--align-list", true);
			myArgsEngine.add("-m", "--min-len", true);
			myArgsEngine.add("-b", "--best", true);
			myArgsEngine.add("-n","--bootstrap", true);
			myArgsEngine.add("-g", "--genome-size", true);
			myArgsEngine.add("-no-dup", "--no-duplicates", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-d", "--debug", false);
			myArgsEngine.add("-dd", "--debug-debug", false);
			myArgsEngine.add("-o", "--out-prefix", true);
			myArgsEngine.parse(args);
		}

		if (!myArgsEngine.getBoolean("-a")&&!myArgsEngine.getBoolean("-aL")) {
			printUsage();
			throw new IllegalArgumentException("Please specify the alignment file(s) using -a or -aL option.");
		}
		
		if (myArgsEngine.getBoolean("-a")&&myArgsEngine.getBoolean("-aL")) {
			printUsage();
			throw new IllegalArgumentException("Options -a and -aL are exclusive.");
		}
		
		if (myArgsEngine.getBoolean("-a")) {
			this.bamList = myArgsEngine.getString("-a").trim().split(":");
		}
		
		if (myArgsEngine.getBoolean("-aL")) {
			try {
				BufferedReader br = Utils.getBufferedReader(myArgsEngine.getString("-aL").trim());
				final List<String> file_list = new ArrayList<String>();
				String line;
				while( (line = br.readLine()) != null) {
					line = line.trim();
					if(line.length()>0) file_list.add(line);
				}
				this.bamList = file_list.toArray(new String[file_list.size()]);
				br.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		}
		
		if (myArgsEngine.getBoolean("-m")) {
			this.minLength = Integer.parseInt(myArgsEngine.getString("-m"));
		}

		if (myArgsEngine.getBoolean("-b")) {
			this.numContigsForInitialEstimate = Integer.parseInt(myArgsEngine.getString("-b"));
		}

		if (myArgsEngine.getBoolean("-n")) {
			this.numIterations = Integer.parseInt(myArgsEngine.getString("-n"));
		}

		if (myArgsEngine.getBoolean("-g")) {
			this.genomeSize = Long.parseLong(myArgsEngine.getString("-g"));
		}

		if (myArgsEngine.getBoolean("-no-dup")) {
			this.bKeepDuplicates = false;
		}

		if (myArgsEngine.getBoolean("-t")) {
			int t = Integer.parseInt(myArgsEngine.getString("-t"));
			if(t<this.num_threads) this.num_threads = t;
			this.THREADS = t;
			Constants.omp_threads = this.num_threads;
			myLogger.info("OMP_THREADS = "+this.num_threads);
		}

		if (myArgsEngine.getBoolean("-d")) {
			this.debug = true;
		}

		if (myArgsEngine.getBoolean("-dd")) {
			this.debug  = true;
			this.ddebug = true;
		}

		if (myArgsEngine.getBoolean("-o")) {
			this.out_prefix = myArgsEngine.getString("-o");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the prefix of output files.");
		}
	}

	final SamReaderFactory factory =
			SamReaderFactory.makeDefault()
			.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
					SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
			.validationStringency(ValidationStringency.SILENT);

	private final static Object lock = new Object();

	private static long totalReads = 0;
	private static long sumReadLength = 0;
	private static Contig[] contigData;

	@Override
	public void run() {
		// TODO Auto-generated method stub
		
		myLogger.info("Reading alignments from "+this.bamList.length+" BAM file"+
				(this.bamList.length>1?"s":"")+":");
		for(String bamfile : this.bamList)
			myLogger.info(bamfile);
		myLogger.info("****");
		
		final SamReader in1 = factory.open(new File(this.bamList[0]));
		final List<SAMSequenceRecord> refseqs = 
				in1.getFileHeader().getSequenceDictionary().getSequences();
		contigData = new Contig[refseqs.size()];
		for(int i=0; i<refseqs.size(); i++)
			contigData[i] = new Contig(refseqs.get(i).getSequenceName(), 
					refseqs.get(i).getSequenceLength());

		this.initial_thread_pool();
		for(String bamfile : this.bamList)
			this.executor.submit(new Runnable() {
				private String bamfile;

				@Override
				public void run() {
					// TODO Auto-generated method stub
					long readCount1 = 0;
					long sumReadLength1 = 0;
					int last_ref_idx = -1;
					int last_pos = -1;
					int ref_idx, pos;
					String ref_name;
					Contig cd;

					try {
						final SamReader in1 = factory.open(new File(bamfile));
						final List<SAMSequenceRecord> refseqs = 
								in1.getFileHeader().getSequenceDictionary().getSequences();
						final Contig[] contigData1 = new Contig[refseqs.size()];
						for(int i=0; i<refseqs.size(); i++)
							contigData1[i] = new Contig(refseqs.get(i).getSequenceName(), 
									refseqs.get(i).getSequenceLength());

						final SAMRecordIterator iter1 = in1.iterator();

						myLogger.info("Reading alignments from "+ this.bamfile);

						SAMRecord alignment;
						while(iter1.hasNext()) {
							alignment = iter1.next();

							if(alignment.getReadUnmappedFlag()) 
								continue;

							ref_idx = alignment.getReferenceIndex();
							ref_name = alignment.getReferenceName();
							pos = alignment.getAlignmentStart();

							if(ref_idx == last_ref_idx && 
									pos == last_pos && 
									!bKeepDuplicates)
								continue;

							cd = contigData1[ref_idx];
							cd.addRead();
							readCount1 += 1;
							sumReadLength1 += alignment.getReadLength();
							if(!cd.name.equals(ref_name)) 
								throw new RuntimeException("!!!");
							last_ref_idx = ref_idx;
							last_pos = pos;
						}

						iter1.close();
						in1.close();

						synchronized(lock) {
							totalReads += readCount1;
							sumReadLength += sumReadLength1;
							for(int i=0; i<contigData1.length; i++) {
								contigData[i].addRead(contigData1[i].n);
							}
						}

					} catch (Exception e) {
						// TODO Auto-generated catch block
						Thread t = Thread.currentThread();
						t.getUncaughtExceptionHandler().uncaughtException(t, e);
						e.printStackTrace();
						executor.shutdown();
						System.exit(1);
					}
				} 

				public Runnable init(String bamfile) {
					// TODO Auto-generated method stub
					this.bamfile = bamfile;
					return this;
				}				
			}.init(bamfile));

		this.waitFor();

		try {
			BufferedWriter bw = Utils.getBufferedWriter(this.out_prefix+".rstat");
			for(Contig cd : contigData) 
				bw.write(cd.toString()+"\n");
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
		}

		
		double avgReadLen = (double)sumReadLength / totalReads;
		Arrays.sort(contigData, new Comparator<Contig>() {

			@Override
			public int compare(Contig arg0, Contig arg1) {
				// TODO Auto-generated method stub
				return arg1.len-arg0.len;
			}

		});

		for(Contig cd : contigData) 
			cd.nlen = (cd.len > avgReadLen ? (cd.len - avgReadLen + 1) : 0);

		long bootstrapLen = 0;
		long bootstrapReads = 0;

		if(genomeSize == 0) {
			for(int i=0; i<Math.min(numContigsForInitialEstimate, contigData.length); i++) {
				Contig cd = contigData[i];
				bootstrapLen += cd.nlen;
				bootstrapReads += cd.n;
			}

			arrivalRate = (double) bootstrapReads / (double) bootstrapLen;
			genomeSize = (long)(totalReads / arrivalRate);

			myLogger.info("Initial arrival rate: " + arrivalRate);
			myLogger.info("Initial genome size estimate: " + genomeSize);

			for(int i=0; i<numIterations; i++) {
				bootstrapLen = 0;
				bootstrapReads = 0;
				for(Contig cd : contigData) {
					cd.astat = computeAStat(arrivalRate, cd.nlen, cd.n);
					cd.bUnique = cd.astat > singleCopyThreshold;

					if(cd.len >= minLength && cd.bUnique) {
						bootstrapLen += cd.nlen;
						bootstrapReads += cd.n;
					}
				}

				arrivalRate = (double) bootstrapReads / (double) bootstrapLen;
				genomeSize = (long)(totalReads / arrivalRate);
				myLogger.info("Iteration "+i+" arrival rate: " + arrivalRate);
				myLogger.info("Iteration "+i+" genome size estimate: " + genomeSize);
			}
		}

		arrivalRate = (double)totalReads / genomeSize;
		myLogger.info("Using genome size: " + genomeSize);
		myLogger.info("Using arrival rate: " + arrivalRate);

		for(Contig cd : contigData) {
			cd.astat = computeAStat(arrivalRate, cd.nlen, cd.n);
			cd.bUnique = cd.astat > singleCopyThreshold;

			if(cd.len >= minLength && cd.bUnique) {
				bootstrapLen += cd.nlen;
				bootstrapReads += cd.n;
			}
		}

		long sumUnique = 0;
		long sumRepeat = 0;
		try {
			BufferedWriter bw = Utils.getBufferedWriter(this.out_prefix+".astat");
			for(Contig cd : contigData) {
				if(cd.len >= minLength && cd.nlen > 0) {
					bw.write(cd.name+"\t"+cd.len+"\t"+cd.nlen+"\t"+cd.n+"\t"+(cd.n / (cd.nlen * arrivalRate))+"\t"+cd.astat+"\n");
					if(cd.bUnique)
						sumUnique += cd.len;
					else
						sumRepeat += cd.len;
				}
			}
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
		}

		myLogger.info("Sum unique bases in contigs >= " +minLength+ "bp in length: "+sumUnique);
		myLogger.info("Sum repeat bases in contigs >= " +minLength+ "bp in length: "+sumRepeat);

	}

	private final static class Contig {
		private final String name;
		private final int len;
		private double nlen;
		private int n = 0;
		private double astat = 0.0;
		private boolean bUnique = true;

		public Contig(String name, int len) {
			this.name = name;
			this.len = len;
			this.nlen = len;
		}

		public void addRead() {
			++this.n;
		}

		public void addRead(int n) {
			this.n += n;
		}
		
		@Override
		public String toString() {
			return this.name+" "+this.len+" "+this.nlen+" "+this.n+" "+this.astat+" "+this.bUnique;
		}
	}

	private static double computeAStat(double arrivalRate, double len, int n) {
		return arrivalRate * len - (n * Math.log(2));
	}
}
