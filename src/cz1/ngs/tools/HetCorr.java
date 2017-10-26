package cz1.ngs.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.RandomStringUtils;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

import java.util.Set;
import java.util.TreeMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class HetCorr extends Executor {

	private String[][] paired_fq = null;
	private String[] single_fq = null;
	private String[] fastq_files = null;
	private String[] bam_list = null;
	private String[] corr_files = null;
	private String out_prefix = null;
	private String target_region = null;
	private int batch_size = 1000000;
	private static enum Task {all, hetcorr, sortq, sortc, corrp, corrs, zzz}
	private Task task_list = Task.zzz;

	private final static String bam_lib = "--b([0-9]+)";
	private final static String corr_lib = "--c([0-9]+)";
	private final static String se_lib = "--s([0-9]+)";
	private final static String pe1_lib = "--pe([0-9]+)-1";
	private final static String pe2_lib = "--pe([0-9]+)-2";

	private final static Pattern bam_pat = Pattern.compile(bam_lib);
	private final static Pattern corr_pat = Pattern.compile(corr_lib);
	private final static Pattern se_pat = Pattern.compile(se_lib);
	private final static Pattern pe1_pat = Pattern.compile(pe1_lib);
	private final static Pattern pe2_pat = Pattern.compile(pe2_lib);

	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		switch(this.task_list) {
		case zzz:
			myLogger.info(
					"\n\nChoose a task from listed:\n"
							+ " hetcorr                 Find heterozygous and sequencing errors.\n"
							+ " sortq                   Sort FASTQ files. The SEQ id of FASTQ records \n"
							+ "                         in a file should not be duplicated.\n"
							+ " sortc                   Sort CORR file generated in \'hetcorr\' step.\n"
							+ " corrs                   Correct single-end FASTQ file.\n"
							+ " corrp                   Correct paired-end FASTQ file.\n"
							+ " all                     Run the whole pipeline (default).\n"
							+ "\n");
			break;
		case hetcorr:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " --b<#>                  List of input BAM files (<#> = 1,2,...).\n"
							+ " -r/--region             Target regions in the form <chr>:<start>-<end> (default: whole genome).\n"
							+ " -t/--threads            Threads to use (default: 1).\n"
							+ " -o/--out-prefix         Prefix of the output CORR files.\n"
							+ "\n");
			break;
		case sortq:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " --s<#>                  Input single-end FASTQ file (<#> = 1,2,...).\n"
							+ " --pe<#>-1               Input  paired-end forward FASTQ file (<#> = 1,2,...).\n"
							+ " --pe<#>-2               Input  paired-end reverse FASTQ file (<#> = 1,2,...).\n"
							+ " -bs/--batch-size        Batch size store in memory (default 1000000). \n"
							+ "                         Reduce this number if run out of memory.\n"
							+ " -t/--threads            Threads to use (default: 1).\n"
							+ " -o/--out                Prefix of the output FASTQ file.\n"
							+ "\n");
			break;
		case sortc:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " --c<#>                  Input CORR file (<#> = 1,2,...).\n"
							+ " -bs/--batch-size        Batch size store in memory (default 1000000). \n"
							+ "                         Reduce this number if run out of memory.\n"
							+ " -t/--threads            Threads to use (default: 1).\n"
							+ " -o/--out                Prefix of the output CORR file.\n"
							+ "\n");
			break;
		case corrs:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " --s<#>                  Input single-end FASTQ file (<#> = 1,2,...).\n"
							+ " --c<#>                  Input CORR file.\n"
							+ " -t/--threads            Threads to use (default: 1).\n"
							+ " -o/--out                Prefix of the output FASTA file.\n"
							+ "\n");
			break;
		case corrp:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " --pe<#>-1               Input  paired-end forward FASTQ file (<#> = 1,2,...).\n"
							+ " --pe<#>-2               Input  paired-end reverse FASTQ file (<#> = 1,2,...).\n"
							+ " --c<#>                  Input CORR file.\n"
							+ " -t/--threads            Threads to use (default: 1).\n"
							+ " -o/--out                Prefix of the output FASTA file.\n"
							+ "\n");
			break;
		case all:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " --b<#>                  List of input BAM files (<#> = 1,2,...).\n"
							+ " --s<#>                  Input single-end FASTQ file (<#> = 1,2,...).\n"
							+ " --pe<#>-1               Input  paired-end forward FASTQ file (<#> = 1,2,...).\n"
							+ " --pe<#>-2               Input  paired-end reverse FASTQ file (<#> = 1,2,...).\n"
							+ " -bs/--batch-size        Batch size store in memory (default 1000000). \n"
							+ "                         Reduce this number if run out of memory.\n"
							+ " -t/--threads            Threads to use (default: 1).\n"
							+ " -o/--out-prefix         Prefix of the output files.\n"
							+ "\n");
			break;
		default:
			throw new RuntimeException("!!!");
		}
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}

		switch(args[0].toUpperCase()) {
		case "HETCORR":
			this.task_list = Task.hetcorr;
			break;
		case "SORTQ":
			this.task_list = Task.sortq;
			break;
		case "SORTC":
			this.task_list = Task.sortc;
			break;
		case "CORRS":
			this.task_list = Task.corrs;
			break;
		case "CORRP":
			this.task_list = Task.corrp;
			break;
		case "ALL":
			this.task_list = Task.all;
			break;
		default:
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");	
		}

		String[] args2 = new String[args.length-1];
		System.arraycopy(args, 1, args2, 0, args2.length);

		switch(this.task_list) {
		case hetcorr:
			if (myArgsEngine == null) {
				myArgsEngine = new ArgsEngine();
				myArgsEngine.add("-r", "--region", true);
				myArgsEngine.add("-t", "--threads", true);
				myArgsEngine.add("-o", "--out-prefix", true);
				myArgsEngine.addWildOptions(bam_lib, true);
				myArgsEngine.parse(args2);
			}

			if (myArgsEngine.getBoolean("-r")) {
				this.target_region = myArgsEngine.getString("-r");
			}

			if (myArgsEngine.getBoolean("-t")) {
				this.THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
			}

			if (myArgsEngine.getBoolean("-o")) {
				this.out_prefix = myArgsEngine.getString("-o");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the output file.");
			}

			this.parseDataLibrary(args2);

			break;
		case sortq:
			if (myArgsEngine == null) {
				myArgsEngine = new ArgsEngine();
				myArgsEngine.add("-bs", "--batch-size", true);
				myArgsEngine.add("-t", "--threads", true);
				myArgsEngine.add("-o", "--out", true);
				myArgsEngine.addWildOptions(se_lib, true);
				myArgsEngine.addWildOptions(pe1_lib, true);
				myArgsEngine.addWildOptions(pe2_lib, true);
				myArgsEngine.parse(args2);
			}

			if (myArgsEngine.getBoolean("-bs")) {
				this.batch_size = Integer.parseInt(myArgsEngine.getString("-bs"));
			}

			if (myArgsEngine.getBoolean("-t")) {
				this.THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
			}

			if (myArgsEngine.getBoolean("-o")) {
				this.out_prefix = myArgsEngine.getString("-o");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the output file.");
			}

			this.parseDataLibrary(args2);

			break;
		case sortc:
			if (myArgsEngine == null) {
				myArgsEngine = new ArgsEngine();
				myArgsEngine.add("-bs", "--batch-size", true);
				myArgsEngine.add("-t", "--threads", true);
				myArgsEngine.add("-o", "--out", true);
				myArgsEngine.addWildOptions(corr_lib, true);
				myArgsEngine.parse(args2);
			}

			if (myArgsEngine.getBoolean("-bs")) {
				this.batch_size = Integer.parseInt(myArgsEngine.getString("-bs"));
			}

			if (myArgsEngine.getBoolean("-t")) {
				this.THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
			}

			if (myArgsEngine.getBoolean("-o")) {
				this.out_prefix = myArgsEngine.getString("-o");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the output file.");
			}

			this.parseDataLibrary(args2);

			break;
		case corrs:
			if (myArgsEngine == null) {
				myArgsEngine = new ArgsEngine();
				myArgsEngine.add("-t", "--threads", true);
				myArgsEngine.add("-o", "--out", true);
				myArgsEngine.addWildOptions(se_lib, true);
				myArgsEngine.addWildOptions(corr_lib, true);
				myArgsEngine.parse(args2);
			}

			if (myArgsEngine.getBoolean("-t")) {
				this.THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
			}

			if (myArgsEngine.getBoolean("-o")) {
				this.out_prefix = myArgsEngine.getString("-o");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the output file.");
			}

			this.parseDataLibrary(args2);

			break;
		case corrp:
			if (myArgsEngine == null) {
				myArgsEngine = new ArgsEngine();
				myArgsEngine.add("-t", "--threads", true);
				myArgsEngine.add("-o", "--out", true);
				myArgsEngine.addWildOptions(pe1_lib, true);
				myArgsEngine.addWildOptions(pe2_lib, true);
				myArgsEngine.addWildOptions(corr_lib, true);
				myArgsEngine.parse(args2);
			}

			if (myArgsEngine.getBoolean("-t")) {
				this.THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
			}

			if (myArgsEngine.getBoolean("-o")) {
				this.out_prefix = myArgsEngine.getString("-o");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the output file.");
			}

			this.parseDataLibrary(args2);

			break;
		case all:
			if (myArgsEngine == null) {
				myArgsEngine = new ArgsEngine();
				myArgsEngine.add("-bs", "--batch-size", true);
				myArgsEngine.add("-t", "--threads", true);
				myArgsEngine.add("-o", "--out-prefix", true);
				myArgsEngine.addWildOptions(se_lib, true);
				myArgsEngine.addWildOptions(pe1_lib, true);
				myArgsEngine.addWildOptions(pe2_lib, true);
				myArgsEngine.addWildOptions(bam_lib, true);
				myArgsEngine.parse(args2);
			}

			if (myArgsEngine.getBoolean("-bs")) {
				this.batch_size = Integer.parseInt(myArgsEngine.getString("-bs"));
			}

			if (myArgsEngine.getBoolean("-t")) {
				this.THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
			}

			if (myArgsEngine.getBoolean("-o")) {
				this.out_prefix = myArgsEngine.getString("-o");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the output file.");
			}

			this.parseDataLibrary(args2);

			break;
		default:
			throw new RuntimeException("!!!");	
		}

	}

	private void parseDataLibrary(String[] args) {
		// TODO Auto-generated method stub

		final Map<Integer, String> bamLib = new HashMap<Integer, String>();
		final Map<Integer, String> corrLib = new HashMap<Integer, String>();
		final Map<Integer, String> seLib = new HashMap<Integer, String>();
		final Map<Integer, String> pe1Lib = new HashMap<Integer, String>();
		final Map<Integer, String> pe2Lib = new HashMap<Integer, String>();

		for(int i=0; i<args.length; i++) {
			String arg = args[i];
			if(arg.matches(bam_lib)) {
				Matcher m = bam_pat.matcher(arg);
				m.find();
				int lib = Integer.parseInt(m.group(1));
				bamLib.put(lib, args[++i]);
			}
			if(arg.matches(corr_lib)) {
				Matcher m = corr_pat.matcher(arg);
				m.find();
				int lib = Integer.parseInt(m.group(1));
				corrLib.put(lib, args[++i]);
			}
			if(arg.matches(se_lib)) {
				Matcher m = se_pat.matcher(arg);
				m.find();
				int lib = Integer.parseInt(m.group(1));
				seLib.put(lib, args[++i]);
			}
			if(arg.matches(pe1_lib)) {
				Matcher m = pe1_pat.matcher(arg);
				m.find();
				int lib = Integer.parseInt(m.group(1));
				pe1Lib.put(lib, args[++i]);
			}
			if(arg.matches(pe2_lib)) {
				Matcher m = pe2_pat.matcher(arg);
				m.find();
				int lib = Integer.parseInt(m.group(1));
				pe2Lib.put(lib, args[++i]);
			}
		}


		int i;
		switch(this.task_list) {
		case hetcorr:
			if(bamLib.isEmpty()) {
				this.printUsage();
				throw new RuntimeException("Please specify the BAM file.");
			}
			i = 0;
			this.bam_list = new String[bamLib.size()];
			for(Map.Entry<Integer, String> entry : bamLib.entrySet()) 
				bam_list[i++] = entry.getValue();
			if(corrLib.size()>0) myLogger.warn("CORR files ignored!!!");
			if(seLib.size()>0) myLogger.warn("Single-end library ignored!!!");
			if(pe1Lib.size()>0||pe2Lib.size()>0) myLogger.warn("Paired-end library ignored!!!");
			break;
		case sortq:
			if(seLib.isEmpty()&&pe1Lib.isEmpty()&&pe2Lib.isEmpty()) {
				this.printUsage();
				throw new RuntimeException("Please specify the FASTQ file.");
			}
			i = 0;
			this.fastq_files = new String[seLib.size()+pe1Lib.size()+pe2Lib.size()];
			for(Map.Entry<Integer, String> entry : seLib.entrySet()) 
				fastq_files[i++] = entry.getValue();
			for(Map.Entry<Integer, String> entry : pe1Lib.entrySet()) 
				fastq_files[i++] = entry.getValue();
			for(Map.Entry<Integer, String> entry : pe2Lib.entrySet()) 
				fastq_files[i++] = entry.getValue();
			if(bamLib.size()>0) myLogger.warn("BAM files ignored!!!");
			if(corrLib.size()>0) myLogger.warn("CORR files ignored!!!");
			break;
		case sortc:
			if(corrLib.isEmpty()) {
				this.printUsage();
				throw new RuntimeException("Please specify the CORR file.");
			}
			i = 0;
			this.corr_files = new String[corrLib.size()];
			for(Map.Entry<Integer, String> entry : corrLib.entrySet()) 
				corr_files[i++] = entry.getValue();
			if(bamLib.size()>0) myLogger.warn("BAM files ignored!!!");
			if(seLib.size()>0) myLogger.warn("Single-end library ignored!!!");
			if(pe1Lib.size()>0||pe2Lib.size()>0) myLogger.warn("Paired-end library ignored!!!");
			break;
		case corrs:
			if(seLib.isEmpty()) {
				this.printUsage();
				throw new RuntimeException("Please specify the single-end FASTQ file.");
			}
			if(corrLib.isEmpty()) {
				this.printUsage();
				throw new RuntimeException("Please specify the CORR file.");
			}
			i = 0;
			this.single_fq = new String[seLib.size()];
			this.corr_files = new String[corrLib.size()];
			if(single_fq.length!=corr_files.length) {
				this.printUsage();
				throw new RuntimeException("Single-end library and CORR file does not match!!!");
			}
			for(Integer key : seLib.keySet()) {
				if(!corrLib.containsKey(key)) {
					this.printUsage();
					throw new RuntimeException("Single-end library and CORR file does not match!!!");
				}
				single_fq[i] = seLib.get(key);
				corr_files[i] = corrLib.get(key);
				++i;
			}
			if(bamLib.size()>0) myLogger.warn("BAM files ignored!!!");
			if(pe1Lib.size()>0||pe2Lib.size()>0) myLogger.warn("Paired-end library ignored!!!");
			break;
		case corrp:
			if(pe1Lib.isEmpty()) {
				this.printUsage();
				throw new RuntimeException("Please specify the paired-end FASTQ file.");
			}
			if(corrLib.isEmpty()) {
				this.printUsage();
				throw new RuntimeException("Please specify the CORR file.");
			}
			i = 0;
			this.paired_fq = new String[pe1Lib.size()][2];
			this.corr_files = new String[corrLib.size()];
			if(paired_fq.length!=corr_files.length || paired_fq.length!=pe2Lib.size()) {
				this.printUsage();
				throw new RuntimeException("Paired-end library and CORR file does not match!!!");
			}
			for(Integer key : pe1Lib.keySet()) {
				if(!corrLib.containsKey(key) || !pe2Lib.containsKey(key)) {
					this.printUsage();
					throw new RuntimeException("Paired-end library and CORR file does not match!!!");
				}
				paired_fq[i][0] = pe1Lib.get(key);
				paired_fq[i][1] = pe2Lib.get(key);
				corr_files[i] = corrLib.get(key);
				++i;
			}
			if(bamLib.size()>0) myLogger.warn("BAM files ignored!!!");
			if(seLib.size()>0) myLogger.warn("Single-end library ignored!!!");
			break;
		case all:
			if(seLib.isEmpty()&&pe1Lib.isEmpty()&&pe2Lib.isEmpty()) {
				this.printUsage();
				throw new IllegalArgumentException("Please specify at least one FASTQ library.");
			}
			if(bamLib.isEmpty()) {
				this.printUsage();
				throw new IllegalArgumentException("Please specify the BAM file.");
			}
			i = 0;
			this.paired_fq = new String[bamLib.size()][2];
			this.single_fq = new String[bamLib.size()];
			this.bam_list = new String[bamLib.size()];

			if(bamLib.size()!=(pe1Lib.size()+seLib.size())) {
				this.printUsage();
				throw new RuntimeException("FASTQ library and CORR file does not match!!!");
			}
			if(pe1Lib.size()!=pe2Lib.size()) {
				this.printUsage();
				throw new RuntimeException("Paired-end file does not match!!!");
			}
			for(Integer key : bamLib.keySet()) {
				bam_list[i] = bamLib.get(key);
				if(seLib.containsKey(key)) {
					single_fq[i] = seLib.get(key);
				} else if(pe1Lib.containsKey(key)) {
					if(!pe2Lib.containsKey(key)) {
						this.printUsage();
						throw new RuntimeException("Paired-end file does not match!!!");
					}
					paired_fq[i][0] = pe1Lib.get(key);
					paired_fq[i][1] = pe2Lib.get(key);
				} else {
					this.printUsage();
					throw new RuntimeException("FASTQ library and CORR file does not match!!!");
				}
				++i;
			}
			if(corrLib.size()>0) myLogger.warn("CORR files ignored!!!");
			break;
		default:
			this.printUsage();
			throw new RuntimeException("!!!");
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub

		String queryRefSeq = null;
		int queryStart = -1, queryEnd = -1, n;
		SAMSequenceDictionary refSeqDict;
		List<SAMSequenceRecord> refSeq;
		int[][] chrSpan, jobs;
		String[] joined_outPref;
		
		switch(this.task_list) {
		
		case zzz:
			myLogger.info("Task list is empty!!!");
			break;
		case hetcorr:
			this.initial_thread_pool();
			if(target_region!=null) {
				String[] s = target_region.split(":");
				queryRefSeq = s[0].trim();
				s = s[1].split("-");
				queryStart = Integer.parseInt(s[0]);
				queryEnd = Integer.parseInt(s[1]);
			}
			refSeqDict = factory.open(new File(this.bam_list[0])).
					getFileHeader().
					getSequenceDictionary();
			refSeq = refSeqDict.getSequences();
			n = target_region==null ? refSeq.size() : 1;
			chrSpan = new int[n][3];
			if(target_region!=null) {
				chrSpan[0] = new int[]{refSeqDict.getSequenceIndex(queryRefSeq), 
						queryStart, queryEnd};
			} else {
				for(int i=0; i!=n; i++) {
					chrSpan[i] = new int[]{i, 0, refSeq.get(i).getSequenceLength()};
				}
			}

			jobs = scheduleJobs(chrSpan);

			joined_outPref = new String[jobs.length];
			for(int i=0; i<jobs.length; i++) {
				String out = this.out_prefix+"_"+
						jobs[i][0]+"_"+
						jobs[i][1]+"_"+
						jobs[i][2]+"_"+
						RandomStringUtils.randomAlphanumeric(20).toUpperCase();
				joined_outPref[i] = out;
				this.executor.submit(new HetCorrecter(jobs[i][0], 
						jobs[i][1], 
						jobs[i][2], 
						refSeq,
						out));
			}
			this.waitFor();

			try {
				for(int i=0; i<this.bam_list.length; i++) {
					BufferedWriter corrOut = Utils.getBufferedWriter(this.out_prefix+
							new File(bam_list[i]).getName()+".corr");
					for(int j=0; j<jobs.length; j++) {
						BufferedReader corrIn = Utils.getBufferedReader(joined_outPref[j]+"_"+
								new File(bam_list[i]).getName()+".corr");
						IOUtils.copy(corrIn, corrOut);
						corrIn.close();
						new File(joined_outPref[j]+"_"+new File(bam_list[i]).getName()+".corr").delete();
					}
					corrOut.close();
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			break;
		case sortq:
			this.initial_thread_pool();
			for(int i=0; i<this.fastq_files.length; i++) {
				String out = this.out_prefix+(new File(this.fastq_files[i]).getName()).
						replaceAll(".fq.gz$", "").
						replaceAll(".fastq.gz$", "").
						replaceAll(".fastq", "").
						replaceAll(".fq$", "")+"_sorted.fastq.gz";
				if(new File(out).exists())
					throw new RuntimeException("Cannot write File!!! Output file "+out+"existed!!!");
				myLogger.info("Writing to "+out);
				this.executor.submit(new FastqFileSorter(this.fastq_files[i], out));
			}
			this.waitFor();
			break;
		case sortc:
			this.initial_thread_pool();
			for(int i=0; i<this.corr_files.length; i++) {
				String out = this.out_prefix+(new File(this.corr_files[i]).getName()).
						replaceAll(".corr$", "")+"_sorted.corr";
				if(new File(out).exists())
					throw new RuntimeException("Cannot write File!!! Output file "+out+"existed!!!");
				myLogger.info("Writing to "+out);
				this.executor.submit(new CorrFileSorter(this.corr_files[i], out));
			}
			this.waitFor();
			break;
		case corrs:
			this.initial_thread_pool();
			for(int i=0; i<this.single_fq.length; i++) {
				String out = this.out_prefix+(new File(this.single_fq[i]).getName()).
						replaceAll(".fq.gz$", "").
						replaceAll(".fastq.gz$", "").
						replaceAll(".fastq", "").
						replaceAll(".fq$", "")+"_corrected.fa.gz";
				if(new File(out).exists())
					throw new RuntimeException("Cannot write File!!! Output file "+out+"existed!!!");
				myLogger.info("Writing to "+out);
				this.executor.submit(new SECorrecter(this.single_fq[i], this.corr_files[i], out));
			}
			this.waitFor();
			break;
		case corrp:
			this.initial_thread_pool();
			for(int i=0; i<this.paired_fq.length; i++) {
				String q1 = this.paired_fq[i][0], q2 = this.paired_fq[i][1];
				int j = 0;
				while(q1.charAt(j)==q2.charAt(j)) j++;
				String out = this.out_prefix+new File(q1.substring(0, j)).getName();
				if(new File(out).exists())
					throw new RuntimeException("Cannot write File!!! Output file "+out+"existed!!!");
				myLogger.info("Writing to "+out);
				this.executor.submit(new PECorrecter(this.paired_fq[i], this.corr_files[i], out));
			}
			this.waitFor();
			break;
		case all:
			
			// sort FASTQ files
			this.initial_thread_pool();
			for(int i=0; i<this.single_fq.length; i++) {
				if(this.single_fq[i]==null) continue;
				String out = this.out_prefix+(new File(this.single_fq[i]).getName()).
						replaceAll(".fq.gz$", "").
						replaceAll(".fastq.gz$", "").
						replaceAll(".fastq", "").
						replaceAll(".fq$", "")+"_sorted.fastq.gz";
				if(new File(out).exists())
					throw new RuntimeException("Cannot write File!!! Output file "+out+"existed!!!");
				myLogger.info("Writing to "+out);
				this.executor.submit(new FastqFileSorter(this.single_fq[i], out));
			}
			for(int i=0; i<this.paired_fq.length; i++) {
				for(int j=0; j<2; j++) {
					if(this.paired_fq[i][j]==null) continue;
					String out = this.out_prefix+(new File(this.paired_fq[i][j]).getName()).
							replaceAll(".fq.gz$", "").
							replaceAll(".fastq.gz$", "").
							replaceAll(".fastq", "").
							replaceAll(".fq$", "")+"_sorted.fastq.gz";
					if(new File(out).exists())
						throw new RuntimeException("Cannot write File!!! Output file "+out+"existed!!!");
					myLogger.info("Writing to "+out);
					this.executor.submit(new FastqFileSorter(this.paired_fq[i][j], out));
				}
			}
			this.waitFor();
			
			
			for(int i=0; i<this.single_fq.length; i++) {
				if(this.single_fq[i]==null) continue;
				String out = this.out_prefix+(new File(this.single_fq[i]).getName()).
						replaceAll(".fq.gz$", "").
						replaceAll(".fastq.gz$", "").
						replaceAll(".fastq", "").
						replaceAll(".fq$", "")+"_sorted.fastq.gz";
				this.single_fq[i] = out;
			}
			for(int i=0; i<this.paired_fq.length; i++) {
				for(int j=0; j<2; j++) {
					if(this.paired_fq[i][j] == null) continue;
					String out = this.out_prefix+(new File(this.paired_fq[i][j]).getName()).
							replaceAll(".fq.gz$", "").
							replaceAll(".fastq.gz$", "").
							replaceAll(".fastq", "").
							replaceAll(".fq$", "")+"_sorted.fastq.gz";
					this.paired_fq[i][j] = out;
				}
			}
			
			// hetcorr
			this.initial_thread_pool();
			if(target_region!=null) {
				String[] s = target_region.split(":");
				queryRefSeq = s[0].trim();
				s = s[1].split("-");
				queryStart = Integer.parseInt(s[0]);
				queryEnd = Integer.parseInt(s[1]);
			}
			refSeqDict = factory.open(new File(this.bam_list[0])).
					getFileHeader().
					getSequenceDictionary();
			refSeq = refSeqDict.getSequences();
			n = target_region==null ? refSeq.size() : 1;
			chrSpan = new int[n][3];
			if(target_region!=null) {
				chrSpan[0] = new int[]{refSeqDict.getSequenceIndex(queryRefSeq), 
						queryStart, queryEnd};
			} else {
				for(int i=0; i!=n; i++) {
					chrSpan[i] = new int[]{i, 0, refSeq.get(i).getSequenceLength()};
				}
			}

			jobs = scheduleJobs(chrSpan);

			joined_outPref = new String[jobs.length];
			for(int i=0; i<jobs.length; i++) {
				String out = this.out_prefix+"_"+
						jobs[i][0]+"_"+
						jobs[i][1]+"_"+
						jobs[i][2]+"_"+
						RandomStringUtils.randomAlphanumeric(20).toUpperCase();
				joined_outPref[i] = out;
				this.executor.submit(new HetCorrecter(jobs[i][0], 
						jobs[i][1], 
						jobs[i][2], 
						refSeq,
						out));
			}
			this.waitFor();
			
			this.corr_files = new String[this.bam_list.length];
			try {
				for(int i=0; i<this.bam_list.length; i++) {
					corr_files[i] = this.out_prefix+new File(bam_list[i]).getName()+".corr";
					BufferedWriter corrOut = Utils.getBufferedWriter(corr_files[i]);
					for(int j=0; j<jobs.length; j++) {
						BufferedReader corrIn = Utils.getBufferedReader(joined_outPref[j]+"_"+
								new File(bam_list[i]).getName()+".corr");
						IOUtils.copy(corrIn, corrOut);
						corrIn.close();
						new File(joined_outPref[j]+"_"+new File(bam_list[i]).getName()+".corr").delete();
					}
					corrOut.close();
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			// sort corr files
			this.initial_thread_pool();
			for(int i=0; i<this.corr_files.length; i++) {
				String out = this.out_prefix+(new File(this.corr_files[i]).getName()).
						replaceAll(".corr$", "")+"_sorted.corr";
				if(new File(out).exists())
					throw new RuntimeException("Cannot write File!!! Output file "+out+"existed!!!");
				myLogger.info("Writing to "+out);
				this.executor.submit(new CorrFileSorter(this.corr_files[i], out));
			}
			this.waitFor();
			for(int i=0; i<this.corr_files.length; i++) {
				String out = this.out_prefix+(new File(this.corr_files[i]).getName()).
						replaceAll(".corr$", "")+"_sorted.corr";
				corr_files[i] = out;
			}
			
			// generate file corr file
			this.initial_thread_pool();
			for(int i=0; i<this.single_fq.length; i++) {
				if(this.single_fq[i]==null) continue;
				String out = this.out_prefix+(new File(this.single_fq[i]).getName()).
						replaceAll(".fq.gz$", "").
						replaceAll(".fastq.gz$", "").
						replaceAll(".fastq", "").
						replaceAll(".fq$", "")+"_corrected.fa.gz";
				if(new File(out).exists())
					throw new RuntimeException("Cannot write File!!! Output file "+out+"existed!!!");
				myLogger.info("Writing to "+out);
				this.executor.submit(new SECorrecter(this.single_fq[i], this.corr_files[i], out));
			}
			
			for(int i=0; i<this.paired_fq.length; i++) {
				if(this.paired_fq[i][0]==null) continue;
				String q1 = this.paired_fq[i][0], q2 = this.paired_fq[i][1];
				int j = 0;
				while(q1.charAt(j)==q2.charAt(j)) j++;
				String out = this.out_prefix+new File(q1.substring(0, j)).getName();
				if(new File(out).exists())
					throw new RuntimeException("Cannot write File!!! Output file "+out+"existed!!!");
				myLogger.info("Writing to "+out);
				this.executor.submit(new PECorrecter(this.paired_fq[i], this.corr_files[i], out));
			}
			this.waitFor();
			break;
		default:
			throw new RuntimeException("!!!");
		}
		return;
	}

	private final static int job_batch_size = 1000000;
	private int[][] scheduleJobs(int[][] chrSpan) {
		// TODO Auto-generated method stub
		List<int[]> jobs = new ArrayList<int[]>();
		for(int i=0; i<chrSpan.length; i++) {
			int[] span = chrSpan[i];
			int j=span[1];
			while( (j+=job_batch_size)<span[2] ) {
				jobs.add(new int[]{span[0], j-job_batch_size, j});
			}
			if(jobs.isEmpty() || jobs.get(jobs.size()-1)[2]<span[2])
				jobs.add(new int[]{span[0], j-job_batch_size, span[2]});
		}
		int[][] jobsAll = new int[jobs.size()][3];
		for(int i=0; i<jobs.size(); i++)
			jobsAll[i] = jobs.get(i);
		return jobsAll;
	}

	private final static SamReaderFactory factory =
			SamReaderFactory.makeDefault()
			.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
					SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
			.validationStringency(ValidationStringency.SILENT);
	private final static char[] nucleotide = new char[]{'A','C','G','T'};
	private final static Map<Character, Character> revcmp = new HashMap<Character, Character>();
	static {
		revcmp.put('A', 'T');
		revcmp.put('C', 'G');
		revcmp.put('G', 'C');
		revcmp.put('T', 'A');
		revcmp.put('N', 'N');
	}
	private final static int max_seq_length = 300; 
	private final static boolean check = false;
	private final static int[] check_stats = new int[10];

	private static void checkErrorCorr(List<CorrRecord> buff) {
		// TODO Auto-generated method stub
		// check substitution and deletion
		boolean[] corr = new boolean[max_seq_length];
		int overhang, oriS;
		CorrOperator operator;
		for(CorrRecord corrRecord : buff) {
			overhang = corrRecord.seq_pos;
			oriS = corrRecord.seq_ori.length();
			operator = corrRecord.operator;
			if(operator==CorrOperator.S || operator==CorrOperator.D) {
				for(int i=0; i!=oriS; i++) {
					if(corr[i+overhang]) 
						throw new RuntimeException("!!!");
					corr[i+overhang] = true;
				}
			}
		}

		// check insertion
		CorrRecord corrRecord = buff.get(0), corrRecord2;
		for(int i=1; i!=buff.size(); i++) {
			corrRecord2 = buff.get(i);
			if( corrRecord.operator==CorrOperator.I &&
					corrRecord2.operator==CorrOperator.I &&
					corrRecord.seq_pos == corrRecord2.seq_pos && 
					Integer.parseInt(corrRecord.seq_ori)==
					Integer.parseInt(corrRecord2.seq_ori) )
				throw new RuntimeException("!!!");
			corrRecord = corrRecord2;
		}
		return;
	}

	private static void bufferCorrRecord(String corr_line,
			List<List<CorrRecord>> buffered_corr) {
		// TODO Auto-generated method stub
		String[] s = corr_line.split(",");
		buffered_corr.get(Integer.parseInt(s[1])).add(
				new CorrRecord(s[0],Integer.parseInt(s[2]),s[3],s[4],s[5]));
	}

	private static String insertion_seq(final SAMRecord record, 
			final int pos, 
			final int[] ht_hc_ins,
			final int[] rl_shift) {
		// TODO Auto-generated method stub
		String dna_seq = record.getReadString();
		int rs = record.getAlignmentStart(), re = record.getAlignmentEnd();
		int a = record.getReadPositionAtReferencePosition(pos-1);
		int shift = 0;
		while(a==0&&pos-shift>rs) {
			shift++;
			a = record.getReadPositionAtReferencePosition(pos-shift); 
		}
		if(rl_shift!=null) rl_shift[0] = shift;
		int b = record.getReadPositionAtReferencePosition(pos);
		shift = 0;
		while(b==0&&pos+shift<re) {
			shift++;
			b = record.getReadPositionAtReferencePosition(pos+shift); 
		}
		if(rl_shift!=null) rl_shift[1] = shift;
		if(ht_hc_ins!=null) {
			ht_hc_ins[0] = a;
			ht_hc_ins[1] = b-1;
		}
		return a==0||b>re?null:dna_seq.substring(a, b-1);
	}

	private static String revcmp_seq(String ins_seq) {
		// TODO Auto-generated method stub
		StringBuilder rev = new StringBuilder();
		for(char b : ins_seq.toCharArray())
			rev.append(revcmp.get(b));
		return rev.reverse().toString();
	}

	private static int seqLength(final SAMRecord record, final int[] ht_hc) {
		// TODO Auto-generated method stub
		CigarElement co;
		Arrays.fill(ht_hc, 0);
		co = record.getCigar().getFirstCigarElement();
		if(co.getOperator()==CigarOperator.H) ht_hc[0] = co.getLength();
		co = record.getCigar().getLastCigarElement();
		if(co.getOperator()==CigarOperator.H) ht_hc[1] = co.getLength();
		return record.getReadLength()+ht_hc[0]+ht_hc[1];
	}

	private static Entry<String, Integer> max_index(final Map<String, Integer> ins_stats) {
		// TODO Auto-generated method stub
		int mx = Integer.MIN_VALUE;
		Entry<String, Integer> mi = null;
		for(Map.Entry<String, Integer> entry : ins_stats.entrySet()) {
			if(entry.getValue()>mx) {
				mx = entry.getValue();
				mi = entry;
			}
		}
		return mi;
	}

	private static int max_index(int[] ints) {
		// TODO Auto-generated method stub
		int mx = Integer.MIN_VALUE;
		int mi = -1;
		for(int i=0; i!=ints.length; i++) {
			if(ints[i]>mx) {
				mx = ints[i];
				mi = i;
			}
		}
		return mi;
	}

	private static void buffer(
			final Set<SAMRecord> record_pool,
			final Set<Integer> ins_pool,
			final SAMRecordIterator[] iter1,
			final SAMRecord[] buffer, 
			final int[] bufferS, 
			final int pos, 
			final int i) {
		// TODO Auto-generated method stub
		SAMRecord tmp_record = buffer[i];
		while( tmp_record!=null && tmp_record.getAlignmentStart()<=pos) {
			// if(tmp_record.getNotPrimaryAlignmentFlag() || 
			//		tmp_record.getSupplementaryAlignmentFlag()&&tmp_record.getMappingQuality()<10 ) {
			//	System.out.println(tmp_record.getSAMString());
			// } else {
			record_pool.add(tmp_record);
			// to decide insertion position which cannot be decided from reference sequence
			if(tmp_record.getCigarString().contains("I")) {
				List<CigarElement> cigars = tmp_record.getCigar().getCigarElements();
				int overhang = tmp_record.getAlignmentStart();
				for(CigarElement cigar : cigars) {
					CigarOperator operator = cigar.getOperator();
					if(operator==CigarOperator.I)
						ins_pool.add(overhang);
					else if(operator!=CigarOperator.D &&
							operator!=CigarOperator.S &&
							operator!=CigarOperator.H)
						overhang += cigar.getLength();
				}
			}
			// }
			tmp_record=iter1[i].hasNext()?iter1[i].next():null;
		}
		buffer[i] = tmp_record;
		bufferS[i] = tmp_record==null?-1:tmp_record.getAlignmentStart();
		return;
	}

	private static void sortOvlFile(String ovlIn, String ovlOut) {
		try {

			String[] ovlIns = ovlIn.trim().replaceAll(";$", "").split(";");
			int nBatch = ovlIns.length;
			BufferedWriter bwOvlOut = new BufferedWriter(new FileWriter(ovlOut));
			BufferedReader[] brOvlOutTmp = new BufferedReader[nBatch];
			final boolean[] reachFileEnd = new boolean[nBatch];
			for(int i=0; i!=nBatch; i++) { 
				brOvlOutTmp[i] = new BufferedReader(new FileReader(ovlIns[i]));
				reachFileEnd[i] = false;
			}
			System.out.println(brOvlOutTmp.length);
			final TreeMap<String, Integer> treeMap = new TreeMap<String, Integer>();
			String key=null;
			for(int i=0; i!=nBatch; i++) {
				key = getKey(brOvlOutTmp[i].readLine());

				while(key!=null&&treeMap.containsKey(key)) 
					key = getKey(brOvlOutTmp[i].readLine());
				if(key!=null) treeMap.put(key, i);
			}
			Entry<String, Integer> firstEntry;
			int bch, nReachFileEnd=0;
			final int zeroPadLen = zeroPad.length()*2;
			while( !treeMap.isEmpty() ) {
				firstEntry = treeMap.pollFirstEntry();
				bwOvlOut.write(firstEntry.getKey().substring(zeroPadLen)+"\n");
				bch = firstEntry.getValue();
				if(!reachFileEnd[bch]) {
					key = getKey(brOvlOutTmp[bch].readLine());
					if(key==null) {
						reachFileEnd[bch] = true;
						nReachFileEnd++;
					} else {
						while(key!=null&&treeMap.containsKey(key)) 
							key = getKey(brOvlOutTmp[bch].readLine());
						if(key!=null) treeMap.put(key, bch);
					}
				}
				if(treeMap.isEmpty()&&nReachFileEnd!=nBatch) {
					for(int i=0; i!=nBatch; i++) 
						if(!reachFileEnd[i]) {
							key = getKey(brOvlOutTmp[i].readLine());
							while(key!=null&&treeMap.containsKey(key)) 
								key = getKey(brOvlOutTmp[i].readLine());
							if(key!=null) treeMap.put(key, i);
						}
				}
			}
			for(int i=0; i!=nBatch; i++) {
				brOvlOutTmp[i].close();
			}
			bwOvlOut.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private final static String zeroPad = "0000000";

	private static String getKey(String line) {
		// TODO Auto-generated method stub
		String key = null;
		if(line!=null) {
			String[] s = line.split("\\s+");
			key = zeroPad.substring(s[0].length())+s[0]+ 
					zeroPad.substring(s[1].length())+s[1]+
					line; 
		}
		return key;
	}

	private static void checkDuplicateRead(List<FastqRecord> batch) {
		// TODO Auto-generated method stub
		String seq_id = batch.get(0).seq_id, seq_id2;
		for(int i=1; i!=batch.size(); i++) {
			seq_id2 = batch.get(i).seq_id;
			if(seq_id.equals(seq_id2)) 
				throw new RuntimeException(seq_id+"!!!\n"
						+ "Duplicated FASTQ record id which is not allowed!!! "
						+ "Please rename your FASTQ file and retry.");
			seq_id = seq_id2;
		}
		return;
	}

	private static enum CorrOperator {S, I, D}

	private final static class CorrRecord implements Comparator<CorrRecord>, Comparable<CorrRecord> {
		private final String seq_id;
		private final int seq_pos;
		private final String seq_ori;
		private final String seq_sub;
		private final CorrOperator operator;

		public CorrRecord(String seq_id,
				int seq_pos,
				String seq_ori,
				String seq_sub,
				String opt) {
			this.seq_id = seq_id;
			this.seq_pos = seq_pos;
			this.seq_ori = seq_ori;
			this.seq_sub = seq_sub;
			switch(opt) {
			case "S":
				this.operator = CorrOperator.S;
				break;
			case "I":
				this.operator = CorrOperator.I;
				break;
			case "D":
				this.operator = CorrOperator.D;
				break;
			default:
				throw new RuntimeException("!!!");
			}
		}

		@Override
		public int compareTo(CorrRecord cr) {
			// TODO Auto-generated method stub
			if(this.seq_pos!=cr.seq_pos) return this.seq_pos-cr.seq_pos;
			final boolean[] is_ins = new boolean[2];
			is_ins[0] = this.operator==CorrOperator.I;
			is_ins[1] = cr.operator==CorrOperator.I;
			if(is_ins[0]&&is_ins[1])
				return Integer.parseInt(this.seq_ori)-Integer.parseInt(cr.seq_ori);
			if(is_ins[0]&&!is_ins[1]) return -1;
			if(!is_ins[0]&&is_ins[1]) return 1;
			//throw new RuntimeException("!!!");
			check_stats[7]++;;
			return 0;
		}

		@Override
		public int compare(CorrRecord cr, CorrRecord cr2) {
			// TODO Auto-generated method stub
			if(cr.seq_pos!=cr2.seq_pos) return cr.seq_pos-cr2.seq_pos;
			final boolean[] is_ins = new boolean[2];
			is_ins[0] = cr.operator==CorrOperator.I;
			is_ins[1] = cr2.operator==CorrOperator.I;
			if(is_ins[0]&&is_ins[1])
				return Integer.parseInt(cr.seq_ori)-Integer.parseInt(cr2.seq_ori);
			if(is_ins[0]&&!is_ins[1]) return -1;
			if(!is_ins[0]&&is_ins[1]) return 1;
			//throw new RuntimeException("!!!");
			check_stats[8]++;
			return 0;
		}
	}

	private final static class FastqRecord implements Comparator<FastqRecord>, Comparable<FastqRecord> {
		private final String seq_id;
		private final String seq_str;
		private final String seq_qual;

		public FastqRecord(String seq_id,
				String seq_str,
				String seq_qual) {
			this.seq_id = seq_id;
			this.seq_str = seq_str;
			this.seq_qual = seq_qual;
		}

		@Override
		public int compare(FastqRecord r0, FastqRecord r1) {
			// TODO Auto-generated method stub
			return r0.seq_id.compareTo(r1.seq_id);
		}

		@Override
		public int compareTo(FastqRecord r1) {
			// TODO Auto-generated method stub
			return this.seq_id.compareTo(r1.seq_id);
		}
	}

	private final class FastqFileSorter implements Runnable {

		private final String fastq_in;
		private final String fastq_out;

		public FastqFileSorter(String fastqIn, String fastqOut) {
			// TODO Auto-generated constructor stub
			this.fastq_in = fastqIn;
			this.fastq_out = fastqOut;
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub
			try {
				BufferedReader brCorrIn = Utils.getBufferedReader(fastq_in);
				final List<FastqRecord> batch = new ArrayList<FastqRecord>();
				int nBatch = 0;
				String seq_id = brCorrIn.readLine(), seq_str, seq_qual;
				while( seq_id!=null ) {
					seq_id = seq_id.split("\\s+")[0];
					seq_str = brCorrIn.readLine();
					brCorrIn.readLine();
					seq_qual = brCorrIn.readLine();
					batch.add(new FastqRecord(seq_id, seq_str, seq_qual));
					while( (seq_id=brCorrIn.readLine())!=null && batch.size()<batch_size ) {
						seq_id = seq_id.split("\\s+")[0];
						seq_str = brCorrIn.readLine();
						brCorrIn.readLine();
						seq_qual = brCorrIn.readLine();
						batch.add(new FastqRecord(seq_id, seq_str, seq_qual));	
					}
					Collections.sort(batch);

					checkDuplicateRead(batch);

					BufferedWriter bwCorrOutTmp = fastq_out.endsWith(".gz") ?
							Utils.getGZIPBufferedWriter(fastq_out+"_"+String.format("%010d", nBatch)+".tmp.gz") : 
								new BufferedWriter(new FileWriter(fastq_out+"_"+String.format("%010d", nBatch)+".tmp"));
							for(FastqRecord b : batch) {
								bwCorrOutTmp.write(b.seq_id+"\n");
								bwCorrOutTmp.write(b.seq_str+"\n");
								bwCorrOutTmp.write("+\n");
								bwCorrOutTmp.write(b.seq_qual+"\n");
							}
							bwCorrOutTmp.close();
							nBatch++;
							batch.clear();				
				}
				brCorrIn.close();

				BufferedWriter bwCorrOut = fastq_out.endsWith(".gz") ? 
						Utils.getGZIPBufferedWriter(fastq_out) : 
							new BufferedWriter(new FileWriter(fastq_out));
						BufferedReader[] brCorrOutTmp = new BufferedReader[nBatch];
						final boolean[] reachFileEnd = new boolean[nBatch];
						for(int i=0; i!=nBatch; i++) { 
							brCorrOutTmp[i] = fastq_out.endsWith(".gz") ? 
									Utils.getBufferedReader(fastq_out+"_"+String.format("%010d", i)+".tmp.gz") :
										new BufferedReader(new FileReader(fastq_out+"_"+String.format("%010d", i)+".tmp"));
									reachFileEnd[i] = false;
						}

						final TreeMap<FastqRecord, Integer> treeMap = new TreeMap<FastqRecord, Integer>();
						for(int i=0; i!=nBatch; i++) {
							seq_id = brCorrOutTmp[i].readLine().split("\\s+")[0];
							seq_str = brCorrOutTmp[i].readLine();
							brCorrOutTmp[i].readLine();
							seq_qual = brCorrOutTmp[i].readLine();
							treeMap.put(new FastqRecord(seq_id, seq_str, seq_qual), i);
						}
						Entry<FastqRecord, Integer> firstEntry;
						int bch, nReachFileEnd=0;
						FastqRecord fq;
						while( !treeMap.isEmpty() ) {
							firstEntry = treeMap.pollFirstEntry();
							fq = firstEntry.getKey();
							bwCorrOut.write(fq.seq_id+"\n");
							bwCorrOut.write(fq.seq_str+"\n");
							bwCorrOut.write("+\n");
							bwCorrOut.write(fq.seq_qual+"\n");
							bch = firstEntry.getValue();
							if(!reachFileEnd[bch]) {
								seq_id = brCorrOutTmp[bch].readLine();
								if(seq_id==null) {
									reachFileEnd[bch] = true;
									nReachFileEnd++;
								} else {
									seq_id = seq_id.split("\\s+")[0];
									seq_str = brCorrOutTmp[bch].readLine();
									brCorrOutTmp[bch].readLine();
									seq_qual = brCorrOutTmp[bch].readLine();
									treeMap.put(new FastqRecord(seq_id, seq_str, seq_qual), bch);
								}
							}
							if(treeMap.isEmpty()&&nReachFileEnd!=nBatch) {
								for(int i=0; i!=nBatch; i++) { 
									if(!reachFileEnd[i]) {
										seq_id = brCorrOutTmp[i].readLine().split("\\s+")[0];
										seq_str = brCorrOutTmp[i].readLine();
										brCorrOutTmp[i].readLine();
										seq_qual = brCorrOutTmp[i].readLine();
										treeMap.put(new FastqRecord(seq_id, seq_str, seq_qual), i);
									}
								}
							}
						}
						for(int i=0; i!=nBatch; i++) {
							brCorrOutTmp[i].close();
							if(fastq_out.endsWith(".gz")) 
								new File(fastq_out+"_"+String.format("%010d", i)+".tmp.gz").delete();
							else
								new File(fastq_out+"_"+String.format("%010d", i)+".tmp").delete();
						}
						bwCorrOut.close();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				Thread t = Thread.currentThread();
				t.getUncaughtExceptionHandler().uncaughtException(t, e);
				e.printStackTrace();
				executor.shutdown();
				System.exit(1);
			}
		}
	}

	private final class CorrFileSorter implements Runnable {
		private final String corr_in;
		private final String corr_out;

		public CorrFileSorter(String corrIn, String corrOut) {
			// TODO Auto-generated constructor stub
			this.corr_in = corrIn;
			this.corr_out = corrOut;
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub
			try {
				BufferedReader brCorrIn = new BufferedReader(new FileReader(corr_in));
				final List<String> batch = new ArrayList<String>();
				int nBatch = 0;
				String line = brCorrIn.readLine();
				while( line!=null ) {
					batch.add(line);
					while( (line=brCorrIn.readLine())!=null && batch.size()<batch_size ) {
						batch.add(line);
					}
					Collections.sort(batch);
					BufferedWriter bwCorrOutTmp = new BufferedWriter(new FileWriter(corr_out+"_"+String.format("%010d", nBatch)+".tmp"));
					for(String b : batch) bwCorrOutTmp.write(b+"\n");
					bwCorrOutTmp.close();
					nBatch++;
					batch.clear();				
				}
				brCorrIn.close();

				BufferedWriter bwCorrOut = new BufferedWriter(new FileWriter(corr_out));
				BufferedReader[] brCorrOutTmp = new BufferedReader[nBatch];
				final boolean[] reachFileEnd = new boolean[nBatch];
				for(int i=0; i!=nBatch; i++) { 
					brCorrOutTmp[i] = new BufferedReader(new FileReader(corr_out+"_"+String.format("%010d", i)+".tmp"));
					reachFileEnd[i] = false;
				}
				System.out.println(brCorrOutTmp.length);
				final TreeMap<String, Integer> treeMap = new TreeMap<String, Integer>();
				for(int i=0; i!=nBatch; i++) {
					line = brCorrOutTmp[i].readLine();
					while(line!=null&&treeMap.containsKey(line)) 
						line = brCorrOutTmp[i].readLine();
					if(line!=null) treeMap.put(line, i);
				}
				Entry<String, Integer> firstEntry;
				int bch, nReachFileEnd=0;
				while( !treeMap.isEmpty() ) {
					firstEntry = treeMap.pollFirstEntry();
					bwCorrOut.write(firstEntry.getKey()+"\n");
					bch = firstEntry.getValue();
					if(!reachFileEnd[bch]) {
						line = brCorrOutTmp[bch].readLine();
						if(line==null) {
							reachFileEnd[bch] = true;
							nReachFileEnd++;
						} else {
							while(line!=null&&treeMap.containsKey(line)) 
								line = brCorrOutTmp[bch].readLine();
							if(line!=null) treeMap.put(line, bch);
						}
					}
					if(treeMap.isEmpty()&&nReachFileEnd!=nBatch) {
						for(int i=0; i!=nBatch; i++) 
							if(!reachFileEnd[i]) {
								line = brCorrOutTmp[i].readLine();
								while(line!=null&&treeMap.containsKey(line)) 
									line = brCorrOutTmp[i].readLine();
								if(line!=null) treeMap.put(line, i);
							}
					}
				}
				for(int i=0; i!=nBatch; i++) {
					brCorrOutTmp[i].close();
					new File(corr_out+"_"+String.format("%010d", i)+".tmp").delete();
				}
				bwCorrOut.close();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				Thread t = Thread.currentThread();
				t.getUncaughtExceptionHandler().uncaughtException(t, e);
				e.printStackTrace();
				executor.shutdown();
				System.exit(1);
			}
		}
	}

	private final class PECorrecter implements Runnable {

		private final String[] sortedFastqIn;
		private final String sortedCorrIn;
		private final String faOutPrefix;

		public PECorrecter(String[] sortedFastqIn, 
				String sortedCorrIn, 
				String faOutPrefix) {
			// TODO Auto-generated constructor stub
			this.sortedFastqIn = sortedFastqIn;
			this.sortedCorrIn = sortedCorrIn;
			this.faOutPrefix = faOutPrefix;
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub
			try {
				final BufferedReader[] brFq = new BufferedReader[2];
				brFq[0] = Utils.getBufferedReader(sortedFastqIn[0]);
				brFq[1] = Utils.getBufferedReader(sortedFastqIn[1]);
				BufferedReader brCorr = Utils.getBufferedReader(sortedCorrIn);
				final BufferedWriter[] bwOut = new BufferedWriter[2];
				bwOut[0] = Utils.getBufferedWriter(faOutPrefix+".R1.fa");
				bwOut[1] = Utils.getBufferedWriter(faOutPrefix+".R2.fa");

				final List<List<CorrRecord>> buffered_corr = new ArrayList<List<CorrRecord>>();
				List<CorrRecord> buff;
				buffered_corr.add(new ArrayList<CorrRecord>()); // first read
				buffered_corr.add(new ArrayList<CorrRecord>()); // second read
				String corr_line = brCorr.readLine(), fq_line;
				String seq_id, seq_id_plus_comma, fq_id, fa_id, seq_ori, seq_sub;
				final String[] seq_pair = new String[2];
				final StringBuilder seq_sb = new StringBuilder();
				int shift, offset, seq_pos, oriS, subS, seqS;
				while( corr_line!=null ) {
					buffered_corr.get(0).clear();
					buffered_corr.get(1).clear();
					seq_id = corr_line.split(",")[0];
					seq_id_plus_comma = seq_id+",";

					bufferCorrRecord(corr_line, buffered_corr);
					while( (corr_line=brCorr.readLine())!=null &&
							corr_line.startsWith(seq_id_plus_comma) ) {
						bufferCorrRecord(corr_line, buffered_corr);
					}

					fq_id = "@"+seq_id;
					fa_id = ">"+seq_id+"\n";
					Arrays.fill(seq_pair, null);
					int skipped_lnno = 0;
					fq_line = brFq[0].readLine();
					while( !fq_line.startsWith(fq_id) ) {
						bwOut[0].write(">"+fq_line.split("\\s+")[0].substring(1)+"\n");
						bwOut[0].write(brFq[0].readLine()+"\n");
						brFq[0].readLine();
						brFq[0].readLine();
						fq_line = brFq[0].readLine();
						skipped_lnno++;
					}
					seq_pair[0] = brFq[0].readLine();
					brFq[0].readLine();
					brFq[0].readLine();

					for(int i=0; i!=skipped_lnno; i++) {
						bwOut[1].write(">"+brFq[1].readLine().split("\\s+")[0].substring(1)+"\n");
						bwOut[1].write(brFq[1].readLine()+"\n");
						brFq[1].readLine();
						brFq[1].readLine();
					}
					brFq[1].readLine();
					seq_pair[1] = brFq[1].readLine();
					brFq[1].readLine();
					brFq[1].readLine();

					if(seq_pair[0]==null||seq_pair[1]==null)
						throw new RuntimeException("!!!");

					for(int i=0; i!=2; i++) {
						buff = buffered_corr.get(i);
						if(buff.isEmpty()) {
							bwOut[i].write(fa_id);
							bwOut[i].write(seq_pair[i]+"\n");
							continue;
						}
						Collections.sort(buff);
						if(check) checkErrorCorr(buff);
						seq_sb.setLength(0);
						seq_sb.append(seq_pair[i]);
						shift = 0;
						seq_pos = 0;
						for(CorrRecord bf : buff) {
							seqS = seq_pos; // record last position updated

							seq_pos = bf.seq_pos;
							seq_ori = bf.seq_ori;
							seq_sub = bf.seq_sub;
							oriS = seq_ori.length();
							subS = seq_sub.length();
							offset = shift+seq_pos;

							switch(bf.operator) {
							case S:
								if(offset<seqS||!seq_sb.substring(offset, offset+oriS).equals(seq_ori)) {
									check_stats[0]++;
									continue;
									// throw new RuntimeException("!!!");
								}
								seq_sb.replace(offset, offset+oriS, seq_sub);
								shift += subS-oriS;
								check_stats[9]++;
								break;
							case D:
								if(offset<seqS||!seq_sb.substring(offset, offset+oriS).equals(seq_ori)) {
									check_stats[1]++;
									continue;
									//throw new RuntimeException("!!!");
								}
								seq_sb.delete(offset, offset+oriS);
								shift -= oriS;
								check_stats[9]++;
								break;
							case I:
								if(offset<seqS) {
									check_stats[2]++;
									continue;
									//throw new RuntimeException("!!!");
								}
								seq_sb.insert(offset+1, seq_sub);
								shift += subS;
								check_stats[9]++;
								break;
							default:
								throw new RuntimeException("!!!");
							}
						}
						seq_sb.append("\n");

						bwOut[i].write(fa_id);
						bwOut[i].write(seq_sb.toString());
					}
				}

				for(int i=0; i!=2; i++) {
					while( (fq_line=brFq[i].readLine())!=null ) {
						bwOut[i].write(">"+fq_line.split("\\s+")[0].substring(1)+"\n");
						bwOut[i].write(brFq[i].readLine()+"\n");
						brFq[i].readLine();
						brFq[i].readLine();
					}
				}

				brFq[0].close();
				brFq[1].close();
				brCorr.close();
				bwOut[0].close();
				bwOut[1].close();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				Thread t = Thread.currentThread();
				t.getUncaughtExceptionHandler().uncaughtException(t, e);
				e.printStackTrace();
				executor.shutdown();
				System.exit(1);
			}
		}

	}

	private final class SECorrecter implements Runnable {

		private final String sorted_fastqIn;
		private final String sorted_corrIn;
		private final String out_prefix;

		public SECorrecter(String sortedFastqIn, String sortedCorrIn, String faOutPrefix) {
			// TODO Auto-generated constructor stub
			this.sorted_corrIn = sortedFastqIn;
			this.sorted_fastqIn = sortedFastqIn;
			this.out_prefix = faOutPrefix;
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub
			try {
				final BufferedReader brFq = Utils.getBufferedReader(sorted_fastqIn);
				BufferedReader brCorr = Utils.getBufferedReader(sorted_corrIn);
				final BufferedWriter bwOut = Utils.getBufferedWriter(out_prefix+".fa");

				final List<CorrRecord> buffered_corr = new ArrayList<CorrRecord>();
				String corr_line = brCorr.readLine(), fq_line;
				String seq_id, seq_str, seq_id_plus_comma, fq_id, fa_id, seq_ori, seq_sub;
				final StringBuilder seq_sb = new StringBuilder();
				int shift, offset, seq_pos, oriS, subS, seqS;
				while( corr_line!=null ) {
					buffered_corr.clear();
					seq_id = corr_line.split(",")[0];
					seq_id_plus_comma = seq_id+",";

					String[] s = corr_line.split(",");
					buffered_corr.add(new CorrRecord(s[0],Integer.parseInt(s[2]),s[3],s[4],s[5]));
					while( (corr_line=brCorr.readLine())!=null &&
							corr_line.startsWith(seq_id_plus_comma) ) {
						s = corr_line.split(",");
						buffered_corr.add(new CorrRecord(s[0],Integer.parseInt(s[2]),s[3],s[4],s[5]));
					}

					fq_id = "@"+seq_id;
					fa_id = ">"+seq_id+"\n";

					fq_line = brFq.readLine();
					while( !fq_line.startsWith(fq_id) ) {
						bwOut.write(">"+fq_line.split("\\s+")[0].substring(1)+"\n");
						bwOut.write(brFq.readLine()+"\n");
						brFq.readLine();
						brFq.readLine();
						fq_line = brFq.readLine();
					}
					seq_str = brFq.readLine();
					brFq.readLine();
					brFq.readLine();

					if(buffered_corr.isEmpty()) {
						bwOut.write(fa_id);
						bwOut.write(seq_str+"\n");
						continue;
					}
					Collections.sort(buffered_corr);
					if(check) checkErrorCorr(buffered_corr);

					seq_sb.setLength(0);
					seq_sb.append(seq_str);
					shift = 0;
					seq_pos = 0;
					for(CorrRecord bf : buffered_corr) {
						seqS = seq_pos; // record last position updated

						seq_pos = bf.seq_pos;
						seq_ori = bf.seq_ori;
						seq_sub = bf.seq_sub;
						oriS = seq_ori.length();
						subS = seq_sub.length();
						offset = shift+seq_pos;

						switch(bf.operator) {
						case S:
							if(offset<seqS||!seq_sb.substring(offset, offset+oriS).equals(seq_ori)) {
								check_stats[0]++;
								continue;
								//throw new RuntimeException("!!!");
							}
							seq_sb.replace(offset, offset+oriS, seq_sub);
							shift += subS-oriS;
							check_stats[9]++;
							break;
						case D:
							if(offset<seqS||!seq_sb.substring(offset, offset+oriS).equals(seq_ori)) {
								check_stats[1]++;
								continue;
								//throw new RuntimeException("!!!");
							}
							seq_sb.delete(offset, offset+oriS);
							shift -= oriS;
							check_stats[9]++;
							break;
						case I:
							if(offset<seqS) {
								check_stats[2]++;
								continue;
								//throw new RuntimeException("!!!");
							}
							seq_sb.insert(offset+1, seq_sub);
							shift += subS;
							check_stats[9]++;
							break;
						default:
							throw new RuntimeException("!!!");
						}
					}
					seq_sb.append("\n");

					bwOut.write(fa_id);
					bwOut.write(seq_sb.toString());
				}

				while( (fq_line=brFq.readLine())!=null ) {
					bwOut.write(">"+fq_line.split("\\s+")[0].substring(1)+"\n");
					bwOut.write(brFq.readLine()+"\n");
					brFq.readLine();
					brFq.readLine();
				}

				brFq.close();
				brCorr.close();
				bwOut.close();

			} catch (Exception e) {
				// TODO Auto-generated catch block
				Thread t = Thread.currentThread();
				t.getUncaughtExceptionHandler().uncaughtException(t, e);
				e.printStackTrace();
				executor.shutdown();
				System.exit(1);
			}
		}
	}


	private final class HetCorrecter implements Runnable {

		private final int chr;
		private final int queryStart;
		private final int queryEnd;
		private final String corrOutPref;
		private final List<SAMSequenceRecord> refSeq;

		public HetCorrecter(int chr, 
				int queryStart, 
				int queryEnd, 
				List<SAMSequenceRecord> refSeq, 
				String out) {
			// TODO Auto-generated constructor stub
			this.chr = chr;
			this.queryStart = queryStart;
			this.queryEnd = queryEnd;
			this.refSeq = refSeq;
			this.corrOutPref = out;
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub

			try {
				// nucleotide x4 + deletion
				// A C G T D
				final int[] allele_stats = new int[5];
				// insertion is processed separately
				final Set<Integer> ins_pool = new HashSet<Integer>();
				// insertion statistics - could be multiple alleles
				// probably contain key " " - means identical to the reference/no insertion
				final Map<String, Integer> ins_stats = new HashMap<String, Integer>();

				Set<SAMRecord> tmp_set;
				final Set<SAMRecord> tmp_removal = new HashSet<SAMRecord>();
				final int[] ht_hc = new int[2]; // head and tail hard clip length
				String dna_seq, ins_seq;
				int pos, tmp_int, sel, len_seq, pair_flag;
				Map.Entry<String, Integer> ins_sel;
				char nucl, nucl2;
				boolean ins = false;

				final SamReader[] in1 = new SamReader[bam_list.length];
				final SAMRecordIterator[] iter1 = new SAMRecordIterator[in1.length];
				final SAMRecord[] buffer = new SAMRecord[iter1.length];
				final int[] bufferS = new int[buffer.length];
				final List<Set<SAMRecord>> record_pool = new ArrayList<Set<SAMRecord>>();
				final BufferedWriter[] corrOut = new BufferedWriter[bam_list.length];
				for(int i=0; i!=corrOut.length; i++) {
					try {
						corrOut[i] = new BufferedWriter(new FileWriter(corrOutPref+"_"+
								new File(bam_list[i]).getName()+".corr"));
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}

				Arrays.fill(in1, null);
				Arrays.fill(iter1, null);
				Arrays.fill(buffer, null);
				Arrays.fill(bufferS, 0);
				record_pool.clear();

				for(int i=0; i!=bam_list.length; i++) {
					in1[i] = factory.open(new File(bam_list[i]));
					iter1[i] = in1[i].queryOverlapping(refSeq.get(chr).getSequenceName(), 
							queryStart, queryEnd);
					buffer[i] = iter1[i].hasNext()?iter1[i].next():null;
					bufferS[i] = buffer[i]==null?-1:buffer[i].getAlignmentStart();
					record_pool.add(new HashSet<SAMRecord>());
				}

				for(int p=queryStart; p!=queryEnd; p++) {

					if(p%1000000==0) System.err.println(refSeq.get(chr).getSequenceName()+":"+p);
					pos = p+1; // BAM format is 1-based coordination

					// if(pos==3000) break outerloop;
					ins = ins_pool.contains(pos);
					if(ins) {
						ins_pool.remove(pos);
						ins_stats.clear();
					}
					Arrays.fill(allele_stats, 0);

					// CAPTURE:
					// allele statistics
					for(int i=0; i!=buffer.length; i++) {
						tmp_set = record_pool.get(i);
						if(bufferS[i]<=pos) 
							buffer(tmp_set, ins_pool, iter1, buffer, bufferS, pos, i);
						if(tmp_set.isEmpty()) continue;

						tmp_removal.clear();
						for(SAMRecord record : tmp_set) {

							if(record.getAlignmentEnd()<pos) {
								tmp_removal.add(record);
								continue;
							}
							dna_seq = record.getReadString();
							tmp_int = record.getReadPositionAtReferencePosition(pos);
							nucl = tmp_int==0?'D':dna_seq.charAt(tmp_int-1);
							switch(nucl) {
							case 'A':
								allele_stats[0]++;
								break;
							case 'C':
								allele_stats[1]++;
								break;
							case 'G':
								allele_stats[2]++;
								break;
							case 'T':
								allele_stats[3]++;
								break;
							case 'D':
								allele_stats[4]++;
								break;
							case 'N':
								// do nothing here
								break;
							default:
								throw new RuntimeException("!!!");
							}
						}
						tmp_set.removeAll(tmp_removal);

						if(ins) {
							for(SAMRecord record : tmp_set) {
								ins_seq = insertion_seq(record, pos, null, null);
								if(ins_seq==null) continue;
								if(!ins_stats.containsKey(ins_seq)) ins_stats.put(ins_seq, 0);
								ins_stats.put(ins_seq, ins_stats.get(ins_seq)+1);
							}
						}
					}

					try {
						// CORR: substitution and deletion
						// allele to select
						sel = max_index(allele_stats);
						if(sel==4) {
							// will choose the deletion
							for(int i=0; i!=buffer.length; i++) {
								tmp_set = record_pool.get(i);
								for(SAMRecord record : tmp_set) {
									len_seq  = seqLength(record, ht_hc);
									dna_seq = record.getReadString();
									tmp_int = record.getReadPositionAtReferencePosition(pos);
									nucl = tmp_int==0?'D':dna_seq.charAt(tmp_int-1);
									pair_flag = record.getReadPairedFlag()?(record.getFirstOfPairFlag()?0:1):0;
									if(nucl!='D') {
										if(record.getReadNegativeStrandFlag()) {
											//System.out.println("File "+i+","+record.getReadName()+","+(len_seq-tmp_int-ht_hc[0]-1)+","+revcmp.get(nucl)+",,D");
											corrOut[i].write(record.getReadName()+","+pair_flag+","+(len_seq-tmp_int-ht_hc[0])+","+revcmp.get(nucl)+",,D\n");
										} else {
											//System.out.println("File "+i+","+record.getReadName()+","+(tmp_int+ht_hc[0])+","+nucl+",,D");
											corrOut[i].write(record.getReadName()+","+pair_flag+","+(tmp_int+ht_hc[0]-1)+","+nucl+",,D\n");
										}
									}
									// if(ht_hc[0]!=0||ht_hc[1]!=0) System.out.println(record.getSAMString());
								}
							}
						} else {
							// will choose a substitution/nucleotide base
							nucl2 = nucleotide[sel];
							int shift;
							for(int i=0; i!=buffer.length; i++) {
								tmp_set = record_pool.get(i);
								for(SAMRecord record : tmp_set) {

									len_seq  = seqLength(record, ht_hc);
									dna_seq = record.getReadString();
									tmp_int = record.getReadPositionAtReferencePosition(pos);
									nucl = tmp_int==0?'D':dna_seq.charAt(tmp_int-1);
									pair_flag = record.getReadPairedFlag()?(record.getFirstOfPairFlag()?0:1):0;
									if(nucl=='D') {
										// find insert position on the read
										shift = 0;
										int rs = record.getAlignmentStart();
										while(tmp_int==0&&pos-shift>rs) {
											shift++;
											tmp_int = record.getReadPositionAtReferencePosition(pos-shift); 
										}
										if(record.getReadNegativeStrandFlag()) {
											//System.out.println("File "+i+","+record.getReadName()+","+(len_seq-tmp_int-ht_hc[0]-2)+","+shift+","+revcmp.get(nucl2)+",I");
											corrOut[i].write(record.getReadName()+","+pair_flag+","+(len_seq-tmp_int-ht_hc[0]-1)+",-"+shift+","+revcmp.get(nucl2)+",I\n");
										} else {
											//System.out.println("File "+i+","+record.getReadName()+","+(tmp_int+ht_hc[0])+","+shift+","+nucl2+",I");
											corrOut[i].write(record.getReadName()+","+pair_flag+","+(tmp_int+ht_hc[0]-1)+","+shift+","+nucl2+",I\n");
										}
										// if(ht_hc[0]!=0||ht_hc[1]!=0) System.out.println(record.getSAMString());
									} else if(nucl!=nucl2) {
										if(record.getReadNegativeStrandFlag()) {
											//System.out.println("File "+i+","+record.getReadName()+","+(len_seq-tmp_int-ht_hc[0]-1)+","+revcmp.get(nucl)+","+revcmp.get(nucl2)+",S");
											corrOut[i].write(record.getReadName()+","+pair_flag+","+(len_seq-tmp_int-ht_hc[0])+","+revcmp.get(nucl)+","+revcmp.get(nucl2)+",S\n");
										} else {
											//System.out.println("File "+i+","+record.getReadName()+","+(tmp_int+ht_hc[0])+","+nucl+","+nucl2+",S");
											corrOut[i].write(record.getReadName()+","+pair_flag+","+(tmp_int+ht_hc[0]-1)+","+nucl+","+nucl2+",S\n");
										}
										// if(ht_hc[0]!=0||ht_hc[1]!=0) System.out.println(record.getSAMString());
									}
								}
							}
						}

						// CORR: insertion
						if(ins) {
							ins_sel = max_index(ins_stats);
							String ins_allele = ins_sel.getKey();
							final int[] ht_hc_ins = new int[2];
							final int[] rl_shift = new int[2];
							for(int i=0; i!=buffer.length; i++) {
								tmp_set = record_pool.get(i);
								for(SAMRecord record : tmp_set) {									
									len_seq  = seqLength(record, ht_hc);
									ins_seq = insertion_seq(record, pos, ht_hc_ins, rl_shift);
									if(ins_seq==null) continue;
									pair_flag = record.getReadPairedFlag()?(record.getFirstOfPairFlag()?0:1):0;
									if(record.getReadNegativeStrandFlag()) {
										if(ins_allele.equals("")) {
											if(!ins_seq.equals("")) {
												//System.out.println("File "+i+","+record.getReadName()+","+(len_seq-ht_hc_ins[1]-ht_hc[0])+","+revcmp_seq(ins_seq)+",,D");
												corrOut[i].write(record.getReadName()+","+pair_flag+","+(len_seq-ht_hc_ins[1]-ht_hc[0])+","+revcmp_seq(ins_seq)+",,D\n");
											}
										} else {
											if(ins_seq.equals("")) {
												//System.out.println("File "+i+","+record.getReadName()+","+(len_seq-ht_hc_ins[1]-ht_hc[0]-2)+","+revcmp_seq(ins_allele)+","+rl_shift[1]+",I");
												corrOut[i].write(record.getReadName()+","+pair_flag+","+(len_seq-ht_hc_ins[1]-ht_hc[0]-1)+","+rl_shift[1]+","+revcmp_seq(ins_allele)+",I\n");
											} else if(!ins_seq.equals(ins_allele)){
												//System.out.println("File "+i+","+record.getReadName()+","+(len_seq-ht_hc_ins[1]-ht_hc[0])+","+revcmp_seq(ins_seq)+","+revcmp_seq(ins_allele)+",S");
												corrOut[i].write(record.getReadName()+","+pair_flag+","+(len_seq-ht_hc_ins[1]-ht_hc[0])+","+revcmp_seq(ins_seq)+","+revcmp_seq(ins_allele)+",S\n");
											}
										}
									} else {
										if(ins_allele.equals("")) {
											if(!ins_seq.equals("")) {
												//System.out.println("File "+i+","+record.getReadName()+","+(ht_hc_ins[0]+ht_hc[0])+","+ins_seq+",,D");
												corrOut[i].write(record.getReadName()+","+pair_flag+","+(ht_hc_ins[0]+ht_hc[0])+","+ins_seq+",,D\n");
											}
										} else {
											if(ins_seq.equals("")) {
												//System.out.println("File "+i+","+record.getReadName()+","+(ht_hc_ins[0]+ht_hc[0])+","+revcmp_seq(ins_allele)+","+rl_shift[0]+",I");
												corrOut[i].write(record.getReadName()+","+pair_flag+","+(ht_hc_ins[0]+ht_hc[0]-1)+","+rl_shift[0]+","+ins_allele+",I\n");
											} else if(!ins_seq.equals(ins_allele)) {
												//System.out.println("File "+i+","+record.getReadName()+","+(ht_hc_ins[0]+ht_hc[0])+","+revcmp_seq(ins_seq)+","+revcmp_seq(ins_allele)+",S");
												corrOut[i].write(record.getReadName()+","+pair_flag+","+(ht_hc_ins[0]+ht_hc[0])+","+ins_seq+","+ins_allele+",S\n");
											}
										}
									}
								}
							}
						}
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}

				try {
					for(int i=0; i!=bam_list.length; i++) {
						iter1[i].close();
						in1[i].close();
					}
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

				try {
					for(int i=0; i!=bam_list.length; i++) 
						corrOut[i].close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
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
	}
}
