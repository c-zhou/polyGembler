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
import java.io.Reader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Map.Entry;
import java.util.Set;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.RejectedExecutionHandler;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.lang3.RandomStringUtils;
import org.apache.commons.lang3.StringUtils;

import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

import cz1.ngs.model.AlignmentSegment;
import cz1.ngs.model.Blast6Segment;
import cz1.ngs.model.Sequence;
import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Consensus extends Executor {

	private final static String s_lib = "--s([0-9]+)";
	private final static String pe_lib = "--pe([0-9]+)";
	private final static String mp_lib = "--mp([0-9]+)";
	private final static String ins_lib = "--f([0-9]+)";
	private final static String w_lib = "--w([0-9]+)";
	
	private final static Pattern s_pat = Pattern.compile(s_lib);
	private final static Pattern pe_pat = Pattern.compile(pe_lib);
	private final static Pattern mp_pat = Pattern.compile(mp_lib);
	private final static Pattern ins_pat = Pattern.compile(ins_lib);
	private final static Pattern w_pat = Pattern.compile(w_lib);
	
	private static enum Task {all, count, parse, zzz}
	private Task task_list = Task.zzz;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		switch(this.task_list) {
		case zzz:
			myLogger.info(
					"\n\nChoose a task from listed:\n"
							+ " count                   Calculate weighted link counts.\n"
							+ " parse                   Parse consensus scaffolds. \n"
							+ " all                     Count links and then parse scaffolds.\n"
							+ "\n");
			break;
		case count:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " --s<#>                  Input BAM file (sorted by read name) for long read library (<#> = 1,2,...).\n"
							+ " --pe<#>                 Input BAM file (sorted by read name) for paired-end library (<#> = 1,2,...).\n"
							+ " --mp<#>                 Input BAM file (sorted by read name) for mate-pair library (<#> = 1,2,...).\n"
							+ " --f<#>                  Fragment/insert size threshold of this library (<#> = 1,2,...).\n"
							+ "                         NOTE: required for paired-end and mate-pair libraries,\n"
							+ "                               and will be ignored if provided for long read libraries.\n"
							+ " --w<#>                  Inverse weight of the links from this library (<#> = 1,2,...).\n"
							+ "                         This parameter defines the minimum number of links in this library are required \n"
							+ "                         to confirm a consensus. Generally larger insert size libraries are less realiable \n"
							+ "                         thus a larger number should be assigned. \n"
							+ "                         NOTE: if multiple libraries (BAM files) are provided, the weights are summed. \n"
							+ "                               If the total weight is no less than 1, the consensus is confirmed.\n"
							+ "                         FOR EXAMPLE: --w1 3 --w2 5 --w3 10. To confirm a consensus, \n"
							+ "                               1. we need 3 links from the library 1 (W=1.0).\n"
							+ "                               Or,\n"
							+ "                               2. we need 2 links from the library 1 and 2 links from the library 2 (W=1.06). \n"
							+ "                               Or,\n"
							+ "                               3. we need 4 links from the library 2 and 2 links from the library 3 (W=1.0).\n"
							+ " -m/--map                Map file generated in the anchor step indicating the order of the contigs. \n"
							+ " -t/--threads            Threads to use (default 1). \n"
							+ "                         The maximum number of threads to use will be the number of BAM files.\n"
							+ " -o/--out                Prefix of the output link count file.\n"
							+ "\n");
			break;
		case parse:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -q/--query              The FASTA file contain query sequences to anchor. \n"
							+ " -m/--map                Map file generated in the link count step indicating the order of the contigs. \n"
							+ " -l/--min-size           Minimum size of scaffolds to output (default 300). \n"
							+ " -o/--out                Output FASTA file name.\n"
							+ "\n");
			break;
		case all:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " --s<#>                  Input BAM file (sorted by read name) for long read library (<#> = 1,2,...).\n"
							+ " --pe<#>                 Input BAM file (sorted by read name) for paired-end library (<#> = 1,2,...).\n"
							+ " --mp<#>                 Input BAM file (sorted by read name) for mate-pair library (<#> = 1,2,...).\n"
							+ " --f<#>                  Fragment/insert size threshold of this library (<#> = 1,2,...).\n"
							+ "                         NOTE: required for paired-end libraries,\n"
							+ "                               and will be ignored if provided for long read libraries.\n"
							+ " --w<#>                  Inverse weight of the links from this library (<#> = 1,2,...).\n"
							+ "                         This parameter defines the minimum number of links in this library are required \n"
							+ "                         to confirm a consensus. Generally larger insert size libraries are less realiable \n"
							+ "                         thus a larger number should be assigned. \n"
							+ "                         NOTE: if multiple libraries (BAM files) are provided, the weights are summed. \n"
							+ "                               If the total weight is no less than 1, the consensus is confirmed.\n"
							+ "                         FOR EXAMPLE: --w1 3 --w2 5 --w3 10. To confirm a consensus, \n"
							+ "                               1. we need 3 links from the library 1 (W=1.0).\n"
							+ "                               Or,\n"
							+ "                               2. we need 2 links from the library 1 and 2 links from the library 2 (W=1.06). \n"
							+ "                               Or,\n"
							+ "                               3. we need 4 links from the library 2 and 2 links from the library 3 (W=1.0).\n"
							+ " -q/--query              The FASTA file contain query sequences to anchor. \n"
							+ " -m/--map                Map file indicate the order of the contigs. \n"
							+ " -t/--threads            Threads to use (default 1). \n"
							+ "                         The maximum number of threads to use will be the number of BAM files.\n"
							+ " -l/--min-size           Minimum size of scaffolds to output (default 300). \n"
							+ " -o/--out                Prefix of the output files.\n"
							+ "\n");
			break;
		default:
			throw new RuntimeException("!!!");
		}
	}

	private static enum Library {s, pe, mp}
	private String sequence_file = null;
	private String map_file = null;
	private int min_size = 300;
	private String out_prefix = null;
	private String[] bam_list = null;
	// specify the library types: s, pe, mp
	private Library[] lib_list= null;
	private int[] ins_thres = null;
	private double[] link_w = null;
	
	private double collinear_shift = 0.5;
	
	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}
		
		switch(args[0].toUpperCase()) {
		case "COUNT":
			this.task_list = Task.count;
			break;
		case "PARSE":
			this.task_list = Task.parse;
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
		
		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-q", "--query", true);
			myArgsEngine.add("-m", "--map", true);
			myArgsEngine.add("-as", "--block-size", true);
			myArgsEngine.add("-bs", "--batch-size", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-l", "--min-size", true);
			myArgsEngine.add("-o", "--out", true);
			myArgsEngine.addWildOptions(s_lib, true);
			myArgsEngine.addWildOptions(pe_lib, true);
			myArgsEngine.addWildOptions(mp_lib, true);
			myArgsEngine.addWildOptions(ins_lib, true);
			myArgsEngine.addWildOptions(w_lib, true);
			myArgsEngine.parse(args2);
		}
		
		switch(this.task_list) {
		case count:

			if (myArgsEngine.getBoolean("-m")) {
				this.map_file = myArgsEngine.getString("-m");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the map file.");
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
		case parse:
			
			if (myArgsEngine.getBoolean("-q")) {
				this.sequence_file = myArgsEngine.getString("-q");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the contig file.");
			}

			if (myArgsEngine.getBoolean("-m")) {
				this.map_file = myArgsEngine.getString("-m");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the map file.");
			}
			
			if (myArgsEngine.getBoolean("-l")) {
				this.min_size = Integer.parseInt(myArgsEngine.getString("-l"));
			}
			
			if (myArgsEngine.getBoolean("-o")) {
				this.out_prefix = myArgsEngine.getString("-o");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the output file.");
			}
			
			break;
		case all:
			
			if (myArgsEngine.getBoolean("-q")) {
				this.sequence_file = myArgsEngine.getString("-q");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the contig file.");
			}

			if (myArgsEngine.getBoolean("-m")) {
				this.map_file = myArgsEngine.getString("-m");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the map file.");
			}
			
			if (myArgsEngine.getBoolean("-t")) {
				this.THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
			}
			
			if (myArgsEngine.getBoolean("-l")) {
				this.min_size = Integer.parseInt(myArgsEngine.getString("-l"));
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
		final Map<Integer, String> sLib = new HashMap<Integer, String>();
		final Map<Integer, String> peLib = new HashMap<Integer, String>();
		final Map<Integer, String> mpLib = new HashMap<Integer, String>();
		final Map<Integer, Integer> insLib = new HashMap<Integer, Integer>();
		final Map<Integer, Integer> wLib = new HashMap<Integer, Integer>();

		int n = 0;
		for(int i=0; i<args.length; i++) {
			String arg = args[i];
			if(arg.matches(s_lib)) {
				Matcher m = s_pat.matcher(arg);
				m.find();
				int lib = Integer.parseInt(m.group(1));
				sLib.put(lib, args[++i]);
				++n;
			}
			if(arg.matches(pe_lib)) {
				Matcher m = pe_pat.matcher(arg);
				m.find();
				int lib = Integer.parseInt(m.group(1));
				peLib.put(lib, args[++i]);
				++n;
			}
			if(arg.matches(mp_lib)) {
				Matcher m = mp_pat.matcher(arg);
				m.find();
				int lib = Integer.parseInt(m.group(1));
				mpLib.put(lib, args[++i]);
				++n;
			}
			if(arg.matches(ins_lib)) {
				Matcher m = ins_pat.matcher(arg);
				m.find();
				int lib = Integer.parseInt(m.group(1));
				insLib.put(lib, Integer.parseInt(args[++i]));
			}
			if(arg.matches(w_lib)) {
				Matcher m = w_pat.matcher(arg);
				m.find();
				int lib = Integer.parseInt(m.group(1));
				wLib.put(lib, Integer.parseInt(args[++i]));
			}
		}

		if(wLib.size()!=n) {
			printUsage();
			throw new IllegalArgumentException("Library parameters do not match!!!");
		}
		
		this.bam_list = new String[n];
		this.lib_list = new Library[n];
		this.ins_thres = new int[n];
		this.link_w = new double[n];
		
		int i=0;
		for(Integer key : sLib.keySet()) {
			if(!wLib.containsKey(key)) {
				printUsage();
				throw new IllegalArgumentException("Library parameters do not match!!!");
			}
			bam_list[i] = sLib.get(key);
			lib_list[i] = Library.s;
			ins_thres[i] = Integer.MAX_VALUE;
			link_w[i] = 1.0/wLib.get(key);
			i++;
		}
		
		for(Integer key : peLib.keySet()) {
			if(!insLib.containsKey(key) || !wLib.containsKey(key)) {
				printUsage();
				throw new IllegalArgumentException("Library parameters do not match!!!");
			}
			bam_list[i] = peLib.get(key);
			lib_list[i] = Library.pe;
			ins_thres[i] = insLib.get(key);
			link_w[i] = 1.0/wLib.get(key);
			i++;
		}
		
		for(Integer key : mpLib.keySet()) {
			if(!insLib.containsKey(key) || !wLib.containsKey(key)) {
				printUsage();
				throw new IllegalArgumentException("Library parameters do not match!!!");
			}
			bam_list[i] = mpLib.get(key);
			lib_list[i] = Library.mp;
			ins_thres[i] = insLib.get(key);
			link_w[i] = 1.0/wLib.get(key);
			i++;
		}
	}
	
	private String[] link_file = null;
	@Override
	public void run() {
		// TODO Auto-generated method stub
	
		String link_out, scaff_out;
		switch(this.task_list) {
		
		case zzz:
			myLogger.info("Task list is empty!!!");
			break;

		case count:
			// STEP 1. count links
			//         merge links
			myLogger.info("STEP 1. count links");
			link_file = new String[this.bam_list.length];
			this.initial_thread_pool();
			for(int i=0; i<this.bam_list.length; i++) {
				String out = this.out_prefix+"_"+
						new File(bam_list[i]).getName().replaceAll(".bam$", "")+".linkC";
				link_file[i] = out;
				switch(lib_list[i]) {
				case s:
					this.executor.submit(new LongReadLinkCounter(this.bam_list[i], this.map_file, out));
					break;
				case pe:
					this.executor.submit(new PairedReadLinkCounter(this.bam_list[i], this.ins_thres[i], this.map_file, true, out));
					break;
				case mp:
					this.executor.submit(new PairedReadLinkCounter(this.bam_list[i], this.ins_thres[i], this.map_file, false, out));
					break;
				default:
					throw new RuntimeException("Undefined library type!!!");	
				}
			}
			this.waitFor();

			myLogger.info("STEP 2. parse links");
			link_out = this.out_prefix+".linkW";
			this.parseLink(link_file, link_out);
			break;
			
		case parse:
			// STEP 3. parse scaffolds
			myLogger.info("STEP 3. parse scaffolds");
			scaff_out = this.out_prefix+"_parsedScaffold.fa";

			this.parseScaffold(this.sequence_file, 
					this.map_file,
					this.min_size, 
					scaff_out);
			break;
			
		case all:
			// STEP 1. count links
			//         merge links
			myLogger.info("STEP 1. count links");
			link_file = new String[this.bam_list.length];
			this.initial_thread_pool();
			for(int i=0; i<this.bam_list.length; i++) {
				String out = this.out_prefix+"_"+
						new File(bam_list[i]).getName().replaceAll(".bam$", "")+".linkC";
				link_file[i] = out;
				switch(lib_list[i]) {
				case s:
					this.executor.submit(new LongReadLinkCounter(this.bam_list[i], this.map_file, out));
					break;
				case pe:
					this.executor.submit(new PairedReadLinkCounter(this.bam_list[i], this.ins_thres[i], this.map_file, true, out));
					break;
				case mp:
					this.executor.submit(new PairedReadLinkCounter(this.bam_list[i], this.ins_thres[i], this.map_file, false, out));
					break;
				default:
					throw new RuntimeException("Undefined library type!!!");	
				}
			}
			this.waitFor();

			myLogger.info("STEP 2. parse links");
			link_out = this.out_prefix+".linkW";
			this.parseLink(link_file, link_out);
			
			myLogger.info("STEP 3. parse scaffolds");
			scaff_out = this.out_prefix+"_parsedScaffold.fa";

			this.parseScaffold(this.sequence_file, 
					link_out,
					this.min_size, 
					scaff_out);
			break;
			
		default:
			throw new RuntimeException("!!!");
		}
		return;
	}

	private final static Object lock = new Object();
	private final static SamReaderFactory factory =
			SamReaderFactory.makeDefault()
			.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
					SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
			.validationStringency(ValidationStringency.SILENT);
	private final static int refMask = Integer.MAX_VALUE; // 31 bits reference index
	private final static int dirMask = 1;                 // 1  bit  direction
	
	
	private final int insz_thres = 2;
	public void buildReadPairGraph(String bam_in,
			String sz_in,
			String out,
			//String contig_file,
			//int max_gap,
			int mapQ_thres,
			int max_NM) {
		try {
			/***
			BufferedReader br_tmp = new BufferedReader(new FileReader(new File(contig_file)));
			int contig_num = 0;
			String line;
			while( (line=br_tmp.readLine())!=null )
				if(line.startsWith(">")) contig_num++;
			br_tmp.close();

			BufferedReader br_contig_file = new BufferedReader(new FileReader(new File(contig_file)));
			String[] contig_id = new String[contig_num];
			String[] contig_str = new String[contig_num];
			int[] contig_sz = new int[contig_num];
			int c = 0;
			line = br_contig_file.readLine();
			while( line!=null ) {
				if(line.startsWith(">")) {
					contig_id[c] = line.replace(">", "").trim();
					StringBuilder bs = new StringBuilder();
					while( (line=br_contig_file.readLine())!=null
							&& !line.startsWith(">")) {
						bs.append(line.trim());
					}
					contig_str[c] = bs.toString();
					contig_sz[c] = bs.length();
					c++;
				}
			}
			br_contig_file.close();
			 **/
			
			String[] bam_files = bam_in.trim().split(";");
			String[] s = sz_in.trim().split(";");
			final Map<String, Integer> szs= new HashMap<String, Integer>();
			for(int i=0; i!=bam_files.length; i++) 
				szs.put(bam_files[i], Integer.parseInt(s[i]));
			
			final Map<Long, Integer> links = new HashMap<Long, Integer>();
			
			this.initial_thread_pool();
			
			for(String bam_file : bam_files) {

				this.executor.submit(new Runnable() {
					private String bam_file;
					
					@Override
					public void run() {
						System.err.println("Process file ... "+bam_file);
						
						SAMRecord tmp_record, sam_record, mate_record;
						String sam_id;
						final List<SAMRecord> record_pool = new ArrayList<SAMRecord>();
						int sam_ref, mate_ref;
						long key_ref;
						
						final SamReader in1 = factory.open(new File(bam_file));
						final List<SAMSequenceRecord> refSeq = in1.getFileHeader().
								getSequenceDictionary().getSequences();
						final int[] refSz = new int[refSeq.size()];
						for(int i=0; i!=refSz.length; i++) 
							refSz[i] = refSeq.get(i).getSequenceLength();
						
						SAMRecordIterator iter1 = in1.iterator();
						tmp_record = iter1.next();
						while( tmp_record!=null ) {
							record_pool.clear();

							record_pool.add(tmp_record);
							sam_id = tmp_record.getReadName();

							while( (tmp_record = iter1.hasNext() ? iter1.next() : null) !=null &&
									tmp_record.getReadName().equals(sam_id) ) 
								record_pool.add(tmp_record);

							sam_record = null;
							mate_record = null;
							for(SAMRecord record : record_pool) {
								if(record.getFirstOfPairFlag() && !record.getNotPrimaryAlignmentFlag())
									sam_record = record;
								if(record.getSecondOfPairFlag() && !record.getNotPrimaryAlignmentFlag())
									mate_record = record;
							}
							if(sam_record==null || mate_record==null)
								throw new RuntimeException(record_pool.get(0).getSAMString()+"!!!");

							sam_ref = sam_record.getReferenceIndex();
							mate_ref = mate_record.getReferenceIndex();
							if( sam_record.getReadUnmappedFlag() || 
									mate_record.getReadUnmappedFlag() || 
									sam_ref==mate_ref || 
									sam_record.getMappingQuality()<mapQ_thres ||
									mate_record.getMappingQuality()<mapQ_thres ||
									sam_record.getIntegerAttribute("NM")>max_NM ||
									mate_record.getIntegerAttribute("NM")>max_NM)
								continue;

							if(sam_ref>mate_ref) {
								int tmp_int = sam_ref;
								sam_ref = mate_ref;
								mate_ref = tmp_int;
								SAMRecord tmp_rec = sam_record;
								sam_record = mate_record;
								mate_record = tmp_rec;
							}
							
							int infSz;
							double maxSz = insz_thres*szs.get(bam_file);
							key_ref = 0L;
							
							/**FF*
							 * ---=>---		---=>--
							 */
							if(!sam_record.getReadNegativeStrandFlag() &&
									!mate_record.getReadNegativeStrandFlag()) {
								infSz = refSz[sam_ref]-sam_record.getAlignmentStart()+
										refSz[mate_ref]-mate_record.getAlignmentStart();
								
								if(infSz<=maxSz) {
									key_ref = 1; // 1 - end
									key_ref <<= 31;
									key_ref += sam_ref;
									key_ref <<= 1;
									key_ref += 1; // 1 - end
									key_ref <<= 31;
									key_ref += mate_ref;
								}
							} 
							
							/**FR*
							 * ---=>---		---<=--
							 */
							else if(!sam_record.getReadNegativeStrandFlag() &&
									mate_record.getReadNegativeStrandFlag()) {
								infSz = refSz[sam_ref]-sam_record.getAlignmentStart()+
										mate_record.getAlignmentEnd();
								
								if(infSz<=maxSz) {
									key_ref = 1; // 1 - end
									key_ref <<= 31;
									key_ref += sam_ref;
									key_ref <<= 1;
									key_ref += 0; // 0 - start
									key_ref <<= 31;
									key_ref += mate_ref;
								}
							}
							
							/**RF*
							 * ---<=---		---=>--
							 */
							else if(sam_record.getReadNegativeStrandFlag() &&
									!mate_record.getReadNegativeStrandFlag()) {
								infSz = sam_record.getAlignmentEnd()+
										refSz[mate_ref]-mate_record.getAlignmentStart();
								
								if(infSz<=maxSz) {
									key_ref = 0; // 0 - start
									key_ref <<= 31;
									key_ref += sam_ref;
									key_ref <<= 1;
									key_ref += 1; // 1 - end
									key_ref <<= 31;
									key_ref += mate_ref;
								}
							}
							
							/**RR*
							 * ---<=---		---<=--
							 */
							else if(sam_record.getReadNegativeStrandFlag() &&
									mate_record.getReadNegativeStrandFlag()) {
								infSz = sam_record.getAlignmentEnd()+
										mate_record.getAlignmentEnd();
								
								if(infSz<=maxSz) {
									key_ref = 0; // 0 - start
									key_ref <<= 31;
									key_ref += sam_ref;
									key_ref <<= 1;
									key_ref += 0; // 0 - start
									key_ref <<= 31;
									key_ref += mate_ref;
								}
							}
							
							if(key_ref==0L) continue;
							
							synchronized(lock) {
								if(links.containsKey(key_ref))
									links.put(key_ref, links.get(key_ref)+1);
								else links.put(key_ref, 1);
							}
						}
						iter1.close();
						try {
							in1.close();
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
						
						System.err.println("Process file ... "+bam_file+" done.");
					}
					
					public Runnable init(final String bam_file) {
				        this.bam_file = bam_file;
				        return(this);
				    }
            	}.init(bam_file));
			}
			
			this.waitFor();
			
			int mate_ref, sam_ref, sam_dir, mate_dir;
			long key_ref;
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(out)));
			for(Map.Entry<Long, Integer> entry : links.entrySet()) {
				key_ref = entry.getKey();
				mate_ref = (int) (key_ref&refMask);
				key_ref >>= 31;
				mate_dir = (int) (key_ref&dirMask);
				key_ref >>= 1;
				
				sam_ref = (int) (key_ref&refMask);
				key_ref >>= 31;
				sam_dir = (int) (key_ref&dirMask);

				bw.write(sam_ref+"\t"+mate_ref+"\t"+sam_dir+"\t"+mate_dir+"\t"+entry.getValue()+"\n");
			}
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	final static long mask_24bits = (long)Math.pow(2,24)-1;
	final static long mask_08bits = (long)Math.pow(2, 8)-1;
	private final static DecimalFormat df3 = new DecimalFormat("#.000");
	
	public void parseLink(final String[] links_in,
			final String link_out) {

		myLogger.info("Insert size threshold: "+strcat(this.ins_thres,","));
		myLogger.info("Link weights: "+strcat(link_w,","));
		
		try {
			int n = links_in.length;
			BufferedReader[] br_maps = new BufferedReader[n];
			for(int i=0; i<n; i++)
				br_maps[i] = new BufferedReader(new FileReader(links_in[i]));
			BufferedWriter link_bw = new BufferedWriter(new FileWriter(link_out));
			
			String line;
			String[] s;

			double weightedLinkCounts;
			int scaff_i = 0;
			String scaff = "Scaffold"+String.format("%08d", ++scaff_i);
			
			while( (line=br_maps[0].readLine())!=null ) {
				
				s = line.split("\\s+");
				
				link_bw.write(s[0]+"\t"+
						s[1]+"\t"+
						s[2]+"\t"+
						s[3]+"\t"+
						s[4]+"\t");
				
				if(line.startsWith("GAP")) {
					weightedLinkCounts = Double.parseDouble(s[8])*link_w[0];
					for(int i=1; i<n; i++) {
						line = br_maps[i].readLine();
						s = line.split("\\s+");
						weightedLinkCounts += Double.parseDouble(s[8])*link_w[i];
					}
					if(weightedLinkCounts<1.0) {
						link_bw.write("N/A"+"\t"+
								s[6]+"\t"+
								s[7]+"\t"+
								df3.format(weightedLinkCounts)+"\t"+
								s[5]+"\n");
						scaff = "Scaffold"+String.format("%08d", ++scaff_i);
					} else {
						link_bw.write(scaff+"\t"+
								s[6]+"\t"+
								s[7]+"\t"+
								df3.format(weightedLinkCounts)+"\t"+
								s[5]+"\n");
					}
				} else {
					for(int i=1; i<n; i++) br_maps[i].readLine();
					link_bw.write(scaff+"\t"+
							s[6]+"\t"+
							s[7]+"\t"+
							"-1.0"+"\t"+
							s[5]+"\n");
				}
			}

			for(int i=0; i<n; i++) br_maps[i].close();
				
			link_bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private String strcat(int[] is, String sep) {
		// TODO Auto-generated method stub
		StringBuilder str = new StringBuilder(""+is[0]);
		for(int i=1; i<is.length; i++) {
			str.append(sep);
			str.append(is[i]);
		}
		return str.toString();
	}

	private String strcat(double[] db, String sep) {
		// TODO Auto-generated method stub
		StringBuilder str = new StringBuilder(""+db[0]);
		for(int i=1; i<db.length; i++) {
			str.append(sep);
			str.append(db[i]);
		}
		return str.toString();
	}

	private Map<String, Sequence> parseContigFromBam() {
		// TODO Auto-generated method stub
		
		final Map<String, Sequence> contig_list = new HashMap<String, Sequence>();
		final SamReader in1 = factory.open(new File(this.bam_list[0]));
		try {
			List<SAMSequenceRecord> seqs = 
					in1.getFileHeader().getSequenceDictionary().getSequences();
			for(SAMSequenceRecord seq : seqs) 
				contig_list.put( seq.getSequenceName(), 
						new Sequence(seq.getSequenceIndex(), seq.getSequenceLength()) );
			in1.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return contig_list;
	}
	
	public Map<String, Sequence> parseContig(final String contig_index) {

		final Map<String, Sequence> contig_list = new HashMap<String, Sequence>();
		try {
			BufferedReader ctg_br = new BufferedReader(new FileReader(contig_index));
			String line;
			String[] s;
			int index = 0;
			while( (line=ctg_br.readLine())!=null ) {
				s = line.split("\\s+");
				contig_list.put( s[1].split(":")[1], 
						new Sequence(index, Integer.parseInt(s[2].split(":")[1])) );
				index++;
			}
			ctg_br.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
		return contig_list;
		
	}
	
	public List<List<Segment>> parseSegment(final String segment) {
		
		final List<List<Segment>> segment_list = new ArrayList<List<Segment>>();
		try {
			BufferedReader seg_br = new BufferedReader(new FileReader(segment));
			String line = seg_br.readLine();
			String[] s = line.split("\\s+");
			String chr = s[5];
			while( line!=null ) {
				final List<Segment> segment_chr = new ArrayList<Segment>();
				segment_chr.add( new Segment(
						s[0],
						Integer.parseInt(s[1]),
						Integer.parseInt(s[2]),
						Integer.parseInt(s[3]),
						s[4],
						s[5],
						Integer.parseInt(s[6]),
						Integer.parseInt(s[7])) );
				while( (line=seg_br.readLine())!=null ) {
					s = line.split("\\s+");
					if(chr.equals(s[5])) {
						segment_chr.add( new Segment(
								s[0],
								Integer.parseInt(s[1]),
								Integer.parseInt(s[2]),
								Integer.parseInt(s[3]),
								s[4],
								s[5],
								Integer.parseInt(s[6]),
								Integer.parseInt(s[7])) );
					} else {
						chr = s[5];
						break;
					}
				}
				segment_list.add(segment_chr);
			}
			seg_br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return segment_list;
	}
	
	private final static int lineWidth = 50;
	
	private final void parseScaffold(final String contig_fa,
			final String agp_map,
			final int minL,
			final String out_fa) {
		final Map<String, Sequence> contig_map = Sequence.parseFastaFileAsMap(contig_fa);
		
		try {
			BufferedReader br_agp = Utils.getBufferedReader(agp_map);
			final List<Sequence> scaffolds = new ArrayList<Sequence>();
			final List<Segment> segments = new ArrayList<Segment>();
			final StringBuilder str_buf = new StringBuilder();
			String line = br_agp.readLine();
			String[] s;
			String seq_sn, ref_sn, anchor_start, anchor_end;
			Sequence contig;
			final Set<String> anchored_seqs = new HashSet<String>();
			while(line!=null) {
				s = line.split("\\s+");
				seq_sn = s[5];
				anchor_start = s[6];
				segments.clear();
				
				segments.add( new Segment(
						s[0],
						Integer.parseInt(s[1]),
						Integer.parseInt(s[2]),
						Integer.parseInt(s[3]),
						s[4],
						s[5],
						Integer.parseInt(s[6]),
						Integer.parseInt(s[7])) );
				
            	while( (line=br_agp.readLine())!=null ) {
            		s = line.split("\\s+");
            		if(s[5].equals(seq_sn))
            			segments.add( new Segment(
        						s[0],
        						Integer.parseInt(s[1]),
        						Integer.parseInt(s[2]),
        						Integer.parseInt(s[3]),
        						s[4],
        						s[5],
        						Integer.parseInt(s[6]),
        						Integer.parseInt(s[7])) );
            		else break;
            	}
            	
            	anchor_end = s[7];
				ref_sn = s.length>9 ? s[9]+":"+anchor_start+"_"+anchor_end : null;
            	
            	str_buf.setLength(0);
            	for(Segment seg : segments) {
            		if(seg.type==MAP_ENUM.GAP) {
            			str_buf.append(Sequence.polyN(seg.seq_ln));
            		} else{
            			contig = contig_map.get(seg.seq_sn);
            			str_buf.append(seg.seq_rev ? 
            					Sequence.revCompSeq(contig.seq_str().substring(seg.seq_start, seg.seq_end)) : 
            					contig.seq_str().substring(seg.seq_start, seg.seq_end));
            			anchored_seqs.add(seg.seq_sn);
            		}
            	}
            	scaffolds.add(new Sequence(ref_sn, str_buf.toString().replaceAll("N{1,}$", "").replaceAll("^N{1,}", "") ));
            }
			br_agp.close();
			
			for(String seq : anchored_seqs) contig_map.remove(seq);
			
			BufferedWriter bw_fa = Utils.getBufferedWriter(out_fa);
			
			// whether to add unplaced contigs ?
			// alternatively can refer to the map file and do it later
			System.err.println("Buffer contigs...");
			for(Map.Entry<String, Sequence> entry : contig_map.entrySet())
				scaffolds.add(entry.getValue());
			
			System.err.println("Sort scaffolds...");
			Collections.sort(scaffolds);
			Collections.reverse(scaffolds);
			
			int n = scaffolds.size(), rm = 0, seq_ln;
			String seq_str;
			
			System.err.println(n+" scaffolds to parse...");
			for(int i=0; i!=n; i++) {
				if(i%10000==0) System.err.println(i+" scaffolds parsed...");
				
				contig = scaffolds.get(i);
				if(contig.seq_ln()<minL || contig.seq_ln()==0) 
					break;
				
				rm++;
				
				//bw_fa.write(contig.seq_str().contains("N") ? ">S" : ">C");
				bw_fa.write(">S");
				bw_fa.write(String.format("%08d", i+1));
				if(contig.seq_sn()!=null) bw_fa.write("|"+contig.seq_sn());
				bw_fa.write("\n");
				
				/***
				// not too slow as string buffer insertion is O(n)
				str_buf.setLength(0);
				str_buf.append(contig.seq_str);
				
				int offset=0, seq_ln=str_buf.length();
				while(offset<seq_ln) {
					offset += lineWidth;
					str_buf.insert( Math.min(offset, seq_ln), "\n" );
					seq_ln++;
					offset++;
				}
				
				bw_fa.write(str_buf.toString());
				**/
				
				seq_str = contig.seq_str();
				
				seq_ln = contig.seq_ln();
				bw_fa.write(seq_str.charAt(0));
				for(int j = 1; j!=seq_ln; j++) {
					if(j%lineWidth==0) bw_fa.write("\n");
					bw_fa.write(seq_str.charAt(j));
				}
				bw_fa.write("\n");
			}
			bw_fa.close();
			
			System.err.println("Calculate statistics...");
			// for(int i=rm; i!=n; i++) scaffolds.remove(i);
			scaffolds.subList(rm, scaffolds.size()).clear();
			int gap_ln = 0;
			seq_ln = 0;
			for(int i=0; i!=rm; i++) {
				seq_ln += scaffolds.get(i).seq_ln();
				gap_ln += StringUtils.countMatches(scaffolds.get(i).seq_str(), "N");
			}
			System.err.println("# SEQ_LEN\t"+seq_ln);
			System.err.println("# GAP_LEN\t"+gap_ln);
			System.err.println("#     N50\t"+calcWeightedMedianSatistic(scaffolds, 0.5));
		
		} catch (IOException e) {
            e.printStackTrace();
        }
	}
	
    public int calcWeightedMedianSatistic (List<Sequence> contigSortedAscending, double p) {
    	
    	if(contigSortedAscending.isEmpty()) return -1;
    	
    	int L = 0;
        for(int i=0; i<contigSortedAscending.size(); i++) {
            L += contigSortedAscending.get(i).seq_ln();
        }

        int N = (int) Math.ceil(p*L);
        int C = 0, B = 0, R = 0;
        while(C<N) {
            R = N-C;
            C+=contigSortedAscending.get(B++).seq_ln();
        }
        B--;

        if(L%2==0 && R==0) {
            return (contigSortedAscending.get(B).seq_ln()+contigSortedAscending.get(B-1).seq_ln())/2;
        } else if (L%2!=0 && R==0) {
            return contigSortedAscending.get(B-1).seq_ln();
        } else {
            return contigSortedAscending.get(B).seq_ln();
        }
    }
	
	private enum MAP_ENUM {CONTIG, GAP}

	private class Segment {
		private final MAP_ENUM type;
		
		private final String seq_sn;
		private final int seq_ln;
		private final int seq_start;
		private final int seq_end;
		private final boolean seq_rev;
		private final String mol_sn;
		private final int mol_start;
		private final int mol_end;
	
		public Segment(String seq_sn, 
				int seq_ln, 
				int seq_start, 
				int seq_end, 
				String seq_rev, 
				String mol_sn,
				int mol_start, 
				int mol_end) {
			// TODO Auto-generated constructor stub
			this.type = seq_sn.equals("GAP") ? MAP_ENUM.GAP : MAP_ENUM.CONTIG;
			this.seq_sn = seq_sn;
			this.seq_ln = seq_ln;
			this.seq_start = seq_start;
			this.seq_end = seq_end;
			this.seq_rev = seq_rev.equals("-");
			this.mol_sn = mol_sn;
			this.mol_start = mol_start;
			this.mol_end = mol_end;
		}
	}
	
	
	private abstract class LinkCounter implements Runnable {

		protected final String bam_in;
		protected final String out_prefix;
		protected final Map<String, RangeSet<Integer>> pseudo_gap = 
				new HashMap<String, RangeSet<Integer>>();
		protected final Map<String, Map<Range<Integer>, Integer>> link_count = 
				new HashMap<String, Map<Range<Integer>, Integer>>();
		
		public LinkCounter(String bam_in, 
				String map_file,
				String out_prefix) {
			// TODO Auto-generated constructor stub
			this.bam_in = bam_in;
			this.out_prefix = out_prefix;
			this.parseMapFile(map_file);
		}
		
		private void parseMapFile(String map_file) {
			// TODO Auto-generated method stub
			try {
				BufferedReader br_map = Utils.getBufferedReader(map_file);
				String line;
				String[] s;
				while( (line=br_map.readLine())!=null ) {
					if(!line.startsWith("GAP")) continue;
					s= line.split("\\s+");
					if(!pseudo_gap.containsKey(s[5])) {
						pseudo_gap.put(s[5], TreeRangeSet.create());
						link_count.put(s[5], new HashMap<Range<Integer>, Integer>());
					}
					Range<Integer> range = Range.closedOpen(Integer.parseInt(s[6]), 
							Integer.parseInt(s[7]));
					pseudo_gap.get(s[5]).add(range);
					link_count.get(s[5]).put(range, 0);
				}
				br_map.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub
			try {
				this.buildReadPairGraph(bam_in, out_prefix);
				this.writeLinkCount();
			} catch (Exception e) {
				// TODO Auto-generated catch block
				Thread t = Thread.currentThread();
				t.getUncaughtExceptionHandler().uncaughtException(t, e);
				e.printStackTrace();
				executor.shutdown();
				System.exit(1);
			}
		}
		
		private void writeLinkCount() {
			// TODO Auto-generated method stub
			try {
				BufferedWriter bw = Utils.getBufferedWriter(this.out_prefix);
				BufferedReader br_map = Utils.getBufferedReader(map_file);
				
				String line;
				String[] s;
				while( (line=br_map.readLine())!=null ) {
					if(!line.startsWith("GAP")) {
						bw.write(line+"\t-1.0\n");
						continue;
					}
					s= line.split("\\s+");
					Range<Integer> range = Range.closedOpen(Integer.parseInt(s[6]), 
							Integer.parseInt(s[7]));
					
					bw.write(line+"\t"+link_count.get(s[5]).get(range)+"\n");
				}
				br_map.close();
				bw.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		public abstract void buildReadPairGraph(final String bam_in,
				final String out);
	}
	
	private final class LongReadLinkCounter extends LinkCounter {
		
		private final Map<String, Set<Long>> map_segs;
		private final List<List<Segment>> list_segs; 
		
		public LongReadLinkCounter(String bam_in, 
				String map_file, 
				String out_prefix) {
			// TODO Auto-generated constructor stub
			super(bam_in, map_file, out_prefix);
			this.list_segs = parseSegment(map_file);
			map_segs = new HashMap<String, Set<Long>>();
			// 32 bits chromosome id + 32 bits positions on chromosome
			String seq_sn;
			
			for(int i=0; i<list_segs.size(); i++) {
				List<Segment> segs = list_segs.get(i);
				
				for(int j=0; j<segs.size(); j++) {
					seq_sn = segs.get(j).seq_sn;
					long key = i;
					key    <<= 32;
					key     += j;
					if(!map_segs.containsKey(seq_sn))
						map_segs.put(seq_sn, new HashSet<Long>());
					map_segs.get(seq_sn).add(key);
				}
			}
		}
		
		@Override
		public void buildReadPairGraph(final String bam_in, 
				final String out) {
			// TODO Auto-generated method stub
		
			try {	
				myLogger.info("Process file ... "+bam_in);
				
				final List<AlignmentSegment> record_list = new ArrayList<AlignmentSegment>();
				final Set<SAMRecord> record_buffer = new HashSet<SAMRecord>();
				AlignmentSegment record, primary_record, secondary_record;
				SAMRecord tmp_record;
				String sam_id;
				
				final SamReader in1 = factory.open(new File(bam_in));
				
				SAMRecordIterator iter1 = in1.iterator();
				long record_count = 1;
				tmp_record = iter1.next();
				int primary_ref;
				
				Range<Integer> span;
				String refSeq;
				RangeSet<Integer> intersect;
				Map<Range<Integer>, Integer> links;
				
				while( tmp_record!=null ) {
					
					record_buffer.clear();
					if( !tmp_record.getReadUnmappedFlag() &&
							!tmp_record.getNotPrimaryAlignmentFlag() ) 
						record_buffer.add(tmp_record);

					sam_id = tmp_record.getReadName();
					refSeq = tmp_record.getReferenceName();
					
					while( (tmp_record = iter1.hasNext() ? iter1.next() : null) !=null &&
							tmp_record.getReadName().equals(sam_id) ) {
						if(++record_count%1000000==0)
							myLogger.info(""+record_count+" ... "+bam_in);

						if(!tmp_record.getReadUnmappedFlag() &&
								!tmp_record.getNotPrimaryAlignmentFlag() ) 
							record_buffer.add(tmp_record);
					}
					
					primary_ref = -1;
					for(SAMRecord r : record_buffer)
						if(!r.getSupplementaryAlignmentFlag()) {
							primary_ref = r.getReferenceIndex();
							break;
						}
					
					record_list.clear();
					for(SAMRecord r : record_buffer)
						if(r.getReferenceIndex()==primary_ref) 
							record_list.add(samRecordToAlignmentRecord(r));
					
					// is empty
					if( record_list.isEmpty() ) continue;
					
					// merge collinear alignment records
					Collections.sort(record_list, new AlignmentSegment.SubjectCoordinationComparator());
					
					for(int i=0; i<record_list.size(); i++) {
						primary_record = record_list.get(i);
						for(int j=i+1; j<record_list.size(); j++) {
							secondary_record = record_list.get(j);
							double max_shift = collinear_shift*
									Math.min(primary_record.qlength(), secondary_record.qlength());
							if( (record=AlignmentSegment.collinear(primary_record, 
									secondary_record, max_shift))!=null ) {
								record_list.set(i, record);
								record_list.remove(j);
								--i;
								break;
							}
						}
					}
					
					span = Range.closed(0, 0);
					for(AlignmentSegment r : record_list) {
						if(r.slength()>span.upperEndpoint()-span.lowerEndpoint()+1)
							if(r.sstart()>r.send())
								span = Range.closed(r.send(), r.sstart());
							else
								span = Range.closed(r.sstart(), r.send());
					}
					
					intersect = pseudo_gap.get(refSeq).subRangeSet(span);
					links = link_count.get(refSeq);
					
					for(Range<Integer> range : intersect.asRanges())
						if(links.containsKey(range))
							links.put(range, links.get(range)+1);
					
				}
				
				iter1.close();
				in1.close();
				
				myLogger.info("Process file ... "+bam_in+" done.");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		}
		
		private AlignmentSegment samRecordToAlignmentRecord(SAMRecord samRecord) {
			// TODO Auto-generated method stub
			return new AlignmentSegment( samRecord.getReadName(),
					samRecord.getReferenceName(),
					samRecord.getReadPositionAtReferencePosition(samRecord.getAlignmentStart())-1,
					samRecord.getReadPositionAtReferencePosition(samRecord.getAlignmentEnd())-1,
					samRecord.getAlignmentStart()-1,
					samRecord.getAlignmentEnd()-1 );
		}

		private int readCoveredBase(SAMRecord samRecord) {
			// TODO Auto-generated method stub
			if(samRecord.getReadUnmappedFlag()) return 0;
			return samRecord.getReadPositionAtReferencePosition(samRecord.getAlignmentEnd())-
					samRecord.getReadPositionAtReferencePosition(samRecord.getAlignmentStart());
		}
	}
	
	private final class PairedReadLinkCounter extends LinkCounter {
		
		private final int ins_thres;
		// paired-end  -> <-
		private final boolean fr;
		// mate-pair   <- ->
		private final boolean rf;
		
		public PairedReadLinkCounter(String bam_in, 
				int ins_thres,
				String map_file,
				boolean fr,
				String out_prefix) {
			// TODO Auto-generated constructor stub
			super(bam_in, map_file, out_prefix);
			this.fr =  fr;
			this.rf = !fr;
			this.ins_thres = ins_thres;
		}
		
		@Override
		public void buildReadPairGraph(final String bam_in,
				final String out) {
			try {	
				myLogger.info("Process file ... "+bam_in);
				
				final SAMRecord[] record_pair = new SAMRecord[2];
				SAMRecord tmp_record;
				String sam_id;
				
				final SamReader in1 = factory.open(new File(bam_in));

				SAMRecordIterator iter1 = in1.iterator();
				long record_count = 1;
				tmp_record = iter1.next();
				
				Range<Integer> span;
				RangeSet<Integer> intersect;
				Map<Range<Integer>, Integer> links;
				String  refSeq;
				
				long badReadPair = 0, goodReadPair = 0;
				
				while( tmp_record!=null ) {

					Arrays.fill(record_pair, null);
					
					// skip supplementary and secondary alignments
					if( !tmp_record.getSupplementaryAlignmentFlag() &&
							!tmp_record.getNotPrimaryAlignmentFlag() ) {
						if(tmp_record.getFirstOfPairFlag()) 
							record_pair[0] = tmp_record;
						else record_pair[1] = tmp_record;
					}

					sam_id = tmp_record.getReadName();

					while( (tmp_record = iter1.hasNext() ? iter1.next() : null) !=null &&
							tmp_record.getReadName().equals(sam_id) ) {
						if(++record_count%1000000==0)
							myLogger.info(""+record_count+" ... "+bam_in);

						if( !tmp_record.getSupplementaryAlignmentFlag() &&
								!tmp_record.getNotPrimaryAlignmentFlag() ) {
							if(tmp_record.getFirstOfPairFlag()) 
								record_pair[0] = tmp_record;
							else record_pair[1] = tmp_record;
						}
					}
					
					if(record_pair[0]==null||record_pair[1]==null) 
						throw new RuntimeException("!!!");
					
					if( (span=parseReadPairSpan(record_pair))==null) {
						++badReadPair;
						continue;
					}
					++goodReadPair;
					refSeq = record_pair[0].getReferenceName();
					intersect = pseudo_gap.get(refSeq).subRangeSet(span);
					links = link_count.get(refSeq);
					
					for(Range<Integer> range : intersect.asRanges())
						if(links.containsKey(range))
							links.put(range, links.get(range)+1);
				}
				iter1.close();
				in1.close();
				
				long total_pair = goodReadPair+badReadPair;
				myLogger.info("Process file ... "+bam_in+" done.");
				myLogger.info("Total read pairs, "+total_pair+"; "
						+ "good read pairs, "+goodReadPair+"/"+total_pair+"("+(double)goodReadPair/(total_pair)+"); "
						+ "bad read pairs, "+badReadPair+"/"+total_pair+"("+(double)badReadPair/(total_pair)+")\n");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		private Range<Integer> parseReadPairSpan(SAMRecord[] record_pair) {
			// TODO Auto-generated method stub
			// ensure both mapped
			if( record_pair[0].getReadUnmappedFlag() || 
					record_pair[1].getReadUnmappedFlag() ) return null;
			// ensure mapped to the same reference chromosome
			if( record_pair[0].getReferenceIndex() !=
					record_pair[1].getReferenceIndex() ) return null;
			// ensure one forward and one reverse mapped
			if( record_pair[0].getReadNegativeStrandFlag() == 
					record_pair[1].getReadNegativeStrandFlag() ) return null;
			// swap read pair if mapping position inversed
			if( record_pair[0].getAlignmentStart()>record_pair[1].getAlignmentStart()) {
				SAMRecord record = record_pair[0];
				record_pair[0] = record_pair[1];
				record_pair[1] = record;
			}
			// ensure fr mapping for paired-end libraries
			// ensure rf mapping for mate-pair libraries
			if( record_pair[0].getReadNegativeStrandFlag() && fr ||
					record_pair[1].getReadNegativeStrandFlag() && rf) 
				return null;
			// ensure insert size no greater that the predefined threshold
			if( Math.abs(record_pair[0].getInferredInsertSize())>ins_thres) 
				return null;
			return Range.closed(record_pair[0].getAlignmentStart()-1, 
					record_pair[1].getAlignmentEnd()-1);
		}
	}
}
