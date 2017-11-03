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

import cz1.ngs.model.Sequence;
import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Consensus extends Executor {

	private final static String bam_lib = "--b([0-9]+)";
	private final static String ins_lib = "--f([0-9]+)";
	private final static String w_lib = "--w([0-9]+)";
	
	private final static Pattern bam_pat = Pattern.compile(bam_lib);
	private final static Pattern ins_pat = Pattern.compile(ins_lib);
	private final static Pattern w_pat = Pattern.compile(w_lib);
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " --b<#>                  Input BAM file for this library (<#> = 1,2,...).\n"
						+ " --f<#>                  Fragment/insert size threshold of this library (<#> = 1,2,...).\n"
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
						+ " -c/--contig             The FASTA file contain contigs. \n"
						+ " -m/--map                Map file indicate the order of the contigs. \n"
						+ " -bs/--batch-size        Batch size store in memory (default 4000000). \n"
						+ "                         Reduce this number if run out of memory.\n"
						+ " -t/--threads            Threads to use (default 1). \n"
						+ "                         The maximum number of threads to use will be the number of BAM files.\n"
						+ " -l/--min-size           Minimum size of scaffolds to output (default 0). \n"
						+ " -o/--out                Prefix of the output FASTQ file.\n"
						+ "\n");
	}

	private String contig_file = null;
	private String map_file = null;
	private int batch_size = 4000000;
	private int min_size = 0;
	private String out_prefix = null;
	private String[] bam_list = null;
	private int[] ins_thres = null;
	private double[] link_w = null;
	
	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-c", "--contig", true);
			myArgsEngine.add("-m", "--map", true);
			myArgsEngine.add("-bs", "--batch-size", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-l", "--min-size", true);
			myArgsEngine.add("-o", "--out", true);
			myArgsEngine.addWildOptions(bam_lib, true);
			myArgsEngine.addWildOptions(ins_lib, true);
			myArgsEngine.addWildOptions(w_lib, true);
			myArgsEngine.parse(args);
		}
		
		if (myArgsEngine.getBoolean("-c")) {
			this.contig_file = myArgsEngine.getString("-c");
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
		
		if (myArgsEngine.getBoolean("-bs")) {
			this.batch_size = Integer.parseInt(myArgsEngine.getString("-bs"));
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
		
		this.parseDataLibrary(args);
	}

	private void parseDataLibrary(String[] args) {
		// TODO Auto-generated method stub
		final Map<Integer, String> bamLib = new HashMap<Integer, String>();
		final Map<Integer, Integer> insLib = new HashMap<Integer, Integer>();
		final Map<Integer, Integer> wLib = new HashMap<Integer, Integer>();

		for(int i=0; i<args.length; i++) {
			String arg = args[i];
			if(arg.matches(bam_lib)) {
				Matcher m = bam_pat.matcher(arg);
				m.find();
				int lib = Integer.parseInt(m.group(1));
				bamLib.put(lib, args[++i]);
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

		int n = bamLib.size();
		if(insLib.size()!=n || wLib.size()!=n) {
			printUsage();
			throw new IllegalArgumentException("Please library parameters do not match!!!");
		}
		
		this.bam_list = new String[n];
		this.ins_thres = new int[n];
		this.link_w = new double[n];
		int i=0;
		for(Integer key : bamLib.keySet()) {
			if(!insLib.containsKey(key) || !wLib.containsKey(key)) {
				printUsage();
				throw new IllegalArgumentException("Please library parameters do not match!!!");
			}
			bam_list[i] = bamLib.get(key);
			ins_thres[i] = insLib.get(key);
			link_w[i] = 1.0/wLib.get(key);
			i++;
		}
	}
	
	private String[] link_file = null;
	@Override
	public void run() {
		// TODO Auto-generated method stub
		// STEP 1. count links
		//         merge links
		myLogger.info("STEP 1. count links");
		link_file = new String[this.bam_list.length];
		this.initial_thread_pool();
		for(int i=0; i<this.bam_list.length; i++) {
			String out = this.out_prefix+"_"+
					new File(bam_list[i]).getName().replaceAll(".bam$", "")+
					RandomStringUtils.randomAlphanumeric(20).toUpperCase()+
					".mergedLink.txt";
			link_file[i] = out;
			this.executor.submit(new LinkCounter(this.bam_list[i], out));
		}
		this.waitFor();
		
		// STEP 2. parse links
		myLogger.info("STEP 2. parse links");
		String link_out = this.out_prefix+"_"+
				RandomStringUtils.randomAlphanumeric(20).toUpperCase()+
				".parsedLink";
		this.parseLink(link_file, this.map_file, link_out);
		
		// STEP 3. parse scaffolds
		myLogger.info("STEP 3. parse scaffolds");
		String scaff_out = this.out_prefix+"_parsedScaffold.fa";
		
		this.parseScaffold(this.contig_file, 
				link_out+".map"+String.format("%04d", this.link_file.length), 
				this.min_size, 
				scaff_out);
	}

	private final static Object lock = new Object();
	private final static SamReaderFactory factory =
			SamReaderFactory.makeDefault()
			.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
					SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
			.validationStringency(ValidationStringency.SILENT);
	private final static int refMask = Integer.MAX_VALUE; // 31 bits reference index
	private final static int dirMask = 1;                 // 1  bit  direction
	
	/***
	public void buildReadPairGraph(String bam_in,
			String out) {
		try {	
			String[] bam_files = bam_in.trim().split(";");
			final Map<Long, Integer> links = new HashMap<Long, Integer>();
			
			this.initial_thread_pool();
			
			for(String bam_file : bam_files) {

				this.executor.submit(new Runnable() {
					private String bam_file;
					
					@Override
					public void run() {
						System.err.println("Process file ... "+bam_file);
						
						final Map<Long, Integer> link = new HashMap<Long, Integer>(); 
						
						final List<SAMRecord> sam_record_list = new ArrayList<SAMRecord>();
						final List<SAMRecord> mat_record_list = new ArrayList<SAMRecord>();
						
						SAMRecord tmp_record;
						String sam_id;
						int sam_ref, mat_ref;
						long key_ref;
						
						final SamReader in1 = factory.open(new File(bam_file));
						
						SAMRecordIterator iter1 = in1.iterator();
						tmp_record = iter1.next();
						while( tmp_record!=null ) {
							
							sam_record_list.clear();
							mat_record_list.clear();

							if(!tmp_record.getReadUnmappedFlag()) {
								if(tmp_record.getFirstOfPairFlag()) 
									sam_record_list.add(tmp_record);
								else mat_record_list.add(tmp_record);
							}
							
							sam_id = tmp_record.getReadName();

							while( (tmp_record = iter1.hasNext() ? iter1.next() : null) !=null &&
									tmp_record.getReadName().equals(sam_id) ) {
								if(!tmp_record.getReadUnmappedFlag()) {
									if(tmp_record.getFirstOfPairFlag()) 
										sam_record_list.add(tmp_record);
									else mat_record_list.add(tmp_record);
								}
							}

							if(sam_record_list.isEmpty()||mat_record_list.isEmpty()) continue;
							
							for(SAMRecord sam_record : sam_record_list) {
								sam_ref = sam_record.getReferenceIndex();
								
								for(SAMRecord mat_record : mat_record_list) {
									mat_ref = mat_record.getReferenceIndex();
									
									if( sam_ref==mat_ref ) continue;

									if(sam_ref>mat_ref) {
										int tmp_int = sam_ref;
										sam_ref = mat_ref;
										mat_ref = tmp_int;
									}
									
									key_ref = sam_ref;
									key_ref <<= 32;
									key_ref += mat_ref;
									
									
									if(link.containsKey(key_ref))
										link.put(key_ref, link.get(key_ref)+1);
									else link.put(key_ref, 1);
								}
							}
						}
						iter1.close();
						
						synchronized(lock) {
							for(Map.Entry<Long, Integer> entry : link.entrySet()) {
								key_ref = entry.getKey();

								if(links.containsKey(key_ref))
									links.put(key_ref, 
											links.get(key_ref)+link.get(key_ref));
								else links.put(key_ref, link.get(key_ref));
							}
						}
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
			
			int mat_ref, sam_ref;
			long key_ref;
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(out)));
			for(Map.Entry<Long, Integer> entry : links.entrySet()) {
				key_ref = entry.getKey();
				mat_ref = (int) (key_ref);
				key_ref >>= 32;
				sam_ref = (int) (key_ref);
				
				bw.write(sam_ref+"\t"+mat_ref+"\t"+entry.getValue()+"\n");
			}
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	**/
	
	
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
	private final static DecimalFormat df2 = new DecimalFormat("#.00");
	
	public void parseLink(final String[] links_in,
			final String segment_in,
			final String link_out) {
		
		Map<String, Sequence> contig_list = parseContigFromBam();
		List<List<Segment>> segment_list = parseSegment(segment_in);
		Map<Integer, Set<Integer>> contig_seg_map = new HashMap<Integer, Set<Integer>>();
		int chr_index = 0, seg_index, seq_index;
		int seg_key; // 8  bits pseudo molecule indices
					 // 24 bits segment indices 
		// List<String> pseudo_mol_sn = new ArrayList<String>();
		List<double[]> pseudo_mol_sf = new ArrayList<double[]>();
		// pseudo_mol_sn.add("");
		for(List<Segment> seg_list : segment_list) {
			// pseudo_mol_sn.add("Chr"+(chr_index<10?"0":"")+chr_index);
			seg_index = 0;
			for(Segment  seg : seg_list) {
				if(seg.type==MAP_ENUM.CONTIG) {
					seg_key = chr_index;
					seg_key <<= 24;
					seg_key += seg_index;
					seq_index = contig_list.get(seg.seq_sn).seq_no();
					if(!contig_seg_map.containsKey(seq_index))
						contig_seg_map.put(seq_index, new HashSet<Integer>());
					contig_seg_map.get(seq_index).add(seg_key);	
				}
				seg_index++;
			}
			pseudo_mol_sf.add(new double[seg_index-1]);
			chr_index++;
		}
		myLogger.info("Segments: "+segment_list.size());
		myLogger.info("SegMap:   "+contig_seg_map.size());
		
		//for(Map.Entry<Integer, Integer> entry : contig_seg_map.entrySet()) 
		//	System.err.println(entry.getKey()+" "+entry.getValue());
			
		Map<Long, Integer> links = new HashMap<Long, Integer>();
		
		myLogger.info("Insert size threshold: "+strcat(this.ins_thres,","));
		myLogger.info("Link weights: "+strcat(link_w,","));
		
		long key;
		int val;
		try {
			for(int i=0; i<links_in.length; i++) {
				links.clear();
				System.gc ();
				System.runFinalization ();
				
				String link_i = links_in[i];
				BufferedReader link_br = new BufferedReader(new FileReader(link_i));
				String line;
				String[] s;
				Set<Integer> C1, C2;
				int z1, z2;
				while( (line=link_br.readLine())!=null ) {
					s = line.split("\\s+");
					val = Integer.parseInt(s[2]);
					
                    z1 = Integer.parseInt(s[0]);
                    z2 = Integer.parseInt(s[1]);
                    if( !contig_seg_map.containsKey(z1) ||
                            !contig_seg_map.containsKey(z2) )
                        continue;
                    C1 = contig_seg_map.get(z1);
                    C2 = contig_seg_map.get(z2);

                    for(Integer c1 : C1) {
                    	for(Integer c2 : C2) {
                    		if(c1<c2) {
                                key = c1;
                                key <<= 32;
                                key += c2;
                            } else {
                            	key = c2;
                                key <<= 32;
                                key += c1;
                            }
                            if(links.containsKey(key)) {
                            	link_br.close();
                            	throw new RuntimeException("!!!");
                            	// links.put(key, val+links.get(key));
                            } else {
                                links.put(key, val);
                            }
                    	}
                    }
				}

				link_br.close();
				
				double insz_ub = ins_thres[i];
				double w = link_w[i];
				List<Segment> seg_list;
				double[] mol_sf;
				
				BufferedWriter link_bw = new BufferedWriter(new FileWriter(link_out+
                            ".link"+String.format("%04d", i+1)));
				int n1, n2, ins;
				Segment seg;
				
				outerloop:
					for(Map.Entry<Long, Integer> entry : links.entrySet()) {
						key = entry.getKey();
						n2 = (int) (key&mask_24bits);
						key >>= 24;
	                    z2 = (int) (key&mask_08bits);
	                    key >>=  8;
						n1 = (int) (key&mask_24bits);
						key >>= 24;
						z1 = (int) key;
						val = entry.getValue();
						link_bw.write(n1+" "+n2+" "+z1+" "+z2+" "+val+"\n");
						
						// if(c1!=c2 || val<link_thres) continue;
						if(z1!=z2) continue;
						// System.out.println(n1+" "+n2+" "+c1+" "+c2+" "+val+"\n");
                        
						seg_list = segment_list.get(z1);
						mol_sf = pseudo_mol_sf.get(z1);
						
                        ins = 0;
						if(n1>n2) throw new RuntimeException("!!!");
						for(int j=n1+1; j<n2; j++) {
							seg = seg_list.get(j);
							if(seg.type==MAP_ENUM.GAP)
								continue;
							ins += seg.seq_ln;
							if(ins>insz_ub) continue outerloop;
						}
						
						for(int j=n1; j<n2; j++) mol_sf[j] += w*val;
					}
				link_bw.close();
				
				BufferedWriter map_bw = new BufferedWriter(new FileWriter(link_out+".map"+
						String.format("%04d", i+1)));
				int scaff_i = 0;
				
				for(int j=0; j!=pseudo_mol_sf.size(); j++) {
					mol_sf = pseudo_mol_sf.get(j);
					seg_list = segment_list.get(j);
					if(mol_sf.length!=seg_list.size()-1) {
						map_bw.close();
						throw new RuntimeException("!!!");
					}
					String scaff = "Scaffold"+String.format("%08d", ++scaff_i);
					seg = seg_list.get(0);
					
					map_bw.write(seg.seq_sn+"\t"+
							seg.seq_ln+"\t"+
							seg.seq_start+"\t"+
							seg.seq_end+"\t"+
							(seg.seq_rev?"-":"+")+"\t"+
							scaff+"\t"+
							seg.mol_start+"\t"+
							seg.mol_end+"\t"+
							df2.format(-1)+"\n" );
					
					int n = seg_list.size();
					for(int k=1; k!=n; k++) {
						seg = seg_list.get(k);
						if(mol_sf[k-1]<1.0) 
							scaff = "Scaffold"+String.format("%08d", ++scaff_i);
						map_bw.write(seg.seq_sn+"\t"+
								seg.seq_ln+"\t"+
								seg.seq_start+"\t"+
								seg.seq_end+"\t"+
								(seg.seq_rev?"-":"+")+"\t"+
								scaff+"\t"+
								seg.mol_start+"\t"+
								seg.mol_end+"\t"+
								df2.format(mol_sf[k-1])+"\n" );
					}
				}
				map_bw.close();
			}
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
			String seq_sn;
			Sequence contig;
			final Set<String> anchored_seqs = new HashSet<String>();
			while(line!=null) {
				s = line.split("\\s+");
				seq_sn = s[5];
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
            	scaffolds.add(new Sequence( str_buf.toString().replaceAll("N{1,}$", "").replaceAll("^N{1,}", "") ));
            }
			br_agp.close();
			
			for(String seq : anchored_seqs) contig_map.remove(seq);
			
			BufferedWriter bw_fa = Utils.getBufferedWriter(out_fa);
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
				
				bw_fa.write(contig.seq_str().contains("N") ? ">S" : ">C");
				bw_fa.write(String.format("%08d\n", i+1));
				
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
	
	private final class LinkCounter implements Runnable {

		private final String bam_in;
		private final String out_prefix;
		
		public LinkCounter(String bam_in, 
				String out_prefix) {
			// TODO Auto-generated constructor stub
			this.bam_in = bam_in;
			this.out_prefix = out_prefix;
		}

		private final Set<String> tmpLinkFile = new HashSet<String>();
		
		@Override
		public void run() {
			// TODO Auto-generated method stub
			this.buildReadPairGraph(bam_in, out_prefix, true, batch_size);
			this.mergeLinks(tmpLinkFile.toArray(new String[tmpLinkFile.size()]), 
					out_prefix);
			for(String tmpf : tmpLinkFile) new File(tmpf).delete();
		}
		
		public void buildReadPairGraph(final String bam_in,
				final String out,
				final boolean ignore_sa, // ignore secondary alignment
				final int batch_size) {
			buildReadPairGraph(bam_in, out, ignore_sa, batch_size, 0, null);
		}
		
		public void buildReadPairGraph(final String bam_in,
				final String out,
				final boolean ignore_sa, // ignore secondary alignment
				final int batch_size,
				final int batch_start,
				final String resume_r_name) {
			try {	
				myLogger.info("Process file ... "+bam_in);
				
				final Map<Long, Integer> links = new HashMap<Long, Integer>();
				final List<SAMRecord> sam_record_list = new ArrayList<SAMRecord>();
				final List<SAMRecord> mat_record_list = new ArrayList<SAMRecord>();
				final List<SAMRecord> sam_sec_sa_list = new ArrayList<SAMRecord>();
				final List<SAMRecord> mat_sec_sa_list = new ArrayList<SAMRecord>();

				SAMRecord tmp_record;
				String sam_id;
				int sam_ref, mat_ref, batch = batch_start;
				long key_ref;
				
				final SamReader in1 = factory.open(new File(bam_in));

				SAMRecordIterator iter1 = in1.iterator();
				long record_count = 1;
				tmp_record = iter1.next();
				
				if(resume_r_name!=null) {
					myLogger.info("Resuming...");
					while(!tmp_record.getReadName().equals(resume_r_name)) {
						tmp_record = iter1.next();
						record_count++;
					}
					myLogger.info("Resuming... "+record_count+"th record: "+tmp_record.getSAMString());
				}
				
				boolean _ignore_sa = !ignore_sa;
				while( tmp_record!=null ) {

					sam_record_list.clear();
					mat_record_list.clear();
					sam_sec_sa_list.clear();
					mat_sec_sa_list.clear();

					if(!tmp_record.getReadUnmappedFlag()) {
						if(tmp_record.getFirstOfPairFlag()) 
							sam_record_list.add(tmp_record);
						else mat_record_list.add(tmp_record);
					}

					sam_id = tmp_record.getReadName();

					while( (tmp_record = iter1.hasNext() ? iter1.next() : null) !=null &&
							tmp_record.getReadName().equals(sam_id) ) {
						if(++record_count%1000000==0)
							myLogger.info(""+record_count+" ... "+bam_in);

						if(!tmp_record.getReadUnmappedFlag()) {
							if(tmp_record.getFirstOfPairFlag()) { 
								if(tmp_record.getNotPrimaryAlignmentFlag()) 
									sam_sec_sa_list.add(tmp_record);
								else sam_record_list.add(tmp_record);
							} else {
								if(tmp_record.getNotPrimaryAlignmentFlag())
									mat_sec_sa_list.add(tmp_record);
								else mat_record_list.add(tmp_record);
							}
						}
					}

					if(_ignore_sa) {
						sam_record_list.addAll(sam_sec_sa_list);
						mat_record_list.addAll(mat_sec_sa_list);
					}
					
					if(sam_record_list.isEmpty()||mat_record_list.isEmpty()) continue;

					for(SAMRecord sam_record : sam_record_list) {
						sam_ref = sam_record.getReferenceIndex();

						for(SAMRecord mat_record : mat_record_list) {
							mat_ref = mat_record.getReferenceIndex();

							if( sam_ref==mat_ref ) continue;

							if(sam_ref>mat_ref) {
								int tmp_int = sam_ref;
								sam_ref = mat_ref;
								mat_ref = tmp_int;
							}

							key_ref = sam_ref;
							key_ref <<= 32;
							key_ref += mat_ref;

							if(links.containsKey(key_ref))
								links.put(key_ref, links.get(key_ref)+1);
							else links.put(key_ref, 1);
						}
					}
					
					if(links.size()>batch_size) {
						String outf = out+".tmp"+String.format("%016d", batch);
						
						BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outf)));
						SortedSet<Long> keys = new TreeSet<Long>(links.keySet());
						for(Long key : keys) 
							bw.write(key+"\t"+links.get(key)+"\n");
						bw.close();
					    myLogger.info(out+".tmp"+String.format("%016d", batch)+" written with "+record_count+"th record: "+tmp_record.getSAMString());
	                    batch++;

	                    tmpLinkFile.add(outf);
	                    
	                    links.clear();
	                    System.gc ();
	    				System.runFinalization ();
	                }
				}
				iter1.close();
				in1.close();

				if(links.size()>0) {
					String outf = out+".tmp"+String.format("%016d", batch);
					
					BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outf)));
					SortedSet<Long> keys = new TreeSet<Long>(links.keySet());
					for(Long key : keys) 
						bw.write(key+"\t"+links.get(key)+"\n");
					bw.close();
					
					tmpLinkFile.add(outf);
				}
				
				myLogger.info("Process file ... "+bam_in+" done.");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		private void mergeLinks(final String[] link_in,
				final String link_out) {
			try {
				
				int nBatch = link_in.length;
				BufferedReader[] brLinkInTmp = new BufferedReader[nBatch];
				final boolean[] reachFileEnd = new boolean[nBatch];
				final int[] link_count = new int[nBatch];
				for(int i=0; i!=nBatch; i++) { 
					brLinkInTmp[i] = new BufferedReader(new FileReader(link_in[i]));
					reachFileEnd[i] = false;
					link_count[i] = 0;
				}
				
				final TreeMap<Long, Integer> treeMap = new TreeMap<Long, Integer>();
				String line;
				int nReachFileEnd = 0;
				
				for(int i=0; i!=nBatch; i++) {
					line = brLinkInTmp[i].readLine();
					while( line!=null&&!fillLinkTreeMap(treeMap, link_count, line, i) ) {
						line = brLinkInTmp[i].readLine();
					}
					if(line==null) {
						reachFileEnd[i] = true;
						nReachFileEnd++;
					}
				}
				Entry<Long, Integer> firstEntry;
				long link_key;
				int link_val;
				BufferedWriter bwLinkOut = new BufferedWriter(new FileWriter(link_out));
				while( !treeMap.isEmpty() ) {
					firstEntry = treeMap.pollFirstEntry();
					link_key = firstEntry.getKey();
					link_val = firstEntry.getValue();
					bwLinkOut.write( (int)(link_key>>32)+" "+(int)link_key+" "+link_count[link_val]+"\n" );
					
					if(!reachFileEnd[link_val]) {
						line = brLinkInTmp[link_val].readLine();
						while( line!=null&&!fillLinkTreeMap(treeMap, link_count, line, link_val) ) {
							line = brLinkInTmp[link_val].readLine();
						}
						if(line==null) {
							reachFileEnd[link_val] = true;
							nReachFileEnd++;
						}
					}
					
					if(treeMap.isEmpty()&&nReachFileEnd!=nBatch) {
						for(int i=0; i!=nBatch; i++) {
							if(!reachFileEnd[i]) {
								line = brLinkInTmp[i].readLine();
								while( line!=null&&!fillLinkTreeMap(treeMap, link_count, line, i) ) {
									line = brLinkInTmp[i].readLine();
								}
								if(line==null) {
									reachFileEnd[i] = true;
									nReachFileEnd++;
								}
							}
						}
					}
				}
				bwLinkOut.close();
				for(int i=0; i!=nBatch; i++) brLinkInTmp[i].close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		private boolean fillLinkTreeMap(final TreeMap<Long, Integer> treeMap, 
				final int[] link_count,
				final String line, 
				final int i) {
			// TODO Auto-generated method stub
			String[] parse_link = line.split("\\s+");
			long link_key = Long.parseLong(parse_link[0]); 
			int link_val = Integer.parseInt(parse_link[1]);
			if(treeMap.containsKey(link_key)) {
				link_count[treeMap.get(link_key)] += link_val;
				return false;
			} else {
				link_count[i] = link_val;
				treeMap.put(link_key, i);
				return true;
			}
		}
	}
}
