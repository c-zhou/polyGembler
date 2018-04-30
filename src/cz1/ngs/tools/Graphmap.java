package cz1.ngs.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.TreeMap;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;
import org.w3c.dom.stylesheets.LinkStyle;

import com.google.common.collect.DiscreteDomain;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

import cz1.ngs.model.AlignmentSegment;
import cz1.ngs.model.GFA;
import cz1.ngs.model.OverlapEdge;
import cz1.ngs.model.SAMSegment;
import cz1.ngs.model.Sequence;
import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Graphmap extends Executor {

	private final static boolean USE_OS_BUFFER = true;
	// 8Mb buffer size
	private final static int buffSize8Mb = 8388608;
	private final static Writer STD_OUT_BUFFER = USE_OS_BUFFER ? 
			new BufferedWriter(new OutputStreamWriter(System.out), buffSize8Mb) : new OutputStreamWriter(System.out);
	private static enum Task {hash, map, zzz};
	private static enum Library {pe, r454, long3};
	private Task task_list = Task.zzz;

	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		switch(this.task_list) {
		case zzz:
			myLogger.info(
					"\n\nChoose a task from listed:\n"
							+ " hash                    Create and save hash table.\n"
							+ " map                     Run graphmap to map sequences to an assembly graph.\n"
							+ "\n");
			break;
		case hash:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -s/--subject            The FASTA file contain subject/reference sequences. \n"
							+ " -k/--kmer-size          K-mer size (no greater than 16, default 12).\n"
							+ " -x/--max-mer-count      Maxmium mer count (no limit).\n"
							+ " -t/--threads            Threads to use (default 1). \n"
							+ "                         The maximum number of threads to use will be the number of BAM files.\n"
							+ " -o/--out                Prefix of the output files.\n"
							+ "\n");
			break;
		case map:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -s/--subject            The FASTA file contain subject/reference sequences. \n"
							+ " -q/--query              The FASTA/FASTQ/BAM file contain query sequences to map. \n"
							+ " -l/--library            The read library type (pe, 454 or long).\n"
							+ " -g/--graph              Assembly graph (GFA) format. Currently, the program only accept \n"
							+ "                         the assembly graph format used by the assembler SPAdes (de-bruijn \n"
							+ "                         graph) or CANU (overlap). For assembly graphs in other formats: \n"
							+ "                         1. for contigs generated with de-bruijn graph, please ignore this \n"
							+ "                         option and provide k-mer size, the program is able to rebuild the \n"
							+ "                         assembly graph from contigs; and, \n"
							+ "                         2. for contigs generated with overlapping algorithms, please convert \n"
							+ "                         assembly graph to the format used by CANU.\n"
							+ "                         NOTE: it is possible to run the program without an assembly graph, \n"
							+ "                         however, the assembly might be less accurate.\n"
							+ " -k/--kmer-size          K-mer size (no greater than 16, default 12).\n"
							+ " -H/--hash-table         Hash table. If provided will load hash table from the file instead of \n"
							+ "                         reconstruct it.\n"
							+ " -x/--max-mer-count      Maxmium mer count (no limit).\n"
							+ " -t/--threads            Threads to use (default 1). \n"
							+ "                         The maximum number of threads to use will be the number of BAM files.\n"
							+ " -d/--debug              Debugging mode will have extra information printed out.\n"
							+ " -dd/--debug-debug       Debugging mode will have more information printed out than -d mode.\n"
							+ " -o/--out                Prefix of the output files.\n"
							+ "\n");
			break;
		default:
			throw new RuntimeException("!!!");
		}
	}

	private String subject_file;
	private String graph_file;
	private String[] query_file;
	private String hash_file = null;
	private String out_prefix;
	private Library library = null;
	private int merK = 12;
	private int maxC = Integer.MAX_VALUE-1;
	private boolean debug = false;
	private boolean ddebug = false;

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}

		switch(args[0].toUpperCase()) {
		case "HASH":
			this.task_list = Task.hash;
			break;
		case "MAP":
			this.task_list = Task.map;
			break;
		default:
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");	
		}

		String[] args2 = new String[args.length-1];
		System.arraycopy(args, 1, args2, 0, args2.length);

		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-s", "--subject", true);
			myArgsEngine.add("-q", "--query", true);
			myArgsEngine.add("-l", "--library", true);
			myArgsEngine.add("-g", "--graph", true);
			myArgsEngine.add("-k", "--kmer-size", true);
			myArgsEngine.add("-H", "--hash-table", true);
			myArgsEngine.add("-x", "--max-mer-count", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-d", "--debug", false);
			myArgsEngine.add("-dd", "--debug-debug", false);
			myArgsEngine.add("-o", "--out", true);
			myArgsEngine.parse(args2);
		}

		if (myArgsEngine.getBoolean("-s")) {
			this.subject_file = myArgsEngine.getString("-s");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the contig file.");
		}

		if (myArgsEngine.getBoolean("-k")) {
			this.merK = Integer.parseInt(myArgsEngine.getString("-k"));
			if(this.merK>16) {
				merK = 16;
				myLogger.warn("Set mer size K=16.");
			}
		}

		if (myArgsEngine.getBoolean("-x")) {
			this.maxC = Integer.parseInt(myArgsEngine.getString("-x"));
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

		if (myArgsEngine.getBoolean("-d")) {
			this.debug = true;
		}

		if (myArgsEngine.getBoolean("-dd")) {
			this.debug  = true;
			this.ddebug = true;
		}

		switch(this.task_list) {
		case hash:

			break;
		case map:

			if (myArgsEngine.getBoolean("-q")) {
				this.query_file = myArgsEngine.getString("-q").split(",");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the contig file.");
			}

			if (myArgsEngine.getBoolean("-l")) {
				switch(myArgsEngine.getString("-l")) {
				case "pe":
					this.library = Library.pe;
					break;
				case "454":
					this.library = Library.r454;
					break;
				case "long":
					this.library = Library.long3;
					break;
				default:
					throw new IllegalArgumentException("Please specify the read library.");
				}
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the read library.");
			}

			if (myArgsEngine.getBoolean("-g")) {
				this.graph_file = myArgsEngine.getString("-g");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the map file.");
			}

			if (myArgsEngine.getBoolean("-H")) {
				this.hash_file = myArgsEngine.getString("-H");
			}
			break;
		default:
			throw new RuntimeException("!!!");	
		}
	}

	private final static int[] char_table = new int[256];
	static {
		char_table['A'] = 0;
		char_table['C'] = 1;
		char_table['G'] = 2;
		char_table['T'] = 3;
		char_table['a'] = 0;
		char_table['c'] = 1;
		char_table['g'] = 2;
		char_table['t'] = 3;
	}

	private final static char[] int_table = new char[4];
	static {
		int_table[0] = 'A';
		int_table[1] = 'C';
		int_table[2] = 'G';
		int_table[3] = 'T';
	}

	private Map<String, Sequence> sub_seqs;
	private GFA gfa;
	private final BidiMap<String, Integer> seq_index = new DualHashBidiMap<String, Integer>();
	private final Map<String, String> symm_seqsn = new HashMap<String, String>();
	// kmer hash table
	// key   :  mer
	// value :  positions of the mer.
	//          32bits sequence index + 32bits sequence position
	private final Map<Integer, Set<Long>> kmer_ht = new HashMap<Integer, Set<Long>>();
	private final static Object lock = new Object();
	private static long cons_progress = 0L, cons_size = 0L;

	@Override
	public void run() {
		// TODO Auto-generated method stub

		sub_seqs = Sequence.parseFastaFileWithRevCmpAsMap(subject_file);
		int index = 0;
		for(String seq : sub_seqs.keySet()) seq_index.put(seq, ++index);
		for(String seq : sub_seqs.keySet()) {
			if(!seq.endsWith("'")) {
				symm_seqsn.put(seq, seq+"'");
				symm_seqsn.put(seq+"'", seq);
			}
		}
		
		switch(this.task_list) {

		case zzz:
			myLogger.info("Task list is empty!!!");
			break;
		case hash:
			this.hash(true);
			break;
		case map:
			switch(this.library) {
			case pe:
				this.map_pe();
				break;
			case r454:
				this.map_r454();
				break;
			case long3:
				this.map_long3();
				break;
			default:
				throw new RuntimeException("!!!");
			}
			break;
		default:
			throw new RuntimeException("!!!");
		}
		return;
	}

	final static SamReaderFactory factory =
			SamReaderFactory.makeDefault()
			.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
					SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
			.validationStringency(ValidationStringency.SILENT);

	private void map_r454() {
		// TODO Auto-generated method stub
		gfa = new GFA(subject_file, graph_file);
		String source, target, seqid;
		Sequence source_seq, target_seq, dualSeq;
		Set<String> dualSeqid = new HashSet<String>();
		try {
			BufferedWriter bw_out = Utils.getBufferedWriter(out_prefix+".fa");
			for(OverlapEdge edge : gfa.edgeSet()) {
				source = gfa.getEdgeSource(edge);
				target = gfa.getEdgeTarget(edge);
				seqid = source+"_"+target;
				if(dualSeqid.contains(seqid)) continue;
				source_seq = sub_seqs.get(source);
				target_seq = sub_seqs.get(target);
				dualSeq = new Sequence(seqid, source_seq.seq_str()+target_seq.seq_str().substring((int)edge.olapR()));
				dualSeqid.add(seqid);
				dualSeqid.add(gfa.getRev(target)+"_"+gfa.getRev(source));
				bw_out.write(dualSeq.formatOutput());
			}
			bw_out.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
	}

	private final static Map<Long, Integer> linkCount = new HashMap<Long, Integer>();
	private final static Map<Long, Double> linkRadius = new HashMap<Long, Double>();
	private final static int m_ins = 5000; // maximum insert size for pe read library
	private final static int m_lnk = 3;    // minimum #link to confirm a link
	private final static int m_qual = 20;  // minimum alignment quality
	private static long exceed_ins = 0, links = 0, contained_single = 0, contained_multi = 0;
			
	private void map_pe() {
		// TODO Auto-generated method stub
		gfa = new GFA(subject_file, graph_file);
		try {
			this.initial_thread_pool();
			for(final String qf : query_file) {
				final SamReader in1 = factory.open(new File(qf));
				final SAMRecordIterator iter1 = in1.iterator();
				SAMRecord[] sam_records = new SAMRecord[2];
				SAMRecord sam_record;
				while(iter1.hasNext()) {
					sam_record = iter1.next();
					if(sam_record.getNotPrimaryAlignmentFlag() || 
							sam_record.getSupplementaryAlignmentFlag())
						continue;
					if(sam_record.getFirstOfPairFlag()) 
						sam_records[0] = sam_record;
					else sam_records[1] = sam_record;
					
					if(sam_records[0]!=null&&sam_records[1]!=null) {
						
						executor.submit(new Runnable() {
							private SAMRecord[] sam_records;
							
							@Override
							public void run() {
								// TODO Auto-generated method stub
								try {
									// so we have a confident read pair aligned to two contigs
									// check return
									if(sam_records[0].getReadUnmappedFlag() ||
											sam_records[1].getReadUnmappedFlag())
										return;
									if(!sam_records[0].getReadName().
											equals(sam_records[1].getReadName()))
										throw new RuntimeException("!!!");
									if(sam_records[0].getReferenceIndex().intValue()==
											sam_records[1].getReferenceIndex().intValue())
										return;
									if(sam_records[0].getMappingQuality()<m_qual || 
											sam_records[1].getMappingQuality()<m_qual)
										return;
									
									boolean rev0 = sam_records[0].getReadNegativeStrandFlag();
									boolean rev1 = sam_records[1].getReadNegativeStrandFlag();
									int reflen0  = sub_seqs.get(sam_records[0].getReferenceName()).seq_ln();
									int reflen1  = sub_seqs.get(sam_records[1].getReferenceName()).seq_ln();
									
									final String[] refstr = new String[2];
									double radius;
									int a, b;
									
									if(rev0&&rev1) {
										//       0                  1
										// ---------------     -------------
										//     <===                <===
										
										// reverse 0
										// ---------------     -------------
										//     ===>                <===
										a = reflen0-sam_records[0].getAlignmentStart()+1+sam_records[1].getAlignmentEnd();
										// reverse 1
										// ---------------     -------------
										//     <===                ===>
										b = sam_records[0].getAlignmentEnd()+reflen1-sam_records[1].getAlignmentStart()+1;
										
										if(ddebug)
											STD_OUT_BUFFER.write(">>>>"+Math.min(a, b)+"\n"+
													sam_records[0].getSAMString()+sam_records[1].getSAMString()+"<<<<\n");
										if(a>m_ins&&b>m_ins) {
											synchronized(lock) {++exceed_ins;}
											return;
										}
										if(a<b) {
											// reverse 0
											refstr[0] = sam_records[0].getReferenceName()+"'";
											refstr[1] = sam_records[1].getReferenceName();
											radius = a;
										} else {
											// reverse 1
											refstr[0] = sam_records[1].getReferenceName()+"'";
											refstr[1] = sam_records[0].getReferenceName();
											radius = b;
										}
									} else if(rev0&&!rev1) {
										//       0                  1
										// ---------------     -------------
										//     <===                ===>
										a = sam_records[0].getAlignmentEnd()+reflen1-sam_records[1].getAlignmentStart()+1;
										
										// reverse 0 & reverse1
										// ---------------     -------------
										//     ===>                <===
										b = reflen0-sam_records[0].getAlignmentStart()+1+sam_records[1].getAlignmentEnd();
										
										if(ddebug)
											STD_OUT_BUFFER.write(">>>>"+Math.min(a, b)+"\n"+
													sam_records[0].getSAMString()+sam_records[1].getSAMString()+"<<<<\n");
										if(a>m_ins&&b>m_ins) {
											synchronized(lock) {++exceed_ins;}
											return;
										}
										
										if(a<b) {
											refstr[0] = sam_records[1].getReferenceName();
											refstr[1] = sam_records[0].getReferenceName();
											radius = a;
										} else {
											refstr[0] = sam_records[0].getReferenceName()+"'";
											refstr[1] = sam_records[1].getReferenceName()+"'";
											radius = b;
										}
										
									} else if(!rev0&&rev1) {
										//       0                  1
										// ---------------     -------------
										//     ===>                <===
										a = reflen0-sam_records[0].getAlignmentStart()+1+sam_records[1].getAlignmentEnd();
										
										// reverse 0 & reverse1
										// ---------------     -------------
										//     <===                ===>
										b = sam_records[0].getAlignmentEnd()+reflen1-sam_records[1].getAlignmentStart()+1;
										
										if(ddebug)
											STD_OUT_BUFFER.write(">>>>"+Math.min(a, b)+"\n"+
													sam_records[0].getSAMString()+sam_records[1].getSAMString()+"<<<<\n");
										if(a>m_ins&&b>m_ins) {
											synchronized(lock) {++exceed_ins;}
											return;
										}
										
										if(a<b) {
											refstr[0] = sam_records[0].getReferenceName();
											refstr[1] = sam_records[1].getReferenceName();
											radius = a;
										} else {
											refstr[0] = sam_records[1].getReferenceName()+"'";
											refstr[1] = sam_records[0].getReferenceName()+"'";
											radius = b;
										}
										
									} else {
										//       0                  1
										// ---------------     -------------
										//     ===>                ===>
										
										// reverse 0
										// ---------------     -------------
										//     <===                ===>
										a = sam_records[0].getAlignmentEnd()+reflen1-sam_records[1].getAlignmentStart()+1;
										
										// reverse 1
										// ---------------     -------------
										//     ===>                <===
										b = reflen0-sam_records[0].getAlignmentStart()+1+sam_records[1].getAlignmentEnd();
										
										if(ddebug)
											STD_OUT_BUFFER.write(">>>>"+Math.min(a, b)+"\n"+
													sam_records[0].getSAMString()+sam_records[1].getSAMString()+"<<<<\n");
										if(a>m_ins&&b>m_ins) {
											synchronized(lock) {++exceed_ins;}
											return;
										}
										
										if(a<b) {
											refstr[0] = sam_records[1].getReferenceName();
											refstr[1] = sam_records[0].getReferenceName()+"'";
											radius = a;
										} else {
											refstr[0] = sam_records[0].getReferenceName();
											refstr[1] = sam_records[1].getReferenceName()+"'";
											radius = b;
										}
									}
									
									long refind;
									refind  = seq_index.get(refstr[0]);
									refind <<= 32;
									refind += seq_index.get(refstr[1]);
									synchronized(lock) {
										if(linkCount.containsKey(refind)) {
											linkCount.put(refind, linkCount.get(refind)+1);
											linkRadius.put(refind, linkRadius.get(refind)+radius);
										} else {
											linkCount.put(refind, 1);
											linkRadius.put(refind, radius);
										}
									}
									
									refind  = seq_index.get(symm_seqsn.get(refstr[1]));
									refind <<= 32;
									refind += seq_index.get(symm_seqsn.get(refstr[0]));
									synchronized(lock) {
										if(linkCount.containsKey(refind)) {
											linkCount.put(refind, linkCount.get(refind)+1);
											linkRadius.put(refind, linkRadius.get(refind)+radius);
										} else {
											linkCount.put(refind, 1);
											linkRadius.put(refind, radius);
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

							public Runnable init(SAMRecord[] sam_records) {
								this.sam_records = sam_records;
								return this;
							}
						}.init(sam_records));
						sam_records = new SAMRecord[2];
					}
				}
				
				iter1.close();
				in1.close();
			}
			this.waitFor();
			
			// now we get link counts
			String source, target;
			long refind;
			double radius;
			this.initial_thread_pool();
			for(Map.Entry<Long, Integer> entry : linkCount.entrySet()) {
				if(entry.getValue()<m_lnk) continue;
				refind = entry.getKey();
				target = seq_index.getKey((int)  refind     );
				source = seq_index.getKey((int) (refind>>32));
				radius = linkRadius.get(refind)/entry.getValue();
				
				++links;
				if(ddebug) STD_OUT_BUFFER.write("#link: "+source+"->"+target+"\n");
				if(gfa.containsEdge(source, target)) {
					++contained_single;
					continue;
				}
				executor.submit(new Runnable() {
					private String source;
					private String target;
					private double radius;
					
					@Override
					public void run() {
						// TODO Auto-generated method stub
						try {
							if(ddebug) {
								synchronized(lock) {
									STD_OUT_BUFFER.write("#search: "+source+"->"+target+","+(m_ins-radius)+"\n");
								}
							}
							TraversalTraceable traceable_target = search(source, target, m_ins-radius);
							if(traceable_target!=null) {
								// so we found a path from source to target
								synchronized(lock) {
									++contained_multi;
								}

							} else {
								// we didn't find a path, so we calculate the overlap between them

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

					public Runnable init(String source, String target, double radius) {
						// TODO Auto-generated method stub
						this.source = source;
						this.target = target;
						this.radius = radius;
						return this;
					}
					
				}.init(source, target, radius));
				
			}
			this.waitFor();
			STD_OUT_BUFFER.flush();
			
			myLogger.info("################################################################");
			myLogger.info("#insert size exceeds "+m_ins+": "+exceed_ins);
			myLogger.info("#good links: "+linkCount.size()/2);
			myLogger.info("#parsed/confident links: "+links/2);
			myLogger.info("#contained single-edges: "+contained_single/2+"/"+links/2);
			myLogger.info("#contained multi-edges : "+contained_multi/2+"/"+links/2);
			myLogger.info("################################################################");
			
		} catch(IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	protected TraversalTraceable search(final String source_str, final String target_str, final double radius) {
		// TODO Auto-generated method stub
		// search a path from source to target within a radius using BFS
		TraversalTraceable source, target, traceable_target;
		final Queue<TraversalTraceable> queue = new LinkedList<TraversalTraceable>();
		final Set<String> visited = new HashSet<String>();
		double distance;
		String reached;
		
		queue.offer(new TraversalTraceable(source_str));
		traceable_target = null;
		
		innerloop:
			while(!queue.isEmpty()) {
				source = queue.poll();
				visited.add(source.getId());
				
				for(final OverlapEdge outEdge : gfa.outgoingEdgesOf(source.getId())) {
					reached = gfa.getEdgeTarget(outEdge);
					
					if(reached.equals(target_str)) {
						// hit
						traceable_target = new TraversalTraceable(reached);
						traceable_target.setTraceBackward(source);
						break innerloop;
					}

					if(!visited.contains(reached) && // not visited yet
							(distance=source.getGap()+sub_seqs.get(reached).seq_ln()-
							outEdge.olapR())<=radius) { // radius limit 
						target = new TraversalTraceable(reached);
						target.setTraceBackward(source);
						target.setGap(distance);
					}
				}
			}
		return traceable_target;
	}

	private void hash() {
		// TODO Auto-generated method stub
		this.hash(false);
	}

	private void hash(boolean writeHashTable) {
		// TODO Auto-generated method stub
		// final long chunkSize = seqTotalSize/this.THREADS/10;
		final long chunkSize = 1000000;

		myLogger.info("++++JVM memory after loading data++++");
		myLogger.info("Total memory : "+totalMemory()+"Mb");
		myLogger.info("Free memory  : "+freeMemory() +"Mb");
		myLogger.info("Used memory  : "+usedMemory() +"Mb");

		// initialise the kmer hash table 
		myLogger.info("Construct initialise "+merK+"-mer hash table using "+this.THREADS+" threads.");
		long elapsed_start = System.nanoTime();


		this.initial_thread_pool();
		List<Sequence> sequences = new ArrayList<Sequence>();
		long seq_sz = 0;
		Iterator<Map.Entry<String, Sequence>> it = sub_seqs.entrySet().iterator();
		while(it.hasNext()) {
			Sequence seq = it.next().getValue();
			sequences.add(seq);
			seq_sz += seq.seq_ln();

			// we process chunkSize chunks for parallelism
			// no need to gain lock frequently compared to parallelise in sequence level
			// however will end up with extra work on copy hash table and extra memory consumption
			if(seq_sz<chunkSize&&it.hasNext()) continue; 

			executor.submit(new Runnable() {
				private List<Sequence> sequences;

				@Override
				public void run() {
					// TODO Auto-generated method stub

					try {
						final Map<Integer, Set<Long>> ht = new HashMap<Integer, Set<Long>>();
						for(Sequence sequence : sequences) {
							String seq_sn = sequence.seq_sn();
							long long_key = seq_index.get(seq_sn);
							long_key <<= 32;
							String seq_str = sequence.seq_str();
							int seq_ln = seq_str.length()-merK+1;
							String kmer;
							for(int i=0; i!=seq_ln; i++) {
								// process each mer
								kmer = seq_str.substring(i, i+merK);
								if(kmer.contains("N")||kmer.contains("n")) 
									continue;
								int kmer_hash = int_hash(kmer);

								if(!ht.containsKey(kmer_hash)) 
									ht.put(kmer_hash, new HashSet<Long>());
								ht.get(kmer_hash).add(long_key+i);
							}
						}

						synchronized(lock) {

							for(Map.Entry<Integer, Set<Long>> entry : ht.entrySet()) {
								if(!kmer_ht.containsKey(entry.getKey())) {
									kmer_ht.put(entry.getKey(), entry.getValue());
								} else {
									kmer_ht.get(entry.getKey()).addAll(entry.getValue());
									// not sure if this is necessary
									// likely will simplify the garbage collection? 
									// although this is not correct
									// ht.remove(entry.getKey());
								}
							}

							cons_progress += this.sequences.size();
							for(Sequence sequence : sequences) cons_size += sequence.seq_ln();
						}

						myLogger.info("#"+cons_progress+"/"+cons_size+"bp sequences processed.");
					} catch (Exception e) {
						// TODO Auto-generated catch block
						Thread t = Thread.currentThread();
						t.getUncaughtExceptionHandler().uncaughtException(t, e);
						e.printStackTrace();
						executor.shutdown();
						System.exit(1);
					}
				}

				public Runnable init(List<Sequence> sequences) {
					this.sequences = sequences;
					return this;
				}
			}.init(sequences));

			seq_sz = 0;
			sequences = new ArrayList<Sequence>();
		}
		this.waitFor();

		// filtering by max mer count
		for(Map.Entry<Integer, Set<Long>> entry : kmer_ht.entrySet()) 
			if(entry.getValue().size()>maxC) kmer_ht.remove(entry.getKey());

		long elapsed_end = System.nanoTime();
		myLogger.info(merK+"-mer hash table construction completed: "+kmer_ht.size()+" "+
				merK+"-mers in "+(elapsed_end-elapsed_start)/1e9+" secondes");
		myLogger.info("++++JVM memory after Kmer hash table construction++++");
		myLogger.info("Total memory : "+totalMemory()+"Mb");
		myLogger.info("Free memory  : "+freeMemory() +"Mb");
		myLogger.info("Used memory  : "+usedMemory() +"Mb");

		if(writeHashTable) this.writeHashTable();
	}

	private void writeHashTable() {
		// TODO Auto-generated method stub
		myLogger.info("Writing "+merK+"-mer hash table to File.");
		long elapsed_start = System.nanoTime();
		try {
			BufferedWriter bw_ht = Utils.getBufferedWriter(this.out_prefix+".h");
			for(Map.Entry<Integer, Set<Long>> entry : kmer_ht.entrySet()) {
				bw_ht.write(entry.getKey().toString());
				for(long val : entry.getValue()) bw_ht.write(" "+val);
				bw_ht.write("\n");
			}
			bw_ht.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		long elapsed_end = System.nanoTime();
		myLogger.info(merK+"-mer hash table writing completed: "+kmer_ht.size()+" "+
				merK+"-mers in "+(elapsed_end-elapsed_start)/1e9+" secondes");
		return;
	}

	private void readHashTable() {
		// TODO Auto-generated method stub
		myLogger.info("Loading "+merK+"-mer hash table from File.");
		long elapsed_start = System.nanoTime();
		long progress = 0;
		try {
			BufferedReader br_ht = Utils.getBufferedReader(this.hash_file);
			String line;
			String[] s;
			while( (line=br_ht.readLine())!=null ) {
				s = line.split(" ");
				if(s.length>maxC+1) continue;
				Set<Long> val = new HashSet<Long>();
				for(int i=1; i<s.length; i++)
					val.add(Long.parseLong(s[i]));
				kmer_ht.put(Integer.parseInt(s[0]), val);
				if(++progress%100000==0) myLogger.info(merK+"-mer processed: "+progress);
			}
			br_ht.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		long elapsed_end = System.nanoTime();
		myLogger.info(merK+"-mer hash table loading completed: "+kmer_ht.size()+" "+
				merK+"-mers in "+(elapsed_end-elapsed_start)/1e9+" secondes");
		myLogger.info("++++JVM memory with Kmer hash table loaded++++");
		myLogger.info("Total memory : "+totalMemory()+"Mb");
		myLogger.info("Free memory  : "+freeMemory() +"Mb");
		myLogger.info("Used memory  : "+usedMemory() +"Mb");
		return;
	}

	private void writeHashTableInParallel() {
		// TODO Auto-generated method stub
		myLogger.info("Writing "+merK+"-mer hash table to File.");
		long elapsed_start = System.nanoTime();
		BufferedWriter bw_ht = Utils.getBufferedWriter(this.out_prefix+".h");
		this.initial_thread_pool();
		for(Map.Entry<Integer, Set<Long>> entry : kmer_ht.entrySet()) {
			executor.submit(new Runnable() {

				@Override
				public void run() {
					// TODO Auto-generated method stub
					try {
						StringBuilder os = new StringBuilder();
						os.append(entry.getKey());
						for(long val : entry.getValue()) {
							os.append(" ");
							os.append(val);
						}
						os.append("\n");
						bw_ht.write(os.toString());
					} catch (Exception e) {
						// TODO Auto-generated catch block
						Thread t = Thread.currentThread();
						t.getUncaughtExceptionHandler().uncaughtException(t, e);
						e.printStackTrace();
						executor.shutdown();
						System.exit(1);
					}
				}
			});
		}
		this.waitFor();
		try {

			bw_ht.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		long elapsed_end = System.nanoTime();
		myLogger.info(merK+"-mer hash table writing completed: "+kmer_ht.size()+" "+
				merK+"-mers in "+(elapsed_end-elapsed_start)/1e9+" secondes");
		return;
	}

	private void readHashTableInParallel() {
		// TODO Auto-generated method stub
		myLogger.info("Loading "+merK+"-mer hash table from File.");
		long elapsed_start = System.nanoTime();
		try {
			BufferedReader br_ht = Utils.getBufferedReader(this.hash_file);
			String line; 
			this.initial_thread_pool();
			while( (line=br_ht.readLine())!=null ) {
				executor.submit(new Runnable() {
					private String line;

					@Override
					public void run() {
						// TODO Auto-generated method stub
						try {
							String[] s = line.split(" ");
							if(s.length>maxC+1) return;
							Set<Long> val = new HashSet<Long>();
							for(int i=1; i<s.length; i++)
								val.add(Long.parseLong(s[i]));
							synchronized(lock) {
								kmer_ht.put(Integer.parseInt(s[0]), val);
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

					public Runnable init(String line) {
						// TODO Auto-generated method stub
						this.line = line;
						return this;
					}
				}.init(line));
			}
			this.waitFor();
			br_ht.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		long elapsed_end = System.nanoTime();
		myLogger.info(merK+"-mer hash table loading completed: "+kmer_ht.size()+" "+
				merK+"-mers in "+(elapsed_end-elapsed_start)/1e9+" secondes");
		myLogger.info("++++JVM memory with Kmer hash table loaded++++");
		myLogger.info("Total memory : "+totalMemory()+"Mb");
		myLogger.info("Free memory  : "+freeMemory() +"Mb");
		myLogger.info("Used memory  : "+usedMemory() +"Mb");
		return;
	}

	private final static double d_max = 100; // max gap size
	private final static int k_min = 10; // min kmer count
	// assume the sequencing error is 0.15
	// for k=12
	// we will have (1-0.15)^12=0.142 mer hit per 12bp
	// which means we will see a mer hit per 12-1+1/0.142=18bp
	// we set it to 24bp for a relaxation
	private final static int k_dst = 24; // at least one kmer in k_dst bp on average
	private final static double olap_min = 0.99; // min overlap fraction for containment
	private final static double collinear_shift = 1.0;
	private final static double m_clip = 0.2; // max clip size (%) to treat an alignment end-to-end

	// ok, let call this a backup
	/***
	private void map() {
		// TODO Auto-generated method stub
		gfa = new GFA(subject_file, graph_file);
		if(this.hash_file==null) {
			this.hash();
		} else {
			this.readHashTable();
		}

		// now map each sequence in the query file to the graph
		try {
			BufferedReader br_qry = Utils.getBufferedReader(query_file);
			String line = br_qry.readLine();
			boolean isFASTQ = true;
			if(line.startsWith(">")) isFASTQ = false;
			String qry_sn, qry_str, kmer;
			int qry_ln, kmer_hash, sub_ind, sub_pos, qry_pos, sub_pos2, qry_pos2, score;

			KMP kmp_prev, kmp_curr;
			List<KMP> hits, sort_hits, segs;

			// a hashmap to hold the kmer hits of the query sequence to each subject sequence
			// a list to hold the kmer hits positions to the subject sequence 
			final Map<Integer, List<KMP>> kmer_hits = new HashMap<Integer, List<KMP>>();

			while( line!=null ) {
				qry_sn = line.split("\\s+")[0].substring(1);
				qry_str = br_qry.readLine();

				// we have query sequence now
				// we need to find shared kmers with the subject/reference sequences
				qry_ln = qry_str.length()-merK+1;
				kmer_hits.clear();

				for(int i=0; i!=qry_ln; i++) {
					// process each mer
					kmer = qry_str.substring(i, i+merK);
					if(kmer.contains("N")||kmer.contains("n")) 
						continue;
					kmer_hash = int_hash(kmer);

					if(kmer_ht.containsKey(kmer_hash)) {
						// we need to check the locations of the kmer on the subject sequences
						Set<Long> hts = kmer_ht.get(kmer_hash);
						for(long j : hts) {
							sub_pos = (int) j;
							sub_ind = (int) (j>>32);
							if(!kmer_hits.containsKey(sub_ind)) 
								kmer_hits.put(sub_ind, new ArrayList<KMP>());
							// TODO: this is not right if a kmer is mapped to two positions 
							//       on the reference sequence
							kmer_hits.get(sub_ind).add(new KMP(i, sub_pos));
						}
					}
				}

				sort_hits = new ArrayList<KMP>();
				for(final int i : kmer_hits.keySet()) 
					sort_hits.add(new KMP(i, kmer_hits.get(i).size()));
				Collections.sort(sort_hits, new Comparator<KMP>() {
					@Override
					public int compare(KMP o1, KMP o2) {
						// TODO Auto-generated method stub
						return o2.b-o1.b;
					}
				});	
				// now we get all kmer hits
				// for each subject sequence we now filter non-collinear hits
				// and calculate the alignment
				for(final KMP k : sort_hits) {
					final int i = k.a;
					hits = kmer_hits.get(i);
					Collections.sort(hits);
					// we find best alignment block on query sequence first
					//            #                                    ##
					// ------o----------o-o--o--o--o--o--o--o-oo-ooo---------o-o---o-
					// we break the alignment from the very sparse regions
					// in the example above, we break in # and ## positions
					// TODO: this is not exactly right as we also need to consider 
					//       the positions on the reference sequence

					// this is a placeholder for segments <position, #kmers>
					segs = new ArrayList<KMP>();
					kmp_prev = hits.get(0);
					int start = 0;
					for(int w=1; w<hits.size(); w++) {
						kmp_curr = hits.get(w);
						if(KMP.distance(kmp_prev, kmp_curr)>d_max) {
							// break here
							segs.add(new KMP(start, w-start));
							start = w;
						}
						kmp_prev = kmp_curr;
					}

					System.out.println();
				}

				if(isFASTQ) {
					// skip two lines if is FASTQ file
					br_qry.readLine();
					br_qry.readLine();
				}
				line = br_qry.readLine();
			}

			br_qry.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	 **/

	private void map_long3() {
		// TODO Auto-generated method stub
		gfa = new GFA(subject_file, graph_file);
		if(this.hash_file==null) {
			this.hash();
		} else {
			this.readHashTable();
		}

		// now map each sequence in the query file to the graph
		try {
			this.initial_thread_pool();

			for(final String qf : query_file) {

				BufferedReader br_qry = Utils.getBufferedReader(qf);
				String line = br_qry.readLine();
				boolean isFASTQ = true;
				if(line.startsWith(">")) isFASTQ = false;

				String qry_sn, qry_sq;
				while( line!=null ) {

					qry_sn = line.split("\\s+")[0].substring(1);
					qry_sq = br_qry.readLine();

					executor.submit(new Runnable() {
						private String qry_sn;
						private String qry_sq;

						@Override
						public void run() {
							// TODO Auto-generated method stub
							final StringBuilder std_out = new StringBuilder();

							synchronized(lock) {
								if(++cons_progress%1000==0)
									std_out.append("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
											+ "Reads processed: "+cons_progress+
											"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
							}

							try {
								std_out.append(">>>>>>>>>>"+qry_sn+"<<<<<<<<<<\n");

								String kmer;
								int qry_ln, sub_ln;
								int kmer_hash, sub_ind, sub_pos;
								KMP kmp_prev, kmp_curr;

								// we have query sequence now
								// we need to find shared kmers with the subject/reference sequences
								qry_ln = qry_sq.length();

								// a hashmap to hold the kmer hits of the query sequence to each subject sequence
								// a list to hold the kmer hits positions to the subject sequence 
								final Map<Integer, List<KMP>> kmer_hits = new HashMap<Integer, List<KMP>>();
								final int N = qry_ln-merK+1;
								for(int i=0; i<N; i++) {
									// process each mer
									kmer = qry_sq.substring(i, i+merK);
									if(kmer.contains("N")||kmer.contains("n")) 
										continue;
									kmer_hash = int_hash(kmer);

									if(kmer_ht.containsKey(kmer_hash)) {
										// we need to check the locations of the kmer on the subject sequences
										Set<Long> hts = kmer_ht.get(kmer_hash);
										for(long j : hts) {
											sub_pos = (int) j;
											sub_ind = (int) (j>>32);
											if(!kmer_hits.containsKey(sub_ind)) 
												kmer_hits.put(sub_ind, new ArrayList<KMP>());
											// TODO: this is not right if a kmer is mapped to two positions 
											//       on the reference sequence
											kmer_hits.get(sub_ind).add(new KMP(i, sub_pos));
										}
									}
								}

								final Set<String> sub_hits = new HashSet<String>();
								for(final int i : kmer_hits.keySet()) {
									if(kmer_hits.get(i).size()>=3)
										// we need at least three mers to confirm a hit
										sub_hits.add(seq_index.getKey(i));
								}

								final List<KMP> sort_hits = new ArrayList<KMP>();
								for(final int i : kmer_hits.keySet()) 
									sort_hits.add(new KMP(i, kmer_hits.get(i).size()));
								Collections.sort(sort_hits, new Comparator<KMP>() {
									@Override
									public int compare(KMP o1, KMP o2) {
										// TODO Auto-generated method stub
										return o2.b-o1.b;
									}
								});	

								int qstart, qend, sstart, send, merCount, a, b;
								double clip;
								// now we get all kmer hits
								// for each subject sequence we now filter non-collinear hits
								// and calculate the alignment
								final List<TraceableAlignmentSegment> alignments = new ArrayList<TraceableAlignmentSegment>();
								final Set<Integer> a_rm = new HashSet<Integer>();
								final Set<Integer> b_rm = new HashSet<Integer>();
								final Set<Integer> a_bucket = new HashSet<Integer>();
								final Set<Integer> b_bucket = new HashSet<Integer>();
								TraceableAlignmentSegment alignment;

								for(final KMP k : sort_hits) {
									final int i = k.a;

									final List<KMP> kmps = kmer_hits.get(i);
									if(kmps.size()<k_min) continue;

									// we find best alignment block on query sequence first
									//    ________________________________
									//   |            1/ o        /     / |
									//   |         2 o/          /     /  |
									//   |          3/o         /     /   |
									//   |        4 /o         /     /    |
									//   |         /       5 x/     /     |
									//   |     6 o/          /     /      |
									//   |       /          /  7 x/       |
									//   |     8/o         /     /        |
									//   |_____/__________/_____/_________|
									//    a   b   c   d   e   f   g   h    
									// (1 2 3 4 6 8)     (5)   (7)
									// in the example above, the intercepts on the x-axis gathers around b
									// (1 2 3 4 6 8) are selected kmers
									// TODO: this is not exactly right as we also need to consider 
									//       the positions on the reference sequence
									// TODO: this is not always correct
									//       especially in the repetitive region
									// test case: 27f73126-832e-48a9-b7b1-794c670ec630_Basecall_1D_template
									//  
									//    ________________________________
									//   |     1/o     /          /     / |
									//   |     2/o    /          /     /  |
									//   |     3/o   /          /     /   |
									//   |     4/o  /          /     /    |
									//   |     5/o /       9 x/     /     |
									//   |     6/o/          /     /      |
									//   |    7 o/          /  a x/       |
									//   |     8/o         /     /        |
									//   |_____/__________/_____/_________|
									//    a   b   c   d   e   f   g   h
									//  (1 2 ... 8)      (5)   (6)
									//    ________________________________
									//   |             /          /     / |
									//   |            /          /     /  |
									//   |           /          /     /   |
									//   |          /          /     /    |
									//   |         /       9 o/     /     |
									//   |      1 /2 3 4 5 6 /  a o/      |
									//   |      o/ o o o o o/o o  /       |
									//   |      /          / 7 8 /        |
									//   |_____/__________/_____/_________|
									//    a   b   c   d   e   f   g   h
									//             (1 2 ... a) 
									// in the example above, put (1 2 ... a) together is very dangerous 
									// this is possible if the query sequence is repetitive
									//
									// this is to fix the repetitive problem

									a_rm.clear();
									b_rm.clear();
									a_bucket.clear();
									b_bucket.clear();

									for(int w=0; w<kmps.size(); w++) {
										kmp_prev = kmps.get(w);
										a = kmp_prev.a;
										if(a_bucket.contains(a)) {
											a_rm.add(a);
										} else {
											a_bucket.add(a);
										}
										b = kmp_prev.b;
										if(b_bucket.contains(b)) {
											b_rm.add(b);
										} else {
											b_bucket.add(b);
										}
									}

									final List<KMP> hits = new ArrayList<KMP>(); 
									for(KMP kmp : kmps) {
										if(!a_rm.contains(kmp.a)&&
												!b_rm.contains(kmp.b)) 
											hits.add(kmp);
									}

									if(hits.isEmpty()) continue;

									// this is a placeholder for intercepts <position, #kmers>
									final List<KMP> intercept = new ArrayList<KMP>();
									for(int w=0; w<hits.size(); w++) {
										KMP p = hits.get(w);
										intercept.add(new KMP(p.a-p.b, w));
									}
									Collections.sort(intercept);

									if(ddebug) {
										std_out.append("+");
										std_out.append(sub_seqs.get(seq_index.getKey(i)).seq_sn());
										std_out.append("\n");
										for(KMP kmp : hits) {
											std_out.append(kmp.a+" "+kmp.b);
											std_out.append("\n");
										}
									}

									// this is a placeholder for segments <position, #kmers>
									final List<KMP> segs = new ArrayList<KMP>();
									kmp_prev = intercept.get(0);
									int start = 0, nk;
									int is = intercept.size();
									for(int w=1; w<is; w++) {
										kmp_curr = intercept.get(w);
										if(kmp_curr.a-kmp_prev.a>d_max) {
											// break here
											if( (nk=w-start)>=k_min ) segs.add(new KMP(start, nk));
											start = w;
										}
										kmp_prev = kmp_curr;
									}
									if( (nk=is-start)>=k_min ) segs.add(new KMP(start, nk));
									if(segs.isEmpty()) continue;

									Collections.sort(segs, new Comparator<KMP>() {

										@Override
										public int compare(KMP kmp, KMP kmp2) {
											// TODO Auto-generated method stub
											return kmp2.b-kmp.b;
										}

									});

									List<TraceableAlignmentSegment> seg_list = new ArrayList<TraceableAlignmentSegment>();
									for(KMP seg : segs) {
										qstart = Integer.MAX_VALUE;
										sstart = Integer.MAX_VALUE;
										qend   = Integer.MIN_VALUE;
										send   = Integer.MIN_VALUE;
										for(int w=seg.a; w<seg.a+seg.b; w++) {
											kmp_prev = hits.get(intercept.get(w).b);
											if(kmp_prev.a<qstart) {
												qstart = kmp_prev.a;
												sstart = kmp_prev.b;
											}
											if(kmp_prev.a>qend) {
												qend   = kmp_prev.a;
												send   = kmp_prev.b;
											}
										}
										qend  += merK;
										send  += merK;

										if(qend<qstart) {
											throw new RuntimeException("!!!");
										}
										// filter by segment size
										if(qend-qstart+1<d_max) continue;
										// filter by kmer density
										// ok will do this after processing collinearity
										//if(qend-qstart+1>k_dst*seg.b) continue;
										seg_list.add(new TraceableAlignmentSegment(qry_sn, sub_seqs.get(seq_index.getKey(i)).seq_sn(), qstart, qend, sstart, send, seg.b));
									}

									if(seg_list.isEmpty()) continue;

									// now need to merge collinear segment
									// we always keep the longest hit
									// and check the collinearity with the remains
									TraceableAlignmentSegment primary_seg, secondary_seg;
									final List<Integer> sels = new ArrayList<Integer>();
									sels.add(0);
									int max_shift;
									for(int w=1; w<seg_list.size();w++) {
										primary_seg = seg_list.get(w);
										for(int z : sels) {
											secondary_seg = seg_list.get(z); 
											max_shift = Math.min(primary_seg.qlength(), secondary_seg.qlength());
											if( Math.abs(primary_seg.qintercept()-secondary_seg.qintercept())<=max_shift||
													TraceableAlignmentSegment.pdistance(primary_seg, secondary_seg)<=max_shift ) {
												sels.add(w);
												break;
											}
										}
									}

									qstart = Integer.MAX_VALUE;
									sstart = Integer.MAX_VALUE;
									qend   = Integer.MIN_VALUE;
									send   = Integer.MIN_VALUE;
									merCount = 0;
									for(int z : sels) {
										primary_seg = seg_list.get(z); 
										if(primary_seg.qstart()<qstart) qstart = primary_seg.qstart();
										if(primary_seg.qend()  >qend)   qend = primary_seg.qend();
										if(primary_seg.sstart()<sstart) sstart = primary_seg.sstart();
										if(primary_seg.send()  >send)   send = primary_seg.send();
										merCount += primary_seg.getMerCount();
									}

									if(qend-qstart+1>k_dst*merCount) continue;
									alignment = new TraceableAlignmentSegment(qry_sn, sub_seqs.get(seq_index.getKey(i)).seq_sn(), qstart, qend, sstart, send, merCount);

									// check if this is end-to-end alignment
									sub_ln = sub_seqs.get(seq_index.getKey(i)).seq_ln();
									clip = Math.max(100, (qend-qstart+1)*m_clip);
									alignment.setEndToEnd(qstart<=clip, qend+clip>=qry_ln, sstart<=clip, send+clip>=sub_ln);
									alignment.setClip(qstart-1, qry_ln-qend, sstart-1, sub_ln-send);
									alignment.calcScore();
									if(alignment.getEndToEnd()) alignments.add(alignment);
								}

								Collections.sort(alignments, new TraceableAlignmentSegment.QLengthComparator());

								// RangeSet<Integer> qry_cov = TreeRangeSet.create();
								// qry_cov.add( Range.closed(1, qry_ln).canonical(DiscreteDomain.integers()) );
								// for(TraceableAlignmentSegment as : alignments) 
								// 	qry_cov.remove( Range.closed(as.qstart(), as.qend()).canonical(DiscreteDomain.integers()) );
								//	int as_ln = 0;
								//	for(Range<Integer> r : qry_cov.asRanges()) as_ln += r.upperEndpoint()-r.lowerEndpoint();
								//	myLogger.info(qry_sn+": alignment fraction "+(1-(double)as_ln/qry_ln));


								if(alignments.size()<=1) {
									std_out.append(qry_sn+": alignment fraction "+ 
											(alignments.isEmpty()?0:((double)alignments.get(0).qlength()/qry_ln))+"\n");
									std_out.append("<<<<<<<<<<"+qry_sn+">>>>>>>>>>\n");
									STD_OUT_BUFFER.write(std_out.toString());
									return;
								}
								Collections.sort(alignments, new TraceableAlignmentSegment.QueryCoordinationComparator());

								if(ddebug) for(TraceableAlignmentSegment as : alignments) {
									std_out.append(as.toString()+"\n");
								}

								// TODO: need to check the collinearity, i.e., cannot be too distant alignments
								//       also if the clip is too large, then we don't want to merge them neither
								// test case: e023806d-6dbd-42d2-bac0-da6734d90e52_Basecall_1D_template
								//            f7e121bb-12dd-4144-8d3d-d8426c686d0b_Basecall_1D_template
								//            f57698c4-3e82-4d9f-b7fb-5f99c872bf9c_Basecall_1D_template
								//            b956bf82-e9b9-47d1-9071-af5ea1143ce2_Basecall_1D_template
								int asz = alignments.size();
								TraceableAlignmentSegment source_as, target_as, tmp_as;
								String source_id, target_id;
								int source_qstart, source_qend, target_qstart, target_qend;
								double source_ln, target_ln;
								final Map<String, List<TraceableAlignmentSegment>> merged_seq = new HashMap<String, List<TraceableAlignmentSegment>>(); 
								for(int w=0; w<asz; w++) {
									source_as = alignments.get(w);
									// source need to be no clip to the end
									if(source_as==null||!source_as.getToEnd()) continue;
									source_id = source_as.sseqid();
									source_qend = source_as.qend();
									source_ln = source_as.qend()-source_as.qstart()+1;
									sub_ln = sub_seqs.get(source_id).seq_ln();

									for(int z=w+1; z<asz; z++) {
										target_as = alignments.get(z);
										// target need to be not clip to the start
										if(target_as==null||!target_as.getToStart()) continue;
										target_id = target_as.sseqid();
										target_ln = target_as.qend()-target_as.qstart()+1;
										// if gap size is too big
										clip = Math.max(100, (source_ln+target_ln)*m_clip);
										if(target_as.qstart()-source_qend>clip) continue;

										if(gfa.containsEdge(source_id, target_id)) {
											if(merged_seq.containsKey(source_id)) {
												List<TraceableAlignmentSegment> new_seq = merged_seq.get(source_id);
												new_seq.add(target_as);
												merged_seq.remove(source_id);
												merged_seq.put(target_id, new_seq);
											} else {
												List<TraceableAlignmentSegment> new_seq = new ArrayList<TraceableAlignmentSegment>();
												new_seq.add(source_as);
												new_seq.add(target_as);
												merged_seq.put(target_id, new_seq);
											}
											alignment = new TraceableAlignmentSegment(qry_sn, target_id, 
													source_as.qstart(), 
													Math.max(source_as.qend(),   target_as.qend()),
													-1, 
													-1);
											alignment.setEndToEnd(source_as.getQueryStartClip()<=clip, target_as.getQueryEndClip()<=clip,
													source_as.getSubjectStartClip()<=clip, target_as.getSubjectEndClip()<=clip);
											alignment.setClip(source_as.getQueryStartClip(), target_as.getQueryEndClip(), 
													source_as.getSubjectStartClip(), target_as.getSubjectEndClip());
											alignment.calcScore();
											alignments.set(w,  alignment);
											alignments.set(z, null);
											--w;
											break;
										}
									}
								}

								if(ddebug) {
									std_out.append("--------------------------------\n");
									for(TraceableAlignmentSegment as : alignments) if(as!=null) {
										std_out.append(as.toString()+"\n");
									}
								}

								// we remove containment
								int source_len, target_len;
								double source_olap, target_olap;
								double olap;
								boolean sEndToEnd, tEndToEnd;
								outerloop:
									for(int w=0; w<asz; w++) {
										source_as = alignments.get(w);
										if(source_as==null) continue;
										source_qstart = source_as.qstart();
										source_qend   = source_as.qend();
										source_len    = source_qend-source_qstart+1;
										sEndToEnd = source_as.getEndToEnd();
										for(int z=w+1; z<asz; z++) {
											target_as = alignments.get(z);
											if(target_as==null) continue;
											target_qstart = target_as.qstart();
											target_qend   = target_as.qend();
											target_len    = target_qend-target_qstart+1;
											tEndToEnd = target_as.getEndToEnd();
											if(!sEndToEnd&&!tEndToEnd) continue;
											// we calculate overlap size
											if(target_qend<=source_qend) {
												if(sEndToEnd) alignments.set(z, null);
												continue;
											}
											if((olap=source_qend-target_qstart)<=0)
												continue outerloop;
											source_olap = olap/source_len;
											target_olap = olap/target_len;
											if(source_olap>=olap_min&&
													target_olap>=olap_min) {
												if(sEndToEnd&&tEndToEnd) {
													if(source_olap<target_olap) {
														alignments.set(z, null);
														continue;
													} else {
														alignments.set(w, null);
														continue outerloop;
													}
												} 
												if(sEndToEnd) {
													alignments.set(z, null);
													continue;
												} 
												if(tEndToEnd) {
													alignments.set(w, null);
													continue outerloop;
												}
											}
											if(source_olap>=olap_min&&tEndToEnd) {
												alignments.set(w, null);
												continue outerloop;
											}
											if(target_olap>=olap_min&&sEndToEnd) {
												alignments.set(z, null);
												continue;
											}
										}
									}
								if(ddebug) {
									std_out.append("--------------------------------\n");
									for(TraceableAlignmentSegment as : alignments) if(as!=null) {
										std_out.append(as.toString()+"\n");
									}
								}

								// we find the best path for coverage
								// get a copy of the alignment to remove nulls
								final List<TraceableAlignmentSegment> selected = new ArrayList<TraceableAlignmentSegment>();
								for(TraceableAlignmentSegment as : alignments) if(as!=null) selected.add(as);
								if(selected.size()<=1) {
									if(debug) {
										std_out.append("--------------------------------\n");
										if(!selected.isEmpty()) {
											TraceableAlignmentSegment as = selected.get(0);
											if(merged_seq.containsKey(as.sseqid())) {
												List<TraceableAlignmentSegment> tas = merged_seq.get(as.sseqid());
												for(TraceableAlignmentSegment ta : tas) std_out.append(ta.toString()+"\n");
											} else {
												std_out.append(as.toString()+"\n");
											}	
										}
									}
									std_out.append(qry_sn+": alignment fraction "+ 
											(selected.isEmpty()?0:((double)selected.get(0).qlength()/qry_ln))+"\n");
									std_out.append("<<<<<<<<<<"+qry_sn+">>>>>>>>>>\n");
									STD_OUT_BUFFER.write(std_out.toString());
									return;
								}

								double objective = Double.NEGATIVE_INFINITY, dobj, opt_obj = 0d, tmp; // record the current best path
								int tmp_qstart, tmp_qend;
								TraceableAlignmentSegment traceback = null;
								asz = selected.size();

								for(int w=0; w<asz; w++) {
									source_as = selected.get(w);
									if(source_as.getTraceBackward()!=null) 
										continue;
									objective = source_as.getScore();
									if(objective<0) continue;
									source_id = source_as.sseqid();
									source_qstart = source_as.qstart();
									source_qend   = source_as.qend();

									// so we start from w
									int z = w+1;
									outerloop:
										while(z<asz) {
											// first find next 
											target_as = selected.get(z);

											if(!target_as.getToSubjectStart()||
													!target_as.getToEnd()) {
												++z;
												continue;
											}

											target_qstart = target_as.qstart();
											target_qend   = target_as.qend();

											// calculate addition to the objective
											dobj = target_as.getScore()+
													(target_qstart>source_qend?ScoreMatrix.gap_open:0)-
													(target_qstart>source_qend?0:(source_qend-target_qstart)*ScoreMatrix.match_score);

											// #a <------------------->
											// #b      <----------------->
											// #c            <--------------->
											// we prefer #c over #b
											for(int v=z+1; v<asz; v++) {
												tmp_as = selected.get(v);
												if(!tmp_as.getToSubjectStart()||
														!tmp_as.getToEnd()) 
													continue;

												tmp_qstart = tmp_as.qstart();
												tmp_qend = tmp_as.qend();
												if(tmp_qstart>source_qend && 
														tmp_qstart>target_qstart)
													break;
												tmp = tmp_as.getScore()+
														(tmp_qstart>source_qend?ScoreMatrix.gap_open:0)-
														(tmp_qstart>source_qend?0:(source_qend-tmp_qstart)*ScoreMatrix.match_score);
												if(tmp>dobj) {
													// update selected
													dobj = tmp;
													target_as = tmp_as;
													target_qstart = tmp_qstart;
													target_qend   = tmp_qend;
													z = v;
												}
											}

											// OK now we found it
											if(dobj>0 && dobj+objective>target_as.getObjective()) {
												// if the objective is better, we choose this path
												objective += dobj;
												target_as.setTraceBackward(source_as);
												target_as.setObjective(objective);
												source_as.setTraceForward(target_as);
												// finally we update source
												source_as = target_as;
												source_id = source_as.sseqid();
												source_qstart = target_qstart;
												source_qend   = target_qend;

												if(source_as.getTraceForward()!=null) {
													// so we trace forward from here to avoid redoing this
													while((target_as=source_as.getTraceForward())!=null) {
														objective = target_as.getObjective()+dobj;
														target_as.setObjective(objective);
														source_as = target_as;
													}
													break outerloop;
												}
											}
											// don't forget this
											++z;
										}

									// now we get a path
									if(objective>opt_obj) {
										opt_obj = objective;
										traceback = source_as;	
									}
								}

								// now we get the final path
								// which can be traced back from 'traceback'
								if(traceback==null) {
									std_out.append(qry_sn+": alignment fraction 0\n");
									std_out.append("<<<<<<<<<<"+qry_sn+">>>>>>>>>>\n");
									STD_OUT_BUFFER.write(std_out.toString());
									return;
								}

								final List<TraceableAlignmentSegment> graph_path = new ArrayList<TraceableAlignmentSegment>();
								while(traceback!=null) {
									graph_path.add(traceback);
									traceback = traceback.getTraceBackward();
								}
								Collections.reverse(graph_path);

								source_as = graph_path.get(0);
								source_qstart = source_as.qstart();
								source_qend   = source_as.qend();
								double dcov = source_qend-source_qstart+1;
								for(int z=1; z<graph_path.size(); z++) {
									target_as = graph_path.get(z);
									target_qstart = target_as.qstart();
									target_qend   = target_as.qend();
									dcov += target_qend-target_qstart+1-
											(target_qstart>source_qend?0:(source_qend-target_qstart));
									source_as = target_as;
									source_qstart = target_qstart;
									source_qend   = target_qend;
								}
								std_out.append(qry_sn+": alignment fraction "+ dcov/qry_ln+"\n");

								if(debug) {
									std_out.append("--------------------------------\n");
									for(TraceableAlignmentSegment as : graph_path) { 
										if(merged_seq.containsKey(as.sseqid())) {
											List<TraceableAlignmentSegment> tas = merged_seq.get(as.sseqid());
											for(TraceableAlignmentSegment ta : tas) std_out.append(ta.toString()+"\n");
										} else {
											std_out.append(as.toString()+"\n");
										}	
									}
								}

								// now contigging and updating the assembly graph


								std_out.append("<<<<<<<<<<"+qry_sn+">>>>>>>>>>\n");
								STD_OUT_BUFFER.write(std_out.toString());
							} catch (Exception e) {
								// TODO Auto-generated catch block
								Thread t = Thread.currentThread();
								t.getUncaughtExceptionHandler().uncaughtException(t, e);
								e.printStackTrace();
								executor.shutdown();
								System.exit(1);
							}
						}

						public Runnable init(String qry_sn, String qry_sq) {
							this.qry_sn = qry_sn;
							this.qry_sq = qry_sq;
							return this;
						}
					}.init(qry_sn, qry_sq));

					if(isFASTQ) {
						// skip two lines if is FASTQ file
						br_qry.readLine();
						br_qry.readLine();
					}
					line = br_qry.readLine();
				}

				br_qry.close();
			}
			this.waitFor();
			STD_OUT_BUFFER.flush();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private int int_hash(String kmer) {
		// TODO Auto-generated method stub
		int hash = 0;
		for(int i=0; i<merK; i++) {
			hash <<= 2;
			hash += char_table[kmer.charAt(i)];
		}
		return hash;
	}

	private String str_hash(int hash) {
		// TODO Auto-generated method stub
		char[] c = new char[merK];
		int hash_copy = hash;
		for(int i=merK-1; i>=0; i--) {
			c[i] = int_table[hash_copy&3];
			hash_copy >>= 2;
		}
		return String.valueOf(c);
	}
}

final class KMP implements Comparable<KMP> {
	// this is only a data pair holder which could be used in many cases
	final int a;
	final int b;

	public KMP(final int a, final int b) {
		this.a = a;
		this.b = b;
	}

	@Override
	public int compareTo(KMP kmp) {
		// TODO Auto-generated method stub
		return this.a-kmp.a;
	}

	public double distance(KMP kmp) {
		return distance(this, kmp);
	}

	public static double distance(KMP kmp, KMP kmp2) {
		// TODO Auto-generated method stub
		double da = kmp.a-kmp2.a, db = kmp.b-kmp2.b;
		return Math.sqrt(da*da+db*db);
	}
}

final class ScoreMatrix {
	final static int match_score = 2;
	final static int mismatch_penalty = -4;
	final static int gap_open = -20;
	final static int gap_extension = -1;
}

final class TraversalTraceable {
	private final String id;
	private TraversalTraceable previous = null;
	private double gap = 0;

	public TraversalTraceable(final String id) {
		this.id = id;
	}

	public String getId() {
		return this.id;
	}

	public void setTraceBackward(TraversalTraceable previous) {
		this.previous = previous;
	}

	public TraversalTraceable getTraceBackward() {
		return this.previous;
	}

	public void setGap(double gap) {
		this.gap = gap;
	}

	public double getGap() {
		return this.gap;
	}

	@Override
	public int hashCode(){
		return this.id.hashCode();
	}

	@Override
	public boolean equals(Object obj){
		if (obj instanceof TraversalTraceable) {
			return this.id.equals(((TraversalTraceable)obj).id);
		}
		return false;
	}
}

final class TraceableAlignmentSegment extends AlignmentSegment {

	private TraceableAlignmentSegment next = null;
	private TraceableAlignmentSegment previous = null;
	private double objective = 0d;
	private int mer_count = 0;
	private boolean to_query_start = false;
	private boolean to_query_end   = false;
	private boolean to_subject_start = false;
	private boolean to_subject_end   = false;
	private boolean to_start = false;
	private boolean to_end   = false;
	private boolean end_to_end = false; 
	private int sstart_clip = 0;
	private int send_clip = 0;
	private int qstart_clip = 0;
	private int qend_clip = 0;
	private int score = 0;

	public TraceableAlignmentSegment(String qseqid, String sseqid, 
			int qstart, int qend, int sstart, int send) {
		// TODO Auto-generated constructor stub
		super(qseqid, sseqid, qstart, qend, sstart, send, true);
	}

	public void setClip(int qstart_clip, int qend_clip, int sstart_clip, int send_clip) {
		// TODO Auto-generated method stub
		this.qstart_clip = qstart_clip;
		this.qend_clip = qend_clip;
		this.sstart_clip = sstart_clip;
		this.send_clip = send_clip;
	}

	public void setQueryStartClip(int qstart_clip) {
		// TODO Auto-generated method stub
		this.qstart_clip = qstart_clip;
	}

	public void setQueryEndClip(int qend_clip) {
		// TODO Auto-generated method stub
		this.qend_clip = qend_clip;
	}

	public int getQueryStartClip() {
		// TODO Auto-generated method stub
		return this.qstart_clip;
	}

	public int getQueryEndClip() {
		// TODO Auto-generated method stub
		return this.qend_clip;
	}

	public void setSubjectStartClip(int sstart_clip) {
		// TODO Auto-generated method stub
		this.sstart_clip = sstart_clip;
	}

	public void setSubjectEndClip(int send_clip) {
		// TODO Auto-generated method stub
		this.send_clip = send_clip;
	}

	public int getSubjectStartClip() {
		// TODO Auto-generated method stub
		return this.sstart_clip;
	}

	public int getSubjectEndClip() {
		// TODO Auto-generated method stub
		return this.send_clip;
	}

	public int getClip() {
		// TODO Auto-generated method stub
		return this.sstart_clip+this.send_clip;
	}

	public void setEndToEnd(boolean to_query_start, boolean to_query_end,
			boolean to_subject_start, boolean to_subject_end) {
		// TODO Auto-generated method stub
		this.to_query_start = to_query_start;
		this.to_query_end = to_query_end;
		this.to_subject_start = to_subject_start;
		this.to_subject_end = to_subject_end;
		this.to_start = to_query_start||to_subject_start;
		this.to_end   = to_query_end||to_subject_end;
		this.end_to_end = to_start&&to_end;
	}

	public void setToQueryStart(boolean to_query_start) {
		// TODO Auto-generated method stub
		this.to_query_start = to_query_start;
	}

	public void setToQueryEnd(boolean to_query_end) {
		// TODO Auto-generated method stub
		this.to_query_end = to_query_end;
	}

	public void setToSubjectStart(boolean to_subject_start) {
		// TODO Auto-generated method stub
		this.to_subject_start = to_subject_start;
	}

	public void setToSubjectEnd(boolean to_subject_end) {
		// TODO Auto-generated method stub
		this.to_subject_end = to_subject_end;
	}

	public boolean getToQueryStart() {
		// TODO Auto-generated method stub
		return this.to_query_start;
	}

	public boolean getToQueryEnd() {
		// TODO Auto-generated method stub
		return this.to_query_end;
	}

	public boolean getToSubjectStart() {
		// TODO Auto-generated method stub
		return this.to_subject_start;
	}

	public boolean getToSubjectEnd() {
		// TODO Auto-generated method stub
		return this.to_subject_end;
	}

	public void setEndToEnd(boolean end_to_end) {
		// TODO Auto-generated method stub
		this.end_to_end = end_to_end;
	}

	public void setToStart(boolean to_start) {
		// TODO Auto-generated method stub
		this.to_start = to_start;
	}

	public void setToEnd(boolean to_end) {
		// TODO Auto-generated method stub
		this.to_end = to_end;
	}

	public boolean getToStart() {
		// TODO Auto-generated method stub
		return this.to_start;
	}


	public boolean getToEnd() {
		// TODO Auto-generated method stub
		return this.to_end;
	}

	public boolean getEndToEnd() {
		// TODO Auto-generated method stub
		return this.end_to_end;
	}

	public TraceableAlignmentSegment(String qseqid, String sseqid, 
			int qstart, int qend, int sstart, int send, int mer_count) {
		// TODO Auto-generated constructor stub
		super(qseqid, sseqid, qstart, qend, sstart, send, true);
		this.mer_count = mer_count;
	}

	public void setTraceForward(TraceableAlignmentSegment next) {
		this.next = next;
	}

	public void setTraceBackward(TraceableAlignmentSegment previous) {
		this.previous = previous;
	}

	public void calcScore() {
		this.score = (this.qend-this.qstart+1)*ScoreMatrix.match_score+
				Math.min(qstart_clip, sstart_clip)*ScoreMatrix.gap_extension+
				Math.min(qend_clip, send_clip)*ScoreMatrix.gap_extension;
	}

	public void setScore(int score) {
		this.score = score;
	}

	public int getScore() {
		return this.score;
	}

	public void setObjective(double objective) {
		this.objective = objective;
	}

	public void addObjective(double dobj) {
		this.objective += dobj;
	}

	public void setMerCount(int mer_count) {
		this.mer_count = mer_count;
	}

	public TraceableAlignmentSegment getTraceForward() {
		return this.next;
	}

	public TraceableAlignmentSegment getTraceBackward() {
		return this.previous;
	}

	public double getObjective() {
		return this.objective;
	}

	public int getMerCount() {
		return this.mer_count;
	}

	public static TraceableAlignmentSegment collinear(final TraceableAlignmentSegment record1, 
			final TraceableAlignmentSegment record2, final double max_shift) {
		// TODO Auto-generated method stub

		// if(TraceableAlignmentSegment.sdistance(record1, record2)>max_shift ||
		// 		TraceableAlignmentSegment.qdistance(record1, record2)>max_shift ||
		//		TraceableAlignmentSegment.pdistance(record1, record2)>max_shift) {
		//	return null;
		// }

		// how about we only check the pdistance?
		if(TraceableAlignmentSegment.pdistance(record1, record2)>max_shift)
			return null;

		// merge collinear alignment segments
		int qstart = Math.min(record1.qstart(), record2.qstart());
		int qend = Math.max(record1.qend(), record2.qend());
		int sstart = Math.min(record1.sstart(), record2.sstart());
		int send = Math.max(record1.send(), record2.send());

		return new TraceableAlignmentSegment(record1.qseqid, record1.sseqid,qstart,qend,sstart,send);
	}
}









