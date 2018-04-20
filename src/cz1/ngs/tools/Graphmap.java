package cz1.ngs.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;

import com.google.common.collect.DiscreteDomain;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

import cz1.ngs.model.AlignmentSegment;
import cz1.ngs.model.GFA;
import cz1.ngs.model.SAMSegment;
import cz1.ngs.model.Sequence;
import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

public class Graphmap extends Executor {

	protected final static BufferedWriter STD_OUT_BUFFER = new BufferedWriter(new OutputStreamWriter(System.out),65536);
	private static enum Task {hash, map, zzz}
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
							+ " -t/--threads            Threads to use (default 1). \n"
							+ "                         The maximum number of threads to use will be the number of BAM files.\n"
							+ " -o/--out                Prefix of the output files.\n"
							+ "\n");
			break;
		case map:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -s/--subject            The FASTA file contain subject/reference sequences. \n"
							+ " -q/--query              The FASTA/FASTQ file contain query sequences to map. \n"
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
	private String query_file;
	private String hash_file = null;
	private String out_prefix;
	private int merK = 12;
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
			myArgsEngine.add("-g", "--graph", true);
			myArgsEngine.add("-k", "--kmer-size", true);
			myArgsEngine.add("-H", "--hash-table", true);
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
				this.query_file = myArgsEngine.getString("-q");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the contig file.");
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
	// kmer hash table
	// key   :  mer
	// value :  positions of the mer.
	//          32bits sequence index + 32bits sequence position
	private final Map<Integer, Set<Long>> kmer_ht = new HashMap<Integer, Set<Long>>();
	private final static Object lock = new Object();
	private final static int match_score = 2;
	private final static int mismatch_penalty = -4;
	private final static int[] gap_open = new int[]{-4, -24};
	private final static int[] gap_extension = new int[]{-2, -1};
	private static long cons_progress = 0L, cons_size = 0L;

	@Override
	public void run() {
		// TODO Auto-generated method stub

		sub_seqs = Sequence.parseFastaFileWithRevCmpAsMap(subject_file);
		int index = 0;
		for(String seq : sub_seqs.keySet()) seq_index.put(seq, ++index);

		switch(this.task_list) {

		case zzz:
			myLogger.info("Task list is empty!!!");
			break;
		case hash:
			this.hash(true);
			break;
		case map:
			this.map();
			break;
		default:
			throw new RuntimeException("!!!");
		}
		return;
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
		try {
			BufferedReader br_ht = Utils.getBufferedReader(this.hash_file);
			String line;
			String[] s;
			while( (line=br_ht.readLine())!=null ) {
				s = line.split(" ");
				Set<Long> val = new HashSet<Long>();
				for(int i=1; i<s.length; i++)
					val.add(Long.parseLong(s[i]));
				kmer_ht.put(Integer.parseInt(s[0]), val);
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

			this.initial_thread_pool();

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
							double qry_ln;
							int kmer_hash, sub_ind, sub_pos;
							KMP kmp_prev, kmp_curr;
							List<KMP> hits, sort_hits, intercept, segs;

							// we have query sequence now
							// we need to find shared kmers with the subject/reference sequences
							qry_ln = qry_sq.length()-merK+1;

							// a hashmap to hold the kmer hits of the query sequence to each subject sequence
							// a list to hold the kmer hits positions to the subject sequence 
							final Map<Integer, List<KMP>> kmer_hits = new HashMap<Integer, List<KMP>>();

							for(int i=0; i!=qry_ln; i++) {
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

							int qstart, qend, sstart, send, merCount;
							// now we get all kmer hits
							// for each subject sequence we now filter non-collinear hits
							// and calculate the alignment
							final List<TraceableAlignmentSegment> alignments = new ArrayList<TraceableAlignmentSegment>();
							for(final KMP k : sort_hits) {
								final int i = k.a;

								hits = kmer_hits.get(i);
								if(hits.size()<k_min) continue;

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

								// this is a placeholder for intercepts <position, #kmers>
								intercept = new ArrayList<KMP>();
								for(int w=0; w<hits.size(); w++) {
									KMP p = hits.get(w);
									intercept.add(new KMP(p.a-p.b, w));
								}
								Collections.sort(intercept);

								if(ddebug) {
									for(KMP kmp : hits) {
										std_out.append(kmp.a+" "+kmp.b);
										std_out.append("\n");
									}
								}
								
								// this is a placeholder for segments <position, #kmers>
								segs = new ArrayList<KMP>();
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
								alignments.add(new TraceableAlignmentSegment(qry_sn, sub_seqs.get(seq_index.getKey(i)).seq_sn(), qstart, qend, sstart, send, merCount));
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
										(alignments.isEmpty()?0:(alignments.get(0).qlength()/qry_ln))+"\n");
								std_out.append("<<<<<<<<<<"+qry_sn+">>>>>>>>>>\n");
								STD_OUT_BUFFER.write(std_out.toString());
								return;
							}
							Collections.sort(alignments, new TraceableAlignmentSegment.QueryCoordinationComparator());

							if(ddebug) for(TraceableAlignmentSegment as : alignments) {
								std_out.append(as.toString()+"\n");
							}

							int asz = alignments.size();
							TraceableAlignmentSegment source_as, target_as, tmp_as;
							String source_id, target_id;
							final Map<String, List<TraceableAlignmentSegment>> merged_seq = new HashMap<String, List<TraceableAlignmentSegment>>(); 
							for(int w=0; w<asz; w++) {
								source_as = alignments.get(w);
								if(source_as==null) continue;
								source_id = source_as.sseqid();
								for(int z=w+1; z<asz; z++) {
									target_as = alignments.get(z);
									if(target_as==null) continue;
									target_id = target_as.sseqid();
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
										alignments.set(w, new TraceableAlignmentSegment(qry_sn, target_id, 
												source_as.qstart(), 
												Math.max(source_as.qend(),   target_as.qend()),
												-1, -1) );
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
							int source_qstart, source_qend, target_qstart, target_qend, source_len, target_len;
							double source_olap, target_olap;
							double olap;
							outerloop:
								for(int w=0; w<asz; w++) {
									source_as = alignments.get(w);
									if(source_as==null) continue;
									source_qstart = source_as.qstart();
									source_qend   = source_as.qend();
									source_len    = source_qend-source_qstart+1;
									for(int z=w+1; z<asz; z++) {
										target_as = alignments.get(z);
										if(target_as==null) continue;
										target_qstart = target_as.qstart();
										target_qend   = target_as.qend();
										target_len    = target_qend-target_qstart+1;
										// we calculate overlap size
										if(target_qend<=source_qend) {
											alignments.set(z, null);
											continue;
										}
										if((olap=source_qend-target_qstart)<=0)
											continue outerloop;
										source_olap = olap/source_len;
										target_olap = olap/target_len;
										if(source_olap>=olap_min&&
												target_olap>=olap_min) {
											if(source_olap<target_olap) {
												alignments.set(z, null);
												continue;
											} else {
												alignments.set(w, null);
												continue outerloop;
											}
										}

										if(source_olap>=olap_min) {
											alignments.set(w, null);
											continue outerloop;
										}
										if(target_olap>=olap_min) {
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
								std_out.append(qry_sn+": alignment fraction "+ 
										(selected.isEmpty()?0:(selected.get(0).qlength()/qry_ln))+"\n");
								std_out.append("<<<<<<<<<<"+qry_sn+">>>>>>>>>>\n");
								STD_OUT_BUFFER.write(std_out.toString());
								return;
							}

							double dcov = 0d, objective = 0d; // record the current best path
							TraceableAlignmentSegment traceback = null;
							asz = selected.size();

							for(int w=0; w<asz; w++) {
								source_as = selected.get(w);
								if(source_as.getTraceBack()!=null || 
										source_as.qstart()>(qry_ln-dcov)) 
									continue;
								source_id = source_as.sseqid();
								source_qstart = source_as.qstart();
								source_qend   = source_as.qend();
								objective = source_qend-source_qstart+1;
								source_as.setObjective(objective); 

								// so we start from w
								int z = w+1;
								while(z<asz) {
									// first find next 
									// #a <------------------->
									// #b      <----------------->
									// #c            <--------------->
									// we prefer #c over #b
									target_as = selected.get(z);
									target_qstart = target_as.qstart();
									target_qend    = target_as.qend();
									for(int v=z+1; v<asz; v++) {
										tmp_as = selected.get(v);
										if(tmp_as.qstart()>source_qend)
											break;
										if(tmp_as.qend()>=target_qend) {
											target_as = tmp_as;
											target_qstart = target_as.qstart();
											target_qend   = target_as.qend();
											z = v;
										}
									}
									// OK now we found it
									// calculate addition do the objective
									objective += target_qend-target_qstart+1-
											(target_qstart>source_qend?0:(source_qend-target_qstart));

									if(objective>target_as.getObjective()) {
										// if the objective is better, we choose this path
										target_as.setTraceBack(source_as);

										// finally we update source
										source_as = target_as;
										source_qstart = target_qstart;
										source_qend   = target_qend;

										// don't forget this
										++z;
									}
								}

								// now we get a path
								if(objective>dcov) {
									dcov = objective;
									traceback = source_as;	
								}
							}

							// now we get the final path
							// which can be traced back from 'traceback'
							std_out.append(qry_sn+": alignment fraction "+ objective/qry_ln+"\n");
							final List<TraceableAlignmentSegment> graph_path = new ArrayList<TraceableAlignmentSegment>();
							while(traceback!=null) {
								graph_path.add(traceback);
								traceback = traceback.getTraceBack();

							}
							Collections.reverse(graph_path);
							if(debug) {
								std_out.append("--------------------------------\n");
								for(TraceableAlignmentSegment as : graph_path) { 
									if(merged_seq.containsKey(as.sseqid())) {
										List<TraceableAlignmentSegment> ma = merged_seq.get(as.sseqid());
										for(TraceableAlignmentSegment a : ma) std_out.append(a.toString()+"\n");
									} else {
										std_out.append(as.toString()+"\n");
									}	
								}
							}

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
			this.waitFor();
			STD_OUT_BUFFER.flush();
			
			br_qry.close();
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


final class TraceableAlignmentSegment extends AlignmentSegment {

	private TraceableAlignmentSegment trace = null;
	private double objective = 0d;
	private int mer_count = 0;

	public TraceableAlignmentSegment(String qseqid, String sseqid, 
			int qstart, int qend, int sstart, int send) {
		// TODO Auto-generated constructor stub
		super(qseqid, sseqid, qstart, qend, sstart, send, true);
	}

	public TraceableAlignmentSegment(String qseqid, String sseqid, 
			int qstart, int qend, int sstart, int send, int mer_count) {
		// TODO Auto-generated constructor stub
		super(qseqid, sseqid, qstart, qend, sstart, send, true);
		this.mer_count = mer_count;
	}

	public void setTraceBack(TraceableAlignmentSegment trace) {
		this.trace = trace;
	}

	public void setObjective(double objective) {
		this.objective = objective;
	}

	public void setMerCount(int mer_count) {
		this.mer_count = mer_count;
	}

	public TraceableAlignmentSegment getTraceBack() {
		return this.trace;
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









