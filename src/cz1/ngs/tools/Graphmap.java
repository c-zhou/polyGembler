package cz1.ngs.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;

import cz1.ngs.model.GFA;
import cz1.ngs.model.Sequence;
import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

public class Graphmap extends Executor {

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
		switch(this.task_list) {

		case zzz:
			myLogger.info("Task list is empty!!!");
			break;
		case hash:
			this.hash(true);
			break;
		case map:
			this.graphmap();
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
		sub_seqs = Sequence.parseFastaFileWithRevCmpAsMap(subject_file);
		int index = 0;
		long seqTotalSize = 0L;
		for(String seq : sub_seqs.keySet()) {
			seq_index.put(seq, ++index);
			seqTotalSize += sub_seqs.get(seq).seq_ln();
		}
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

	private void graphmap() {
		// TODO Auto-generated method stub
		gfa = new GFA(subject_file, graph_file);
		sub_seqs = gfa.getSequenceMap();
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

			// a hashmap to hold the kmer hits of the query sequence to each subject sequence
			// a red-black tree to hold the kmer hits positions to the subject sequence 
			final Map<Integer, TreeMap<Integer, Integer>> kmer_hits = new HashMap<Integer, TreeMap<Integer, Integer>>();
			// a red-black tree to compare the kmer hits to subject sequences
			// we calculate a score for each subject sequence
			final TreeMap<Integer, Integer> hits_scorer = new TreeMap<Integer, Integer>();

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
						Set<Long> hits = kmer_ht.get(kmer_hash);
						for(long j : hits) {
							sub_pos = (int) j;
							sub_ind = (int) (j>>32);
							if(!kmer_hits.containsKey(sub_ind)) 
								kmer_hits.put(sub_ind, new TreeMap<Integer, Integer>());
							// TODO: this is not right if a kmer is mapped to two positions 
							//       on the reference sequence
							kmer_hits.get(sub_ind).put(i, sub_pos);
						}
					}
				}

				// now we get all kmer hits
				// we calculate a score for each hit
				hits_scorer.clear();
				for(int i : kmer_hits.keySet()) {
					TreeMap<Integer, Integer> hits = kmer_hits.get(i);
					Map.Entry<Integer, Integer> entry = hits.firstEntry();
					qry_pos = entry.getKey();
					sub_pos = entry.getValue();
					score = merK*match_score;
					while( (entry=hits.ceilingEntry(qry_pos))!=null ) {
						qry_pos2 = entry.getKey();
						sub_pos2 = entry.getValue();

					}
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












