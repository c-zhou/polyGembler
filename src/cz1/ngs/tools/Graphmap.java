package cz1.ngs.tools;

import java.io.BufferedReader;
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

	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
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
						+ " -t/--threads            Threads to use (default 1). \n"
						+ "                         The maximum number of threads to use will be the number of BAM files.\n"
						+ " -o/--out                Prefix of the output files.\n"
						+ "\n");
	}
	
	private String subject_file;
	private String graph_file;
	private String query_file;
	private String out_prefix;
	private int merK = 12;

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}

		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-s", "--subject", true);
			myArgsEngine.add("-q", "--query", true);
			myArgsEngine.add("-g", "--graph", true);
			myArgsEngine.add("-k", "--kmer-size", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-o", "--out", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-s")) {
			this.subject_file = myArgsEngine.getString("-s");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the contig file.");
		}
		
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
	final BidiMap<String, Integer> seq_index = new DualHashBidiMap<String, Integer>();
	// kmer hash table
	// key   :  mer
	// value :  positions of the mer.
	//          32bits sequence index + 32bits sequence position
	final Map<Integer, Set<Long>> kmer_ht = new HashMap<Integer, Set<Long>>();
	private final static Object lock = new Object();
	private static long cons_progress = 0L, cons_size = 0L;
	
	@Override
	public void run() {
		// TODO Auto-generated method stub

		gfa = new GFA(subject_file, graph_file);
		sub_seqs = gfa.getSequenceMap();
		int index = 0;
		long seqTotalSize = 0L;
		for(String seq : sub_seqs.keySet()) {
			seq_index.put(seq, ++index);
			seqTotalSize += sub_seqs.get(seq).seq_ln();
		}
		final long chunkSize = seqTotalSize/this.THREADS/10;
		
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
									ht.remove(entry.getKey());
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
	
		// now map each sequence in the query file to the graph
		try {
			BufferedReader br_qry = Utils.getBufferedReader(query_file);
			String line = br_qry.readLine();
			boolean isFASTQ = true;
			if(line.startsWith(">")) isFASTQ = false;
			String qry_sn, qry_str, kmer;
			int qry_ln, kmer_hash, sub_ind, sub_pos;
			Set<Long> hits;
			
			// a hashmap to hold the kmer hits of the query sequence to each subject sequence
			// a red-black tree to hold the kmer hits positions to the subject sequence 
			final Map<Integer, TreeMap<Integer, Integer>> kmer_hits = new HashMap<Integer, TreeMap<Integer, Integer>>();
			// a red-black tree to compare the kmer hits to subject sequences
			// we calculate a score for each subject sequence
			final TreeMap<Integer, Double> hits_scorer = new TreeMap<Integer, Double>();
			
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
						hits = kmer_ht.get(kmer_hash);
						for(long j : hits) {
							sub_pos = (int) j;
							sub_ind = (int) (j>>32);
							if(!kmer_hits.containsKey(sub_ind)) 
								kmer_hits.put(sub_ind, new TreeMap<Integer, Integer>());
							kmer_hits.get(sub_ind).put(i, sub_pos);
						}
					}
				}
				
				// now we get all kmer hits
				// we calculate a score for each hit
				hits_scorer.clear();
				
				
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












