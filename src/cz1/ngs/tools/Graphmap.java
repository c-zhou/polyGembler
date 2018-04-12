package cz1.ngs.tools;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;

import cz1.ngs.model.GFA;
import cz1.ngs.model.Sequence;
import cz1.util.ArgsEngine;
import cz1.util.Executor;

public class Graphmap extends Executor {

	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -s/--subject            The FASTA file contain subject/reference sequences. \n"
						+ " -q/--query              The FASTA file contain query sequences to map. \n"
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
						+ " -k/--kmer-size          K-mer size (no greater than 16, default 12)."
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
	private static long cons_progress = 0;
	
	@Override
	public void run() {
		// TODO Auto-generated method stub

		gfa = new GFA(subject_file, graph_file);
		sub_seqs = gfa.getSequenceMap();
		int index = 0;
		for(String seq : sub_seqs.keySet()) seq_index.put(seq, ++index);
		
		// initialise the kmer hash table 
		myLogger.info("Construct initialise "+merK+"-mer hash table using "+this.THREADS+" threads.");
		long elapsed_start = System.nanoTime();
		
		this.initial_thread_pool();
		for(Map.Entry<String, Sequence> entry : sub_seqs.entrySet()) {
			executor.submit(new Runnable() {
				private Sequence sequence;
				
				@Override
				public void run() {
					// TODO Auto-generated method stub
					String seq_sn = sequence.seq_sn();
					long long_key = seq_index.get(seq_sn);
					long_key <<= 32;
					String seq_str = sequence.seq_str();
					int seq_ln = seq_str.length()-merK;
					String kmer;
					final Map<Integer, Set<Long>> ht = new HashMap<Integer, Set<Long>>();
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
					synchronized(lock) {
						for(Map.Entry<Integer, Set<Long>> entry : ht.entrySet()) {
							if(!kmer_ht.containsKey(entry.getKey())) {
								kmer_ht.put(entry.getKey(), entry.getValue());
							} else {
								kmer_ht.get(entry.getKey()).addAll(entry.getValue());
							}
						}
						if(++cons_progress%10000==0) 
							myLogger.info(cons_progress+" sequences processed.");
					}
				}
				
				public Runnable init(Sequence sequence) {
					this.sequence = sequence;
					return this;
				}
			}.init(entry.getValue()));
		}
		this.waitFor();
		long elapsed_end = System.nanoTime();
		myLogger.info(merK+"-mer hash table construction completed: "+kmer_ht.size()+" "+
				merK+"-mers in "+(elapsed_end-elapsed_start)/1e9+" secondes");
		
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












