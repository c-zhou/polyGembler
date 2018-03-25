package cz1.ngs.model;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.RejectedExecutionHandler;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import org.apache.log4j.Logger;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.jgrapht.GraphPath;

import cz1.util.Constants;
import cz1.util.Utils;
import edu.umd.marbl.mhap.impl.MinHashSearch;
import edu.umd.marbl.mhap.impl.SequenceId;
import edu.umd.marbl.mhap.impl.SequenceSketchStreamer;

public class GFA {

	protected final static Logger myLogger = Logger.getLogger(GFA.class);

	private enum Graph_type {de_bruijn, olc, overlap};
	private Graph_type type = null;
	private final DirectedWeightedOverlapPseudograph<String> gfa = 
			new DirectedWeightedOverlapPseudograph<String>(OverlapEdge.class);
	private final Map<String, Sequence> seg = new HashMap<String, Sequence>();
	private final Map<String, String>   rev = new HashMap<String, String>();
	private BFSShortestPath<String, OverlapEdge> bfs;
	
	public GFA(String seq_file, String graph_file) {
		this.readSequence(seq_file);
		if(graph_file==null) 
			this.createOverlapGraph(seq_file);
		else this.readGraph(graph_file);
	}

	private static final double DEFAULT_OVERLAP_ACCEPT_SCORE = 0.78;
	private static final double DEFAULT_REPEAT_WEIGHT= 0.9;
	private static final int DEFAULT_KMER_SIZE = 16;
	private static final double DEFAULT_MAX_SHIFT_PERCENT = 0.2;
	private static final int DEFAULT_MIN_STORE_LENGTH = 0;
	private static final int DEFAULT_MIN_OVL_LENGTH = 30;
	private static final int DEFAULT_NUM_MIN_MATCHES = 30;
	private static final int DEFAULT_NUM_THREADS = Constants.omp_threads;
	private static final int DEFAULT_NUM_WORDS = 512;
	private static final int DEFAULT_ORDERED_KMER_SIZE = 12;
	private static final int DEFAULT_ORDERED_SKETCH_SIZE = 1536;
	
	private static final int default_flank = 30;
	private static final Object lock = new Object();
	
	private static final String tmp_dir = System.getProperty("java.io.tmpdir");
	private final boolean mhap_on_the_fly = false;
	
	private void createOverlapGraph(String seq_file) {
		// TODO Auto-generated method stub
		if(!mhap_on_the_fly)
			throw new RuntimeException("Run MHAP on the fly has yet to be implemented. "
					+ "Please run MHAP manually and provide the file!!!");
		
		// detect overlaps with MinHash (MHAP)
		SequenceId.STORE_FULL_ID = true;

		// read and index the kmers
		int seqNumberProcessed = 0;	
		//create search object
		SequenceSketchStreamer seqStreamer;
		try {
			seqStreamer = new SequenceSketchStreamer(
					seq_file, 
					DEFAULT_MIN_OVL_LENGTH, 
					DEFAULT_KMER_SIZE, 
					DEFAULT_NUM_WORDS,
					DEFAULT_ORDERED_KMER_SIZE, 
					DEFAULT_ORDERED_SKETCH_SIZE, 
					null, 
					true, 
					DEFAULT_REPEAT_WEIGHT, 
					seqNumberProcessed);
			/***
			 * write overlaps to file
			 **/ 
			MinHashSearch hashSearch = new MinHashSearch(seqStreamer, 
					DEFAULT_NUM_WORDS, 
					DEFAULT_NUM_MIN_MATCHES, 
					DEFAULT_NUM_THREADS, 
					false, // do not store results, write to file
					DEFAULT_MIN_STORE_LENGTH, 
					DEFAULT_MAX_SHIFT_PERCENT, 
					DEFAULT_OVERLAP_ACCEPT_SCORE, 
					true);

			/**
			System.setOut(new PrintStream(new BufferedOutputStream(new FileOutputStream("C:/users/chenxi.zhou/desktop/file.txt"))));
			hashSearch.findMatches();
			System.out.flush();
			System.setOut(System.out);
			**/
			
			// TODO: make this work
			
			// we redirect System.out to a input stream and read it
			/***
			MinHashSearch hashSearch = new MinHashSearch(seqStreamer, 
					DEFAULT_NUM_WORDS, 
					DEFAULT_NUM_MIN_MATCHES, 
					DEFAULT_NUM_THREADS, 
					false, // do not store results
					DEFAULT_MIN_STORE_LENGTH, 
					DEFAULT_MAX_SHIFT_PERCENT, 
					DEFAULT_OVERLAP_ACCEPT_SCORE, 
					true);

			PipedInputStream pipe_in = new PipedInputStream(1024);
			PipedOutputStream pipe_out = new PipedOutputStream(pipe_in);
			PrintStream ps = new PrintStream(pipe_out);
			PrintStream sys_os = System.out;
			System.setOut(ps);

			Thread findMatches = new Thread(new Runnable() {
				@Override
				public void run() {
					// TODO Auto-generated method stub
					hashSearch.findMatches();
					System.out.flush();
				}
			});

			Thread streamReader = new Thread(new Runnable() {
				@Override
				public void run() {
					// TODO Auto-generated method stub
					try {
						BufferedReader br = new BufferedReader(new InputStreamReader(pipe_in));
						String olap;
						while( (olap=br.readLine())!=null ) {

						}
						br.close();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			});
			findMatches.start();
			streamReader.start();
			try {
				findMatches.join();
				streamReader.join();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			System.setOut(sys_os);
			**/
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	private int olap_loading_progress = 0;
	private long elapsed_time = 0;
	
	private void readOverlapGraph(String olap_file) {
		// TODO Auto-generated method stub
		myLogger.info("Loading assembly graph from file.");
		// 1. read file to initialise assembly graph
		// remember to make it symmetric 
		myLogger.info("    intialising assembly graph.");
		
		BufferedReader br_olap = Utils.getBufferedReader(olap_file);
		String olap;
		
		this.initial_thread_pool();

		try {
			
			while( (olap = br_olap.readLine())!=null ) {
				
				this.executor.submit(new Runnable() {
					private String olap;
					
					@Override
					public void run() {
						// TODO Auto-generated method stub
						
						
						synchronized(lock) {
							olap_loading_progress++;
							if(olap_loading_progress%1000==0) 
								myLogger.info("      olap_loading_progress: "+olap_loading_progress+"\n"
										     +"           loading_elapsed : "+elapsed_time);
						}
						
						String[] olap_str;
						String fromId, toId, fromSeq, toSeq, fromSeq_flank, toSeq_flank;
						boolean fromFw, toFw, olapFw;
						double score, rawScore;
						int a1, a2, b1, b2, fromLen, toLen, a1_flank, a2_flank, b1_flank, b2_flank, aLen;
						SequencePair<DNASequence, NucleotideCompound> seqPair;
						OverlapEdge edge; 
						
						olap_str =   olap.split("\\s+");
						fromId   =   olap_str[0];
						toId     =   olap_str[1];
						score    =   Double.parseDouble(olap_str[2]);
						fromFw   =   olap_str[4].equals("0");
						a1       =   Integer.parseInt(olap_str[5]);
						a2       =   Integer.parseInt(olap_str[6]);
						fromLen  =   Integer.parseInt(olap_str[7]);
						toFw     =   olap_str[8].equals("0");
						b1       =   Integer.parseInt(olap_str[9]);
						b2       =   Integer.parseInt(olap_str[10]);
						toLen    =   Integer.parseInt(olap_str[11]);

						fromSeq   = seg.get(fromId).seq_str();
						a1_flank  = Math.max(0, a1-default_flank);
						a2_flank  = Math.min(fromLen, a2+default_flank);
						toSeq     = seg.get(toId).seq_str();
						b1_flank  = Math.max(0, b1-default_flank);
						b2_flank  = Math.min(toLen, b2+default_flank);
						fromSeq_flank = fromFw ? fromSeq.substring(a1_flank, a2_flank) : Sequence.revCompSeq(fromSeq.substring(a1_flank, a2_flank));
						toSeq_flank   = toFw   ? toSeq.substring  (b1_flank, b2_flank) : Sequence.revCompSeq(toSeq.substring  (b1_flank, b2_flank));
						
						long ss = System.nanoTime();
						seqPair = this.pairMatcher(fromSeq_flank, toSeq_flank);
						synchronized(lock) {
							elapsed_time += System.nanoTime()-ss;
						}
						
						// alignment coordinates on the flanked sequences
						aLen = seqPair.getLength();
						a1 = seqPair.getIndexInQueryAt(1);
						a2 = seqPair.getIndexInQueryAt(aLen);
						b1 = seqPair.getIndexInTargetAt(1);
						b2 = seqPair.getIndexInTargetAt(aLen);

						// convert them back to the original sequence coordinates
						if(fromFw) {
							a1 += a1_flank;
							a2 += a1_flank;
						} else {
							a1 += fromLen-a2_flank;
							a2 += fromLen-a2_flank;
						}
						if(toFw) {
							b1 += b1_flank;
							b2 += b1_flank;
						} else {
							b1 += toLen-b2_flank;
							b2 += toLen-b2_flank;
						}

						// use overlap length as raw score
						rawScore = (double) aLen;

						// decide the overlap direction
						// true : a->b
						// false: b->a
						olapFw = a1+toLen-b2>b1+fromLen-a2;

						// add vertices
						fromId = fromFw ? fromId : rev.get(fromId);
						toId   = toFw   ? toId   : rev.get(toId  );
						
						synchronized(lock) {
							if(!gfa.containsVertex(fromId))  gfa.addVertex(fromId);
							if(!gfa.containsVertex(toId)  )  gfa.addVertex(toId  );

							if(olapFw) {
								// a -> b
								edge = gfa.addEdge(fromId, toId);
								gfa.setEdgeWeight(edge, 1.0);
								gfa.setEdgeOverlapInfo(edge, new OverlapResult(fromId, toId, score, rawScore, a1, a2, fromLen, b1, b2, toLen));
							} else {
								// b -> a
								edge = gfa.addEdge(toId, fromId);
								gfa.setEdgeWeight(edge, 1.0);
								gfa.setEdgeOverlapInfo(edge, new OverlapResult(toId, fromId, score, rawScore, b1, b2, toLen, a1, a2, fromLen));
							}

							// symmetrizing 
							fromId = rev.get(fromId);
							toId   = rev.get(toId  );

							if(!gfa.containsVertex(fromId))  gfa.addVertex(fromId);
							if(!gfa.containsVertex(toId)  )  gfa.addVertex(toId  );
							if(olapFw) {
								// a -> b
								edge = gfa.addEdge(toId, fromId);
								gfa.setEdgeWeight(edge, 1.0);
								gfa.setEdgeOverlapInfo(edge, new OverlapResult(toId, fromId, score, rawScore, toLen-b2+1, toLen-b1+1, toLen, fromLen-a2+1, fromLen-a2+1, fromLen));
							} else {
								// b -> a
								edge = gfa.addEdge(fromId, toId);
								gfa.setEdgeWeight(edge, 1.0);
								gfa.setEdgeOverlapInfo(edge, new OverlapResult(fromId, toId, score, rawScore, fromLen-a2+1, fromLen-a1+1, fromLen, toLen-b2+1, toLen-b1+1, toLen));
							}
						}
					}
					
					public Runnable init(String olap) {
						// TODO Auto-generated method stub
						this.olap = olap;
						return this;
					}

					private SequencePair<DNASequence, NucleotideCompound> pairMatcher(String targetSeq, String querySeq) {
						try {			
							DNASequence target = new DNASequence(targetSeq, AmbiguityDNACompoundSet.getDNACompoundSet());
							DNASequence query  = new DNASequence(querySeq, AmbiguityDNACompoundSet.getDNACompoundSet());
							
							PairwiseSequenceAligner<DNASequence, NucleotideCompound> aligner = Alignments.getPairwiseAligner(
									target,
									query,
									PairwiseSequenceAlignerType.LOCAL,
									new SimpleGapPenalty(6,1), 
									SubstitutionMatrixHelper.getNuc4_4());
							SequencePair<DNASequence, NucleotideCompound> seqAlnPair = aligner.getPair();
							return seqAlnPair;
						} catch (CompoundNotFoundException e1) {
							// TODO Auto-generated catch block
							e1.printStackTrace();
						}
						return null;
					}
					
				}.init(olap));
			}
			
			bfs = new BFSShortestPath<String, OverlapEdge>(gfa);
		} catch(IOException e) {
			e.printStackTrace();
		}
		
		this.waitFor();
		myLogger.info("    finished.");
		
		return;
	}

	
	public GraphPath<String, OverlapEdge> getPath(String source, String sink) {
		return bfs.getPath(source, sink);
	}
	
	private ExecutorService executor;
	private BlockingQueue<Runnable> tasks = null;

	private void initial_thread_pool(int thread_num) {
		if(this.tasks!=null&&tasks.size()>0)
			throw new RuntimeException("Thread pool initialisation "
					+ "exception!!! Task queue is not emputy!!!");
		if(this.executor!=null&&!executor.isShutdown())
			throw new RuntimeException("Thread pool initialisation "
					+ "exception!!! Executor is on the fly!!!");
		
		tasks = new ArrayBlockingQueue<Runnable>(thread_num);
		executor = new ThreadPoolExecutor(thread_num, 
				thread_num, 
				1, 
				TimeUnit.SECONDS, 
				tasks, 
				new RejectedExecutionHandler(){
			@Override
			public void rejectedExecution(Runnable task,
					ThreadPoolExecutor arg1) {
				// TODO Auto-generated method stub
				try {
					tasks.put(task);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		});
	}
	
	private void initial_thread_pool() {
		this.initial_thread_pool(GFA.DEFAULT_NUM_THREADS);
	}
	
	private void waitFor() {
		// TODO Auto-generated method stub
		try {
			executor.shutdown();
			executor.awaitTermination(365, TimeUnit.DAYS);
		} catch (InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
	}

	public GFA(String seq_file, String graph_file, int k) {
		this.readSequence(seq_file);
		if(graph_file==null) 
			this.createDeBruijnGraph(seq_file, k);
		else this.readGraph(graph_file, k);
	}

	private void createDeBruijnGraph(String seq_file, int k) {
		// TODO Auto-generated method stub
		// a pseudo directed weighted graph
		// the nodes are contigs
		// edge weights are all 1.0
		// this is used to 
		// 1. calculate the (pseudo) distance between each pair of contigs
		// 2. as well as to extract scaffolds
		// NOTE: adjacent contigs will have a distance 1.0
		// (index     , contig    )
		// (index+size, rev_contig)

		List<Sequence> qry_seqs = Sequence.parseFastaFileWithRevCmpAsList(seq_file);
		Map<String, Set<String>> pref_mer = new HashMap<String, Set<String>>();
		Map<String, Set<String>> suff_mer = new HashMap<String, Set<String>>();
		String seq_str, kmer;
		int seq_ln;
		Sequence seq;
		for(int i=0; i<qry_seqs.size(); i++) {
			seq = qry_seqs.get(i);
			gfa.addVertex(seq.seq_sn());

			seq_str = seq.seq_str();
			seq_ln = seq.seq_ln();
			kmer = seq_str.substring(0, k);
			if(!pref_mer.containsKey(kmer))
				pref_mer.put(kmer, new HashSet<String>());
			pref_mer.get(kmer).add(seq.seq_sn());
			kmer = seq_str.substring(seq_ln-k, seq_ln);
			if(!suff_mer.containsKey(kmer))
				suff_mer.put(kmer, new HashSet<String>());
			suff_mer.get(kmer).add(seq.seq_sn());
		}
		Set<String> pref_set = pref_mer.keySet();
		Set<String> suff_set = suff_mer.keySet();
		for(String mer : suff_set) 
			// from a suff_mer to a pref_mer
			if(pref_set.contains(mer)) 
				for(String i : suff_mer.get(mer)) 
					for(String j : pref_mer.get(mer)) 
						gfa.setEdgeWeight(gfa.addEdge(i, j), 1.0);

		myLogger.info("Assembly graph "+ gfa.vertexSet().size()+" vetices and "+gfa.edgeSet().size()+" edges.");

		return;		
	}

	private void readGraph(String graph_file) {
		this.type = this.graphType(graph_file);

		switch (this.type) {
		case de_bruijn:
			throw new RuntimeException("need to provide k-mer size for reading de bruijn assembly graph!!!");
		case olc:
			this.readOLC(graph_file); 
			break;
		case overlap:
			this.readOverlapGraph(graph_file);
			break;
		default:
			throw new RuntimeException("errors occured reading the assembly graph file!!!");
		}
	}

	private void readGraph(String graph_file, int k) {
		this.type = this.graphType(graph_file);

		switch (this.type) {
		case de_bruijn:
			this.readDeBruijnGraph(graph_file, k);
			break;
		case olc:
			this.readOLC(graph_file);
			break;
		case overlap:
			this.readOverlapGraph(graph_file);
			break;
		default:
			throw new RuntimeException("errors occured reading the assembly graph file!!!");
		}	
	}

	private DirectedWeightedOverlapPseudograph<String> readDeBruijnGraph(String graph_file, int k) {
		// TODO Auto-generated method stub
		throw new RuntimeException("de bruijn assembly graph reader yet to implement!!!");
		// return null;
	}

	private void readOLC(String graph_file) {
		// TODO Auto-generated method stub
		// overlap-layout-consensus (OLC) graphical fragment assembly (GFA)
		myLogger.info("Loading assembly graph from file.");
		// 1. read file to initialise assembly graph
		// remember to make it symmetric 
		myLogger.info("    intialising assembly graph.");

		String source, target, overlap;
		try {
			BufferedReader br_gfa = Utils.getBufferedReader(graph_file);
			String line;
			String[] s;

			while( (line=br_gfa.readLine())!=null ) {
				s = line.split("\\s+");
				switch(s[0]) {
				case "#":
					// this is comment
					// print and skip this
					myLogger.info(line);
					break;
				case "H":
					// this is header
					// we check GFA version and expect to see V1
					if(s.length>1 && s[1].startsWith("VN:")) {
						String version = s[1].split(":")[2];
						if(!version.startsWith("1."))
							throw new RuntimeException("cannot handle GFA version "+version);
						else
							myLogger.info("      file is GFA version "+version);
					}
					break;
				case "S":
					// this is segment identifier
					// add this to graph vertex set
					gfa.addVertex(s[1]    );
					gfa.addVertex(s[1]+"'");
					rev.put(s[1], s[1]+"'");
					rev.put(s[1]+"'", s[1]);
					break;
				case "L":
					// this is link
					// add this to graph edge set
					source = s[1]+(s[2].equals("+")?"":"'");
					target = s[3]+(s[4].equals("+")?"":"'");
					overlap = s[5];
					if(overlap.equals("*")) {
						// overlap cigar not available?
						// TODO deal with this
						// currently throw a runtime exception 
						// need a overlap algorithm to detect end-to-end overlaps here
						throw new RuntimeException("need a CIGAR string for this overlap, this is yet to implement!!!");
					} else {
						// overlap cigar available
						// calculate overlap length from cigar string
						OverlapEdge e = gfa.addEdge(source, target);
						gfa.setEdgeWeight(e, 1.0);
						gfa.setEdgeCigar(e, overlap);
						gfa.setEdgeOverlap(e, Constants.getOlapFromCigar(overlap));
					}
					break;
				case "C":
					// this is containment
					// throw warning and skip this
					myLogger.warn("Containment was skipped.");
					break;
				case "P":
					// this is path
					// throw warning and skip this
					myLogger.warn("Path was skipped.");
					break;
				default:
					throw new RuntimeException("not a standard GFA 1 format!!!");
				}
			}
			br_gfa.close();
		} catch(Exception e) {
			e.printStackTrace();
			myLogger.error("Loading assembly graph errors!!!");
			System.exit(1);
		}
		myLogger.info("    finished.");

		// 2. make the graph symmetric and flooding 
		myLogger.info("    symmetrizing and flooding assembly graph.");
		int link_added = 0;
		int round = 0;
		myLogger.info("    Edge number "+gfa.edgeSet().size());

		while(true) {
			int link_new = 0;
			// need a copy of the edges to avoid concurrent modification exception
			// flooding
			myLogger.info("      Skip flooding. Yet to implement. ");

			Set<OverlapEdge> symmE = new HashSet<OverlapEdge>();
			symmE.addAll(gfa.edgeSet());
			// symmetrizing
			for(OverlapEdge olapE : symmE) {
				source = rev.get(gfa.getEdgeTarget(olapE));
				target = rev.get(gfa.getEdgeSource(olapE)); 
				if(!gfa.containsEdge(source, target)) {
					// mark a symmetric edge
					OverlapEdge e = gfa.addEdge(source, target);
					gfa.setEdgeWeight(e, 1.0);
					overlap = Constants.cgRevCmp(olapE.cigar);
					gfa.setEdgeCigar(e, overlap);
					gfa.setEdgeOverlap(e, Constants.getOlapFromCigar(overlap));
					++link_new;
				}
			}

			myLogger.info("      Round "+(++round)+": "+link_new+" new links were added to graph after symmetrizing and flooding.");

			if(link_new==0 || round>9) break;
			link_added += link_new;
		}
		myLogger.info("      "+link_added+" new links were added to graph after symmetrizing and flooding.");
		myLogger.info("    Edge number "+gfa.edgeSet().size());
		myLogger.info("    finished.");

		return;
	}
	
	private Graph_type graphType(String asm_graph) {
		// TODO Auto-generated method stub

		// TODO distinguish between the de-bruijn, overlap-layout-consensus and overlap assembly graph
		return Graph_type.overlap;
	}

	public String graph_type() {
		if(this.type==null) return null;
		switch(this.type) {
		case de_bruijn:
			return "de_bruijn";
		case olc:
			return "overlap_layout_consensus";
		case overlap:
			return "overlap";
		default:
			return null;
		}
	}

	private void readSequence(String seq_file) {
		List<Sequence> seqs = Sequence.parseFastaFileAsList(seq_file);
		for(Sequence seq : seqs) {
			seg.put(seq.seq_sn(),                         seq );
			seg.put(seq.seq_sn()+"'", Sequence.revCompSeq(seq));
			rev.put(seq.seq_sn(), seq.seq_sn()+"'");
			rev.put(seq.seq_sn()+"'", seq.seq_sn());
		}
		return;
	}

	public Set<String> vertexSet() {
		return this.gfa.vertexSet();
	}

	public Set<OverlapEdge> edgeSet() {
		return this.gfa.edgeSet();
	}

	public Set<String> sequenceSet() {
		return this.seg.keySet();
	}

	public Sequence getSequence(String sq_sn) {
		return this.seg.get(sq_sn);
	}

	public boolean containsVertex(String v) {
		return this.gfa.containsVertex(v);
	}

	public boolean containsEdge(String source, String target) {
		return this.gfa.containsEdge(source, target);
	}

	public Set<OverlapEdge> edgesOf(String v) {
		return this.gfa.edgesOf(v);
	}

	public Set<OverlapEdge> getAllEdges(String source, String target) {
		return this.gfa.getAllEdges(source, target);
	}

	public OverlapEdge getEdge(String source, String target) {
		return this.gfa.getEdge(source, target);
	}

	public String getEdgeCigar(String source, String target) {
		return this.gfa.getEdgeCigar(this.gfa.getEdge(source, target));
	}

	public String getEdgeCigar(OverlapEdge e) {
		return this.gfa.getEdgeCigar(e);
	}

	public double getEdgeWeight(OverlapEdge e) {
		return this.gfa.getEdgeWeight(e);
	}

	public Set<OverlapEdge> incomingEdgesOf(String v) {
		return this.gfa.incomingEdgesOf(v);
	}

	public Set<OverlapEdge> outgoingEdgesOf(String v) {
		return this.gfa.incomingEdgesOf(v);
	}

	public Set<String> getSeqSet() {
		return this.seg.keySet();
	}

	public Sequence getSeq(String seqid) {
		return this.seg.get(seqid);
	}

	public String getRev(String seqid) {
		return this.rev.get(seqid);
	}

	public Sequence getRevSeq(String seqid) {
		return this.seg.get(this.rev.get(seqid));
	}

	public Map<String, Sequence> getSequenceMap() {
		return this.seg;
	}
}
