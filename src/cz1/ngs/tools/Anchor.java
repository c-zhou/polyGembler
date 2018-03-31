package cz1.ngs.tools;

import java.io.File;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.JFrame;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.TreeBidiMap;
import org.jgraph.JGraph;
import org.jgrapht.Graph;
import org.jgrapht.ext.JGraphModelAdapter;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.jgrapht.traverse.BreadthFirstIterator;

import com.google.common.collect.DiscreteDomain;
import com.google.common.collect.Range;
import com.google.common.collect.TreeRangeSet;

import cz1.ngs.model.AlignmentSegment;
import cz1.ngs.model.GFA;
import cz1.ngs.model.OverlapEdge;
import cz1.ngs.model.SAMSegment;
import cz1.ngs.model.Sequence;
import cz1.ngs.model.TraceableDirectedWeightedPseudograph;
import cz1.ngs.model.TraceableEdge;
import cz1.ngs.model.TraceableVertex;
import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Executor;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Anchor extends Executor {
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -s/--subject            Subject/reference sequences file in FASTA format.\n"
						+ " -q/--query              Query sequences to anchor in FASTA format.\n"
						+ " -a/--align              Alignment file of query sequences to the subject sequences. \n"
						+ "                         IMPORTANT: The alignment file need to be grouped by the query sequences,\n"
						+ "                         i.e., the program assume all the alignment records for a query sequence \n"
						+ "                         can be continuously read from the alignment file.\n"
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
						+ " -k/--kmer-size          Kmer size used for assembly graph construction. Which is supposed \n"
						+ "                         to be the same to the k-mer size used for contig construction. \n"
						+ " -i/--min-identity       Minimum identity between the query and subject sequences \n"
						+ "                         for an alignment record to consider (default 0.90).\n"
						+ " -f/--min-fraction       Minimum alignment fraction of the query sequence (default 0.5).\n"
						+ " -di/--diff-identity     Threshold of the difference of the identity between the primary and secondary \n"
						+ "                         alignments. If the difference is smaller than this value, the query \n"
						+ "                         sequence will be regarded as duplications. Otherwise, the secondary \n"
						+ "                         alignments will be discared (default 0.01).\n"
						+ " -df/--diff-fraction     Threshold of the difference of the alignment fraction between the primary and \n"
						+ "                         secondary alignments. If the difference is smaller than this value, the query \n"
						+ "                         sequence will be regarded as duplications. Otherwise, the secondary \n"
						+ "                         alignments will be discared (default 0.05).\n"
						+ " -t/--threads            Number of threads to use (default 16).\n"
						+ " -d/--debug              Debugging mode will have extra information printed out.\n"
						+ " -dd/--debug-debug       Debugging mode will have more information printed out than -d mode.\n"
						+ " -o/--out-prefix         Prefix of the output files.\n"
						+ "\n");	
	}

	private String subject_file = null;
	private String query_file = null;
	private String align_file = null;
	private String asm_graph = null;
	private double min_ident = 90;        // ignore alignment records with identity smaller than this
	private double min_frac = 0.5;        // ignore alignment records with completeness smaller than this unless 
	                                      // the segment is greater than min_alen
	private double min_alen = 300;        // keep the alignment record if the segment is no smaller than this
	                                      // without considering the completeness
	private double diff_ident = 0.05;     // keep the secondary alignment if the identity score difference is no greater than this
	private double diff_frac = 0.05;      // keep the secondary alignment if the completeness difference is no greater than this
	private int min_overlap = 10;         // minimum overlap length
	private double collinear_shift = 0.5; // maximum shift distance for two collinear alignment segments - 50% of the smaller segment size
	private int kmer_size = -1;           // kmer size to construct the de bruijn assembly graph

	private int num_threads = Runtime.getRuntime().availableProcessors();
	private String out_prefix = null;
	private boolean debug = false;
	private boolean ddebug = false;
	
	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-s", "--subject", true);
			myArgsEngine.add("-q", "--query", true);
			myArgsEngine.add("-a", "--align", true);
			myArgsEngine.add("-g","--graph", true);
			myArgsEngine.add("-k", "--kmer-size", true);
			myArgsEngine.add("-i", "--min-identity", true);
			myArgsEngine.add("-f", "--min-fraction", true);
			myArgsEngine.add("-di", "--diff-identity", true);
			myArgsEngine.add("-df", "--diff-fraction", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-d", "--debug", false);
			myArgsEngine.add("-dd", "--debug-debug", false);
			myArgsEngine.add("-o", "--out-prefix", true);
			myArgsEngine.parse(args);
		}
		
		if (myArgsEngine.getBoolean("-s")) {
			this.subject_file = myArgsEngine.getString("-s");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the subject/reference file.");
		}

		if (myArgsEngine.getBoolean("-q")) {
			this.query_file = myArgsEngine.getString("-q");
			if(!this.query_file.endsWith(".fasta")) 
				throw new RuntimeException("Query file need to be in FASTA format!!!");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the query file.");
		}
		
		if (myArgsEngine.getBoolean("-a")) {
			this.align_file = myArgsEngine.getString("-a");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the alignment file of query sequences to subject sequences.");
		}
				
		if (myArgsEngine.getBoolean("-g")) {
			this.asm_graph = myArgsEngine.getString("-g");
		}
		
		if (myArgsEngine.getBoolean("-k")) {
			this.kmer_size = Integer.parseInt(myArgsEngine.getString("-k"));
		}
		
		if (myArgsEngine.getBoolean("-o")) {
			this.out_prefix = myArgsEngine.getString("-o");
			if(new File(out_prefix+".fa").exists() || 
					new File(out_prefix+".map").exists() || 
					new File(out_prefix+"_unplaced.map").exists()) {
				throw new RuntimeException("Output files exist. Please specify a different name.");
			}
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the prefix of output files.");
		}
		
		if (myArgsEngine.getBoolean("-i")) {
			this.min_ident = 100*Math.min(1.0, Double.parseDouble(myArgsEngine.getString("-i")));
		}
		
		if (myArgsEngine.getBoolean("-f")) {
			this.min_frac = Double.parseDouble(myArgsEngine.getString("-f"));
		}
		
		if (myArgsEngine.getBoolean("-di")) {
			this.diff_ident = Double.parseDouble(myArgsEngine.getString("-di"));
		}
		
		if (myArgsEngine.getBoolean("-t")) {
			int t = Integer.parseInt(myArgsEngine.getString("-t"));
			if(t<this.num_threads) this.num_threads = t;
			Constants.omp_threads = this.num_threads;
			myLogger.info("OMP_THREADS = "+this.num_threads);
		}
		
		if (myArgsEngine.getBoolean("-df")) {
			this.diff_frac = Double.parseDouble(myArgsEngine.getString("-df"));
		}
		
		if (myArgsEngine.getBoolean("-d")) {
			this.debug = true;
		}
		
		if (myArgsEngine.getBoolean("-dd")) {
			this.debug  = true;
			this.ddebug = true;
		}
	}

	private final static int max_clip = 30; // maximum clip of query alignment allowed
	private final static int gap_buff = 30; // buffer size for subject/reference sequences gap clips
	private final static int min_gap  = 10; // take this if estimated gap size is smaller than this
	private final static int max_gap  = 100; // take this if estimated gap size is larger than this
	private Map<String, Sequence> qry_seqs;
	private Map<String, Sequence> sub_seqs;
	private Map<String, TreeRangeSet<Integer>> sub_gaps;
	
	private final static int match_score  = 1;
	private final static int clip_penalty = 1;
	private final static int hc_gap  = 100000;
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		
		// read assembly graph file
		final GFA gfa = new GFA(query_file, asm_graph);
		qry_seqs = gfa.getSequenceMap();
		sub_seqs = Sequence.parseFastaFileAsMap(subject_file);
		
		myLogger.info("  GFA vertices: "+gfa.vertexSet().size());
		myLogger.info("  GFA edges   : "+gfa.edgeSet().size()  );
		// myLogger.info("  GFA edges --- ");
		// for(OverlapEdge olap : gfa.edgeSet()) 
		//	 myLogger.info(olap.olapInfo().toString());
		
		// find 'N/n's in subject/reference sequences
		// which could have impact on parsing the blast records
		sub_gaps = new HashMap<String, TreeRangeSet<Integer>>();

		for(Map.Entry<String, Sequence> entry : sub_seqs.entrySet()) {
			String seq_sn = entry.getKey();
			String seq_str = entry.getValue().seq_str();

			final TreeRangeSet<Integer> tmp_rangeSet = TreeRangeSet.create();
			for(int j=0; j<seq_str.length(); j++) {
				if(seq_str.charAt(j)=='N'||seq_str.charAt(j)=='n')
					// blast record is 1-based closed coordination
					tmp_rangeSet.add(Range.closed(j+1, j+1).
							canonical(DiscreteDomain.integers()));
			}
			int seq_ln = seq_str.length();
			final TreeRangeSet<Integer> range_set = TreeRangeSet.create();
			for(Range<Integer> range : tmp_rangeSet.asRanges()) {
				int lowerend = range.hasLowerBound() ? Math.max(0, range.lowerEndpoint()-gap_buff) : 0;
				int upperend = range.hasUpperBound() ? Math.min(seq_ln, range.upperEndpoint()+gap_buff-1) : seq_ln;
				range_set.add( Range.closed(lowerend, upperend).canonical(DiscreteDomain.integers()) );
			}
			sub_gaps.put(seq_sn, range_set);
		}
		
		// read alignment file and place the query sequences
		final Map<String, Set<SAMSegment>> initPlace = new HashMap<String, Set<SAMSegment>>();
		final Map<String, List<SAMSegment>> initPseudoAssembly = new HashMap<String, List<SAMSegment>>();
		for(String sub_seq : sub_seqs.keySet()) initPseudoAssembly.put(sub_seq, new ArrayList<SAMSegment>());
		
		try {
			final SamReaderFactory factory =
					SamReaderFactory.makeDefault()
					.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
							SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
					.validationStringency(ValidationStringency.SILENT);
			final SamReader in1 = factory.open(new File(align_file));
			final SAMRecordIterator iter1 = in1.iterator();
			
			
			String qry;
			int qry_ln;
			double min_aln;
			final List<SAMSegment> buff = new ArrayList<SAMSegment>();
			SAMRecord rc = iter1.next();
			
			while(rc!=null) {
				qry = rc.getReadName();
				qry_ln = qry_seqs.get(qry).seq_ln();			
				
				buff.clear();
                if(!rc.getReadUnmappedFlag())
                    buff.add(SAMSegment.samRecord(rc, true, qry_ln));

				
				while( (rc=iter1.next())!=null
						&&
						rc.getReadName().equals(qry) ) {
					buff.add(SAMSegment.samRecord(rc, true, qry_ln));
				}
			    
                if(buff.isEmpty()) continue;

				min_aln = 0.9*buff.get(0).qlength();
				
				// keep alignment fragment that has qual>0
				Set<SAMSegment> init_f = new HashSet<SAMSegment>();
				Set<SAMSegment> init_r = new HashSet<SAMSegment>();
				for(SAMSegment record : buff) {
					if(record.qual()==0&&record.qlength()<min_aln) 
						continue;
					if(record.qseqid().equals(qry)) 
						init_f.add(record);
					else init_r.add(record);
					initPseudoAssembly.get(record.sseqid()).add(record);
				}
				if(!init_f.isEmpty()) initPlace.put(qry,     init_f);
				if(!init_r.isEmpty()) initPlace.put(qry+"'", init_r);
			}
			iter1.close();
			in1.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Collections.sort(initPseudoAssembly.get("1_pilon"), new AlignmentSegment.SubjectCoordinationComparator());
		
		// if(debug) {
		//	for(SAMSegment record : initPseudoAssembly.get("1_pilon")) {
		//		System.out.println(record.qseqid()+":"+record.sstart()+"-"+record.send());
		//	}
		// }
		
		final Set<SAMSegment> contained = new HashSet<SAMSegment>();
		final Set<SAMSegment> placed    = new HashSet<SAMSegment>();
		final int flank_size = 10000;
		int distance;

		for(String sub_seq : sub_seqs.keySet()) {

            // sub_seq = "Chr10";
            if(sub_seq.equals("Chr00")) continue;
            
            myLogger.info(">>>>>>>>>>>>>"+sub_seq+"<<<<<<<<<<<<<<<<");
            
            final List<SAMSegment> seq_by_sub = initPseudoAssembly.get(sub_seq);
		    Collections.sort(seq_by_sub, new AlignmentSegment.SubjectCoordinationComparator());
        
            placed.clear();
            
			int nSeq = seq_by_sub.size();
			double edge_penalty, edge_score;
			SAMSegment root_seq, source_seq, target_seq;
			Set<SAMSegment> target_seqs;
			Set<OverlapEdge> outgoing;
			TraceableEdge edge;
			String root_seqid, source_seqid, target_seqid;
			TraceableVertex<String> root_vertex, source_vertex, target_vertex;
			Deque<SAMSegment> deque = new ArrayDeque<SAMSegment>();
			final List<TraceableVertex<String>> traceable = new ArrayList<TraceableVertex<String>>();

			for(int i=0; i<nSeq; i++) {
				
				root_seq = seq_by_sub.get(i);
				root_seqid = root_seq.qseqid();
				
				if(placed.contains(root_seq)) 
					continue;
				
				final TraceableDirectedWeightedPseudograph<String> razor = 
						new TraceableDirectedWeightedPseudograph<String>(TraceableEdge.class);
				

				// final ListenableDirectedWeightedGraph<TraceableVertex<String>, DefaultWeightedEdge> razor = 
				//		new ListenableDirectedWeightedGraph<TraceableVertex<String>, DefaultWeightedEdge>(DefaultWeightedEdge.class);
				
				// JGraphModelAdapter<TraceableVertex<String>, DefaultWeightedEdge> jgAdapter = 
				//		 new JGraphModelAdapter<TraceableVertex<String>, DefaultWeightedEdge>(razor);

				// JGraph jgraph = new JGraph(jgAdapter);
				
				deque.clear();
				deque.push(root_seq);
				contained.clear();
				
				while(!deque.isEmpty()) {
					
					source_seq = deque.pop();
					source_seqid = source_seq.qseqid();
					
					if(contained.contains(source_seq))
						continue;
					contained.add(source_seq);

					source_vertex = new TraceableVertex<String>(source_seqid);
					source_vertex.setSAMSegment(source_seq);
					
					if(!razor.containsVertex(source_vertex)) 
						razor.addVertex(source_vertex);

					outgoing = gfa.outgoingEdgesOf(source_seqid);
					
					for(OverlapEdge out : outgoing) {
						target_seqid = gfa.getEdgeTarget(out);
						if(!initPlace.containsKey(target_seqid))
							continue;
						target_seqs = initPlace.get(target_seqid);
						
						distance = Integer.MAX_VALUE;
						target_seq = null;
						for(SAMSegment seq : target_seqs) {
							int d = AlignmentSegment.sdistance(source_seq, seq);	
							if(d<distance) {
								distance = d;
								target_seq = seq;
							}
						}
						if(distance<=flank_size) {
							target_vertex = new TraceableVertex<String>(target_seqid);
							target_vertex.setSAMSegment(target_seq);
							
							if(!razor.containsVertex(target_vertex))
								razor.addVertex(target_vertex);
							
							if(razor.containsEdge(source_vertex, target_vertex))
								continue;
							
							edge = razor.addEdge(source_vertex, target_vertex);
							// calculate edge weight
							// higher weight edges are those,
							
							/****
							//       1.  large/long alignment segments vertices
							// TODO: 2*. small gaps on the reference
							edge_weight = qry_seqs.get(source_seqid).seq_ln()+
									qry_seqs.get(target_seqid).seq_ln()-
									gfa.getEdge(source_seqid, target_seqid).olap();
							**/

							// TODO: 1*. large/long alignment segments vertices
							//       2.  small gaps on the reference
							edge_penalty = AlignmentSegment.sdistance(source_seq, target_seq);
							edge.setPenalty(edge_penalty);
							
							edge_score = qry_seqs.get(source_seqid).seq_ln()+
									qry_seqs.get(target_seqid).seq_ln()-
									gfa.getEdge(source_seqid, target_seqid).olap();
							edge.setScore(edge_score);
							
							deque.push(target_seq);
						}
					}
				}
				if(ddebug) myLogger.info(root_seqid+" "+razor.vertexSet().size()+" "+razor.edgeSet().size()+" done");

				// JFrame frame = new JFrame();
			    // frame.getContentPane().add(jgraph);
			    // frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			    // frame.pack();
			    // frame.setVisible(true);

			    // "pseudo"-DFS to find the route with the highest score
                final Map<String, TraceableVertex<String>> razv_map  =new HashMap<String, TraceableVertex<String>>();
                for(TraceableVertex<String> v : razor.vertexSet()) razv_map.put(v.getId(), v);
                
                // we use a bidirectional hashmap to simulate the deque
                // this is because we may need to do deletions
		        // Deque<TraceableVertex<String>> queue = new ArrayDeque<TraceableVertex<String>>();
		        final TreeBidiMap<Long, TraceableVertex<String>> bidiQ = new TreeBidiMap<Long, TraceableVertex<String>>();

                root_vertex = razv_map.get(root_seqid);
		        root_vertex.setSAMSegment(root_seq);
		        root_vertex.setScore(qry_seqs.get(root_seqid).seq_ln());
		        root_vertex.setPenalty(0);
                root_vertex.setStatus(true);
		        
                bidiQ.put(0L, root_vertex);
                double max_ws = Double.NEGATIVE_INFINITY,
                		source_penalty, target_penalty, source_score, target_score, 
                		penalty, score, target_ws, source_ws, ws;
				int source_ln;
		        Set<TraceableEdge> out_edges;
				TraceableVertex<String> opt_vertex = null;
				long sizeQ;
				boolean isLeaf;
				
				if(ddebug)
					for(TraceableEdge e : razor.edgeSet()) 
						myLogger.info(e.toString()+"("+razor.getEdgeSource(e).getSAMSegment().toString()+"|"+razor.getEdgeTarget(e).getSAMSegment().toString()+"|"+e.getScore()+"-"+e.getPenalty()+")");

		        while(!bidiQ.isEmpty()) {
		        	sizeQ = bidiQ.lastKey();
		        	source_vertex = bidiQ.get(sizeQ);
		        	bidiQ.remove(sizeQ);
		        	
		        	source_ln = qry_seqs.get(source_vertex.getId()).seq_ln();
		        	source_score = source_vertex.getScore()-source_ln;
		        	source_penalty = source_vertex.getPenalty();
		        	source_ws = source_score-source_penalty;
		        	
		        	isLeaf = true;
		        	out_edges = razor.outgoingEdgesOf(source_vertex);
		        	for(TraceableEdge out : out_edges) {
		        		// this is not right because graph edges are immutable?
		        		// target_vertex = razor.getEdgeTarget(out);
		        		target_vertex = razv_map.get(razor.getEdgeTarget(out).getId());
		        		target_score = target_vertex.getScore();
		        		target_penalty = target_vertex.getPenalty();
		        		target_ws = target_score-target_penalty;
		        		
		        		edge_penalty = out.getPenalty();
		        		penalty = source_penalty+edge_penalty;
		        		edge_score = out.getScore();
		        		score = source_score+edge_score;
		        		ws = score-penalty;
		        		
		        		if( edge_penalty>flank_size || 
		        				target_vertex.getStatus() && 
		        				(ws<=target_ws ||
		        				isLoopback(razor, source_vertex, target_vertex)) ) 
		        			continue;

		        		isLeaf = false;
		        		target_vertex.setBackTrace(source_vertex);
		        		target_vertex.setScore(score);
		        		target_vertex.setPenalty(penalty);
		        		target_vertex.setStatus(true);
		        		
                        bidiQ.put(sizeQ++, target_vertex);
                    }
		        	
		        	if(isLeaf && source_ws>max_ws) {
		        		penalty = source_vertex.getPenalty();
		        		score   = source_vertex.getScore();
		        		max_ws = source_ws;
		        		opt_vertex  = source_vertex;
		        		
		        		if(ddebug) {
		        			
		        			String trace = opt_vertex.toString()+":"+
		        					opt_vertex.getSAMSegment().sstart()+"-"+
		        					opt_vertex.getSAMSegment().send()+"("+
		        					opt_vertex.getScore()+"-"+
		        					opt_vertex.getPenalty()+")";
		        			
		        			TraceableVertex<String> optx = opt_vertex;
		        			while( (optx = optx.getBackTrace())!=null ) {
		        				trace += ","+optx.toString()+":"+
		        						optx.getSAMSegment().sstart()+"-"+
		        						optx.getSAMSegment().send()+"("+
		        						optx.getScore()+"-"+
		        						optx.getPenalty()+")";
		        			}
		        			myLogger.info("trace back ["+score+", "+penalty+"]: "+trace);
		        		}
		        	}
		        }
		        
		        traceable.add(opt_vertex);
		        
		        Set<TraceableVertex<String>> optx = new HashSet<TraceableVertex<String>>();
		        optx.add(opt_vertex);
		        while( (opt_vertex = opt_vertex.getBackTrace())!=null ) optx.add(opt_vertex);
		    
		        for(TraceableVertex<String> v : optx) placed.add(v.getSAMSegment());
			}
			
			// sort traceable by size
			Collections.sort(traceable, new Comparator<TraceableVertex<String>>() {

				@Override
				public int compare(TraceableVertex<String> t0, TraceableVertex<String> t1) {
					// TODO Auto-generated method stub
					return Double.compare(t1.getScore(), t0.getScore());
				}
			});

			if(debug) {
				for(TraceableVertex<String> opt_vertex : traceable) {

					double score = opt_vertex.getScore();
					double penalty = opt_vertex.getPenalty(); 
					
					String trace = opt_vertex.toString()+":"+
        					opt_vertex.getSAMSegment().sstart()+"-"+
        					opt_vertex.getSAMSegment().send()+"("+
        					opt_vertex.getScore()+"-"+
        					opt_vertex.getPenalty()+")";

					while( (opt_vertex = opt_vertex.getBackTrace())!=null ) {
						trace += ","+opt_vertex.toString()+":"+
								opt_vertex.getSAMSegment().sstart()+"-"+
								opt_vertex.getSAMSegment().send()+"("+
								opt_vertex.getScore()+"-"+
								opt_vertex.getPenalty()+")";
					}
					myLogger.info("trace back ["+score+", "+penalty+"]: "+trace);
				}
			}

			// we generate a compound alignment record for each traceable
			for(TraceableVertex<String> opt_vertex : traceable) {
				
			}
		}
	}

	private boolean isLoopback(DirectedWeightedPseudograph<TraceableVertex<String>, TraceableEdge> graph,
			TraceableVertex<String> source,
			TraceableVertex<String> target) {
		// TODO Auto-generated method stub
		if(source.equals(target)) return true;
		while( (source = source.getBackTrace())!=null )
        	if(source.equals(target)) return true;
		return false;
	}
}















