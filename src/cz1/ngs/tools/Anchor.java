package cz1.ngs.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;
import org.jgrapht.Graph;
import org.jgrapht.ListenableGraph;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.alg.CycleDetector;
import org.jgrapht.ext.JGraphXAdapter;
import org.jgrapht.graph.AsSubgraph;
import org.jgrapht.graph.DefaultListenableGraph;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.DirectedWeightedPseudograph;

import com.google.common.collect.DiscreteDomain;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;
import com.mxgraph.layout.mxOrganicLayout;
import com.mxgraph.swing.mxGraphComponent;

import cz1.ngs.model.AlignmentSegment;
import cz1.ngs.model.AlignmentSegmentDirectedWeightedPseudograph;
import cz1.ngs.model.CompoundAlignmentSegment;
import cz1.ngs.model.DirectedWeightedOverlapPseudograph;
import cz1.ngs.model.GFA;
import cz1.ngs.model.OverlapEdge;
import cz1.ngs.model.OverlapResult;
import cz1.ngs.model.SAMSegment;
import cz1.ngs.model.Sequence;
import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Executor;
import cz1.util.Utils;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Anchor extends Executor implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 8177690774607523148L;

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
						+ " -1                      BAM file for the alignment records of the 1st reads in read pairs.\n"
						+ " -2                      BAM file for the alignment records of the 2nd reads in read pairs.\n"
						+ " -12                     Configuration file for BAM files.\n"
						+ " -mp/--mate-pair         Mate-pair library.\n"
						+ " -k/--kmer-size          Kmer size used for assembly graph construction. Which is supposed \n"
						+ "                         to be the same to the k-mer size used for contig construction. \n"
						+ " -i/--min-identity       Minimum identity between the query and subject sequences \n"
						+ "                         for an alignment record to consider (default 0.90).\n"
						+ " -f/--min-fraction       Minimum alignment fraction of the query sequence (default 0.5).\n"
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
	private double min_frac = 0.8;        // ignore alignment records with completeness smaller than this unless 
	// the segment is greater than min_alen
	private double min_alen = 0;          // keep the alignment record if the segment is no smaller than this
	// without considering the completeness
	private int min_overlap = 10;         // minimum overlap length
	private double collinear_shift = 1.0; // maximum shift distance for two collinear alignment segments - 50% of the smaller segment size
	private int kmer_size = -1;           // kmer size to construct the de bruijn assembly graph
	private int minQual = 1;
	
	private int num_threads = Runtime.getRuntime().availableProcessors();
	private String out_prefix = null;
	private boolean debug  = false;
	private boolean ddebug = false;
	private int[] inst;
	private int[] maxInst;
	private String[][] bamList;
	private boolean[] matePair; // mate-pair library
	
	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-s", "--subject", true);
			myArgsEngine.add("-q", "--query", true);
			myArgsEngine.add("-a", "--align", true);
			myArgsEngine.add("-g","--graph", true);
			myArgsEngine.add("-1", null, true);
			myArgsEngine.add("-2", null, true);
			myArgsEngine.add("-12", null, true);
			myArgsEngine.add("-mp", "--mate-pair", false);
			myArgsEngine.add("-k", "--kmer-size", true);
			myArgsEngine.add("-i", "--min-identity", true);
			myArgsEngine.add("-f", "--min-fraction", true);
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
		
		if (!(myArgsEngine.getBoolean("-1")&&myArgsEngine.getBoolean("-2"))&&!myArgsEngine.getBoolean("-12")) {
			printUsage();
			throw new IllegalArgumentException("Please specify the alignment file(s) using (-1, -2) or -12 option.");
		}
		
		if ((myArgsEngine.getBoolean("-1")||myArgsEngine.getBoolean("-2"))&&myArgsEngine.getBoolean("-12")) {
			printUsage();
			throw new IllegalArgumentException("Options (-1, -2) and -12 are exclusive.");
		}
		
		if (myArgsEngine.getBoolean("-1")!=myArgsEngine.getBoolean("-2")) {
			printUsage();
			throw new IllegalArgumentException("Options (-1, -2) need to be used in pair.");
		}
		
		if (myArgsEngine.getBoolean("-1")&&myArgsEngine.getBoolean("-2")) {
			this.bamList = new String[1][2];
			this.bamList[0][0] = myArgsEngine.getString("-1");
			this.bamList[0][1] = myArgsEngine.getString("-2");
			this.matePair = new boolean[1];
			if(myArgsEngine.getBoolean("-mp")) {
				this.matePair[0] = true;
				myLogger.info("a mate-pair library.");
			}
			this.inst = new int[1];
			this.maxInst = new int[1];
			if (myArgsEngine.getBoolean("-f")) {
				this.inst[0] = Integer.parseInt(myArgsEngine.getString("-f"));
				this.maxInst[0] = 2*this.inst[0];
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify insert size of the library (-f).");
			}
		} else {
			if(myArgsEngine.getBoolean("-mp"))
				myLogger.info("mate-pair (-mp) indicator ignored.");
			if(myArgsEngine.getBoolean("-f"))
				myLogger.info("insert size (-f) option ignored.");
		}
		
		if (myArgsEngine.getBoolean("-12")) {
			try {
				BufferedReader br = Utils.getBufferedReader(myArgsEngine.getString("-12").trim());
				final List<String> f1 = new ArrayList<String>();
				final List<String> f2 = new ArrayList<String>();
				final List<Boolean> rev = new ArrayList<Boolean>();
				final List<Integer> ins = new ArrayList<Integer>();
				String line;
				String[] s;
				while( (line = br.readLine()) != null) {
					line = line.trim();
					s = line.split("\\s+");
					if(s.length>0) {
						f1.add(s[0].trim());
						f2.add(s[1].trim());
						rev.add(s[2].trim().equals("0")?false:true);
						ins.add(Integer.parseInt(s[3].trim()));
					}
				}
				br.close();
				
				final int n = f1.size();
				this.bamList = new String[n][2];
				this.matePair = new boolean[n];
				this.inst = new int[n];
				this.maxInst = new int[n];
				for(int i=0; i<n; i++) {
					this.bamList[i][0] = f1.get(i);
					this.bamList[i][1] = f2.get(i);
					this.matePair[i] = rev.get(i);
					this.inst[i] = ins.get(i);
					this.maxInst[i] = 2*this.inst[i];
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				myLogger.error("configuration file with invalid format.");
				e.printStackTrace();
			}
			
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

		if (myArgsEngine.getBoolean("-t")) {
			int t = Integer.parseInt(myArgsEngine.getString("-t"));
			if(t<this.num_threads) this.num_threads = t;
			this.THREADS = t;
			Constants.omp_threads = this.num_threads;
			myLogger.info("OMP_THREADS = "+this.num_threads);
		}

		if (myArgsEngine.getBoolean("-d")) {
			this.debug = true;
		}

		if (myArgsEngine.getBoolean("-dd")) {
			this.debug  = true;
			this.ddebug = true;
		}
	}

	private Map<String, Sequence> qry_seqs;
	private Map<String, Sequence> sub_seqs;
	private final static int hc_gap  = 100000;
	private final static int as_gap  = 100000;

	@Override
	public void run() {
		// this.blastReader();
		// this.run_link();
		
		// this.run1();
		// this.run2();
		// this.run3();
		this.run4();
	}
	
	
	GFA gfa;
	final Map<String, List<Set<SAMSegment>>> initPlace = new HashMap<String, List<Set<SAMSegment>>>();
	final Map<String, List<SAMSegment>> highQualSeg = new HashMap<String, List<SAMSegment>>();
	SAMSequenceDictionary dict = null;
	final long[] elapsed = new long[10];
	
	private void run4() {
		// TODO Auto-generated method stub
		elapsed[0]  = System.nanoTime();
				
		// read assembly graph file
		gfa = new GFA(query_file, asm_graph);
		myLogger.info("  GFA vertices: "+gfa.vertexSet().size());
		myLogger.info("  GFA edges   : "+gfa.edgeSet().size()  );

		qry_seqs = gfa.getSequenceMap();
		sub_seqs = Sequence.parseFastaFileAsMap(subject_file);

		// read alignment file and place the query sequences
		
		try {
			final SamReader in1 = factory.open(new File(align_file));
			final SAMRecordIterator iter1 = in1.iterator();
			dict  =  in1.getFileHeader().getSequenceDictionary();
			final int nSubSeq = dict.getSequences().size();
			if(nSubSeq!=sub_seqs.size()) throw new RuntimeException("!!!");

			myLogger.info("--"+nSubSeq+" Subject sequences");
			for(int i=0; i<nSubSeq; i++) 
				myLogger.info("  "+dict.getSequences().get(i).getSequenceName());
			myLogger.info("----------------");

			for(int i=0; i<nSubSeq; i++) 
				highQualSeg.put(dict.getSequences().get(i).getSequenceName(), new ArrayList<>());
			
			String qry;
			final List<SAMSegment> buff = new ArrayList<SAMSegment>();
			SAMRecord rc = iter1.next();

			int nq = 0;

			while(rc!=null) {
				if(++nq%100000==0) myLogger.info("#query sequences processed "+nq);

				qry = rc.getReadName();	

				buff.clear();
				if(!rc.getReadUnmappedFlag())
					buff.add(SAMSegment.samRecord(rc, true));

				while( (rc=iter1.next())!=null
						&&
						rc.getReadName().equals(qry) ) {
					if(!rc.getReadUnmappedFlag())
						buff.add(SAMSegment.samRecord(rc, true));
				}

				if(buff.isEmpty()) continue;

				for(SAMSegment seg : buff) 
					if(seg.qual()>=this.minQual)
						highQualSeg.get(seg.sseqid()).add(seg);
				
				Set<SAMSegment> init_f = new HashSet<SAMSegment>();
				Set<SAMSegment> init_r = new HashSet<SAMSegment>();
				for(SAMSegment record : buff) {
					if(record.qlength()<min_alen) 
						continue;
					if(record.qseqid().equals(qry)) 
						init_f.add(record);
					else init_r.add(record);
				}
				if(!initPlace.containsKey(qry)) {
					final List<Set<SAMSegment>> list_f = new ArrayList<Set<SAMSegment>>(nSubSeq);
					for(int i=0; i<nSubSeq; i++) list_f.add(new HashSet<SAMSegment>());
					initPlace.put(qry,     list_f);
					final List<Set<SAMSegment>> list_r = new ArrayList<Set<SAMSegment>>(nSubSeq);
					for(int i=0; i<nSubSeq; i++) list_r.add(new HashSet<SAMSegment>());
					initPlace.put(qry+"'", list_r);
				}

				if(!init_f.isEmpty()) {
					final List<Set<SAMSegment>> list_f = initPlace.get(qry);
					for(SAMSegment seg : init_f) 
						list_f.get(dict.getSequenceIndex(seg.sseqid())).add(seg);
				}
				if(!init_r.isEmpty()) {
					final List<Set<SAMSegment>> list_r = initPlace.get(qry+"'");
					for(SAMSegment seg : init_r) 
						list_r.get(dict.getSequenceIndex(seg.sseqid())).add(seg);
				}
			}
			iter1.close();
			in1.close();

			for(List<SAMSegment> segs : highQualSeg.values()) 
				Collections.sort(segs, new AlignmentSegment.SubjectCoordinationComparator());
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		int nQrySeq = 0, totSAMSeg = 0;
		long totQryLen = 0;
		final Set<String> ignFRQry = new HashSet<String>();
		for(final String qseqid : initPlace.keySet()) {
			final List<Set<SAMSegment>> list = initPlace.get(qseqid);
			int n = 0;
			for(final Set<SAMSegment> segs : list) n += segs.size();
			if(n>0) {
				++nQrySeq;
				totSAMSeg+=n;
				totQryLen+=qry_seqs.get(qseqid).seq_ln();
				ignFRQry.add(qseqid.substring(0,11));
			}
		}
		myLogger.info("#Query sequences with SAM segment "+nQrySeq+"/"+ignFRQry.size()+" summed to "+totQryLen+"bp");
		myLogger.info("#Total SAM segments "+totSAMSeg);

		/***
		
		final DirectedWeightedOverlapPseudograph<String> subgraph = new DirectedWeightedOverlapPseudograph<>(OverlapEdge.class);
		final Set<OverlapEdge> edgesToRetain = new HashSet<>();
		int slen, tlen;
		
		int ne = 0;
		for(final OverlapEdge edge : gfa.edgeSet()) {
			
			if(++ne%1000000==0) myLogger.info("#edge processed "+ne);
			
			Set<SAMSegment> source_segs, target_segs;
			Range<Integer> source_uniq, target_uniq, source_cov, target_cov, source_olap, target_olap; 

			String source = gfa.getEdgeSource(edge);
			String target = gfa.getEdgeTarget(edge);

			if(!initPlace.containsKey(source)||!initPlace.containsKey(target))
				continue;

			slen = qry_seqs.get(source).seq_ln();
			tlen = qry_seqs.get(target).seq_ln();
			
			source_uniq = Range.closed(0,  Math.max(0, slen-(int) edge.olapF()));
			target_uniq = Range.closed(Math.min(tlen, (int) edge.olapR()), tlen);

			if(source_uniq.upperEndpoint()-source_uniq.lowerEndpoint()<min_uniq_olap || 
					target_uniq.upperEndpoint()-target_uniq.lowerEndpoint()<min_uniq_olap) 
				// unique region is too small
				continue;

			loop_on_subseqs:
				for(int j=0; j<sub_seqs.size(); j++) {

					final List<SAMSegment> source_sel = new ArrayList<>();
					source_segs = initPlace.get(source).get(j);
					for(final SAMSegment source_seg : source_segs) {
						source_cov = Range.closed(source_seg.qstart(), source_seg.qend());
						if(!source_cov.isConnected(source_uniq)) continue;
						source_olap = source_cov.intersection(source_uniq);
						if(source_olap.upperEndpoint()-source_olap.lowerEndpoint()>=min_uniq_olap)
							source_sel.add(source_seg);
					}

					if(source_sel.isEmpty()) continue;

					final List<SAMSegment> target_sel = new ArrayList<>();
					target_segs = initPlace.get(target).get(j);
					for(final SAMSegment target_seg : target_segs) {
						target_cov = Range.closed(target_seg.qstart(), target_seg.qend());
						if(!target_cov.isConnected(target_uniq)) continue;
						target_olap = target_cov.intersection(target_uniq);
						if(target_olap.upperEndpoint()-target_olap.lowerEndpoint()>=min_uniq_olap)
							target_sel.add(target_seg);
					}

					if(target_sel.isEmpty()) continue;

					for(SAMSegment source_seg : source_sel) {
						for(SAMSegment target_seg : target_sel) {
							if(AlignmentSegment.sdistance(source_seg, target_seg)<hc_gap) {
								edgesToRetain.add(edge);
								break loop_on_subseqs;
							}
						}
					}
				}
		}
		
		myLogger.info("####Edge to retain "+edgesToRetain.size());
		
		for(final String v : gfa.vertexSet()) subgraph.addVertex(v);
		String source, target;
		for(final OverlapEdge edge : edgesToRetain) {
			source = gfa.getEdgeSource(edge);
			target = gfa.getEdgeTarget(edge);
			final OverlapEdge e = subgraph.addEdge(source, target);
			e.setOlap(edge.olap());
			e.setOlapF(edge.olapF());
			e.setOlapR(edge.olapR());
			e.setOlapInfo(edge.olapInfo());
			e.setCigar(edge.cigar());
			e.setRealigned(edge.isRealigned());
		}
		myLogger.info("####subgraph #V "+subgraph.vertexSet().size()+", #E "+subgraph.edgeSet().size());
		
		final List<Set<String>> conns = 
				new ArrayList<>(new ConnectivityInspector<String, OverlapEdge>(subgraph).connectedSets());
		int[] nV = new int[conns.size()];
		for(int i=0; i<conns.size(); i++) nV[i] = conns.get(i).size();
		Arrays.sort(nV);
		myLogger.info("####subgraph "+nV.length);
		myLogger.info("####max #V "+nV[nV.length-1]);
	
		final Set<String> conv = new HashSet<String>();
		for(final String v : subgraph.vertexSet()) {
			if(!subgraph.edgesOf(v).isEmpty())
				conv.add(v.substring(0, 11));
		}
		long sL = 0L;
		for(final String v : conv) sL += qry_seqs.get(v).seq_ln();
		myLogger.info("####query seqs "+conv.size());
		myLogger.info("####total length "+sL);
		***/
		
		this.pileup();
	}
	
	private final int niche = 10;
	
	private void pileup() {
		// TODO Auto-generated method stub
		elapsed[1]  = System.nanoTime();
		
		final String subSeq = "Chr07";
		final int subSeqIndex = dict.getSequenceIndex(subSeq);
		
		final List<Traceable> segBySubAll = new ArrayList<>();
		final List<SAMSegment> segs  = new ArrayList<SAMSegment>();
		int search_radius;
		for(final String qry_seq : gfa.vertexSet()) {
            segs.clear();
            
            if(!initPlace.containsKey(qry_seq)) continue;
			for(SAMSegment seg : initPlace.get(qry_seq).get(subSeqIndex))
                segs.add(seg);
			if(segs.isEmpty()) continue;

			Collections.sort(segs, new AlignmentSegment.SubjectCoordinationComparator());
			int n = segs.size();
			int start = 0, end = 0;
			final int qry_ln = qry_seqs.get(qry_seq).seq_ln();
			search_radius = qry_ln*niche;
			Range<Integer> qryRange, subRange;
			while(start<n) {
				end = start+1;
				while(end<n&&AlignmentSegment.sdistance(segs.get(end-1), segs.get(end))<=search_radius)
					++end;
				
				RangeSet<Integer> subCov = TreeRangeSet.create();
				RangeSet<Integer> qryCov = TreeRangeSet.create();
				for(int j=start; j<end; j++) {
					SAMSegment seg = segs.get(j);
					subRange = Range.closed(seg.sstart(), seg.send()).canonical(DiscreteDomain.integers());
					qryRange = Range.closed(seg.qstart(), seg.qend()).canonical(DiscreteDomain.integers());
					if(qryRange.upperEndpoint()-qryRange.lowerEndpoint()>0.7*qry_ln && 
							extension(qryCov, qryRange)<0.3*qry_ln) {
						end = j;
						break;
					}
					subCov.add(subRange);
					qryCov.add(qryRange);
				}

				Range<Integer> subSpan = subCov.span();
				int sstart = subSpan.lowerEndpoint();
				int send   = subSpan.upperEndpoint();
				int ls = 0;
				for(final Range<Integer> r : subCov.asRanges())
					ls += r.upperEndpoint()-r.lowerEndpoint();

				Range<Integer> qrySpan = qryCov.span();
				int qstart = qrySpan.lowerEndpoint();
				int qend   = qrySpan.upperEndpoint();
				int lq = 0;
				for(final Range<Integer> r : qryCov.asRanges())
					lq += r.upperEndpoint()-r.lowerEndpoint();
				int slen = ls;
				int qlen = lq;

				segBySubAll.add(new Traceable(qry_seq, subSeq, 
						qstart, qend, sstart, send, subCov, qryCov, slen, qlen));

				start = end;
			}
		}

		Collections.sort(segBySubAll, new AlignmentSegment.SubjectCoordinationComparator());
		int nSeg = segBySubAll.size();
		myLogger.info("####Subject Sequence "+subSeq+", #SamSegment "+nSeg);
		
		elapsed[2]  = System.nanoTime();
		
		// now extract subgraphs from the alignments
		
		final Set<Integer> processed = new HashSet<>();
		final Set<Integer> contained = new HashSet<>();
		int source_segix, target_segix;
		Traceable source_seg, target_seg;
		String source_segid, target_segid;
		final Deque<Integer> deque = new ArrayDeque<Integer>();
		final List<DirectedWeightedPseudograph<Integer, DefaultWeightedEdge>> subgraphs = new ArrayList<>();
		
		for(int i=0; i<nSeg; i++) {

			if(i%100000==0) myLogger.info("#SAM segments processed "+i+", #subgraph "+subgraphs.size());
			
			if(processed.contains(i)) continue;
			
			deque.clear();
			deque.add(i);
			contained.clear();
			contained.add(i);
			
			final DirectedWeightedPseudograph<Integer, DefaultWeightedEdge> subgraph = 
					new DirectedWeightedPseudograph<>(DefaultWeightedEdge.class);

			while(!deque.isEmpty()) {
				
				source_segix  = deque.pop();
				source_seg    = segBySubAll.get(source_segix);
				source_segid  = source_seg.qseqid();
				search_radius = qry_seqs.get(source_segid).seq_ln()*niche;
				
				if(!subgraph.containsVertex(source_segix))
					subgraph.addVertex(source_segix);
				
				for(int j=source_segix-1; j>0; j--) {
					target_segix = j;
					target_seg   = segBySubAll.get(target_segix);
					target_segid = target_seg.qseqid();

					if(source_seg.sstart()-target_seg.sstart()>search_radius) break;
					
					if(gfa.containsEdge(source_segid, target_segid)) {
						if(!subgraph.containsVertex(target_segix))
							subgraph.addVertex(target_segix);
						if(!subgraph.containsEdge(source_segix, target_segix))
							subgraph.addEdge(source_segix, target_segix);
						if(!contained.contains(target_segix)) {
							deque.add(target_segix);
							contained.add(target_segix);
						}
					}
					
					if(gfa.containsEdge(target_segid, source_segid)) {
						if(!subgraph.containsVertex(target_segix))
							subgraph.addVertex(target_segix);
						if(!subgraph.containsEdge(target_segix, source_segix))
							subgraph.addEdge(target_segix, source_segix);
						if(!contained.contains(target_segix)) {
							deque.add(target_segix);
							contained.add(target_segix);
						}
					}
				}

				for(int j=source_segix+1; j<nSeg; j++) {
					target_segix = j;
					target_seg   = segBySubAll.get(target_segix);
					target_segid = target_seg.qseqid();

					if(target_seg.sstart()-source_seg.sstart()>search_radius) break;

					if(gfa.containsEdge(source_segid, target_segid)) {
						if(!subgraph.containsVertex(target_segix))
							subgraph.addVertex(target_segix);
						if(!subgraph.containsEdge(source_segix, target_segix))
							subgraph.addEdge(source_segix, target_segix);
						if(!contained.contains(target_segix)) {
							deque.add(target_segix);
							contained.add(target_segix);
						}
					}
					
					if(gfa.containsEdge(target_segid, source_segid)) {
						if(!subgraph.containsVertex(target_segix))
							subgraph.addVertex(target_segix);
						if(!subgraph.containsEdge(target_segix, source_segix))
							subgraph.addEdge(target_segix, source_segix);
						if(!contained.contains(target_segix)) {
							deque.add(target_segix);
							contained.add(target_segix);
						}
					}
				}

				processed.add(source_segix);
			}
			subgraphs.add(subgraph);
		}
		
		if(debug) {
			int numV = 0, numE = 0, maxV = -1, minV = Integer.MAX_VALUE;
			int[] sizes = new int[10];
			for(final DirectedWeightedPseudograph<Integer, DefaultWeightedEdge> subgraph : subgraphs) {
				int nV = subgraph.vertexSet().size();
				if(nV<2)        ++sizes[0];
				else if(nV<10)  ++sizes[1];
				else if(nV<50)  ++sizes[2];
				else if(nV<100) ++sizes[3];
				else            ++sizes[4];
				if(maxV<nV) maxV = nV;
				if(minV>nV) minV = nV;
				numV += nV;
				numE += subgraph.edgeSet().size();
			}
			int avgV = numV/subgraphs.size();
			myLogger.info("#### "+subSeq+" ####\n"+
					"#segments  : "+segBySubAll.size()+"\n"+
					"#subgraphs : "+subgraphs.size()+"\n"+
					"#E         : "+numE+"\n"+
					"#V         : "+numV+"\n"+
					"#max V     : "+maxV+"\n"+
					"#min V     : "+minV+"\n"+
					"#avg V     : "+avgV+"\n"+
					"# size     : \n"+
					"     1     : "+sizes[0]+"\n"+
					"   <10     : "+sizes[1]+"\n"+
					"   <50     : "+sizes[2]+"\n"+
					"  <100     : "+sizes[3]+"\n"+
					"  >100     : "+sizes[4]+"\n");
		}
		
		elapsed[3]  = System.nanoTime();
		myLogger.info("#elpased time 0 "+(elapsed[1]-elapsed[0])/1e9+"s");
		myLogger.info("#elpased time 1 "+(elapsed[2]-elapsed[1])/1e9+"s");
		myLogger.info("#elpased time 2 "+(elapsed[3]-elapsed[2])/1e9+"s");
	}

	private int extension(RangeSet<Integer> cov, Range<Integer> range) {
		// TODO Auto-generated method stub
		
		RangeSet<Integer> ins = TreeRangeSet.create();
		for(Range<Integer> r : cov.complement().asRanges())
			if(range.isConnected(r))
				ins.add(range.intersection(r));
		int ln = 0;
		for(Range<Integer> r : ins.asRanges())
			ln += r.upperEndpoint()-r.lowerEndpoint();
		return ln;
	}

	final SamReaderFactory factory =
			SamReaderFactory.makeDefault()
			.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
					SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
			.validationStringency(ValidationStringency.SILENT);

	private final static Object lock = new Object();
	
	private class BAMBarcodeIterator {
		private final SamReader[] samReader = new SamReader[2];
		private final SAMRecordIterator[] iter = new SAMRecordIterator[2];
		private final SAMRecord[] samRecord = new SAMRecord[2];

		public BAMBarcodeIterator(final String f1, final String f2) {
			this.samReader[0] = factory.open(new File(f1));
			this.samReader[1] = factory.open(new File(f2));
			this.iter[0] = this.samReader[0].iterator();
			this.iter[1] = this.samReader[1].iterator();
			this.samRecord[0] = iter[0].hasNext() ? iter[0].next() : null;
			this.samRecord[1] = iter[1].hasNext() ? iter[1].next() : null;
		}

		public boolean hasNext() {
			return samRecord[0]!=null && samRecord[1]!=null;
		}

		public List<Set<SAMRecord>> next() {

			if(!this.hasNext()) throw new RuntimeException("!!!");

			List<Set<SAMRecord>> records = new ArrayList<Set<SAMRecord>>();
			String rn = samRecord[0].getReadName();
			if(!rn.equals(samRecord[1].getReadName()))
				throw new RuntimeException("!!!");

			for(int i=0; i<2; i++) {
				final Set<SAMRecord> r = new HashSet<SAMRecord>();
				while( samRecord[i]!=null && rn.equals(samRecord[i].getReadName()) ) {
					r.add(samRecord[i]);
					samRecord[i] = iter[i].next();
				}
				records.add(r);
			}
			
			return records;
		}

		public void close() {
			try {
				for(int i=0; i<2; i++) {
					this.iter[i].close();
					this.samReader[i].close();
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	private class SAMPair {
		private final SAMRecord  first;
		private final SAMRecord second;
		private final String  firstr;
		private final String secondr;
		private final int ins;
		
		public SAMPair(final SAMRecord first,
				final SAMRecord second,
				final String  firstr,
				final String secondr,
				final int ins) {
			this.first   =   first;
			this.second  =  second;
			this.firstr  =  firstr;
			this.secondr = secondr;
			this.ins = ins;
		}
	}
	
	private BidiMap<String, Integer> seq_index;
	private Map<String, String> symm_seqsn;
	private final static Map<Long, Double> linkCount = new HashMap<>(); // link count
	private static long exceed_ins = 0;
	private final static double olap_ext = 0.3;
	private static long readCount = 0;
	
	private void run_link() {
		// TODO Auto-generated method stub
		final int niof = this.bamList.length;
		myLogger.info("Reading alignments from "+niof+" BAM file"+
				(this.bamList.length>1?"s":"")+":");
		for(int i=0; i<niof; i++)
			myLogger.info(this.bamList[i][0]+" paired with "+this.bamList[i][1]+" "+
					(this.matePair[i]?"mate-pair":"paired-end")+" library.");
		myLogger.info("****");
		
		final GFA gfa = new GFA(this.subject_file, this.asm_graph);
		myLogger.info("Loading assembly graph done.");
		
		sub_seqs = Sequence.parseFastaFileWithRevCmpAsMap(subject_file);
		seq_index = new DualHashBidiMap<String, Integer>();
		symm_seqsn = new HashMap<String, String>();
		
		int index = 0;
		for(String seq : sub_seqs.keySet()) seq_index.put(seq, ++index);
		for(String seq : sub_seqs.keySet()) {
			if(!seq.endsWith("'")) {
				symm_seqsn.put(seq, seq+"'");
				symm_seqsn.put(seq+"'", seq);
			}
		}
		
		for(int i=0; i<niof; i++) {
			myLogger.info("####Processing "+this.bamList[i][0]+" paired with "+this.bamList[i][1]);
			this.initial_thread_pool();
			final BAMBarcodeIterator iter = new BAMBarcodeIterator(this.bamList[i][0], this.bamList[i][1]);
			while(iter.hasNext()) {
				this.executor.submit(new Runnable() {
					private Set<SAMRecord> sr1;
					private Set<SAMRecord> sr2;
					private boolean matePair;
					private int maxInst;
					
					@Override
					public void run() {
						// TODO Auto-generated method stub
						if(sr1.isEmpty()||sr2.isEmpty()) return;
						
						final List<SAMPair> sampairs = new ArrayList<SAMPair>();
						SAMRecord tmp;
						int inst, s1, s2, e1, e2, olap;
						final boolean[] rev = new boolean[2];
						final int[] reflen = new int[2];
						final String[] ref = new String[2];
						final String[] refstr = new String[2];
						
                        synchronized(lock) {
                            ++readCount;
                            if(readCount%1000000==0)
                                myLogger.info("# "+readCount+" reads processed.");
                        }

						try {
							if(ddebug) myLogger.info(sr1.iterator().next().getReadName()+": "+sr1.size()+"x"+sr2.size());
							for(SAMRecord r1 : sr1) {
								// check return
								if(r1.getReadUnmappedFlag() ||
										r1.getMappingQuality()<minQual) 
									continue;
								s1 = r1.getAlignmentStart();
								e1 = r1.getAlignmentEnd();
								rev[0] = r1.getReadNegativeStrandFlag();
								ref[0] = r1.getReferenceName();
								
								for(SAMRecord r2 : sr2) {
									// check return
									if(r2.getReadUnmappedFlag() ||
											r2.getMappingQuality()<minQual) 
										continue;
									s2 = r2.getAlignmentStart();
									e2 = r2.getAlignmentEnd();
									rev[1] = r2.getReadNegativeStrandFlag();
									ref[1] = r2.getReferenceName();
									
									if(!r1.getReadName().
											equals(r2.getReadName()))
										throw new RuntimeException("!!!");
									
									inst = Integer.MAX_VALUE;

									if(r1.getReferenceIndex().intValue()==r2.getReferenceIndex().intValue()) {
										// mapped to the same contig
										if(rev[0]!=rev[1]) {
											// two reads need to be mapped in reverse direction
											if(rev[1]) {
												tmp = r1;
												r1 = r2;
												r2 = tmp;
											}
											if(s1<s2&&matePair) {
												inst = e2-s1;
											} else if(s1>s2&&!matePair) {
												inst = e1-s2;
											}
											if(inst<=maxInst) return;
										}
									} else {
										// mapped to different contigs
										reflen[0]  = sub_seqs.get(r1.getReferenceName()).seq_ln();
										reflen[1]  = sub_seqs.get(r2.getReferenceName()).seq_ln();

										if(matePair) {
											// this is a mate-pair
											if(rev[0]&&rev[1]) {
												//       0                  1
												// ---------------     -------------
												//   <===                 <===

												// reverse 0
												// ---------------     -------------
												//          ===>          <===
												// reverse 1
												// ---------------     -------------
												//   <===                    ===>

												// reverse 0 and 1 are symmetric
												refstr[0] = r1.getReferenceName();
												refstr[1] = r2.getReferenceName()+"'";
												olap = (int) gfa.getEdge(refstr[0], refstr[1]).olap();
												if( reflen[0]-s1+1-olap<olap_ext*(r1.getAlignmentEnd()-r1.getAlignmentStart()) ||
														reflen[1]-s2+1-olap<olap_ext*(r2.getAlignmentEnd()-r2.getAlignmentStart()) )
													continue;
												inst = reflen[0]-s1+1+reflen[1]-s2+1-olap;
												
											} else if(rev[0]&&!rev[1]) {
												//       0                  1
												// ---------------     -------------
												//   <===                    ===>

												// reverse 0 & reverse 1
												// ---------------     -------------
												//          ===>          <===
												refstr[0] = r1.getReferenceName();
												refstr[1] = r2.getReferenceName();
												olap = (int) gfa.getEdge(refstr[0], refstr[1]).olap();
												if( reflen[0]-s1+1-olap<olap_ext*(r1.getAlignmentEnd()-r1.getAlignmentStart()) ||
														e2-olap<olap_ext*(r2.getAlignmentEnd()-r2.getAlignmentStart()) )
													continue;
												inst = reflen[0]-s1+1+e2-olap;
												
											} else if(!rev[0]&&rev[1]) {
												//       0                  1
												// ---------------     -------------
												//          ===>          <===

												// reverse 0 & reverse 1
												// ---------------     -------------
												//   <===                    ===>
												refstr[0] = r2.getReferenceName();
												refstr[1] = r1.getReferenceName();
												olap = (int) gfa.getEdge(refstr[0], refstr[1]).olap();
												if( e1-olap<olap_ext*(r1.getAlignmentEnd()-r1.getAlignmentStart()) ||
														reflen[1]-s2+1-olap<olap_ext*(r2.getAlignmentEnd()-r2.getAlignmentStart()) )
													continue;
												inst = e1+reflen[1]-s2+1-olap;

											} else {
												//       0                  1
												// ---------------     -------------
												//         ===>               ===>

												// reverse 0
												// ---------------     -------------
												//    <===                    ===>
												// reverse 1
												// ---------------     -------------
												//         ===>          <===
												refstr[0] = r1.getReferenceName()+"'";
												refstr[1] = r2.getReferenceName();
												olap = (int) gfa.getEdge(refstr[0], refstr[1]).olap();
												if( e1-olap<olap_ext*(r1.getAlignmentEnd()-r1.getAlignmentStart()) ||
														e2-olap<olap_ext*(r2.getAlignmentEnd()-r2.getAlignmentStart()) )
													continue;
												inst = e1+e2-olap;

											}
										} else {
											// this is a paired-end
											if(rev[0]&&rev[1]) {
												//       0                  1
												// ---------------     -------------
												//   <===                 <===

												// reverse 0
												// ---------------     -------------
												//          ===>          <===
												// reverse 1
												// ---------------     -------------
												//   <===                    ===>
												
												// reverse 0 and 1 are symmetric
												refstr[0] = r1.getReferenceName()+"'";
												refstr[1] = r2.getReferenceName();
												olap = (int) gfa.getEdge(refstr[0], refstr[1]).olap();
												if( e1-olap<olap_ext*(r1.getAlignmentEnd()-r1.getAlignmentStart()) ||
														e2-olap<olap_ext*(r2.getAlignmentEnd()-r2.getAlignmentStart()) )
													continue;
												inst = e1+e2-olap;

											} else if(rev[0]&&!rev[1]) {
												//       0                  1
												// ---------------     -------------
												//   <===                    ===>

												// reverse 0 & reverse 1
												// ---------------     -------------
												//          ===>          <===

												refstr[0] = r2.getReferenceName();
												refstr[1] = r1.getReferenceName();
												olap = (int) gfa.getEdge(refstr[0], refstr[1]).olap();
												if( e1-olap<olap_ext*(r1.getAlignmentEnd()-r1.getAlignmentStart()) ||
														reflen[1]-s2+1-olap<olap_ext*(r2.getAlignmentEnd()-r2.getAlignmentStart()) )
													continue;
												inst = e1+reflen[1]-s2+1-olap;
												
											} else if(!rev[0]&&rev[1]) {
												//       0                  1
												// ---------------     -------------
												//          ===>          <===

												// reverse 0 & reverse 1
												// ---------------     -------------
												//   <===                    ===>

												refstr[0] = r1.getReferenceName();
												refstr[1] = r2.getReferenceName();
												olap = (int) gfa.getEdge(refstr[0], refstr[1]).olap();
												if( reflen[0]-s1+1-olap<olap_ext*(r1.getAlignmentEnd()-r1.getAlignmentStart()) ||
														e2-olap<olap_ext*(r2.getAlignmentEnd()-r2.getAlignmentStart()) )
													continue;
												inst = reflen[0]-s1+1+e2-olap;

											} else {
												//       0                  1
												// ---------------     -------------
												//         ===>               ===>

												// reverse 0
												// ---------------     -------------
												//    <===                    ===>
												// reverse 1
												// ---------------     -------------
												//         ===>          <===

												refstr[0] = r1.getReferenceName();
												refstr[1] = r2.getReferenceName()+"'";
												olap = (int) gfa.getEdge(refstr[0], refstr[1]).olap();
												if( reflen[0]-s1+1-olap<olap_ext*(r1.getAlignmentEnd()-r1.getAlignmentStart()) ||
														reflen[1]-s2+1-olap<olap_ext*(r2.getAlignmentEnd()-r2.getAlignmentStart()) )
													continue;
												inst = reflen[0]-s1+1+reflen[1]-s2+1-olap;
											}
										}
										
										if(inst>maxInst) continue;
										sampairs.add(new SAMPair(r1, r2, refstr[0], refstr[1], inst));
									}
								}
							}
							
							if(sampairs.isEmpty()) return;
							double w = sampairs.size();
							
							long refind;
							for(SAMPair sampair : sampairs) {
								
								refind  = seq_index.get(sampair.firstr);
								refind <<= 32;
								refind += seq_index.get(sampair.secondr);

								synchronized(lock) {
									if(linkCount.containsKey(refind)) 
										linkCount.put(refind, linkCount.get(refind)+w);
									else
										linkCount.put(refind, w);
								}

								refind  = seq_index.get(symm_seqsn.get(sampair.secondr));
								refind <<= 32;
								refind += seq_index.get(symm_seqsn.get(sampair.firstr));

								synchronized(lock) {
									if(linkCount.containsKey(refind)) 
										linkCount.put(refind, linkCount.get(refind)+w);
									else
										linkCount.put(refind, w);
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

					public Runnable init(final List<Set<SAMRecord>> records, final boolean matePair, final int maxInst) {
						// TODO Auto-generated method stub
						this.sr1 = records.get(0);
						this.sr2 = records.get(1);
						this.matePair = matePair;
						this.maxInst = maxInst;
						return this;
					}
					
				}.init(iter.next(), this.matePair[i], this.maxInst[i]));
			}
			iter.close();
			this.waitFor();
		}
	}

	final int max_dist = 10;
	final int min_iden = 95;
	final int min_olen =  0;
	private void blastReader() {
		// TODO Auto-generated method stub
		qry_seqs = Sequence.parseFastaFileWithRevCmpAsMap(query_file);
		try {
			BufferedReader br = Utils.getBufferedReader(align_file);
			BufferedWriter bw1 = Utils.getBufferedWriter(out_prefix+".txt1");
			BufferedWriter bw2 = Utils.getBufferedWriter(out_prefix+".txt2");
			BufferedWriter bw3 = Utils.getBufferedWriter(out_prefix+".gfa" );
			bw3.write("H\tVN:Z:1.0\n");
			for(final Map.Entry<String, Sequence> entry : qry_seqs.entrySet()) {
				if(entry.getKey().endsWith("'")) continue;
				bw3.write("S\t"+entry.getKey()+"\t*\tLN:i:"+entry.getValue().seq_ln()+"\n");
			}
			
			final Set<String> edges = new HashSet<String>();
			
			String line;
			String[] s;
			int sstart, qstart, send, qend;
			int bqstart = max_dist, bqend, bsstart = max_dist, bsend, D, I, qlen, slen;
			boolean q_from_start, q_to_end, s_from_start, s_to_end;
			while( (line=br.readLine())!=null ) {
				s = line.split("\\s+");
				if(s[0].equals(s[1]) ||
						Double.parseDouble(s[2])<min_iden ||
						Integer.parseInt(s[3])<min_olen) 
					continue;
				
				qlen = qry_seqs.get(s[0]).seq_ln();
				slen = qry_seqs.get(s[1]).seq_ln();
				
				bqend = qlen-max_dist;
				bsend = slen-max_dist;
				
				qstart = Integer.parseInt(s[6]);
				qend   = Integer.parseInt(s[7]);
				sstart = Integer.parseInt(s[8]);
				send   = Integer.parseInt(s[9]);
				
				q_from_start = qstart<bqstart;
				q_to_end     = qend>bqend;
				
				if(send>sstart) {
					s_from_start = sstart<bsstart;
					s_to_end     = send>bsend;
					
					if( q_from_start&&s_to_end || q_to_end&&s_from_start ) 
						bw1.write(line+"\n");
					else if(q_from_start&&q_to_end || s_from_start&&s_to_end)
						bw2.write(line+"\n");
					if(q_from_start&&s_to_end) {
						D = qstart-1;
						I = slen-send;
						if(!edges.contains(s[1]+"+"+s[0]+"+")) {
							bw3.write("L\t"+s[1]+"\t+\t"+s[0]+"\t+\t"+(D>0?(D+"D"):"")+(send-sstart+1)+"M"+(I>0?(I+"I"):"")+"\n");
							edges.add(s[1]+"+"+s[0]+"+");
						}
					}
					if(q_to_end&&s_from_start) {
						D = sstart-1;
						I = qlen-qend;
						if(!edges.contains(s[0]+"+"+s[1]+"+")) {
							bw3.write("L\t"+s[0]+"\t+\t"+s[1]+"\t+\t"+(D>0?(D+"D"):"")+(qend-qstart+1)+"M"+(I>0?(I+"I"):"")+"\n");
							edges.add(s[0]+"+"+s[1]+"+");
						}
					}
				} else {
					s_from_start = send<bsstart;
					s_to_end     = sstart>bsend;
					
					if( q_from_start&&s_from_start || q_to_end&&s_to_end ) 
						bw1.write(line+"\n");
					else if(q_from_start&&q_to_end || s_from_start&&s_to_end)
						bw2.write(line+"\n");
					if(q_from_start&&s_from_start) {
						D = qstart-1;
						I = send-1;
						if(!edges.contains(s[1]+"-"+s[0]+"+")) {
							bw3.write("L\t"+s[1]+"\t-\t"+s[0]+"\t+\t"+(D>0?(D+"D"):"")+(sstart-send+1)+"M"+(I>0?(I+"I"):"")+"\n");
							edges.add(s[1]+"-"+s[0]+"+");
						}
					}
					if(q_to_end&&s_to_end) {
						D = slen-sstart;
						I = qlen-qend;
						if(!edges.contains(s[0]+"+"+s[1]+"-")) {
							bw3.write("L\t"+s[0]+"\t+\t"+s[1]+"\t-\t"+(D>0?(D+"D"):"")+(qend-qstart+1)+"M"+(I>0?(I+"I"):"")+"\n");
							edges.add(s[0]+"+"+s[1]+"-");
						}
					}
				}
			}
			bw1.close();
			bw2.close();
			bw3.close();
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}



	private void run3() {
		// TODO Auto-generated method stub

		// read assembly graph file
		final GFA gfa = new GFA(query_file, asm_graph);
		myLogger.info("  GFA vertices: "+gfa.vertexSet().size());
		myLogger.info("  GFA edges   : "+gfa.edgeSet().size()  );

		qry_seqs = gfa.getSequenceMap();
		sub_seqs = Sequence.parseFastaFileAsMap(subject_file);

		// read alignment file and place the query sequences
		final Map<String, List<Set<SAMSegment>>> initPlace = new HashMap<String, List<Set<SAMSegment>>>();
		SAMSequenceDictionary dict = null;
		try {
			final SamReaderFactory factory =
					SamReaderFactory.makeDefault()
					.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
							SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
					.validationStringency(ValidationStringency.SILENT);
			final SamReader in1 = factory.open(new File(align_file));
			final SAMRecordIterator iter1 = in1.iterator();
			dict  =  in1.getFileHeader().getSequenceDictionary();
			final int nSubSeq = dict.getSequences().size();
			if(nSubSeq!=sub_seqs.size()) throw new RuntimeException("!!!");

			myLogger.info("--"+nSubSeq+" Subject sequences");
			for(int i=0; i<nSubSeq; i++) 
				myLogger.info("  "+dict.getSequences().get(i).getSequenceName());
			myLogger.info("----------------");
			
			String qry;
			final List<SAMSegment> buff = new ArrayList<SAMSegment>();
			SAMRecord rc = iter1.next();
			
			int nq = 0;
			
			while(rc!=null) {
				if(++nq%100000==0) myLogger.info("#query sequences processed "+nq);
				
				qry = rc.getReadName();	

				buff.clear();
				if(!rc.getReadUnmappedFlag()&&rc.getMappingQuality()>0)
					buff.add(SAMSegment.samRecord(rc, true));

				while( (rc=iter1.next())!=null
						&&
						rc.getReadName().equals(qry) ) {
					if(!rc.getReadUnmappedFlag()&&rc.getMappingQuality()>0)
						buff.add(SAMSegment.samRecord(rc, true));
				}

				if(buff.isEmpty()) continue;

				Set<SAMSegment> init_f = new HashSet<SAMSegment>();
				Set<SAMSegment> init_r = new HashSet<SAMSegment>();
				for(SAMSegment record : buff) {
					if(record.qlength()<min_alen) 
						continue;
					if(record.qseqid().equals(qry)) 
						init_f.add(record);
					else init_r.add(record);
				}
				if(!initPlace.containsKey(qry)) {
					final List<Set<SAMSegment>> list_f = new ArrayList<Set<SAMSegment>>(nSubSeq);
					for(int i=0; i<nSubSeq; i++) list_f.add(new HashSet<SAMSegment>());
					initPlace.put(qry,     list_f);
					final List<Set<SAMSegment>> list_r = new ArrayList<Set<SAMSegment>>(nSubSeq);
					for(int i=0; i<nSubSeq; i++) list_r.add(new HashSet<SAMSegment>());
					initPlace.put(qry+"'", list_r);
				}

				if(!init_f.isEmpty()) {
					final List<Set<SAMSegment>> list_f = initPlace.get(qry);
					for(SAMSegment seg : init_f) 
						list_f.get(dict.getSequenceIndex(seg.sseqid())).add(seg);
				}
				if(!init_r.isEmpty()) {
					final List<Set<SAMSegment>> list_r = initPlace.get(qry+"'");
					for(SAMSegment seg : init_r) 
						list_r.get(dict.getSequenceIndex(seg.sseqid())).add(seg);
				}
			}
			iter1.close();
			in1.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		int nQrySeq = 0, totSAMSeg = 0;
		long totQryLen = 0;
		final Set<String> ignFRQry = new HashSet<String>();
		for(final String qseqid : initPlace.keySet()) {
			final List<Set<SAMSegment>> list = initPlace.get(qseqid);
			int n = 0;
			for(final Set<SAMSegment> segs : list) n += segs.size();
			if(n>0) {
				++nQrySeq;
				totSAMSeg+=n;
				totQryLen+=qry_seqs.get(qseqid).seq_ln();
				ignFRQry.add(qseqid.substring(0,11));
			}
		}
		myLogger.info("#Query sequences with SAM segment "+nQrySeq+"/"+ignFRQry.size()+" summed to "+totQryLen+"bp");
		myLogger.info("#Total SAM segments "+totSAMSeg);

		final DirectedWeightedOverlapPseudograph<String> graph = gfa.gfa();
		final Set<OverlapEdge> edgesToRetain = new HashSet<>();
		int slen, tlen;
		
		int ne = 0;
		for(final OverlapEdge edge : gfa.edgeSet()) {
			
			if(++ne%1000000==0) myLogger.info("#edge processed "+ne);
			
			Set<SAMSegment> source_segs, target_segs;
			Range<Integer> source_uniq, target_uniq, source_cov, target_cov, source_olap, target_olap; 

			String source = gfa.getEdgeSource(edge);
			String target = gfa.getEdgeTarget(edge);

			if(!initPlace.containsKey(source)||!initPlace.containsKey(target))
				continue;

			slen = qry_seqs.get(source).seq_ln();
			tlen = qry_seqs.get(target).seq_ln();
			
			source_uniq = Range.closed(0,  Math.max(0, slen-(int) edge.olapF()));
			target_uniq = Range.closed(Math.min(tlen, (int) edge.olapR()), tlen);

			if(source_uniq.upperEndpoint()-source_uniq.lowerEndpoint()<min_uniq_olap || 
					target_uniq.upperEndpoint()-target_uniq.lowerEndpoint()<min_uniq_olap) 
				// unique region is too small
				continue;

			loop_on_subseqs:
				for(int j=0; j<sub_seqs.size(); j++) {

					final List<SAMSegment> source_sel = new ArrayList<>();
					source_segs = initPlace.get(source).get(j);
					for(final SAMSegment source_seg : source_segs) {
						source_cov = Range.closed(source_seg.qstart(), source_seg.qend());
						if(!source_cov.isConnected(source_uniq)) continue;
						source_olap = source_cov.intersection(source_uniq);
						if(source_olap.upperEndpoint()-source_olap.lowerEndpoint()>=min_uniq_olap)
							source_sel.add(source_seg);
					}

					if(source_sel.isEmpty()) continue;

					final List<SAMSegment> target_sel = new ArrayList<>();
					target_segs = initPlace.get(target).get(j);
					for(final SAMSegment target_seg : target_segs) {
						target_cov = Range.closed(target_seg.qstart(), target_seg.qend());
						if(!target_cov.isConnected(target_uniq)) continue;
						target_olap = target_cov.intersection(target_uniq);
						if(target_olap.upperEndpoint()-target_olap.lowerEndpoint()>=min_uniq_olap)
							continue;
						target_sel.add(target_seg);
					}

					if(target_sel.isEmpty()) continue;

					for(SAMSegment source_seg : source_sel) {
						for(SAMSegment target_seg : source_sel) {
							if(AlignmentSegment.sdistance(source_seg, target_seg)<hc_gap) {
								edgesToRetain.add(edge);
								break loop_on_subseqs;
							}
						}
					}
				}

		}
		
		myLogger.info("####Edge to retain "+edgesToRetain.size());
		
		final Set<OverlapEdge> edgesToRemove = new HashSet<>(gfa.edgeSet());
		edgesToRemove.removeAll(edgesToRetain);
		gfa.removeAllEdges(edgesToRemove);
		
		myLogger.info("####Edge removed "+edgesToRemove.size());
		
		final List<Set<String>> conns = 
				new ArrayList<>(new ConnectivityInspector<String, OverlapEdge>(graph).connectedSets());
		int[] nV = new int[conns.size()];
		for(int i=0; i<conns.size(); i++) nV[i] = conns.get(i).size();
		Arrays.sort(nV);
		myLogger.info("####subgraph "+nV.length);
		myLogger.info("####max #V "+nV[nV.length-1]);
	
		final Set<String> conv = new HashSet<String>();
		for(final String v : gfa.vertexSet()) {
			if(!gfa.edgesOf(v).isEmpty())
				conv.add(v.substring(0, 11));
		}
		long sL = 0L;
		for(final String v : conv) sL += qry_seqs.get(v).seq_ln();
		myLogger.info("####query seqs "+conv.size());
		myLogger.info("####total length "+sL);
	}
	
	final int min_uniq_olap = 30;

	private void run1() {
		// TODO Auto-generated method stub

		// read assembly graph file
		final GFA gfa = new GFA(query_file, asm_graph);
		myLogger.info("  GFA vertices: "+gfa.vertexSet().size());
		myLogger.info("  GFA edges   : "+gfa.edgeSet().size()  );

		qry_seqs = gfa.getSequenceMap();
		sub_seqs = Sequence.parseFastaFileAsMap(subject_file);

		/***
		int min_uniq = Integer.MAX_VALUE;
		for(final OverlapEdge edge : gfa.edgeSet()) {
			String source = gfa.getEdgeSource(edge);
			String target = gfa.getEdgeTarget(edge);
			int source_uniq = qry_seqs.get(source).seq_ln()-(int)edge.olapF();
			int target_uniq = qry_seqs.get(target).seq_ln()-(int)edge.olapR();
			if(source_uniq<0||target_uniq<0) {
				myLogger.info("olap<0");
			}
			if(source_uniq<min_uniq) min_uniq = source_uniq;
			if(target_uniq<min_uniq) min_uniq = target_uniq;
		}
		myLogger.info("#Minimum unique "+min_uniq+"bp");
		 ***/

		// read alignment file and place the query sequences
		final Map<String, List<Set<SAMSegment>>> initPlace = new HashMap<String, List<Set<SAMSegment>>>();
		SAMSequenceDictionary dict = null;
		try {
			final SamReaderFactory factory =
					SamReaderFactory.makeDefault()
					.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
							SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
					.validationStringency(ValidationStringency.SILENT);
			final SamReader in1 = factory.open(new File(align_file));
			final SAMRecordIterator iter1 = in1.iterator();
			dict  =  in1.getFileHeader().getSequenceDictionary();
			final int nSubSeq = dict.getSequences().size();
			if(nSubSeq!=sub_seqs.size()) throw new RuntimeException("!!!");

			myLogger.info("--"+nSubSeq+"Subject sequences");
			for(int i=0; i<nSubSeq; i++) 
				myLogger.info("  "+dict.getSequences().get(i).getSequenceName());
			myLogger.info("--");
			
			String qry;
			final List<SAMSegment> buff = new ArrayList<SAMSegment>();
			SAMRecord rc = iter1.next();

			while(rc!=null) {
				qry = rc.getReadName();	

				buff.clear();
				if(!rc.getReadUnmappedFlag()&&rc.getMappingQuality()>0)
					buff.add(SAMSegment.samRecord(rc, true));

				while( (rc=iter1.next())!=null
						&&
						rc.getReadName().equals(qry) ) {
					if(!rc.getReadUnmappedFlag()&&rc.getMappingQuality()>0)
						buff.add(SAMSegment.samRecord(rc, true));
				}

				if(buff.isEmpty()) continue;

				Set<SAMSegment> init_f = new HashSet<SAMSegment>();
				Set<SAMSegment> init_r = new HashSet<SAMSegment>();
				for(SAMSegment record : buff) {
					if(record.qlength()<min_alen) 
						continue;
					if(record.qseqid().equals(qry)) 
						init_f.add(record);
					else init_r.add(record);
				}
				if(!initPlace.containsKey(qry)) {
					final List<Set<SAMSegment>> list_f = new ArrayList<Set<SAMSegment>>(nSubSeq);
					for(int i=0; i<nSubSeq; i++) list_f.add(new HashSet<SAMSegment>());
					initPlace.put(qry,     list_f);
					final List<Set<SAMSegment>> list_r = new ArrayList<Set<SAMSegment>>(nSubSeq);
					for(int i=0; i<nSubSeq; i++) list_r.add(new HashSet<SAMSegment>());
					initPlace.put(qry+"'", list_r);
				}

				if(!init_f.isEmpty()) {
					final List<Set<SAMSegment>> list_f = initPlace.get(qry);
					for(SAMSegment seg : init_f) 
						list_f.get(dict.getSequenceIndex(seg.sseqid())).add(seg);
				}
				if(!init_r.isEmpty()) {
					final List<Set<SAMSegment>> list_r = initPlace.get(qry+"'");
					for(SAMSegment seg : init_r) 
						list_r.get(dict.getSequenceIndex(seg.sseqid())).add(seg);
				}
			}
			iter1.close();
			in1.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		int nQrySeq = 0, totSAMSeg = 0;
		long totQryLen = 0;
		final Set<String> ignFRQry = new HashSet<String>();
		for(final String qseqid : initPlace.keySet()) {
			final List<Set<SAMSegment>> list = initPlace.get(qseqid);
			int n = 0;
			for(final Set<SAMSegment> segs : list) n += segs.size();
			if(n>0) {
				++nQrySeq;
				totSAMSeg+=n;
				totQryLen+=qry_seqs.get(qseqid).seq_ln();
				ignFRQry.add(qseqid.substring(0,11));
			}
		}
		myLogger.info("#Query sequences with SAM segment "+nQrySeq+"/"+ignFRQry.size()+" summed to "+totQryLen+"bp");
		myLogger.info("#Total SAM segments "+totSAMSeg);

		final DirectedWeightedOverlapPseudograph<String> graph = gfa.gfa();
		final Set<OverlapEdge> edgesToRemove = new HashSet<>();
		final Set<OverlapEdge> edgesToRetain = new HashSet<>();
		
		for(final String qseqid : gfa.vertexSet()) {
			
			int n;
			String source, target;
			Set<SAMSegment> source_segs, target_segs;
			Range<Integer> source_uniq, target_uniq, source_cov, target_cov, source_olap, target_olap; 
			
			for(int z=0; z<2; z++) {
				final List<OverlapEdge> edges = (z==0?
						new ArrayList<>(gfa.outgoingEdgesOf(qseqid)):
							new ArrayList<>(gfa.incomingEdgesOf(qseqid)));
				
				n = edges.size();
				//if( (n = edges.size())>1 ) {
					final boolean[] linkevi = new boolean[n];

					for(int i=0; i<n; i++) {
						final OverlapEdge edge = edges.get(i);
						source = gfa.getEdgeSource(edge);
						target = gfa.getEdgeTarget(edge);
						
						if(!initPlace.containsKey(source)||!initPlace.containsKey(target))
							continue;
						
						source_uniq = Range.closed(0, qry_seqs.get(source).seq_ln()-(int) edge.olapF());
						target_uniq = Range.closed((int) edge.olapR(), qry_seqs.get(target).seq_ln());

						if(source_uniq.upperEndpoint()-source_uniq.lowerEndpoint()<min_uniq_olap || 
								target_uniq.upperEndpoint()-target_uniq.lowerEndpoint()<min_uniq_olap) 
							// unique region is too small
							continue;

						loop_on_subseqs:
							for(int j=0; j<sub_seqs.size(); j++) {

								final List<SAMSegment> source_sel = new ArrayList<>();
								source_segs = initPlace.get(source).get(j);
								for(final SAMSegment source_seg : source_segs) {
									source_cov = Range.closed(source_seg.qstart(), source_seg.qend());
									if(!source_cov.isConnected(source_uniq)) continue;
									source_olap = source_cov.intersection(source_uniq);
									if(source_olap.upperEndpoint()-source_olap.lowerEndpoint()>=min_uniq_olap)
										source_sel.add(source_seg);
								}

								if(source_sel.isEmpty()) continue;

								final List<SAMSegment> target_sel = new ArrayList<>();
								target_segs = initPlace.get(target).get(j);
								for(final SAMSegment target_seg : target_segs) {
									target_cov = Range.closed(target_seg.qstart(), target_seg.qend());
									if(!target_cov.isConnected(target_uniq)) continue;
									target_olap = target_cov.intersection(target_uniq);
									if(target_olap.upperEndpoint()-target_olap.lowerEndpoint()>=min_uniq_olap)
										continue;
									target_sel.add(target_seg);
								}

								if(target_sel.isEmpty()) continue;

								for(SAMSegment source_seg : source_sel) {
									for(SAMSegment target_seg : source_sel) {
										if(AlignmentSegment.sdistance(source_seg, target_seg)<hc_gap) {
											linkevi[i] = true;
											break loop_on_subseqs;
										}
									}
								}
							}
					}

					for(int i=0; i<n; i++) {
						if(linkevi[i]) {
							edgesToRetain.add(edges.get(i));
						} else {
							edgesToRemove.add(edges.get(i));
						}
					}
				//}
			}
		}
		
		myLogger.info("####Edge to remove "+edgesToRemove.size());
		myLogger.info("####Edge to retain "+edgesToRetain.size());
		
		edgesToRemove.removeAll(edgesToRetain);
		gfa.removeAllEdges(edgesToRemove);
		
		myLogger.info("####Edge removed "+edgesToRemove.size());
		
		final List<Set<String>> conns = 
				new ArrayList<>(new ConnectivityInspector<String, OverlapEdge>(graph).connectedSets());
		final List<DirectedWeightedOverlapPseudograph<String>> subgraphs = new ArrayList<>();
		final Map<String, Integer> hashv = new HashMap<>();
		final int nG = conns.size();
		for(int i=0; i<nG; i++) {
			final Set<String> conn = conns.get(i);
			final DirectedWeightedOverlapPseudograph<String> subgraph = 
					new DirectedWeightedOverlapPseudograph<>(OverlapEdge.class);
			for(final String v : conn) {
				subgraph.addVertex(v);
				hashv.put(v, i);
			}
			subgraphs.add(subgraph);
		}

		String source, target;
		for(final OverlapEdge edge : graph.edgeSet()) {
			source = gfa.getEdgeSource(edge);
			target = gfa.getEdgeTarget(edge);
			if(hashv.get(source).intValue()!=hashv.get(target).intValue())
				throw new RuntimeException("!!!");
			final DirectedWeightedOverlapPseudograph<String> subgraph =
					subgraphs.get(hashv.get(source));
			final OverlapEdge e = subgraph.addEdge(source, target);
			e.setOlap(edge.olap());
			e.setOlapF(edge.olapF());
			e.setOlapR(edge.olapR());
			e.setOlapInfo(edge.olapInfo());
			e.setCigar(edge.cigar());
			e.setRealigned(e.isRealigned());
		}

		int cycled = 0, nV = 0, nE = 0, maxV = 0;
		for(final DirectedWeightedOverlapPseudograph<String> subgraph : subgraphs) {
			nV += subgraph.vertexSet().size();
			nE += subgraph.edgeSet().size();
			maxV = Math.max(maxV, subgraph.vertexSet().size());
			final CycleDetector<String, OverlapEdge> cycleDetector = new CycleDetector<>(subgraph);
			if(cycleDetector.detectCycles()) {
				++cycled;
				int cL = 0;
				for(String v : subgraph.vertexSet()) {
					cL += qry_seqs.get(v).seq_ln();
				}
				myLogger.info("@ "+subgraph.vertexSet().size()+" "+cL);
				// if(subgraph.vertexSet().size()<100)
				//	 plot(subgraph);
			}
		}

		if(nV!=graph.vertexSet().size()||nE!=graph.edgeSet().size())
			throw new RuntimeException("!!!");

		myLogger.info("#G "+nG+", maxG "+maxV+", #O "+cycled);
		
		for(final DirectedWeightedOverlapPseudograph<String> subgraph : subgraphs) {
			
			final CycleDetector<String, OverlapEdge> cycleDetector = new CycleDetector<>(subgraph);
			if(cycleDetector.detectCycles()) continue;
			if(subgraph.vertexSet().size()==1) continue;
			
			final List<AlignmentSegmentDirectedWeightedPseudograph<String>> gs = new ArrayList<>();
			
			for(int i=0; i<this.sub_seqs.size(); i++) {

				final List<Traceable> seqBySubAll = new ArrayList<>();
				final List<SAMSegment> segs  = new ArrayList<SAMSegment>();
				for(final String qry_seq : subgraph.vertexSet()) {
					segs.clear();
					for(SAMSegment seg : initPlace.get(qry_seq).get(i)) 
						segs.add(seg);
					if(segs.isEmpty()) continue;
					Collections.sort(segs, new AlignmentSegment.SubjectCoordinationComparator());
					int n = segs.size();
					int start = 0, end = 0;
					while(start<n) {
						end = start+1;
						while(end<n&&AlignmentSegment.sdistance(segs.get(end-1), segs.get(end))<=as_gap)
							++end;

						RangeSet<Integer> subCov = TreeRangeSet.create();
						RangeSet<Integer> qryCov = TreeRangeSet.create();
						for(int j=start; j<end; j++) {
							SAMSegment seg = segs.get(j);
							subCov.add(Range.closed(seg.sstart(), seg.send()).canonical(DiscreteDomain.integers()));
							qryCov.add(Range.closed(seg.qstart(), seg.qend()).canonical(DiscreteDomain.integers()));
						}

						Range<Integer> subSpan = subCov.span();
						int sstart = subSpan.lowerEndpoint();
						int send   = subSpan.upperEndpoint();
						int ls = 0;
						for(final Range<Integer> r : subCov.asRanges())
							ls += r.upperEndpoint()-r.lowerEndpoint();

						Range<Integer> qrySpan = qryCov.span();
						int qstart = qrySpan.lowerEndpoint();
						int qend   = qrySpan.upperEndpoint();
						int lq = 0;
						for(final Range<Integer> r : qryCov.asRanges())
							lq += r.upperEndpoint()-r.lowerEndpoint();
						int slen = ls;
						int qlen = lq;

						seqBySubAll.add(new Traceable(qry_seq, dict.getSequences().get(i).getSequenceName(), 
								qstart, qend, sstart, send, subCov, qryCov, slen, qlen));

						start = end;
					}
				}

				Collections.sort(seqBySubAll, new AlignmentSegment.SubjectCoordinationComparator());
				int nSeg = seqBySubAll.size();
				
				final Set<Integer> processed = new HashSet<>();
				final Set<Integer> contained = new HashSet<>();
				int root_segix, source_segix, target_segix;
				Traceable root_seg, source_seg, target_seg;
				String root_segid, source_segid, target_segid;
				final Deque<Integer> deque = new ArrayDeque<Integer>();

				for(int j=0; j<nSeg; j++) {

					if(processed.contains(j)) continue;
					deque.clear();
					deque.add(j);
					contained.clear();
					contained.add(j);

					final AlignmentSegmentDirectedWeightedPseudograph<String> g = 
							new AlignmentSegmentDirectedWeightedPseudograph<>(DefaultWeightedEdge.class);
					root_segix = j;
					root_seg   = seqBySubAll.get(root_segix);
					root_segid = root_seg.qseqid();
					g.addVertex(root_segid, root_seg);

					while(!deque.isEmpty()) {
						source_segix = deque.pop();
						source_seg   = seqBySubAll.get(source_segix);
						source_segid = source_seg.qseqid();

						for(int k=source_segix-1; k>0; k--) {
							target_segix = k;
							target_seg   = seqBySubAll.get(target_segix);
							target_segid = target_seg.qseqid();

							if(source_seg.sstart()-target_seg.sstart()>hc_gap) break;

							if(subgraph.containsEdge(source_segid, target_segid)) {
								if(g.containsVertex(target_segid)) {
									if(!contained.contains(target_segix)) {
										myLogger.info("########\n"+
												"#"+target_seg.toString()+"\n"+
												"#"+g.getAlignmentSegmentByVertex(target_segid).toString()+"\n"+
												"########");
										//throw new RuntimeException("!!!");	
									}
								} else {
									g.addVertex(target_segid, target_seg);
								}
								if(!g.containsEdge(source_segid, target_segid))
									g.addEdge(source_segid, target_segid);
							} else if(subgraph.containsEdge(target_segid, source_segid)) {
								if(g.containsVertex(target_segid)) {
									if(!contained.contains(target_segix)) {
										myLogger.info("########\n"+
												"#"+target_seg.toString()+"\n"+
												"#"+g.getAlignmentSegmentByVertex(target_segid).toString()+"\n"+
												"########");
										//throw new RuntimeException("!!!");	
									}
								} else {
									g.addVertex(target_segid, target_seg);
								}
								if(!g.containsEdge(target_segid, source_segid)) 
									g.addEdge(target_segid, source_segid);
							} else continue;

							if(!contained.contains(target_segix)) {
								deque.add(target_segix);
								contained.add(target_segix);
							}
						}

						for(int k=source_segix+1; k<nSeg; k++) {
							target_segix = k;
							target_seg   = seqBySubAll.get(target_segix);
							target_segid = target_seg.qseqid();

							if(target_seg.sstart()-source_seg.sstart()>hc_gap) break;

							if(subgraph.containsEdge(source_segid, target_segid)) {
								if(g.containsVertex(target_segid)) {
									if(!contained.contains(target_segix)) {
										myLogger.info("########\n"+
												"#"+target_seg.toString()+"\n"+
												"#"+g.getAlignmentSegmentByVertex(target_segid).toString()+"\n"+
												"########");
										//throw new RuntimeException("!!!");	
									}
								} else {
									g.addVertex(target_segid, target_seg);
								}
								if(!g.containsEdge(source_segid, target_segid))
									g.addEdge(source_segid, target_segid);
							} else if(subgraph.containsEdge(target_segid, source_segid)) {
								if(g.containsVertex(target_segid)) {
									if(!contained.contains(target_segix)) {
										myLogger.info("\n"+
												"############\n"+
												"##"+target_seg.toString()+"\n"+
												"##"+g.getAlignmentSegmentByVertex(target_segid).toString()+"\n"+
												"############");
										//throw new RuntimeException("!!!");	
									}
								} else {
									g.addVertex(target_segid, target_seg);
								}
								if(!g.containsEdge(target_segid, source_segid)) 
									g.addEdge(target_segid, source_segid);
							} else continue;

							if(!contained.contains(target_segix)) {
								deque.add(target_segix);
								contained.add(target_segix);
							}
						}

						processed.add(source_segix);
					}

					if(contained.size()!=g.vertexSet().size()) {
						myLogger.info("@@@@@@@@");
						for(final int k : contained) {
							myLogger.info("@@"+seqBySubAll.get(k).toString());
						}
						myLogger.info("@@@@@@@@");
					}

					gs.add(g);
				}
			}
			
			if(ddebug) {
				myLogger.info("########");
				myLogger.info("#graphs constructed "+gs.size()+" from subgraph");
				myLogger.info("#subgraph #V "+subgraph.vertexSet().size());
				
				Collections.sort(gs, new Comparator<AlignmentSegmentDirectedWeightedPseudograph<String>>() {

					@Override
					public int compare(AlignmentSegmentDirectedWeightedPseudograph<String> g0,
							AlignmentSegmentDirectedWeightedPseudograph<String> g1) {
						// TODO Auto-generated method stub
						return g1.vertexSet().size()-g0.vertexSet().size();
					}
					
				});
				
				myLogger.info("#       largest #V "+gs.get(0).vertexSet().size());
				if(gs.size()>1)
					myLogger.info("#second largest #V "+gs.get(1).vertexSet().size());
				myLogger.info("########");
			}
		}
	}

	@SuppressWarnings({ "rawtypes", "unchecked" })
	private void plot(Graph graph) {
		// TODO Auto-generated method stub
		final ListenableGraph listenable = new DefaultListenableGraph<>(graph);
		SwingUtilities.invokeLater(new Runnable() {
			public void run(){

				JGraphXAdapter graphAdapter = 
						new JGraphXAdapter<>(listenable);
				mxOrganicLayout layout = new mxOrganicLayout(graphAdapter);

				layout.setOptimizeBorderLine(true);
				layout.setOptimizeEdgeCrossing(true);
				layout.setOptimizeEdgeDistance(true);
				layout.setOptimizeEdgeLength(true);
				layout.setOptimizeNodeDistribution(true);
				layout.setFineTuning(true);
				layout.setMaxIterations(10000);
				layout.execute(graphAdapter.getDefaultParent());

				JFrame frame = new JFrame("Graphmap");
				frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
				frame.add(new mxGraphComponent(graphAdapter));
				frame.pack();
				frame.setSize(800, 800);
				frame.setVisible(true);
			}
		});
	}



	public void run2() {
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

		// read alignment file and place the query sequences
		final Map<String, Set<SAMSegment>> initPlace = new HashMap<String, Set<SAMSegment>>();
		try {
			final SamReaderFactory factory =
					SamReaderFactory.makeDefault()
					.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
							SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
					.validationStringency(ValidationStringency.SILENT);
			final SamReader in1 = factory.open(new File(align_file));
			final SAMRecordIterator iter1 = in1.iterator();


			String qry;
			final List<SAMSegment> buff = new ArrayList<SAMSegment>();
			SAMRecord rc = iter1.next();

			while(rc!=null) {
				qry = rc.getReadName();	

				buff.clear();
				if(!rc.getReadUnmappedFlag()&&rc.getMappingQuality()>0)
					buff.add(SAMSegment.samRecord(rc, true));

				while( (rc=iter1.next())!=null
						&&
						rc.getReadName().equals(qry) ) {
					if(!rc.getReadUnmappedFlag()&&rc.getMappingQuality()>0)
						buff.add(SAMSegment.samRecord(rc, true));
				}

				if(buff.isEmpty()) continue;

				Set<SAMSegment> init_f = new HashSet<SAMSegment>();
				Set<SAMSegment> init_r = new HashSet<SAMSegment>();
				for(SAMSegment record : buff) {
					if(record.qlength()<min_alen) 
						continue;
					if(record.qseqid().equals(qry)) 
						init_f.add(record);
					else init_r.add(record);
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

		this.initial_thread_pool();
		for(String sub_seq : sub_seqs.keySet()) {

			this.executor.submit(new Runnable() {

				String sub_seq = null;

				@Override
				public void run() {
					// TODO Auto-generated method stub
					try {
						if(sub_seq.equals("Chr00")) return;
						myLogger.info(">>>>>>>>>>>>>"+sub_seq+"<<<<<<<<<<<<<<<<");

						final List<Traceable> seqBySubAll = new ArrayList<>();
						final List<SAMSegment> segs  = new ArrayList<SAMSegment>();
						for(final String qry_seq : initPlace.keySet()) {
							segs.clear();
							for(SAMSegment seg : initPlace.get(qry_seq)) 
								if(seg.sseqid().equals(sub_seq))
									segs.add(seg);
							if(segs.isEmpty()) continue;
							Collections.sort(segs, new AlignmentSegment.SubjectCoordinationComparator());
							int n = segs.size();
							int start = 0, end = 0;
							while(start<n) {
								end = start+1;
								while(end<n&&AlignmentSegment.sdistance(segs.get(end-1), segs.get(end))<=as_gap)
									++end;

								//if(end-start>1)
								//	myLogger.info("debug point");

								RangeSet<Integer> subCov = TreeRangeSet.create();
								RangeSet<Integer> qryCov = TreeRangeSet.create();
								for(int i=start; i<end; i++) {
									SAMSegment seg = segs.get(i);
									subCov.add(Range.closed(seg.sstart(), seg.send()).canonical(DiscreteDomain.integers()));
									qryCov.add(Range.closed(seg.qstart(), seg.qend()).canonical(DiscreteDomain.integers()));
								}

								Range<Integer> subSpan = subCov.span();
								int sstart = subSpan.lowerEndpoint();
								int send   = subSpan.upperEndpoint();
								int ls = 0;
								for(final Range<Integer> r : subCov.asRanges())
									ls += r.upperEndpoint()-r.lowerEndpoint();

								Range<Integer> qrySpan = qryCov.span();
								int qstart = qrySpan.lowerEndpoint();
								int qend   = qrySpan.upperEndpoint();
								int lq = 0;
								for(final Range<Integer> r : qryCov.asRanges())
									lq += r.upperEndpoint()-r.lowerEndpoint();
								int slen = ls;
								int qlen = lq;

								seqBySubAll.add(new Traceable(qry_seq, sub_seq, qstart, qend, sstart, send, subCov, qryCov, slen, qlen));

								start = end;
							}
						}


						Collections.sort(seqBySubAll, new AlignmentSegment.SubjectCoordinationComparator());
						int nSeg = seqBySubAll.size();

						myLogger.info("#segments "+seqBySubAll.size());

						final Set<Integer> processed = new HashSet<>();
						final Set<Integer> contained = new HashSet<>();
						int root_segix, source_segix, target_segix;
						Traceable root_seg, source_seg, target_seg;
						String root_segid, source_segid, target_segid;
						final Deque<Integer> deque = new ArrayDeque<Integer>();

						final List<AlignmentSegmentDirectedWeightedPseudograph<String>> subgraphs = new ArrayList<>();
						for(int i=0; i<nSeg; i++) {

							if(ddebug&&i%100000==0) myLogger.info("####"+this.sub_seq+" #processed "+i);

							if(processed.contains(i)) continue;
							deque.clear();
							deque.add(i);
							contained.clear();
							contained.add(i);

							final AlignmentSegmentDirectedWeightedPseudograph<String> subgraph = 
									new AlignmentSegmentDirectedWeightedPseudograph<>(DefaultWeightedEdge.class);
							root_segix = i;
							root_seg   = seqBySubAll.get(root_segix);
							root_segid = root_seg.qseqid();
							subgraph.addVertex(root_segid, root_seg);

							while(!deque.isEmpty()) {
								source_segix = deque.pop();
								source_seg   = seqBySubAll.get(source_segix);
								source_segid = source_seg.qseqid();

								for(int j=source_segix-1; j>0; j--) {
									target_segix = j;
									target_seg   = seqBySubAll.get(target_segix);
									target_segid = target_seg.qseqid();

									if(source_seg.sstart()-target_seg.sstart()>hc_gap) break;

									if(gfa.containsEdge(source_segid, target_segid)) {
										if(subgraph.containsVertex(target_segid)) {
											if(!contained.contains(target_segix)) {
												myLogger.info("########\n"+
														"#"+target_seg.toString()+"\n"+
														"#"+subgraph.getAlignmentSegmentByVertex(target_segid).toString()+"\n"+
														"########");
												//throw new RuntimeException("!!!");	
											}
										} else {
											subgraph.addVertex(target_segid, target_seg);
										}
										if(!subgraph.containsEdge(source_segid, target_segid))
											subgraph.addEdge(source_segid, target_segid);
									} else if(gfa.containsEdge(target_segid, source_segid)) {
										if(subgraph.containsVertex(target_segid)) {
											if(!contained.contains(target_segix)) {
												myLogger.info("########\n"+
														"#"+target_seg.toString()+"\n"+
														"#"+subgraph.getAlignmentSegmentByVertex(target_segid).toString()+"\n"+
														"########");
												//throw new RuntimeException("!!!");	
											}
										} else {
											subgraph.addVertex(target_segid, target_seg);
										}
										if(!subgraph.containsEdge(target_segid, source_segid)) 
											subgraph.addEdge(target_segid, source_segid);
									} else continue;

									if(!contained.contains(target_segix)) {
										deque.add(target_segix);
										contained.add(target_segix);
									}
								}

								for(int j=source_segix+1; j<nSeg; j++) {
									target_segix = j;
									target_seg   = seqBySubAll.get(target_segix);
									target_segid = target_seg.qseqid();

									if(target_seg.sstart()-source_seg.sstart()>hc_gap) break;

									if(gfa.containsEdge(source_segid, target_segid)) {
										if(subgraph.containsVertex(target_segid)) {
											if(!contained.contains(target_segix)) {
												myLogger.info("########\n"+
														"#"+target_seg.toString()+"\n"+
														"#"+subgraph.getAlignmentSegmentByVertex(target_segid).toString()+"\n"+
														"########");
												//throw new RuntimeException("!!!");	
											}
										} else {
											subgraph.addVertex(target_segid, target_seg);
										}
										if(!subgraph.containsEdge(source_segid, target_segid))
											subgraph.addEdge(source_segid, target_segid);
									} else if(gfa.containsEdge(target_segid, source_segid)) {
										if(subgraph.containsVertex(target_segid)) {
											if(!contained.contains(target_segix)) {
												myLogger.info("\n"+
														"############\n"+
														"##"+target_seg.toString()+"\n"+
														"##"+subgraph.getAlignmentSegmentByVertex(target_segid).toString()+"\n"+
														"############");
												//throw new RuntimeException("!!!");	
											}
										} else {
											subgraph.addVertex(target_segid, target_seg);
										}
										if(!subgraph.containsEdge(target_segid, source_segid)) 
											subgraph.addEdge(target_segid, source_segid);
									} else continue;

									if(!contained.contains(target_segix)) {
										deque.add(target_segix);
										contained.add(target_segix);
									}
								}

								processed.add(source_segix);
							}

							if(contained.size()!=subgraph.vertexSet().size()) {
								myLogger.info("@@@@@@@@");
								for(final int k : contained) {
									myLogger.info("@@"+seqBySubAll.get(k).toString());
								}
								myLogger.info("@@@@@@@@");
							}

							subgraphs.add(subgraph);
						}

						if(debug) {
							int numSegs = 0, maxV = -1, minV = Integer.MAX_VALUE, cycled = 0;
							int[] sizes = new int[10];
							for(final AlignmentSegmentDirectedWeightedPseudograph<String> subgraph : subgraphs) {
								int nV = subgraph.vertexSet().size();
								if(nV<2)        ++sizes[0];
								else if(nV<10)  ++sizes[1];
								else if(nV<50)  ++sizes[2];
								else if(nV<100) ++sizes[3];
								else            ++sizes[4];
								if(maxV<nV) maxV = nV;
								if(minV>nV) minV = nV;
								numSegs += nV;

								final CycleDetector<String, DefaultWeightedEdge> cycleDetector = new CycleDetector<>(subgraph);
								if(cycleDetector.detectCycles()) {
									myLogger.info("#cycle V "+subgraph.vertexSet().size()+" "+cycleDetector.findCycles().size());
									++cycled;
								}
							}
							int avgV = numSegs/subgraphs.size();
							myLogger.info("#### "+this.sub_seq+" ####\n"+
									"#segments "+seqBySubAll.size()+"\n"+
									"#placed    : "+numSegs+"\n"+
									"#subgraphs : "+subgraphs.size()+"\n"+
									"#max V     : "+maxV+"\n"+
									"#min V     : "+minV+"\n"+
									"#avg V     : "+avgV+"\n"+
									"#cycled    : "+cycled+"\n"+
									"# size     : \n"+
									"     1     : "+sizes[0]+"\n"+
									"   <10     : "+sizes[1]+"\n"+
									"   <50     : "+sizes[2]+"\n"+
									"  <100     : "+sizes[3]+"\n"+
									"  >100     : "+sizes[4]+"\n");
						}

						for(final AlignmentSegmentDirectedWeightedPseudograph<String> subgraph : subgraphs) {
							if(subgraph.vertexSet().size()==1) continue;

							final CycleDetector<String, DefaultWeightedEdge> cycleDetector = new CycleDetector<>(subgraph);
							// skip subgraph with cycles
							if(cycleDetector.detectCycles()) continue;

							// find all vertices with no incoming edges
							final List<String> starts = new ArrayList<String>();
							for(final String v : subgraph.vertexSet())
								if(subgraph.incomingEdgesOf(v).isEmpty())
									starts.add(v);

							// now for each 'start' find a best path
							// the final path will be the best of the bests
							for(final String start : starts) {

							}
						}
					} catch (Exception e) {
						Thread t = Thread.currentThread();
						t.getUncaughtExceptionHandler().uncaughtException(t, e);
						e.printStackTrace();
						executor.shutdown();
						System.exit(1);
					}
				}

				public Runnable init(String sub_seq) {
					// TODO Auto-generated method stub
					this.sub_seq = sub_seq;
					return this;
				}		
			}.init(sub_seq));
		}
		this.waitFor();
	}
	
	private final class Traceable extends CompoundAlignmentSegment {
		private Traceable prev, next;
		private double score;

		public Traceable(final String qseqid,   // query (e.g., gene) sequence id
				final String sseqid,   // subject (e.g., reference genome) sequence id
				final int qstart,      // start of alignment in query
				final int qend,        // end of alignment in query
				final int sstart,      // start of alignment in subject
				final int send,         // end of alignment in subject
				final RangeSet<Integer> subCov,
				final RangeSet<Integer> qryCov,
				final int slen,
				final int qlen) {
			super(qseqid, sseqid, qstart,qend,sstart,send, subCov, qryCov, slen, qlen);
		}

		public void setPrev(final Traceable prev) {
			this.prev = prev;
		}

		public void setNext(final Traceable next) {
			this.next = next;
		}

		public void setScore(final double score) {
			this.score = score;
		}

		public Traceable getPrev() {
			return this.prev;
		}

		public Traceable getNext() {
			return this.next;
		}

		public double getScore() {
			return this.score;
		}
	}
}















