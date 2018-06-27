package cz1.ngs.tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.apache.commons.collections4.bidimap.TreeBidiMap;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.DirectedPseudograph;
import org.jgrapht.graph.DirectedWeightedPseudograph;

import com.google.common.collect.DiscreteDomain;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

import cz1.ngs.model.AlignmentSegment;
import cz1.ngs.model.GFA;
import cz1.ngs.model.OverlapEdge;
import cz1.ngs.model.SAMSegment;
import cz1.ngs.model.Sequence;
import cz1.ngs.model.TraceableAlignmentSegment;
import cz1.ngs.model.TraceableDirectedWeightedPseudograph;
import cz1.ngs.model.TraceableEdge;
import cz1.ngs.model.TraceableVertex;
import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Executor;
import cz1.util.Utils;
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
	private double min_alen = 100;        // keep the alignment record if the segment is no smaller than this
	// without considering the completeness
	private int min_overlap = 10;         // minimum overlap length
	private double collinear_shift = 1.0; // maximum shift distance for two collinear alignment segments - 50% of the smaller segment size
	private int kmer_size = -1;           // kmer size to construct the de bruijn assembly graph

	private int num_threads = Runtime.getRuntime().availableProcessors();
	private String out_prefix = null;
	private boolean debug  = false;
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

	private final static int max_clip = 100;  // maximum clip of query alignment allowed
	private final static int gap_buff = 30;  // buffer size for subject/reference sequences gap clips
	private final static int min_len  = 100; // minimum alignment length
	private final static int min_gap  = 10;  // take this if estimated gap size is smaller than this
	private final static int max_gap  = 100; // take this if estimated gap size is larger than this
	private final static int min_ext  = 30;  // minimum extension for contigging
	private final static double m_clip = 0.1d; // max clip size (%) to treat an alignment end-to-end
	private final static double olap_min = 0.99d; // min overlap fraction for containment
	
	private Map<String, Sequence> qry_seqs;
	private Map<String, Sequence> sub_seqs;
	private Map<String, TreeRangeSet<Integer>> sub_gaps;
	private static List<SAMSegment> sam_records;

	private final static int max_dist  = 100000;

	private final static int match_score  = 1;
	private final static int clip_penalty = 1;
	private final static int hc_gap  = 100000;
	private final static int max_cov = 127;
	private final static int aln_flank = 50;
	
	private final static Object lock = new Object();

	@Override
	public void run() {
		// this.run1();
		// this.run2();
		this.run3();
	}

	private final static int min_olap  = 50;
	private final static int max_clip2 = 30;
	
	public void run3() {
		sub_seqs = Sequence.parseFastaFileAsMap(subject_file);
		qry_seqs = Sequence.parseFastaFileWithRevCmpAsMap(query_file);

		// find 'N/n's in subject/reference sequences
		// which could have impact on parsing the blast records
		final Map<String, TreeRangeSet<Integer>> sub_gaps = new HashMap<String, TreeRangeSet<Integer>>();
		final Map<String, List<SAMSegment>> anchored_records = new HashMap<String, List<SAMSegment>>();

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
			// initialise an array list for each reference chromosome
			anchored_records.put(seq_sn, new ArrayList<SAMSegment>());
		}

		// parse SAMRecords
		// SAMRecords buffer
		final List<SAMSegment> buff = new ArrayList<SAMSegment>();
		// selected SAMRecords
		final List<SAMSegment> sel_recs = new ArrayList<SAMSegment>();
		// temp list
		final List<SAMSegment> tmp_records = new ArrayList<SAMSegment>();

		try {
			final SamReaderFactory factory =
					SamReaderFactory.makeDefault()
					.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
							SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
					.validationStringency(ValidationStringency.SILENT);
			final SamReader in1 = factory.open(new File(align_file));
			final SAMRecordIterator iter1 = in1.iterator();

			SAMRecord rc = iter1.next();
			SAMSegment primary_record, secondary_record;
			String qry;
			int qry_ln;
			double thres_ln;

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

				sel_recs.clear();
				// merge collinear records
				for(String sub : sub_seqs.keySet()) {
					// get all records for subject/reference sequence sub_seq
					tmp_records.clear();
					for(SAMSegment record : buff)
						if(record.sseqid().equals(sub))
							tmp_records.add(record);

					if(tmp_records.isEmpty()) continue;

					// find alignment segments that can be deleted
					// those that are subsets of larger alignment segments
					// (Sstart, (sstart, send), Send) and (Qstart, (qstart, qend), Qend)
					Collections.sort(tmp_records, new SAMSegment.SegmentSizeComparator());
					final Set<int[][]> ranges = new HashSet<int[][]>();
					outerloop:
						for(int i=0; i<tmp_records.size(); i++) {
							primary_record = tmp_records.get(i);
							int[][] range = new int[2][2];
							if(primary_record.sstart()<primary_record.send()) {
								range[0][0] = primary_record.sstart();
								range[0][1] = primary_record.send();
							} else {
								range[0][0] = primary_record.send();
								range[0][1] = primary_record.sstart();
							}
							if(primary_record.qstart()<primary_record.qend()) {
								range[1][0] = primary_record.qstart();
								range[1][1] = primary_record.qend();
							} else {
								range[1][0] = primary_record.qend();
								range[1][1] = primary_record.qstart();
							}
							for(int[][] r : ranges) {
								if(r[0][0]<=range[0][0]&&
										r[0][1]>=range[0][1]&&
										r[1][0]<=range[1][0]&&
										r[1][1]>=range[1][1]) {
									tmp_records.remove(i--);
									continue outerloop;
								}
							}
							ranges.add(range);
						}

					// find collinear alignment segments that can be merged
					Collections.sort(tmp_records, new AlignmentSegment.SubjectCoordinationComparator());

					SAMSegment record;
					int pqstart, pqend, sqstart, sqend, olap;
					for(int i=0; i<tmp_records.size(); i++) {
						primary_record = tmp_records.get(i);
						pqstart = primary_record.true_qstart();
						pqend   = primary_record.true_qend();

						for(int j=i+1; j<tmp_records.size(); j++) {
							secondary_record = tmp_records.get(j);
							// we need to check if they are two copies that are close to each other
							// TODO: how?
							// we check the region of alignments on the query sequence
							// if the overlap is greater than 90%, we regard them as different copies
							sqstart = secondary_record.true_qstart();
							sqend   = secondary_record.true_qend();

							olap = 0;
							if(pqstart<=sqstart) olap = Math.max(0, pqend-sqstart);
							if(pqstart> sqstart) olap = Math.max(0, sqend-pqstart);

							if( (double)olap/(pqend-pqstart+1)>=0.9||(double)olap/(sqend-sqstart+1)>=0.9 ) continue;

							// now we check if they two are collinear
							double max_shift = collinear_shift*
									Math.min(primary_record.qlength(), secondary_record.qlength());
							if( (record=SAMSegment.collinear(primary_record, secondary_record, max_shift))!=null ) {
								tmp_records.set(i, record);
								tmp_records.remove(j);
								--i;
								break;
							}
						}
					}

					// process blast records that clipped by gaps
					// (sstart, send)---(start2, send2)
					// (sstart  ...  ---  ...    send2)
					TreeRangeSet<Integer> sub_gap = sub_gaps.get(sub);
					Collections.sort(tmp_records, new AlignmentSegment.SubjectCoordinationComparator());

					for(int i=0; i<tmp_records.size(); i++) {
						primary_record = tmp_records.get(i);
						if( sub_gap.contains(primary_record.true_send()) ) {
							secondary_record = null;
							int sec_j = -1;
							for(int j=i+1; j<tmp_records.size(); j++) {
								if( tmp_records.get(j).true_sstart()>=
										primary_record.true_send() ) {
									secondary_record = tmp_records.get(j);
									sec_j = j;
									break;
								}
							}
							if(secondary_record==null || 
									AlignmentSegment.reverse(primary_record, secondary_record)) {
								// no clipping
								// reverse alignment segments
								continue;
							}

							if( sub_gap.contains(secondary_record.true_sstart()) &&
									sub_gap.rangeContaining(primary_record.true_send()).
									equals(sub_gap.rangeContaining(secondary_record.true_sstart()))) {
								// clipping
								// merge two alignment segments
								int qstart = Math.min(primary_record.true_qstart(), secondary_record.true_qstart());
								int qend = Math.max(primary_record.true_qend(), secondary_record.true_qend());
								int sstart = Math.min(primary_record.true_sstart(), secondary_record.true_sstart());
								int send = Math.max(primary_record.true_send(), secondary_record.true_send());

								// replace primary record with merged record
								// delete secondary record
								SAMSegment merged_record = primary_record.forward()?
										new SAMSegment(qry,sub,qstart,qend,sstart,send):
											new SAMSegment(qry,sub,qstart,qend,send,sstart);
										tmp_records.set(i, merged_record);
										tmp_records.remove(sec_j);

										// the merged records need to be processed
										--i;
							}
						}
					}

					// add to sel_recs
					sel_recs.addAll(tmp_records);
				}

				if(sel_recs.isEmpty()) continue;

				Collections.sort(sel_recs, new SAMSegment.SegmentSizeComparator());
				primary_record = sel_recs.get(0);
				thres_ln = primary_record.qlength()*this.min_frac;
				anchored_records.get(primary_record.sseqid()).add(primary_record);

				for(int i=1; i<sel_recs.size(); i++) {
					secondary_record = sel_recs.get(i);
					if(secondary_record.qlength()>=thres_ln)
						anchored_records.get(secondary_record.sseqid()).add(secondary_record);
				}
			}

			iter1.close();
			in1.close();

			for(Map.Entry<String, List<SAMSegment>> entry : anchored_records.entrySet()) {
				myLogger.info("initial placement #entry "+entry.getKey()+": "+entry.getValue().size());
				if(this.ddebug) {
					List<SAMSegment> segs = entry.getValue();
					Collections.sort(segs, new SAMSegment.SubjectCoordinationComparator());
					for(SAMSegment seg : segs)
						myLogger.info(seg.toString());
				}
			}

			myLogger.info("####process each reference chromosome");

			/***
			final BufferedWriter bw_map = Utils.getBufferedWriter(out_prefix+".map");
			final BufferedWriter bw_fa = Utils.getBufferedWriter(out_prefix+".fa");
			***/
			final Set<String> anchored_contigs = new HashSet<String>();
			
			this.initial_thread_pool();
			for(String sub_seq : sub_seqs.keySet()) {

				this.executor.submit(new Runnable() {

					String sub_seq = null;

					@Override
					public void run() {
						// TODO Auto-generated method stub
						try {
							if(sub_seq.equals("Chr00")) return;
							final int sub_ln = sub_seqs.get(sub_seq).seq_ln();
							
							// make traceable alignment segments from SAM segments
							final List<TraceableAlignmentSegment> alignments = new ArrayList<TraceableAlignmentSegment>();
							TraceableAlignmentSegment alignment;
							int qry_ln, qstart, qend, sstart, send;
							double clip;
							for(final SAMSegment record : anchored_records.get(sub_seq)) {
								qstart = record.qstart();
								qend   = record.qend();
								sstart = record.sstart();
								send   = record.send();
								alignment = new TraceableAlignmentSegment(record.qseqid(), record.sseqid(), 
										qstart, qend, sstart, send);
								// check if this is end-to-end alignment
								qry_ln = qry_seqs.get(record.qseqid()).seq_ln();
								clip = Math.min(max_clip, qry_ln*m_clip);
								alignment.setEndToEnd(qstart<=clip, qend+clip>=qry_ln, false, false);
								alignment.setClip(qstart-1, qry_ln-qend, sstart-1, sub_ln-send);
								alignment.calcScore();
								alignments.add(alignment);	
							}
							
							Collections.sort(alignments, new TraceableAlignmentSegment.SubjectCoordinationComparator());
							
							// we remove containment
							TraceableAlignmentSegment source_as, target_as, tmp_as;
							int source_sstart, source_send, target_sstart, target_send;
							int asz = alignments.size();
							boolean sEndToEnd, qEndToEnd;
							outerloop:
								for(int w=0; w<asz; w++) {
									source_as = alignments.get(w);
									if(source_as==null) continue;
									source_sstart = source_as.sstart();
									source_send   = source_as.send();
									sEndToEnd     = source_as.getToQueryStart()&&source_as.getToQueryEnd();
									innerloop:
										for(int z=w+1; z<asz; z++) {
											target_as = alignments.get(z);
											if(target_as==null) continue;
											target_sstart = target_as.sstart();
											target_send   = target_as.send();
											qEndToEnd     = target_as.getToQueryStart()&&target_as.getToQueryEnd();

											if(target_sstart>source_send) continue outerloop;
											
											if(sEndToEnd&&source_sstart<=target_sstart&&source_send>=target_send) {
												// target is contained
												alignments.set(z, null);
												continue innerloop;
											} 
											
											if(qEndToEnd&&target_sstart<=source_sstart&&target_send>=source_send) {
												// source is contained
												alignments.set(w, null);
												continue outerloop;
											}
									}
								}

							final List<TraceableAlignmentSegment> selected = new ArrayList<TraceableAlignmentSegment>();
							for(TraceableAlignmentSegment as : alignments) if(as!=null) selected.add(as);
							
							if(ddebug) {
								synchronized(lock) {
									myLogger.info("containment removed placement #entry "+sub_seq+": "+selected.size());
									for(TraceableAlignmentSegment seg : selected) myLogger.info(seg.toString());
								}
							}
														
							double objective = Double.NEGATIVE_INFINITY, dobj, opt_obj = Double.NEGATIVE_INFINITY, tmp; // record the current best path
							int tmp_sstart, tmp_send;
							// int source_qstart, source_qend, target_qstart, target_qend;
							TraceableAlignmentSegment traceback = null;
							asz = selected.size();

							// find the best path - dynamic programming
							for(int w=0; w<asz; w++) {
								source_as = selected.get(w);
								if(source_as.getTraceBackward()!=null) 
									continue;
								objective     = source_as.getScore();
								source_sstart = source_as.sstart();
								source_send   = source_as.send();

								// so we start from w+1
								int z = w+1;
								outerloop:
									while(z<asz) {
										// first find next 
										target_as     = selected.get(z);
										target_sstart = target_as.sstart();
										target_send   = target_as.send();

										if(source_send>=target_send) {
											++z;
											continue;
										}
										
										// calculate the additive to the objective
										dobj = target_as.getScore()+
												// penalty: gap on reference
												//(target_sstart>source_send?(target_sstart-source_send)*TraceableAlignmentSegment.gap_extension:0)+
												// subtract match score for overlap
												(target_sstart>source_send?0:(target_sstart-source_send)*TraceableAlignmentSegment.match_score);

										// #a <------------------->
										// #b      <----------------->
										// #c            <--------------->
										// we prefer #c over #b
										for(int v=z+1; v<asz; v++) {
											tmp_as     = selected.get(v);
											tmp_sstart = tmp_as.sstart();
											tmp_send   = tmp_as.send();
											
											if(tmp_sstart>source_send && 
													tmp_sstart>target_sstart)
												break;
											
											tmp = tmp_as.getScore()+
													// penalty: gap on reference
													//(tmp_sstart>source_send?(tmp_sstart-source_send)*TraceableAlignmentSegment.gap_extension:0)+
													// subtract match score for overlap
													(tmp_sstart>source_send?0:(tmp_sstart-source_send)*TraceableAlignmentSegment.match_score);
											
											if(tmp>dobj) {
												// update selected
												dobj          = tmp;
												target_as     = tmp_as;
												target_sstart = tmp_sstart;
												target_send   = tmp_send;
												z = v;
											}
										}

										// OK now we found it
										if(dobj+objective>target_as.getObjective()) {
											// if the objective is better, we choose this path
											objective += dobj;
											target_as.setTraceBackward(source_as);
											target_as.setObjective(objective);
											source_as.setTraceForward(target_as);
											
											// finally we update source
											source_as     = target_as;
											source_sstart = target_sstart;
											source_send   = target_send;

											if(source_as.getTraceForward()!=null) {
												// so we trace forward from here to avoid redoing this
												while((target_as=source_as.getTraceForward())!=null) {
													target_sstart = target_as.sstart();
													source_send   = source_as.send();
													dobj = target_as.getScore()+
															// penalty: gap on reference
															//(target_sstart>source_send?(target_sstart-source_send)*TraceableAlignmentSegment.gap_extension:0)+
															// subtract match score for overlap
															(target_sstart>source_send?0:(target_sstart-source_send)*TraceableAlignmentSegment.match_score);
													objective += dobj;
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
							if(traceback==null) throw new RuntimeException("!!!");
							
							final List<TraceableAlignmentSegment> graph_path = new ArrayList<TraceableAlignmentSegment>();
							while(traceback!=null) {
								graph_path.add(traceback);
								traceback = traceback.getTraceBackward();
							}
							Collections.reverse(graph_path);

							RangeSet<Integer> covs = TreeRangeSet.create();
							for(final TraceableAlignmentSegment seg : graph_path) 
								covs.add(Range.closed(seg.sstart(), seg.send()).canonical(DiscreteDomain.integers()));
							
							int cov = countIntervalCoverage(covs);
							myLogger.info("####reference chromosome covered "+sub_seq+": "+(double)cov/sub_ln+"("+cov+"/"+sub_ln+")");
							
							if(ddebug) {
								synchronized(lock) {
									myLogger.info("####final placement #entry "+sub_seq+": "+graph_path.size()+
											"["+(double)cov/sub_ln+"("+cov+"/"+sub_ln+")]");
									for(TraceableAlignmentSegment seg : graph_path) myLogger.info(seg.toString());
								}
							}
							
							// TODO: need to figure out why this is working so slow
							// right now we write the link file out and 
							// run the consensus step with an extra script
							BufferedWriter bw = Utils.getBufferedWriter(out_prefix+"_"+this.sub_seq+".link");
							bw.write("####final placement #entry "+sub_seq+": "+graph_path.size()+
											"["+(double)cov/sub_ln+"("+cov+"/"+sub_ln+")]\n");
							for(TraceableAlignmentSegment seg : graph_path) bw.write(seg.toString()+"\n");
							bw.close();
							synchronized(lock) {
								for(final TraceableAlignmentSegment seg : graph_path) { 
									String qseqid = seg.qseqid();
									qseqid = qseqid.replaceAll("'$", "");
									anchored_contigs.add(qseqid);
								}
							}
							
							/***
							myLogger.info("####consensus "+this.sub_seq+" initiated.");
							
							// now we join the neighbouring placement 
							final Set<String> anchored_seq = new HashSet<String>();
							final StringBuilder pseudo_chr = new StringBuilder();
							final List<String> agpmap_str  = new ArrayList<String>();
							alignment = graph_path.get(0);
							String qseqid = alignment.qseqid();
							pseudo_chr.append(qry_seqs.get(qseqid).seq_str());
							int chrL = pseudo_chr.length();
							if(qseqid.endsWith("'")) {
								// reverse complement
								qseqid = qseqid.replaceAll("'$", "");
								anchored_seq.add(qseqid);
								agpmap_str.add(qseqid+"\t"+chrL+"\t0\t"+chrL+"\t-\t"+sub_seq+"\t0\t"+chrL);
							} else {
								anchored_seq.add(qseqid);
								agpmap_str.add(qseqid+"\t"+chrL+"\t0\t"+chrL+"\t+\t"+sub_seq+"\t0\t"+chrL);
							}
							
							String target_str, source_aln, target_aln;
							int olap, gap, source_len, target_len, source_qend, target_qstart;
							int aLen, a2, b1, b2;
							SequencePair<DNASequence, NucleotideCompound> seqPair;
							final int n = graph_path.size();
							
							for(int i=1; i<graph_path.size(); i++) {
								myLogger.info("####consensus "+this.sub_seq+": "+i+"/"+n);
								
								source_as = graph_path.get(i-1);
								target_as = graph_path.get(i);
								
								source_qend = source_as.qend();
								source_send = source_as.send();
								source_len  = qry_seqs.get(source_as.qseqid()).seq_ln();
								
								target_qstart = target_as.qstart();
								target_sstart = target_as.sstart();
								
								qseqid = target_as.qseqid();
								target_str = qry_seqs.get(qseqid).seq_str();
								target_len = target_str.length();
								
								// see if there is overlap between source tail and target head
								olap = Math.max(0, source_send-target_sstart);
								
								source_aln = pseudo_chr.substring(Math.max(0, 
										chrL-(source_len-source_qend+olap+aln_flank)));
								target_aln = target_str.substring(0, Math.min(target_len,
										target_qstart+olap+aln_flank));
								seqPair = Constants.localPairMatcher(source_aln, target_aln);
								
								aLen = seqPair.getLength();
								a2 = seqPair.getIndexInQueryAt(aLen);
								b1 = seqPair.getIndexInTargetAt(1);
								b2 = seqPair.getIndexInTargetAt(aLen);
								
								if(aLen>=min_olap&&a2+max_clip2>=source_aln.length()&&max_clip2>=b1) {
									// overlap found

									pseudo_chr.append(target_str.substring(b2));
									
									if(qseqid.endsWith("'")) {
										// reverse complement
										qseqid = qseqid.replaceAll("'$", "");
										anchored_seq.add(qseqid);
										agpmap_str.add(qseqid+"\t"+(target_len-b2)+"\t"+0+"\t"+(target_len-b2)+"\t-\t"+sub_seq+"\t"+chrL+"\t"+pseudo_chr.length());
									} else {
										anchored_seq.add(qseqid);
										agpmap_str.add(qseqid+"\t"+(target_len-b2)+"\t"+b2+"\t"+target_len+"\t+\t"+sub_seq+"\t"+chrL+"\t"+pseudo_chr.length());
									}
									
								} else {
									// no overlap found
									gap = Math.max(min_gap, target_sstart-source_send);
									pseudo_chr.append(Sequence.polyN(gap));
									agpmap_str.add("GAP\t"+gap+"\t"+0+"\t"+gap+"\t+\t"+sub_seq+"\t"+chrL+"\t"+pseudo_chr.length());
									chrL = pseudo_chr.length();
									
									pseudo_chr.append(target_str);
									if(qseqid.endsWith("'")) {
										// reverse complement
										qseqid = qseqid.replaceAll("'$", "");
										anchored_seq.add(qseqid);
										agpmap_str.add(qseqid+"\t"+target_len+"\t"+0+"\t"+target_len+"\t-\t"+sub_seq+"\t"+chrL+"\t"+pseudo_chr.length());
									} else {
										anchored_seq.add(qseqid);
										agpmap_str.add(qseqid+"\t"+target_len+"\t"+0+"\t"+target_len+"\t+\t"+sub_seq+"\t"+chrL+"\t"+pseudo_chr.length());
									}
								}
								chrL = pseudo_chr.length();
							}
							
							myLogger.info("####consensus "+this.sub_seq+" finished.");
							
							synchronized(lock) {
								bw_fa.write(Sequence.formatOutput(this.sub_seq, pseudo_chr.toString()));
								for(final String agpstr : agpmap_str) bw_map.write(agpstr+"\n");
								anchored_contigs.addAll(anchored_seq);
							}
							***/
							
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

			/***
			bw_fa.close();
			bw_map.close();
			***/
			
			final BufferedWriter bw_ufa = Utils.getBufferedWriter(out_prefix+"_unplaced.fa");
			for(String seq : qry_seqs.keySet()) 
				if(!anchored_contigs.contains(seq)&&!seq.endsWith("'"))
					bw_ufa.write(qry_seqs.get(seq).formatOutput());
			bw_ufa.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
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

				// keep alignment fragment that has qual>0
				Set<SAMSegment> init_f = new HashSet<SAMSegment>();
				Set<SAMSegment> init_r = new HashSet<SAMSegment>();
				for(SAMSegment record : buff) {
					if(record.qlength()<min_len) 
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

		String source, target;
		Set<SAMSegment> sourcePlacement, targetPlacement;
		boolean removal;
		String pairkey;
		int slen, olap;
		final Set<OverlapEdge> edgeToRemove = new HashSet<OverlapEdge>();

		for(final OverlapEdge edge : gfa.edgeSet()) {
			source = gfa.getEdgeSource(edge);
			target = gfa.getEdgeTarget(edge);
			sourcePlacement = initPlace.containsKey(source)?initPlace.get(source):null;
			targetPlacement = initPlace.containsKey(target)?initPlace.get(target):null;
			pairkey = source+"->"+target;
			if(sourcePlacement==null||targetPlacement==null) {
				if(ddebug) myLogger.info(pairkey+": "+true+ " | "+ 
						(sourcePlacement==null?source+"=null":"") +" "+
						(targetPlacement==null?target+"=null":"") );
				edgeToRemove.add(edge);
				continue;
			}
			removal = true;
			olap = (int) edge.olap();
			slen = qry_seqs.get(source).seq_ln();

			String a = null, b = null;
			int d = -1;
			for(final SAMSegment s : sourcePlacement) {
				if(s.qlength()-Math.max(s.qend()-(slen-olap), 0)<min_len)
					continue;

				for(final SAMSegment t : targetPlacement) {
					if(t.qlength()-Math.max(olap-t.qstart(), 0)<min_len)
						continue;

					if(s.sseqid().equals(t.sseqid()) &&
							(d=AlignmentSegment.sdistance(s, t))<=max_dist) {
						a = s.sseqid()+"_"+s.sstart()+"-"+s.send();
						b = t.sseqid()+"_"+t.sstart()+"-"+t.send();
						removal = false;
					}
				}
			}
			if(removal) edgeToRemove.add(edge);
			if(ddebug) myLogger.info(pairkey+": "+removal+ " | "+d+" | "+a+","+b);
		}

		gfa.removeAllEdges(edgeToRemove);
		gfa.writeGFA(this.out_prefix+".gfa2");

		final Set<String> linkPlace  = new HashSet<String>();
		final Set<String> contigging = new HashSet<String>();
		final List<Contig> linkSeqStr = new ArrayList<Contig>();

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

						final List<SAMSegment> seqBySubAll = initPseudoAssembly.get(sub_seq);


						final int sub_ln = sub_seqs.get(sub_seq).seq_ln();
						// this is to calculate the coverage on the subject sequence
						// we are going to skip certain regions if the coverage is too high
						// say the average coverage is greater than max_cov
						int[] sub_cvg = new int[sub_ln]; 

						// we calculate the coverage across the reference chromosome
						for(SAMSegment sams : seqBySubAll) {
							int a = sams.sstart()-1, b = sams.send();
							for(int w = a; w<b; w++) sub_cvg[w]++;
						}

						int lowCvg = 0;
						for(int w=0; w<sub_ln; w++) {
							if(sub_cvg[w]<=max_cov) ++lowCvg; 
						}

						// we filter out high-coverage/highly-repetitive regions
						// >max_cov x
						final List<SAMSegment> seqBySubLowCov = new ArrayList<SAMSegment>();
						final Set<String> seqPool = new HashSet<String>();
						for(SAMSegment sams : seqBySubAll) {
							int a = sams.sstart()-1, b = sams.send();
							double cov = 0d;
							for(int w=a; w<b; w++) cov += sub_cvg[w];
							if(cov/(b-a)<=max_cov) {
								seqBySubLowCov.add(sams);
								seqPool.add(sams.qseqid());
							}
						}

						myLogger.info(sub_seq+" highly-repetitive regions: "+(sub_ln-lowCvg)+"/"+sub_ln+"bp,"+
								(seqBySubAll.size()-seqBySubLowCov.size())+"/"+seqBySubAll.size()+
								" alignment records filtered out due to high coverage(>"+max_cov+")");

						final Set<SAMSegment> contained = new HashSet<SAMSegment>();
						final Set<SAMSegment> placed    = new HashSet<SAMSegment>();

						double edge_penalty, edge_score;
						SAMSegment root_seq, source_seq, target_seq;
						Set<SAMSegment> target_seqs;
						Set<OverlapEdge> outgoing;
						TraceableEdge edge;
						String root_seqid, source_seqid, target_seqid;
						TraceableVertex<String> root_vertex, source_vertex, target_vertex;
						Deque<SAMSegment> deque = new ArrayDeque<SAMSegment>();
						final List<TraceableVertex<String>> traceable = new ArrayList<TraceableVertex<String>>();

						int distance;

						Collections.sort(seqBySubLowCov, new AlignmentSegment.SubjectCoordinationComparator());
						int nSeq = seqBySubLowCov.size();
						for(int i=0; i<nSeq; i++) {

							root_seq = seqBySubLowCov.get(i);
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

									if(!seqPool.contains(target_seqid)) continue;

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
									if(distance<=hc_gap) {
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
							root_vertex.setSInterval(root_seq.sstart(), root_seq.send());

							bidiQ.put(0L, root_vertex);
							double max_ws = Double.NEGATIVE_INFINITY,
									source_penalty, target_penalty, source_score, target_score, 
									penalty, score, target_ws, source_ws, ws;
							int source_ln, cvg;
							Set<TraceableEdge> out_edges;
							TraceableVertex<String> opt_vertex = null;
							RangeSet<Integer> target_sinterval;
							SAMSegment target_samseg;
							Range<Integer> target_range;
							long sizeQ;
							boolean isLeaf;

							if(ddebug)
								for(TraceableEdge e : razor.edgeSet()) 
									myLogger.info(e.toString()+"("+razor.getEdgeSource(e).getSAMSegment().toString()+"|"+
											razor.getEdgeTarget(e).getSAMSegment().toString()+"|"+
											e.getScore()+"-"+e.getPenalty()+")");

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
									//**** this is not right because graph edges are immutable?
									//**** target_vertex = razor.getEdgeTarget(out);
									target_vertex = razv_map.get(razor.getEdgeTarget(out).getId());

									target_samseg = target_vertex.getSAMSegment();

									// in order to avoid recursive placement in repetitive regions 
									// we need the new contig to expand the contigging on the reference genome
									target_range = Range.closed(target_samseg.sstart(), 
											target_samseg.send()).canonical(DiscreteDomain.integers());
									target_sinterval = source_vertex.getSIntervalCopy();
									target_sinterval.add(target_range);
									cvg = countIntervalCoverage(target_sinterval);

									if(cvg<countIntervalCoverage(source_vertex.getSInterval())+min_ext ||
											cvg<=countIntervalCoverage(target_vertex.getSInterval()))
										// if no significant improvement on coverage
										// if target vertex has been visited
										// and the reference covered is greater
										continue;

									target_score = target_vertex.getScore();
									target_penalty = target_vertex.getPenalty();
									target_ws = target_score-target_penalty;

									edge_penalty = out.getPenalty();
									penalty = source_penalty+edge_penalty;
									edge_score = out.getScore();
									score = source_score+edge_score;
									ws = score-penalty;

									if( edge_penalty>hc_gap || 
											target_vertex.getStatus() && 
											(ws<=target_ws ||
											isLoopback(razor, source_vertex, target_vertex)) ) 
										continue;

									isLeaf = false;
									target_vertex.setBackTrace(source_vertex);
									target_vertex.setScore(score);
									target_vertex.setPenalty(penalty);
									target_vertex.setStatus(true);
									target_vertex.setSInterval(target_sinterval);

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
							myLogger.info("<<<<<<<<<<<<<"+sub_seq+">>>>>>>>>>>>>>>>");

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
								myLogger.info("final trace back ["+score+", "+penalty+"]: "+trace);
							}
						}

						final StringBuilder linkSeq = new StringBuilder();
						final StringBuilder linkContigging = new StringBuilder();
						String contigStr;
						int olap;

						synchronized(lock) {
							for(TraceableVertex<String> opt_vertex : traceable) {

								List<TraceableVertex<String>> opt_vertices = new ArrayList<TraceableVertex<String>>();
								opt_vertices.add(opt_vertex);
								while( (opt_vertex = opt_vertex.getBackTrace())!=null ) opt_vertices.add(opt_vertex);

								if(opt_vertices.size()==1)
									// this is a singleton
									continue;

								Collections.reverse(opt_vertices);

								opt_vertex = opt_vertices.get(0);
								source_seqid = opt_vertex.getId();
								linkSeq.setLength(0);
								linkSeq.append(qry_seqs.get(source_seqid).seq_str());
								linkContigging.setLength(0);
								linkContigging.append(source_seqid);
								linkPlace.add(source_seqid);

								for(int k=1; k<opt_vertices.size(); k++) {

									opt_vertex = opt_vertices.get(k);
									target_seqid = opt_vertex.getId();

									OverlapEdge e = gfa.getEdge(source_seqid, target_seqid);

									if(!e.isRealigned()) gfa.realign(e);

									olap = (int)e.olapR();
									if(olap<0) {
										linkSeq.append(Constants.scaffold_gap_fill);
										linkSeq.append(qry_seqs.get(target_seqid).seq_str());
									} else {
										linkSeq.append(qry_seqs.get(target_seqid).seq_str().substring(olap));
									}
									linkContigging.append("->"+target_seqid);
									linkPlace.add(target_seqid);
									source_seqid = target_seqid;
								}

								contigStr = linkContigging.toString();
								if(contigging.contains(contigStr)) continue;
								linkSeqStr.add(new Contig(linkSeq.toString(),contigStr));
								contigging.add(contigStr);
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

		try {
			BufferedWriter bw_fa = Utils.getBufferedWriter(this.out_prefix+".fa2");
			BufferedWriter bw_map = Utils.getBufferedWriter(this.out_prefix+".map2");

			Collections.sort(linkSeqStr, new Comparator<Contig>() {

				@Override
				public int compare(Contig c1, Contig c2) {
					// TODO Auto-generated method stub
					return c2.seqStr.length()-c1.seqStr.length();
				}
			});

			int scaf = 1;
			for(Contig seq : linkSeqStr) {
				bw_fa.write(Sequence.formatOutput("contig"+String.format("%08d", scaf), seq.seqStr, 80));
				bw_map.write(String.format("%16s", "contig"+String.format("%08d", scaf))+"\t"+seq.components+"\n");
				++scaf;
			}
			for(String seqid : qry_seqs.keySet()) {
				if(seqid.endsWith("'")) continue;
				if(!linkPlace.contains(seqid)) {
					bw_fa.write(qry_seqs.get(seqid).formatOutput(80));
					bw_map.write(String.format("%16s", seqid)+"\t"+seqid+"\n");
				}
			}
			bw_fa.close();
			bw_map.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
	}

	protected int countIntervalCoverage(RangeSet<Integer> interval) {
		// TODO Auto-generated method stub
		if(interval==null) return 0;
		int cvg = 0;
		for(Range<Integer> range : interval.asRanges()) 
			cvg += range.upperEndpoint()-range.lowerEndpoint()+1;
		return cvg;
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

	private class Contig {
		final String components;
		final String seqStr;

		public Contig(String seqStr, String components) {
			this.seqStr = seqStr;
			this.components = components;
		}
	}

	private void run1() {
		// TODO Auto-generated method stub

		// read assembly graph file
		final GFA gfa = new GFA(query_file, asm_graph);
		qry_seqs = gfa.getSequenceMap();
		sub_seqs = Sequence.parseFastaFileAsMap(subject_file);

		myLogger.info("  GFA vertices: "+gfa.vertexSet().size());
		myLogger.info("  GFA edges   : "+gfa.edgeSet().size()  );

		// read alignment file and place the query sequences
		final Map<String, Set<SAMSegment>> initPlace = new HashMap<String, Set<SAMSegment>>();
		final Map<String, TreeMap<Integer, Set<SAMSegment>>> initPseudoAssembly = 
				new HashMap<String, TreeMap<Integer, Set<SAMSegment>>>();
		for(String sub_seq : sub_seqs.keySet()) {
			if(!sub_seq.equals("Chr00"))
				initPseudoAssembly.put(sub_seq, new TreeMap<Integer, Set<SAMSegment>>());
		}

		try {
			final SamReaderFactory factory =
					SamReaderFactory.makeDefault()
					.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
							SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
					.validationStringency(ValidationStringency.SILENT);
			final SamReader in1 = factory.open(new File(align_file));
			final SAMRecordIterator iter1 = in1.iterator();

			String qry;
			int qry_ln, sstart;
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

				Set<SAMSegment> init_f = new HashSet<SAMSegment>();
				Set<SAMSegment> init_r = new HashSet<SAMSegment>();
				for(SAMSegment record : buff) {
					if("Chr00".equals(record.sseqid())) continue;
					if(record.qlength()<min_alen) 
						continue;
					if(record.qseqid().equals(qry)) 
						init_f.add(record);
					else init_r.add(record);
					sstart = record.sstart();
					if(!initPseudoAssembly.get(record.sseqid()).containsKey(sstart)) 
						initPseudoAssembly.get(record.sseqid()).put(sstart, new HashSet<SAMSegment>());
					initPseudoAssembly.get(record.sseqid()).get(sstart).add(record);
				}
				if(!init_f.isEmpty()) initPlace.put(qry,     init_f);
				if(!init_r.isEmpty()) initPlace.put(qry+"'", init_r);
			}
			iter1.close();
			in1.close();

			String source, target;
			Set<SAMSegment> sourcePlacement, targetPlacement;
			boolean removal;
			int slen, olap;
			final Set<OverlapEdge> edgeToRemove = new HashSet<OverlapEdge>();
			final Map<String, Map<String, Map<Integer, Integer>>> edgeSegPair = 
					new HashMap<String, Map<String, Map<Integer, Integer>>>();
			for(String sub_seq : sub_seqs.keySet()) {
				if(!sub_seq.equals("Chr00"))
					edgeSegPair.put(sub_seq, new HashMap<String, Map<Integer, Integer>>());
			}
			String pairkey, sseqid;
			Map<String, Map<Integer, Integer>> segPair;

			for(final OverlapEdge edge : gfa.edgeSet()) {
				source = gfa.getEdgeSource(edge);
				target = gfa.getEdgeTarget(edge);
				sourcePlacement = initPlace.containsKey(source)?initPlace.get(source):null;
				targetPlacement = initPlace.containsKey(target)?initPlace.get(target):null;
				pairkey = source+"->"+target;
				if(sourcePlacement==null||targetPlacement==null) {
					if(ddebug) myLogger.info(pairkey+": "+true+ " | "+ 
							(sourcePlacement==null?source+"=null":"") +" "+
							(targetPlacement==null?target+"=null":"") );
					edgeToRemove.add(edge);
					continue;
				}
				removal = true;
				olap = (int) edge.olap();
				slen = qry_seqs.get(source).seq_ln();

				String a = null, b = null;
				int d = -1;
				for(final SAMSegment s : sourcePlacement) {
					if(s.qlength()-Math.max(s.qend()-(slen-olap), 0)<min_len)
						continue;

					for(final SAMSegment t : targetPlacement) {
						if(t.qlength()-Math.max(olap-t.qstart(), 0)<min_len)
							continue;

						if(s.sseqid().equals(t.sseqid()) &&
								(d=AlignmentSegment.sdistance(s, t))<=max_dist) {

							sseqid = s.sseqid();
							segPair = edgeSegPair.get(sseqid);
							if(!segPair.containsKey(pairkey)) 
								segPair.put(pairkey, new HashMap<Integer, Integer>());
							segPair.get(pairkey).put(s.sstart(), t.sstart());
							a = s.sseqid()+"_"+s.sstart()+"-"+s.send();
							b = t.sseqid()+"_"+t.sstart()+"-"+t.send();
							removal = false;
						}
					}
				}
				if(removal) edgeToRemove.add(edge);
				if(ddebug) myLogger.info(pairkey+": "+removal+ " | "+d+" | "+a+","+b);
			}

			gfa.removeAllEdges(edgeToRemove);
			gfa.writeGFA(this.out_prefix+".gfa2");

			final List<GFAPath> paths = new ArrayList<GFAPath>();
			int rpos;
			Set<SAMSegment> sseg, tseg;
			Map<Integer, Integer> pairs;
			int source_sstart, target_sstart;
			SAMSegment source_seg, target_seg;

			for(final String refSeqStr : initPseudoAssembly.keySet()) {
				segPair = edgeSegPair.get(refSeqStr);
				final TreeMap<Integer, Set<SAMSegment>> initAssembly = initPseudoAssembly.get(refSeqStr);

				for(final Map.Entry<Integer, Set<SAMSegment>> entry : initAssembly.entrySet()) {
					rpos = entry.getKey();
					sseg = entry.getValue();

					for(final SAMSegment seg : sseg) {

						// create a path tree
						final LinkedList<SAMSegment> leaves  = new LinkedList<SAMSegment>();
						final DirectedPseudograph<SAMSegment, DefaultEdge> pathTree = 
								new DirectedPseudograph<SAMSegment, DefaultEdge>(DefaultEdge.class);
						pathTree.addVertex(seg);
						leaves.add(seg);

						while(!leaves.isEmpty()) {
							source_seg = leaves.poll();
							source = source_seg.qseqid();
							source_sstart = source_seg.sstart();

							for(OverlapEdge edge : gfa.outgoingEdgesOf(source)) {
								target = gfa.getEdgeTarget(edge);
								pairkey = source+"->"+target;
								if(!segPair.containsKey(pairkey)) continue;

								pairs = segPair.get(pairkey);
								if(!pairs.containsKey(source_sstart)) continue;

								target_sstart = pairs.get(source_sstart);

								target_seg = null;

								for(final SAMSegment s : initAssembly.get(target_sstart)) {
									if(s.qseqid().equals(target)) {
										target_seg = s;
										break;
									}
								}

								if(!pathTree.containsVertex(target_seg)) {
									pathTree.addVertex(target_seg);
									leaves.offer(target_seg);
								}

								pathTree.addEdge(source_seg, target_seg);
							}
						}

						// now find best path in this region
						if(pathTree.vertexSet().size()==1) continue;

						if(ddebug) {
							myLogger.info("#########################################");
							myLogger.info(seg.toString());
							myLogger.info("#vertices "+pathTree.vertexSet().size());
							myLogger.info("#edge     "+pathTree.edgeSet().size());
							for(DefaultEdge edge : pathTree.edgeSet()) 
								myLogger.info(pathTree.getEdgeSource(edge).qseqid()+"->"+pathTree.getEdgeTarget(edge).qseqid());
						}

					}
				}
			}

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	final private class GFAPath {
		final List<String> path;
		final int sstart;
		final int send;
		final int gap;

		public GFAPath(final List<String> path,
				final int sstart,
				final int send,
				final int gap) {
			this.path   = path;
			this.sstart = sstart;
			this.send   = send;
			this.gap    = gap;
		}
	}
}















