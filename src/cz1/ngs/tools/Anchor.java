package cz1.ngs.tools;

import java.io.BufferedReader;
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
import cz1.ngs.model.Blast6Segment;
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
	private double collinear_shift = 10; // maximum shift distance for two collinear alignment segments - 50% of the smaller segment size
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
			this.THREADS = t;
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
	private static List<SAMSegment> sam_records;

	private final static int match_score  = 1;
	private final static int clip_penalty = 1;
	private final static int hc_gap  = 10000;
	private final static int max_cov = 31;
	
	private final static Object lock = new Object();

	@Override
	public void run() {
		this.run1();
		// this.run2();
	}
	
	public void run2() {
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
		
		// parse blast records
		// blast records buffer
		final List<SAMSegment> buff = new ArrayList<SAMSegment>();
		// selected blast records
		final List<SAMSegment> sel_recs = new ArrayList<SAMSegment>();
		// temp list
		final List<SAMSegment> tmp_records = new ArrayList<SAMSegment>();
		// collinear merged record list
		// final List<Blast6Record> collinear_merged = new ArrayList<Blast6Record>();

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
			double qry_ln, aln_frac;

			while(rc!=null) {
				
				qry = rc.getReadName();
				qry_ln = qry_seqs.get(qry).seq_ln();			

				buff.clear();
				if(!rc.getReadUnmappedFlag())
					buff.add(SAMSegment.samRecord(rc));

				while( (rc=iter1.next())!=null
						&&
						rc.getReadName().equals(qry) ) {
					buff.add(SAMSegment.samRecord(rc));
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
					// more accurate but slower
					Collections.sort(tmp_records, new AlignmentSegment.SubjectCoordinationComparator());

					SAMSegment record;
					for(int i=0; i<tmp_records.size(); i++) {
						primary_record = tmp_records.get(i);
						for(int j=i+1; j<tmp_records.size(); j++) {
							secondary_record = tmp_records.get(j);
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

				// filter by alignment fraction		
				buff.clear();
				buff.addAll(sel_recs);
				sel_recs.clear();
				for(SAMSegment record : buff) {
					if(record.qlength()/qry_ln>=this.min_frac)
						sel_recs.add(record);
				}

				if(sel_recs.isEmpty()) {
					// unplaced query sequences
					// continue
					continue;
				}
				// process primary alignment
				primary_record = sel_recs.get(0);
				anchored_records.get(primary_record.sseqid()).add(primary_record);
				aln_frac = primary_record.qlength()/qry_ln;
				// compare secondary alignments to primary alignment
				// and process
				for(int i=1; i<sel_recs.size(); i++) {
					secondary_record = sel_recs.get(i);
					if(secondary_record.qlength()/qry_ln+this.diff_frac<aln_frac) {
						break;
					} else {
						anchored_records.get(secondary_record.sseqid()).add(secondary_record);
					}
				}
			}

			iter1.close();
			in1.close();
			
			for(Map.Entry<String, List<SAMSegment>> entry : anchored_records.entrySet()) {
				System.out.println(entry.getKey()+": "+entry.getValue().size());
			}

			final BufferedWriter bw_map = Utils.getBufferedWriter(out_prefix+".map");
			final BufferedWriter bw_fa = Utils.getBufferedWriter(out_prefix+".fa");
			final Set<String> anchored_seqs = new HashSet<String>();
			final List<String> sub_list = Sequence.parseSeqList(subject_file);

			for(String sub_sn : sub_list) {

				sam_records = anchored_records.get(sub_sn);
				int nV = sam_records.size(), count = 0;

				// sort blast records
				Collections.sort(sam_records, new SAMSegment.SubjectCoordinationComparator());
				// consensus
				int posUpto = 0, send_clip = 0;
				int sstart, send, qstart, qend, qlen, tmp_int, qstart_clip, qend_clip, gap_size;
				// distance to last 'N', start from next position to find longest common suffix-prefix
				int prev_n = Integer.MAX_VALUE;
				int mol_len = 0;
				String qseq;
				int nS, nQ;

				// first step: construct super scaffold
				// will form a tree graph indicating the path through the contigs

				// convert contig names to integer indices
				// one contig could end up with multiple indices due to repeats
				Map<Integer, String> ss_coordinate = new HashMap<Integer, String>();
				int index = 0;
				for(SAMSegment record : sam_records)
					ss_coordinate.put(index++, record.qseqid()+(record.forward()?"":"'"));

				StringBuilder seq_str = new StringBuilder();
				for(int v=0; v<nV-1; v++) {
					if(++count%10000==0) myLogger.info(sub_sn+" "+count+"/"+nV+" done.");

					SAMSegment record = sam_records.get(v);
					qlen = qry_seqs.get(record.qseqid()).seq_ln();
					sstart = record.sstart();
					send = record.send();
					qstart = record.qstart();
					qend = record.qend();
					if(sstart>send) {
						// make sure sstart<send
						tmp_int = sstart;
						sstart = send;
						send = tmp_int;
						tmp_int = qstart;
						qstart = qend;
						qend = tmp_int;
					}

					if(qstart>qend) {
						qstart_clip = qlen-qstart;
						qend_clip = qend-1;
						qseq = Sequence.revCompSeq(qry_seqs.get(record.qseqid()).seq_str());
					} else {
						qstart_clip = qstart-1;
						qend_clip = qlen-qend;
						qseq = qry_seqs.get(record.qseqid()).seq_str();
					}

					if(send<posUpto||sstart<posUpto&&qstart_clip>max_clip) {
						// skip if it is redundant
						//     ====================
						//      /-------\
						//        /----\
						// skip if contradiction observed
						//     ====================
						//      /-------\
						//       --/----\--
						// TODO process
						continue;
					}

					// find longest suffix-prefix
					nS = seq_str.length();
					nQ = qseq.length();
					int nO = Math.min(prev_n, Math.min(nS, nQ));
					outerloop:
						for(; nO>=min_overlap; nO--) {
							int nS_i = nS-nO;
							for(int i=0; i<nO; i++) {
								if(seq_str.charAt(nS_i+i)!=qseq.charAt(i))
									continue outerloop;
							}
							break outerloop;
						}

					if(nO<min_overlap) {
						// no overlap found

						// simply extend
						//     ====================
						//      /-------\
						//               /----\

						// will insert a GAP anyway
						// if sstart<=posUpto will insert a small GAP min_gap
						// otherwise will insert a large GAP max(pseduo_distance, max_gap)
						if(posUpto>0) {
							if(sstart<=posUpto) {
								// if too much overlap, then treat it as a contradiction
								if(posUpto-sstart>min_overlap) {
									//discard
									continue;
								} else {
									// insert a min_gap
									gap_size = min_gap;
								}
							} else {
								// estimate gap size
								gap_size = (sstart-posUpto)-(send_clip+qstart_clip);
								if(gap_size<max_gap) gap_size = max_gap;	
							}

							bw_map.write("GAP\t"+
									gap_size+
									"\t0\t"+
									gap_size+
									"\t+\t"+
									sub_sn+
									"\t"+
									mol_len+
									"\t"+
									(mol_len+gap_size)+
									"\n");
							seq_str.append( Sequence.polyN(gap_size) );
							mol_len += gap_size;
						}
						bw_map.write(record.qseqid()+
								"\t"+
								qlen+
								"\t");
						if(qstart>qend) {
							// reverse
							bw_map.write("0\t"+qlen+"\t-\t");
						} else {
							// forward
							bw_map.write("0\t"+qlen+"\t+\t");
						}
						seq_str.append(qseq);
						bw_map.write(sub_sn+
								"\t"+
								mol_len+
								"\t"+
								(mol_len+qlen)+
								"\n");
						mol_len += qlen;
						prev_n = qlen;

						anchored_seqs.add(record.qseqid());

					} else {
						// overlap found
						// will not insert gap
						//     ====================
						//      /-------\
						//            /----\
						// calculate overlaps
						// process overlap

						qstart = nO;
						if(qstart==qlen) continue;
						bw_map.write(record.qseqid()+
								"\t"+
								(qlen-qstart)+
								"\t");

						if(qstart>qend) {
							// reverse
							bw_map.write( 0+
									"\t"+
									(qlen-qstart)+
									"\t-\t" );
						} else {
							// forward
							bw_map.write( qstart+
									"\t"+
									qlen+
									"\t+\t" );
						}
						bw_map.write(sub_sn+
								"\t"+
								mol_len+
								"\t"+
								(mol_len+qlen-qstart)+
								"\n");
						mol_len += qlen-qstart;
						prev_n += qlen-qstart;
						seq_str.append( qseq.substring(qstart) );

						anchored_seqs.add(record.qseqid());
					}

					posUpto = send;
					send_clip = qend_clip;

				}

				if(seq_str.length()>0) 
					bw_fa.write(Sequence.formatOutput(sub_sn, seq_str.toString()));
			}

			bw_fa.close();
			bw_map.close();

			final BufferedWriter bw_ufa = Utils.getBufferedWriter(out_prefix+"_unplaced.fa");
			for(String seq : qry_seqs.keySet()) 
				if(!anchored_seqs.contains(seq))
					bw_ufa.write(qry_seqs.get(seq).formatOutput());
			bw_ufa.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	public void run1() {
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

		final Set<String> linkPlace  = new HashSet<String>();
		final List<String> linkSeqStr = new ArrayList<String>();
		
		this.initial_thread_pool();
		for(String sub_seq : sub_seqs.keySet()) {

			this.executor.submit(
					new Runnable() {

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
								for(SAMSegment sams : seqBySubAll) {
									int a = sams.sstart()-1, b = sams.send();
									double cov = 0d;
									for(int w=a; w<b; w++) cov += sub_cvg[w];
									if(cov/(b-a)<=max_cov) seqBySubLowCov.add(sams);
								}
								
								myLogger.info(sub_seq+" highly-repetitive regions: "+(sub_ln-lowCvg)+"/"+sub_ln+"bp,"+
										(seqBySubAll.size()-seqBySubLowCov.size())+"/"+seqBySubAll.size()+
										"alignment records filtered out due to high coverage(>"+max_cov+")");
								
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

								synchronized(lock) {
									for(TraceableVertex<String> opt_vertex : traceable) {
										
										List<TraceableVertex<String>> opt_vertices = new ArrayList<TraceableVertex<String>>();
										opt_vertices.add(opt_vertex);
										while( (opt_vertex = opt_vertex.getBackTrace())!=null ) opt_vertices.add(opt_vertex);
										Collections.reverse(opt_vertices);
										
										opt_vertex = opt_vertices.get(0);
										source_seqid = opt_vertex.getId();
										linkSeq.setLength(0);
										linkSeq.append(qry_seqs.get(source_seqid).seq_str());
										linkPlace.add(source_seqid);
										
										for(int k=1; k<opt_vertices.size(); k++) {
											
											opt_vertex = opt_vertices.get(k);
											target_seqid = opt_vertex.getId();
												
											OverlapEdge e = gfa.getEdge(source_seqid, target_seqid);
											
											if(!e.olapInfo().isRealigned()) gfa.realign(e);
											
											linkSeq.append(qry_seqs.get(target_seqid).seq_str().substring((int)e.olapR()));
											linkPlace.add(target_seqid);
											source_seqid = target_seqid;
										}
										linkSeqStr.add(linkSeq.toString());
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
			BufferedWriter bw = Utils.getBufferedWriter(this.out_prefix+".fa2");
			
			Collections.sort(linkSeqStr, new Comparator<String>() {

				@Override
				public int compare(String s1, String s2) {
					// TODO Auto-generated method stub
					return s2.length()-s1.length();
				}
			});
			
			int scaf = 1;
			for(String seq : linkSeqStr) 
				bw.write(Sequence.formatOutput("Scaffold"+String.format("%08d", scaf++), seq, 80));
			for(String seqid : qry_seqs.keySet()) {
				if(seqid.endsWith("'")) continue;
				if(!linkPlace.contains(seqid))
					bw.write(qry_seqs.get(seqid).formatOutput(80));
			}
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
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















