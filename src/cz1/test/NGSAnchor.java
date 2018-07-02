package cz1.test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;

import com.google.common.collect.DiscreteDomain;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

import cz1.ngs.model.AlignmentSegment;
import cz1.ngs.model.SAMSegment;
import cz1.ngs.model.Sequence;
import cz1.ngs.model.TraceableAlignmentSegment;
import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Executor;
import cz1.util.Utils;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class NGSAnchor extends Executor {
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
	private final static int min_olap  = 50;
	private final static int max_clip2 = 30;

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
		//this.run_stringent();
		this.run_greedy();
	}

	public void run_greedy() {
		// this will trim the original contigs as long as the are significant overlap found
		
		sub_seqs = Sequence.parseFastaFileAsMap(subject_file);
		qry_seqs = Sequence.parseFastaFileWithRevCmpAsMap(query_file);
		
		try {
			final BufferedWriter bw_fa = Utils.getBufferedWriter(out_prefix+".fa");

			this.initial_thread_pool();
			for(String sub_seq : sub_seqs.keySet()) {

				this.executor.submit(new Runnable() {

					String sub_seq = null;

					@Override
					public void run() {
						// TODO Auto-generated method stub
						try {
							if(sub_seq.equals("Chr00")) return;

							final List<AlignmentSegment> graph_path = new ArrayList<AlignmentSegment>();
							BufferedReader br = Utils.getBufferedReader(asm_graph);
							myLogger.info(br.readLine());
							String line;
							String[] s;
							while( (line=br.readLine())!=null ) {
								s = line.split("\\s+");
								if(s[1].equals(this.sub_seq))
									graph_path.add(new AlignmentSegment(s[0],s[1],Integer.parseInt(s[2]),
											Integer.parseInt(s[3]),Integer.parseInt(s[4]),Integer.parseInt(s[5])));
							}
							br.close();
							
							if(graph_path.isEmpty()) {
								myLogger.info(this.sub_seq+" done");
								return;
							}
							
							// now we join the neighbouring placement 
							final StringBuilder pseudo_chr = new StringBuilder();
							AlignmentSegment alignment = graph_path.get(0);
							
							String qseqid = alignment.qseqid();
							pseudo_chr.append(qry_seqs.get(qseqid).seq_str());
							
							String target_str;
							int olap, gap, source_trim, target_trim, source_len;
							AlignmentSegment source_as, target_as;
							int source_qend, target_qstart, source_send, target_sstart;

							for(int i=1; i<graph_path.size(); i++) {
								
								source_as = graph_path.get(i-1);
								target_as = graph_path.get(i);
								
								source_qend   = source_as.qend();
								source_send   = source_as.send();
								source_len    = qry_seqs.get(source_as.qseqid()).seq_ln();

								target_qstart = target_as.qstart();
								target_sstart = target_as.sstart();

								if(target_as.send()<=source_send) {
									graph_path.remove(i);
									--i;
									continue;
								}
								
								qseqid = target_as.qseqid();
								target_str = qry_seqs.get(qseqid).seq_str();

								// see if there is overlap between source tail and target head
								olap = source_send-target_sstart+1;
								if(olap<min_olap) {
									// so no overlap
									// estimate gap and directly join them
									gap = Math.max(min_gap, target_sstart-source_send);
									pseudo_chr.append(Sequence.polyN(gap));
									pseudo_chr.append(target_str);
									
								} else {
									// so overlap found
									// trim the both contigs if necessary and join them
									// trim source str first
									source_trim = source_len-source_qend;
									// trim target str then
									target_trim = target_qstart-1+olap;
									
									if(target_trim>=target_str.length()) {
										graph_path.remove(i);
										--i;
										continue;
									}
									
									pseudo_chr.setLength(pseudo_chr.length()-source_trim);
									pseudo_chr.append(target_str.substring(target_trim));
									myLogger.info(source_as.qseqid()+" end trimmed "+source_trim+", "+
											target_as.qseqid()+" start trimmed "+target_trim);
								}
							}

							synchronized(lock) {
								bw_fa.write(Sequence.formatOutput(this.sub_seq, pseudo_chr.toString()));
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

			bw_fa.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public void run_stringent() {
		// this will not trim too much the original contigs
		
		sub_seqs = Sequence.parseFastaFileAsMap(subject_file);
		qry_seqs = Sequence.parseFastaFileWithRevCmpAsMap(query_file);
		
		try {
			final BufferedWriter bw_map = Utils.getBufferedWriter(out_prefix+".map");
			final BufferedWriter bw_fa = Utils.getBufferedWriter(out_prefix+".fa");

			this.initial_thread_pool();
			for(String sub_seq : sub_seqs.keySet()) {

				this.executor.submit(new Runnable() {

					String sub_seq = null;

					@Override
					public void run() {
						// TODO Auto-generated method stub
						try {
							if(sub_seq.equals("Chr00")) return;

							final List<AlignmentSegment> graph_path = new ArrayList<AlignmentSegment>();
							BufferedReader br = Utils.getBufferedReader(asm_graph);
							myLogger.info(br.readLine());
							String line;
							String[] s;
							while( (line=br.readLine())!=null ) {
								s = line.split("\\s+");
								if(s[1].equals(this.sub_seq))
									graph_path.add(new AlignmentSegment(s[0],s[1],Integer.parseInt(s[2]),
											Integer.parseInt(s[3]),Integer.parseInt(s[4]),Integer.parseInt(s[5])));
							}
							br.close();
							
							if(graph_path.isEmpty()) {
								myLogger.info(this.sub_seq+" done");
								return;
							}
							
							// now we join the neighbouring placement 
							final Set<String> anchored_seq = new HashSet<String>();
							final StringBuilder pseudo_chr = new StringBuilder();
							final List<String> agpmap_str  = new ArrayList<String>();
							AlignmentSegment alignment = graph_path.get(0);
							myLogger.info(alignment.toString());
							
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
							AlignmentSegment source_as, target_as;
							int source_send, target_sstart;

							for(int i=1; i<graph_path.size(); i++) {
								
								source_as = graph_path.get(i-1);
								target_as = graph_path.get(i);

								myLogger.info(target_as.toString());
								
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

							if(ddebug) {
								myLogger.info("####agp map");
								for(final String agpstr : agpmap_str) myLogger.info(agpstr);
							}

							synchronized(lock) {
								bw_fa.write(Sequence.formatOutput(this.sub_seq, pseudo_chr.toString()));
								for(final String agpstr : agpmap_str) bw_map.write(agpstr+"\n");
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

			bw_fa.close();
			bw_map.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		NGSAnchor anchor = new NGSAnchor();
		anchor.setParameters(args);
		anchor.run();
	}
}
