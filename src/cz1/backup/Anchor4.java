package cz1.backup;

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

public class Anchor4 extends Executor {
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
			myArgsEngine.add("-g", "--graph", true);
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
		sub_seqs = Sequence.parseFastaFileAsMap(subject_file);
		qry_seqs = Sequence.parseFastaFileWithRevCmpAsMap(query_file);

		final Map<String, List<SAMSegment>> anchored_recs = new HashMap<String, List<SAMSegment>>();
		for(String seq_sn : sub_seqs.keySet()) anchored_recs.put(seq_sn, new ArrayList<SAMSegment>());
		
		// parse SAMRecords
		// SAMRecords buffer
		final List<SAMSegment> buff = new ArrayList<SAMSegment>();
		// selected SAMRecords
		final List<SAMSegment> sel_recs = new ArrayList<SAMSegment>();
		// temp list
		final List<SAMSegment> tmp_recs = new ArrayList<SAMSegment>();

		try {
			final SamReaderFactory factory =
					SamReaderFactory.makeDefault()
					.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
							SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
					.validationStringency(ValidationStringency.SILENT);
			final SamReader in1 = factory.open(new File(align_file));
			final SAMRecordIterator iter1 = in1.iterator();
			// find clipping sites
			
			iter1.close();
			in1.close();
			
			
			
			
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
					tmp_recs.clear();
					for(SAMSegment record : buff)
						if(record.sseqid().equals(sub))
							tmp_recs.add(record);
					if(tmp_recs.isEmpty()) continue;

					// find collinear alignment segments that can be merged
					Collections.sort(tmp_recs, new AlignmentSegment.SubjectCoordinationComparator());

					SAMSegment record;
					int pqstart, pqend, sqstart, sqend, olap;
					for(int i=0; i<tmp_recs.size(); i++) {
						primary_record = tmp_recs.get(i);
						pqstart = primary_record.qstart();
						pqend   = primary_record.qend();

						for(int j=i+1; j<tmp_recs.size(); j++) {
							secondary_record = tmp_recs.get(j);
							// we need to check if they are two copies that are close to each other
							// TODO: how?
							// we check the region of alignments on the query sequence
							// if the overlap is greater than 90%, we regard them as different copies
							sqstart = secondary_record.qstart();
							sqend   = secondary_record.qend();

							olap = 0;
							if(pqstart<=sqstart) olap = Math.max(0, pqend-sqstart);
							if(pqstart> sqstart) olap = Math.max(0, sqend-pqstart);

							if( (double)olap/(pqend-pqstart+1)>=0.9||(double)olap/(sqend-sqstart+1)>=0.9 ) continue;

							// now we check if they two are collinear
							double max_shift = collinear_shift*
									Math.min(primary_record.qlength(), secondary_record.qlength());
							if( (record=SAMSegment.collinear(primary_record, secondary_record, max_shift))!=null ) {
								tmp_recs.set(i, record);
								tmp_recs.remove(j);
								--i;
								break;
							}
						}
					}
				}
			}

			

			for(Map.Entry<String, List<SAMSegment>> entry : anchored_recs.entrySet()) {
				myLogger.info("initial placement #entry "+entry.getKey()+": "+entry.getValue().size());
				if(this.ddebug) {
					List<SAMSegment> segs = entry.getValue();
					Collections.sort(segs, new SAMSegment.SubjectCoordinationComparator());
					for(SAMSegment seg : segs)
						myLogger.info(seg.toString());
				}
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
