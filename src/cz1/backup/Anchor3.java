package cz1.backup;

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

import com.google.common.collect.DiscreteDomain;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

import cz1.ngs.model.Blast6Segment;
import cz1.ngs.model.GFA;
import cz1.ngs.model.OverlapEdge;
import cz1.ngs.model.AlignmentSegment;
import cz1.ngs.model.Sequence;
import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

public class Anchor3 extends Executor {
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -s/--subject            Subject/reference sequences file in FASTA format.\n"
						+ " -q/--query              Query sequences to anchor in FASTA format.\n"
						+ " -g/--graph              Assembly graph file (GFA format).\n"
						+ " -b/--blast              BLAST output (-outfmt 6) of query sequences against subject sequences.\n"
						+ " -i/--min-identity       Minimum identity between the query and subject sequences \n"
						+ "                         for an alignment record to consider (default 0.9).\n"
						+ " -c/--min-coverage       Minimum alignment coverage of the query sequence (default 0.95).\n"
						+ " -di/--diff-identity     Threshold of the difference of the identity between the primary and secondary \n"
						+ "                         alignments. If the difference is smaller than this value, the query \n"
						+ "                         sequence will be regarded as duplications. Otherwise, the secondary \n"
						+ "                         alignments will be discared (default 0.05).\n"
						+ " -dc/--diff-coverage     Threshold of the difference of the alignment coverage between the primary and \n"
						+ "                         secondary alignments. If the difference is smaller than this value, the query \n"
						+ "                         sequence will be regarded as duplications. Otherwise, the secondary \n"
						+ "                         alignments will be discared (default 0.05).\n"
						+ " -o/--out-prefix         Prefix of the output FASTQ file.\n"
						+ "\n");	
	}

	private String subject_file = null;
	private String query_file = null;
	private String blast_out = null;
	private double min_ident = 90;
	private double min_cov = 0.90;
	private double diff_ident = 5;
	private double diff_cov = 0.05;
	private String asm_graph; // assembly graph (GFA) format
	// maximum shift distance for two collinear alignment segments
	// 50% of the smaller segment size
	private String out_prefix = null;
	
	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-s", "--subject", true);
			myArgsEngine.add("-q", "--query", true);
			myArgsEngine.add("-g", "--graph", true);
			myArgsEngine.add("-b", "--blast", true);
			myArgsEngine.add("-i", "--min-identity", true);
			myArgsEngine.add("-c", "--min-fraction", true);
			myArgsEngine.add("-di", "--diff-identity", true);
			myArgsEngine.add("-dc", "--diff-coverage", true);
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
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the query file.");
		}
		
		if (myArgsEngine.getBoolean("-g")) {
			this.asm_graph = myArgsEngine.getString("-g");
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify the graph file.");
		}
		
		if (myArgsEngine.getBoolean("-b")) {
			this.blast_out = myArgsEngine.getString("-b");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the BLAST file.");
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
		
		if (myArgsEngine.getBoolean("-c")) {
			this.min_cov = Double.parseDouble(myArgsEngine.getString("-c"));
		}
		
		if (myArgsEngine.getBoolean("-di")) {
			this.diff_ident = 100*Double.parseDouble(myArgsEngine.getString("-di"));
		}
		
		if (myArgsEngine.getBoolean("-dc")) {
			this.diff_cov = Double.parseDouble(myArgsEngine.getString("-dc"));
		}
	}
	
	private final static int gap_buff = 30; // buffer size for subject/reference sequences gap clips
	private final static int gap_min  = 10000;
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		final Map<String, Sequence> sub_seqs = Sequence.parseFastaFileAsMap(subject_file);
		final Map<String, Sequence> qry_seqs = Sequence.parseFastaFileWithRevCmpAsMap(query_file);
		
		//final GFA gfa = new GFA(this.query_file, this.asm_graph);
		myLogger.info("Loading assembly graph done.");
		
		// find 'N/n's in subject/reference sequences
		// which could have impact on parsing the blast records
		final Map<String, TreeRangeSet<Integer>> sub_gaps = new HashMap<String, TreeRangeSet<Integer>>();
		final Map<String, List<Blast6Segment>> anchored_records = new HashMap<String, List<Blast6Segment>>();
		/**
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
			anchored_records.put(seq_sn, new ArrayList<Blast6Segment>());
		}
		**/

		// parse blast records
		// blast records buffer
		final List<Blast6Segment> buff = new ArrayList<Blast6Segment>();
		// selected blast records
		final List<Blast6Segment> sel_recs = new ArrayList<Blast6Segment>();
		// temp list
		final List<Blast6Segment> tmp_records = new ArrayList<Blast6Segment>();
		// collinear merged record list
		// final List<Blast6Record> collinear_merged = new ArrayList<Blast6Record>();
		
		try {
			/***
			BufferedWriter bw_blast = Utils.getBufferedWriter(this.out_prefix+".out6");
			BufferedReader br_blast = Utils.getBufferedReader(blast_out);
			
			Blast6Segment tmp_record = Blast6Segment.blast6Segment(br_blast.readLine(), true);
			Blast6Segment primary_record, secondary_record;
			String qry;
			double qry_ln, cov, max_cov;
			long rec_count = 1;
			/***
			while(tmp_record!=null) {
				qry = tmp_record.qseqid();
				if(qry.endsWith("'")) qry = qry.substring(0,qry.length()-1);
				
				qry_ln = qry_seqs.get(qry).seq_ln();
				buff.clear();
				buff.add(tmp_record);
				max_cov = Math.min(1.0, tmp_record.length()/qry_ln);
				while( (tmp_record = Blast6Segment.blast6Segment(br_blast.readLine(), true))!=null
						&&
						tmp_record.qseqid().startsWith(qry) ) {
					++rec_count;
					if(rec_count%1000000==0) myLogger.info(rec_count/1000000+" millon records loaded.");
					// filter by identity
					buff.add(tmp_record);
					cov = Math.min(1.0, tmp_record.length()/qry_ln);
					if(cov>max_cov) max_cov = cov;
				}
				if(max_cov<min_cov) continue;
				for(Blast6Segment seg : buff) 
					if( (cov=Math.min(1.0, seg.length()/qry_ln)+diff_cov)>=max_cov)
						//sel_recs.add(seg);
						bw_blast.write(seg.toString()+"\t"+Math.min(1.0,cov)+"\n");
			}
			br_blast.close();
			bw_blast.close();
			**/
			
			Blast6Segment tmp_record;
			BufferedReader br_blast = Utils.getBufferedReader(blast_out);
			final Map<String, TreeRangeSet<Integer>> sub_cov = new HashMap<>();
			for(final String sub : sub_seqs.keySet()) sub_cov.put(sub, TreeRangeSet.create());
			while( (tmp_record = Blast6Segment.blast6Segment(br_blast.readLine(), true)) !=null ) {
				sub_cov.get(tmp_record.sseqid()).add(Range.closed(tmp_record.sstart(), 
						tmp_record.send()).canonical(DiscreteDomain.integers()));
			}
			br_blast.close();
			
			BufferedWriter bw = Utils.getBufferedWriter(this.out_prefix+".bed");
			for(final String sub : sub_seqs.keySet()) {
				int cov = 0;
				for(Range<Integer> r : sub_cov.get(sub).asRanges()) {
					cov += r.upperEndpoint()-r.lowerEndpoint();
					bw.write(sub+"\t"+r.lowerEndpoint()+"\t"+r.upperEndpoint()+"\n");
				}
				myLogger.info("##subject chromosome "+sub+" covered "+(double)cov/sub_seqs.get(sub).seq_ln());
			}
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	
	
	/***
	sel_recs.clear();
	// merge collinear records
	for(int t=0; t<2; t++) {
		
		if(t==1) qry += "'";
		for(String sub : sub_seqs.keySet()) {
			// get all records for subject/reference sequence sub_seq
			tmp_records.clear();
			for(Blast6Segment record : buff)
				if(record.sseqid().equals(sub)&&record.qseqid().equals(qry))
					tmp_records.add(record);

			if(tmp_records.isEmpty()) continue;

			// find alignment segments that can be deleted
			// those that are subsets of larger alignment segments
			// (Sstart, (sstart, send), Send) and (Qstart, (qstart, qend), Qend)
			Collections.sort(tmp_records, new Blast6Segment.SegmentSizeComparator());
			final Set<int[][]> ranges = new HashSet<int[][]>();
			outerloop:
				for(int i=0; i<tmp_records.size(); i++) {
					primary_record = tmp_records.get(i);
					int[][] range = new int[2][2];
					range[0][0] = primary_record.sstart();
					range[0][1] = primary_record.send();
					range[1][0] = primary_record.qstart();
					range[1][1] = primary_record.qend();
					
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

			Collections.sort(tmp_records, new AlignmentSegment.SubjectCoordinationComparator());
			for(int i=1; i<tmp_records.size(); i++) {
				primary_record   = tmp_records.get(i-1);
				secondary_record = tmp_records.get( i );
				if(secondary_record.sstart()-primary_record.send()<gap_min) {
					
				}
			}
			
			// add to sel_recs
			sel_recs.addAll(tmp_records);
		}

		// filter by alignment coverage		
		buff.clear();
		buff.addAll(sel_recs);
		sel_recs.clear();
		for(Blast6Segment record : buff) {
			if(record.length()/qry_ln>=this.min_cov)
				sel_recs.add(record);
		}

		if(sel_recs.isEmpty()) {
			// unplaced query sequences
			// continue
			continue;
		}
		// filter blast records
		Collections.sort(sel_recs, new Blast6Segment.MatchIndentityComparator());
		// process primary alignment
		primary_record = sel_recs.get(0);
		anchored_records.get(primary_record.sseqid()).add(primary_record);
		aln_cov = primary_record.length()/qry_ln;
		// compare secondary alignments to primary alignment
		// and process
		for(int i=1; i<sel_recs.size(); i++) {
			secondary_record = sel_recs.get(i);
			if(secondary_record.pident()+this.diff_ident<primary_record.pident() ||
					secondary_record.length()/qry_ln+this.diff_cov<aln_cov) {
				break;
			} else {
				anchored_records.get(secondary_record.sseqid()).add(secondary_record);
			}
		}
	}
	 ***/
}

