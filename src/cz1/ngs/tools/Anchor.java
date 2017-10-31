package cz1.ngs.tools;

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

import cz1.ngs.model.Blast6Record;
import cz1.ngs.model.BlastRecord;
import cz1.ngs.model.Sequence;
import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

public class Anchor extends Executor {
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -s/--subject            Subject/reference sequences file in FASTA format.\n"
						+ " -q/--query              Query sequences to anchor in FASTA format.\n"
						+ " -b/--blast              BLAST output (-outfmt 6) of query sequences against subject sequences.\n"
						+ " -i/--min-identity       Minimum identity between the query and subject sequences \n"
						+ "                         for an alignment record to consider (default 0.9).\n"
						+ " -f/--min-fraction       Minimum alignment fraction of the query sequence (default 0.5).\n"
						+ " -di/--diff-identity     Threshold of the difference of the identity between the primary and secondary \n"
						+ "                         alignments. If the difference is smaller than this value, the query \n"
						+ "                         sequence will be regarded as duplications. Otherwise, the secondary \n"
						+ "                         alignments will be discared (default 0.01).\n"
						+ " -df/--diff-fraction     Threshold of the difference of the alignment fraction between the primary and \n"
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
	private double min_frac = 0.50;
	private double diff_ident = 0.01;
	private double diff_frac = 0.05;
	private double collinear_shift = 0.5;
	private String out_prefix = null;
	
	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-s", "--subject", true);
			myArgsEngine.add("-q", "--query", true);
			myArgsEngine.add("-b", "--blast", true);
			myArgsEngine.add("-i", "--min-identity", true);
			myArgsEngine.add("-f", "--min-fraction", true);
			myArgsEngine.add("-di", "--diff-identity", true);
			myArgsEngine.add("-df", "--diff-fraction", true);
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
		
		if (myArgsEngine.getBoolean("-f")) {
			this.min_frac = Double.parseDouble(myArgsEngine.getString("-f"));
		}
		
		if (myArgsEngine.getBoolean("-di")) {
			this.diff_ident = Double.parseDouble(myArgsEngine.getString("-di"));
		}
		
		if (myArgsEngine.getBoolean("-df")) {
			this.diff_frac = Double.parseDouble(myArgsEngine.getString("-df"));
		}
	}

	private final static int max_clip = 10; // maximum clip of query alignment allowed
	private final static int gap_buff = 10; // buffer size for subject/reference sequences gap clips
	private final static int min_gap =  10; // take this if estimated gap size is smaller than this
	private final static int max_gap = 100; // take this if estimated gap size is larger than this
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		final Map<String, Sequence> sub_seqs = Sequence.parseFastaFileAsMap(subject_file);
		final Map<String, Sequence> qry_seqs = Sequence.parseFastaFileAsMap(query_file);
		
		// find 'N/n's in subject/reference sequences
		// which could have impact on parsing the blast records
		final Map<String, RangeSet<Integer>> sub_gaps = new HashMap<String, RangeSet<Integer>>();
		final Map<String, List<Blast6Record>> anchored_records = new HashMap<String, List<Blast6Record>>();
		
		for(Map.Entry<String, Sequence> entry : sub_seqs.entrySet()) {
			String seq_sn = entry.getKey();
			String seq_str = entry.getValue().seq_str();
			
			final RangeSet<Integer> range_set = TreeRangeSet.create();
			for(int j=0; j<seq_str.length(); j++) {
				if(seq_str.charAt(j)=='N'||seq_str.charAt(j)=='n')
					// blast record is 1-based closed coordination
					range_set.add(Range.closed(j+1, j+1).
							canonical(DiscreteDomain.integers()));
			}
			sub_gaps.put(seq_sn, range_set);
			// initialise an array list for each reference chromosome
			anchored_records.put(seq_sn, new ArrayList<Blast6Record>());
		}

		// parse blast records
		// blast records buffer
		final List<Blast6Record> buff = new ArrayList<Blast6Record>();
		// selected blast records
		final List<Blast6Record> sel_recs = new ArrayList<Blast6Record>();
		// temp list
		final List<Blast6Record> tmp_records = new ArrayList<Blast6Record>();
		try {
			BufferedReader br_blast = Utils.getBufferedReader(blast_out);
			Blast6Record tmp_record = Blast6Record.blast6Record(br_blast.readLine());
			Blast6Record primary_record, secondary_record;
			String qry;
			double qry_ln, aln_frac;
			while(tmp_record!=null) {
				qry = tmp_record.qseqid();
				qry_ln = qry_seqs.get(qry).seq_ln();
				buff.clear();
				buff.add(tmp_record);
				while( (tmp_record=Blast6Record.blast6Record(br_blast.readLine()))!=null
						&&
						tmp_record.qseqid().equals(qry) ) {
					// filter by identity
					if(tmp_record.pident()>=this.min_ident)
						buff.add(tmp_record);
				}
				
				sel_recs.clear();
				// merge collinear records
				for(String sub : sub_seqs.keySet()) {
					// get all records for subject/reference sequence sub_seq
					tmp_records.clear();
					for(Blast6Record record : buff)
						if(record.sseqid().equals(sub))
							tmp_records.add(record);
					
					// find collinear alignment segments that can be merged
					Collections.sort(tmp_records, new BlastRecord.SInterceptComparator());
					
					for(int i=0; i<tmp_records.size(); ) {
						Blast6Record record = tmp_records.get(i);
						double sintercept = record.sintercept();
						int qlen = record.length();
						double max_shift;
						final List<Blast6Record> temp = new ArrayList<Blast6Record>();
						temp.add(record);
						
						while( (++i)<tmp_records.size() ) {
							record = tmp_records.get(i);
							max_shift = collinear_shift*Math.min(qlen, record.length());
							if(record.sintercept()-sintercept>max_shift) {
								break;
							} else {
								temp.add(record);
								sintercept = record.sintercept();
								qlen = record.length();
							}
						}
						
						int qstart = Integer.MAX_VALUE;
						int qend = Integer.MIN_VALUE;
						int sstart = Integer.MAX_VALUE;
						int send = Integer.MIN_VALUE;
						double pident = 0;
						int length = 0;
						
						for(int j=0; j<temp.size(); j++) {
							record = temp.get(j);
							if(record.qstart()<qstart) {
								qstart = record.qstart();
								sstart = record.sstart();
							}
							if(record.qend()>qend) {
								qend = record.qend();
								send = record.send();
							}
							if(record.pident()>pident)
								pident = record.pident();
							length += record.length();
						}
						
						sel_recs.add(new Blast6Record(qry,sub,pident,length,-1,-1,qstart,qend,sstart,send,-1,-1));
					}
				}
				
				// filter by alignment fraction		
				buff.clear();
				buff.addAll(sel_recs);
				sel_recs.clear();
				for(Blast6Record record : buff) {
					if(record.length()/qry_ln>=this.min_frac) {
						sel_recs.add(record);
					} else {
						// TODO blast records that clipped by 
						// 1. gaps: can be bridged
						// or
						// 2. ends: can be extended
						// on subject sequences
						
					}
				}
				
				sel_recs.clear();
				sel_recs.addAll(buff);
				if(sel_recs.isEmpty()) {
					// unplaced query sequences
					// continue
					continue;
				}
				// filter blast records
				Collections.sort(sel_recs, new Blast6Record.MatchIndentityComparator());
				// process primary alignment
				primary_record = sel_recs.get(0);
				anchored_records.get(primary_record.sseqid()).add(primary_record);
				aln_frac = primary_record.length()/qry_ln;
				// compare secondary alignments to primary alignment
				// and process
				for(int i=1; i<sel_recs.size(); i++) {
					secondary_record = sel_recs.get(i);
					if(secondary_record.pident()+this.diff_ident<primary_record.pident() ||
							secondary_record.length()/qry_ln+this.diff_frac<aln_frac) {
						break;
					} else {
						anchored_records.get(secondary_record.sseqid()).add(secondary_record);
					}
				}
			}
			br_blast.close();
		
			final BufferedWriter bw_map = Utils.getBufferedWriter(out_prefix+".map");
			final BufferedWriter bw_fa = Utils.getBufferedWriter(out_prefix+".fa");
			final Set<String> anchored_seqs = new HashSet<String>();
			
			for(Map.Entry<String, List<Blast6Record>> entry : anchored_records.entrySet()) {
				List<Blast6Record> blast6_records = entry.getValue();
				// sort blast records
				Collections.sort(blast6_records, 
						new Blast6Record.SubjectCoordinationComparator());
				// consensus
				String sub_sn = entry.getKey();
				int posUpto = 0, send_clip = 0;
				int sstart, send, qstart, qend, qlen, tmp_int, qstart_clip, qend_clip, gap_size, overlap;
				int slen = sub_seqs.get(sub_sn).seq_ln();
				int mol_len = 0;
				
				// will form a pseudomolecule
				final List<Sequence> sequences = new ArrayList<Sequence>();
				for(Blast6Record record : blast6_records) {
					
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
					} else {
						qstart_clip = qstart-1;
						qend_clip = qlen-qend;
					}

					if(send<posUpto) {
						// skip
						//     ====================
						//      /-------\
						//        /----\
						// TODO process
						continue;
					} else if(sstart<posUpto) {
						// overlap found
						//     ====================
						//      /-------\
						//            /----\
						// calculate overlaps
						if(send_clip>max_clip || qstart_clip>max_clip) {
							// to much clipping
							// cannot be treated as overlap
							// insert max_gap
							bw_map.write("GAP\t"+
									max_gap+
									"\t0\t"+
									max_gap+
									"\t+\t"+
									sub_sn+
									"\t"+
									mol_len+
									"\t"+
									(mol_len+max_gap)+
									"\n");
							sequences.add(Sequence.polyN(max_gap));
							mol_len += max_gap;
							bw_map.write(record.qseqid()+
									"\t"+
									qlen+
									"\t");
							if(qstart>qend) {
								// reverse
								bw_map.write("0\t"+qlen+"\t-\t");
								sequences.add( new Sequence( record.qseqid(), 
										qry_seqs.get(record.qseqid()).revCompSeq() ) );
							} else {
								// forward
								bw_map.write("0\t"+qlen+"\t+\t");
								sequences.add(qry_seqs.get(record.qseqid()));
							}
							bw_map.write(sub_sn+
									"\t"+
									mol_len+
									"\t"+
									(mol_len+qlen)+
									"\n");
							mol_len += qlen;
						} else {
							// process overlap
							overlap = posUpto-sstart+1;
							// no gap inserted
							
							if(qstart>qend&&qstart<=overlap || 
									qstart<=qend&&qstart+overlap>qlen) {
								continue;
							}
							
							// get sequence clipped from overlap
							if(qstart>qend) {
								// reverse
								qstart -= overlap;
								bw_map.write(record.qseqid()+
										"\t"+
										qstart+
										"\t");
								bw_map.write( 0+"\t"+qstart+"\t-\t" );
								
								bw_map.write(sub_sn+
										"\t"+
										mol_len+
										"\t"+
										(mol_len+qstart)+
										"\n");
								mol_len += qstart;
								sequences.add( new Sequence(record.qseqid(), 
										Sequence.revCompSeq(qry_seqs.get(record.qseqid()).seq_str().substring(0, qstart)) ) );
							} else {
								// forward
								qstart += overlap;
								bw_map.write(record.qseqid()+
										"\t"+
										(qlen-qstart+1)+
										"\t");
								bw_map.write( (qstart-1)+"\t"+qlen+"\t+\t" );
								
								bw_map.write(sub_sn+
										"\t"+
										mol_len+
										"\t"+
										(mol_len+qlen-qstart+1)+
										"\n");
								mol_len += qlen-qstart+1;
								sequences.add( new Sequence(record.qseqid(), 
										qry_seqs.get(record.qseqid()).seq_str().substring(qstart-1)) );
							}
						}	
					} else if(sstart>=posUpto) {
						// simply extend
						//     ====================
						//      /-------\
						//               /----\

						if(posUpto>0) {
							// estimate gap size
							gap_size = (sstart-posUpto)-(send_clip+qstart_clip);
							if(gap_size<min_gap) gap_size = min_gap;
							if(gap_size>max_gap) gap_size = max_gap;

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
							sequences.add(Sequence.polyN(gap_size));
							mol_len += gap_size;
						}
						bw_map.write(record.qseqid()+
								"\t"+
								qlen+
								"\t");
						if(qstart>qend) {
							// reverse
							bw_map.write("0\t"+qlen+"\t-\t");
							sequences.add( new Sequence( record.qseqid(), 
									qry_seqs.get(record.qseqid()).revCompSeq() ) );
						} else {
							// forward
							bw_map.write("0\t"+qlen+"\t+\t");
							sequences.add(qry_seqs.get(record.qseqid()));
						}
						bw_map.write(sub_sn+
								"\t"+
								mol_len+
								"\t"+
								(mol_len+qlen)+
								"\n");
						mol_len += qlen;
					} else {
						throw new RuntimeException("!!!");
					}
					
					posUpto = send;
					send_clip = qend_clip;
				}
				
				StringBuilder seq_str = new StringBuilder();
				for(Sequence seq : sequences) {
					seq_str.append(seq.seq_str());
					anchored_seqs.add(seq.seq_sn());
				}
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
}










