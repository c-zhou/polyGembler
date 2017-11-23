package cz1.tenx.tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class BAMstats extends Executor {
	
	private final static int int_min = Integer.MIN_VALUE;
	private final static int int_max = Integer.MAX_VALUE;
	private final static float float_min = Float.MIN_VALUE;
	private final static float float_max = Float.MAX_VALUE;
	
	private final static int min_seqLen = 35;
	
	// mapped to 1, 2, and intersect
	// alignment score difference 0
	private final static long record_mapped_as0[] = new long[3];
	// alignment score difference infinite
	private final static long record_mapped_asInf[] = new long[3]; // mapped to 1, 2, and intersect
	// BWA-MEM read alignment score 
	private final static Map<Float, Integer> alignment_score_c1 = new HashMap<Float, Integer>();
	static {
		alignment_score_c1.put(float_min, 0);
		alignment_score_c1.put(float_max, 0);
	};
	private final static Map<Float, Integer> alignment_score_c2 = new HashMap<Float, Integer>();
	static {
		alignment_score_c2.put(float_min, 0);
		alignment_score_c2.put(float_max, 0);
	};
	// BWA-MEM read alignment score difference between 1 and 2 (1 minus 2)
	private final static Map<Float, Integer> alignment_score_cdiff = new HashMap<Float, Integer>();
	static {
		for(float f=-1000; f<=1000; f+=0.5) 
			alignment_score_cdiff.put(f, 0);
		alignment_score_cdiff.put(float_min, 0);
		alignment_score_cdiff.put(float_max, 0);
	};
	// insert size and counts on chromosomes
	// none of the two reads mapped to the unanchored scaffolds - Chr00
	private final static Map<Integer, Integer> insert_size_c1 = new HashMap<Integer, Integer>();
	static {
		insert_size_c1.put(int_min, 0);
		insert_size_c1.put(int_max, 0);	
	};
	private final static Map<Integer, Integer> insert_size_c2 = new HashMap<Integer, Integer>();
	static {
		insert_size_c2.put(int_min, 0);
		insert_size_c2.put(int_max, 0);
	};
	// insert size difference between 1 and 2 (1 minus 2), for those both mapped in pair and on chromosomes
	private final static Map<Integer, Integer> insert_size_cdiff = new HashMap<Integer, Integer>();
	
	private static long record_counter = 0;
	
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i1/--in-bam1       First BAM file.\n"
						+ " -i2/--in-bam2       Second BAM file. \n"
						+ " -o/--prefix         Output prefix.\n\n");	
	}

	private String bam_in1 = null;
	private String bam_in2 = null;
	private String out = null;
	
	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-i1", "--in-bam1", true);
			myArgsEngine.add("-i2", "--in-bam2", true);
			myArgsEngine.add("-o", "--out", true);
			myArgsEngine.parse(args);
		}
		if (myArgsEngine.getBoolean("-i1")) {
			this.bam_in1 = myArgsEngine.getString("-i1");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the input BAM file.");
		}
		
		if (myArgsEngine.getBoolean("-i2")) {
			this.bam_in2 = myArgsEngine.getString("-i2");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the input BAM file.");
		}
		
		if (myArgsEngine.getBoolean("-o")) {
			this.out = myArgsEngine.getString("-o");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the output file.");
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		final SamReaderFactory factory =
		          SamReaderFactory.makeDefault()
		              .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
		            		  SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
		              .validationStringency(ValidationStringency.SILENT);
		
		try {
			final SamReader in1 = factory.open(new File(bam_in1));
			final SamReader in2 = factory.open(new File(bam_in2));
			SAMRecordIterator iter1 = in1.iterator();
			SAMRecordIterator iter2 = in2.iterator();
			
			SAMRecord tmp_record1 = iter1.next(), tmp_record2 = iter2.next(); 
			final SAMRecord[] sam_record1 = new SAMRecord[2],
					sam_record2 = new SAMRecord[2];
			final float[] sam_as1 = new float[2], sam_as2 = new float[2];
			final boolean[] sam_mapped1 = new boolean[2],
					sam_mapped2 = new boolean[2];
			String record_id1 = tmp_record1.getReadName(),
					record_id2 = tmp_record2.getReadName();
			int i_sz;
			float i_as;

			boolean b1, b2;
			
			while(true) {
				
				if(record_counter%4000000==0) myLogger.info(record_counter+" processed.");
				
				while(record_id1!=null&&record_id2!=null&&!record_id1.equals(record_id2)) {
					
					while(record_id1!=null&&record_id1.compareTo(record_id2)<0) {
						tmp_record1 = bufferRecord(iter1, 
								tmp_record1, 
								sam_record1, 
								sam_as1, 
								sam_mapped1, 
								record_id1);

						processBuffer(sam_record1,
								sam_as1,
								sam_mapped1,
								0);
						
						record_id1 = tmp_record1==null?null:tmp_record1.getReadName();
					}

					if(record_id1==null) break; 
					
					while(record_id2!=null&&record_id2.compareTo(record_id1)<0) {
						
						tmp_record2 = bufferRecord(iter2, 
								tmp_record2, 
								sam_record2, 
								sam_as2, 
								sam_mapped2, 
								record_id2);
						
						processBuffer(sam_record2,
								sam_as2,
								sam_mapped2,
								1);
						
						record_id2 = tmp_record2==null?null:tmp_record2.getReadName();
					}
				}
				
				if(record_id1==null || record_id2==null) break;
				
				tmp_record1 = bufferRecord(iter1, 
						tmp_record1, 
						sam_record1, 
						sam_as1, 
						sam_mapped1, 
						record_id1);
				tmp_record2 = bufferRecord(iter2, 
						tmp_record2, 
						sam_record2, 
						sam_as2, 
						sam_mapped2, 
						record_id2);
				
				record_id1 = tmp_record1==null?null:tmp_record1.getReadName();
				record_id2 = tmp_record2==null?null:tmp_record2.getReadName();
				
				if( (sam_record2[0]==null || 
						sam_record2[1]==null) && 
						(sam_record1[0]==null || 
						sam_record1[1]==null) )
					continue;
				
				if( (sam_record1[0]==null || 
						sam_record1[1]==null) && 
						sam_record2[0]!=null && 
						sam_record2[1]!=null ) {
					processBuffer(sam_record2,
							sam_as2,
							sam_mapped2,
							1);
					continue;
				}
				
				if( (sam_record2[0]==null || 
						sam_record2[1]==null) && 
						sam_record1[0]!=null && 
						sam_record1[1]!=null ) {
					processBuffer(sam_record1,
							sam_as1,
							sam_mapped1,
							0);
					continue;
				}
				
				if(sam_record1[0].getDuplicateReadFlag() ||
						sam_record1[1].getDuplicateReadFlag() ||
						sam_record2[0].getDuplicateReadFlag() ||
						sam_record2[1].getDuplicateReadFlag() ||
						sam_record1[0].getReadString().replaceAll("N+$", "").replaceAll("^N+", "").length()<min_seqLen ||
						sam_record1[1].getReadString().replaceAll("N+$", "").replaceAll("^N+", "").length()<min_seqLen ||
						sam_record2[0].getReadString().replaceAll("N+$", "").replaceAll("^N+", "").length()<min_seqLen ||
						sam_record2[1].getReadString().replaceAll("N+$", "").replaceAll("^N+", "").length()<min_seqLen ) 
					// we need both reads longer than 36
					continue;
				
				
				record_counter += 2;
				
				if( sam_mapped1[0] ) record_mapped_asInf[0]++;
				if( sam_mapped1[1] ) record_mapped_asInf[0]++;
				if( sam_mapped2[0] ) record_mapped_asInf[1]++;
				if( sam_mapped2[1] ) record_mapped_asInf[1]++;
				if( sam_mapped1[0] && sam_mapped2[0] ) record_mapped_asInf[2]++;
				if( sam_mapped1[1] && sam_mapped2[1] ) record_mapped_asInf[2]++;
				
				if( b1=(sam_mapped1[0] && sam_mapped1[1] &&
					!sam_record1[0].getReferenceName().equals("Chr00") && 
					!sam_record1[1].getReferenceName().equals("Chr00")) ) {
						i_sz = Math.abs(sam_record1[0].getInferredInsertSize());
						if(!insert_size_c1.containsKey(i_sz))
							insert_size_c1.put(i_sz, 1);
						else
							insert_size_c1.put(i_sz, insert_size_c1.get(i_sz)+1);
					}
				
				if( b2=(sam_mapped2[0] && sam_mapped2[1] &&
					!sam_record2[0].getReferenceName().equals("Chr00") && 
					!sam_record2[1].getReferenceName().equals("Chr00")) ) {
						i_sz = Math.abs(sam_record2[0].getInferredInsertSize());
						if(!insert_size_c2.containsKey(i_sz))
							insert_size_c2.put(i_sz, 1);
						else
							insert_size_c2.put(i_sz, insert_size_c2.get(i_sz)+1);
					}

				if(!(b1&b2)) continue;
				
				i_sz = Math.abs(sam_record1[0].getInferredInsertSize())-
						Math.abs(sam_record2[0].getInferredInsertSize());
				if(!insert_size_cdiff.containsKey(i_sz))
					insert_size_cdiff.put(i_sz, 1);
				else
					insert_size_cdiff.put(i_sz, insert_size_cdiff.get(i_sz)+1);
				
				// alignment score for pairs always the same
				i_as = sam_as1[0]-sam_as2[0];
				if(!alignment_score_cdiff.containsKey(i_as)) 
					alignment_score_cdiff.put(i_as, 1);
				else
					alignment_score_cdiff.put(i_as, alignment_score_cdiff.get(i_as)+1);

				int i = i_as > 0 ? 0 : (i_as < 0 ? 1 : 2);
				if( sam_mapped1[0] && !sam_mapped2[0]) record_mapped_as0[0]++;
				if( !sam_mapped1[0] && sam_mapped2[0]) record_mapped_as0[1]++;
				if( sam_mapped1[0] && sam_mapped2[0] ) record_mapped_as0[i]++;
				
				if( sam_mapped1[1] && !sam_mapped2[1]) record_mapped_as0[0]++;
				if( !sam_mapped1[1] && sam_mapped2[1]) record_mapped_as0[1]++;
				if( sam_mapped1[1] && sam_mapped2[1] ) record_mapped_as0[i]++;
				
			}
			
			
			if(record_id1==null) {
				
				while(record_id2!=null) {
					tmp_record2 = bufferRecord(iter2, 
							tmp_record2, 
							sam_record2, 
							sam_as2, 
							sam_mapped2, 
							record_id2);
					
					processBuffer(sam_record2,
							sam_as2,
							sam_mapped2,
							1);
					
					record_id2 = tmp_record2==null?null:tmp_record2.getReadName();
				}
				
			}
			
			if(record_id2==null) {
				
				while(record_id1!=null) {
					tmp_record1 = bufferRecord(iter1, 
							tmp_record1, 
							sam_record1, 
							sam_as1, 
							sam_mapped1, 
							record_id1);

					processBuffer(sam_record1,
							sam_as1,
							sam_mapped1,
							0);

					record_id1 = tmp_record1==null?null:tmp_record1.getReadName();
				}
			}
			
			iter1.close();
			iter2.close();
			in1.close();
			in2.close();
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(out));
			bw.write("###record_total\n"+record_counter+"\n");
			bw.write("###record_mapped ( 0 )\n"+record_mapped_as0[0]+", "+record_mapped_as0[1]+", "+record_mapped_as0[2]+"\n");
			bw.write("###record_mapped (Inf)\n"+record_mapped_asInf[0]+", "+record_mapped_asInf[1]+", "+record_mapped_asInf[2]+"\n");
			bw.write("###alignment_score_c1\n");
			for(float i : alignment_score_c1.keySet()) bw.write(i+"\t"+alignment_score_c1.get(i)+"\n");
			bw.write("###alignment_score_c2\n");
			for(float i : alignment_score_c2.keySet()) bw.write(i+"\t"+alignment_score_c2.get(i)+"\n");
			bw.write("###alignment_score_cdiff\n");
			for(float i : alignment_score_cdiff.keySet()) bw.write(i+"\t"+alignment_score_cdiff.get(i)+"\n");
			bw.write("###insert_size_c1\n");
			for(int i : insert_size_c1.keySet()) bw.write(i+"\t"+insert_size_c1.get(i)+"\n");
			bw.write("###insert_size_c2\n");
			for(int i : insert_size_c2.keySet()) bw.write(i+"\t"+insert_size_c2.get(i)+"\n");
			bw.write("###insert_size_cdiff\n");
			for(int i : insert_size_cdiff.keySet()) bw.write(i+"\t"+insert_size_cdiff.get(i)+"\n");
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static void processBuffer(final SAMRecord[] sam_record,
			final float[] sam_as, 
			final boolean[] sam_mapped,
			final int i) {
		// TODO Auto-generated method stub
		
		if(sam_record[0]==null || 
				sam_record[1]==null || 
				sam_record[0].getDuplicateReadFlag() || 
				sam_record[1].getDuplicateReadFlag() ||
				sam_record[0].getReadString().replaceAll("N+$", "").replaceAll("^N+", "").length()<min_seqLen ||
				sam_record[1].getReadString().replaceAll("N+$", "").replaceAll("^N+", "").length()<min_seqLen ) 
			// we need both reads longer than 36
			return;

		record_counter += 2;
		
		if(sam_mapped[0]) {
			++record_mapped_as0[i];
			++record_mapped_asInf[i];
		}
		if(sam_mapped[1]) {
			++record_mapped_as0[i];
			++record_mapped_asInf[i];
		}
		
		if(sam_mapped[0]&&sam_mapped[1]&&
				!sam_record[0].getReferenceName().equals("Chr00") &&
				!sam_record[1].getReferenceName().equals("Chr00") ) {
			
			if(i==0) {
				if(!alignment_score_c1.containsKey(sam_as[0]))
					alignment_score_c1.put(sam_as[0], 1);
				else
					alignment_score_c1.put(sam_as[0], 
							alignment_score_c1.get(sam_as[0])+1);
				alignment_score_cdiff.put(float_max, 
						alignment_score_cdiff.get(float_max)+1);
			} else {
				if(!alignment_score_c2.containsKey(sam_as[0]))
					alignment_score_c2.put(sam_as[0], 1);
				else
					alignment_score_c2.put(sam_as[0], 
							alignment_score_c2.get(sam_as[0])+1);
				alignment_score_cdiff.put(float_min, 
						alignment_score_cdiff.get(float_min)+1);
			}
			
			if(sam_record[0].getReferenceIndex()==
					sam_record[1].getReferenceIndex()) {
				int i_sz = sam_record[0].getInferredInsertSize();
				if(i_sz<0) i_sz = -i_sz;
				if(i==0) {
					if(!insert_size_c1.containsKey(i_sz)) insert_size_c1.put(i_sz, 1);
					else insert_size_c1.put(i_sz, insert_size_c1.get(i_sz)+1);
				} else {
					if(!insert_size_c2.containsKey(i_sz)) insert_size_c2.put(i_sz, 1);
					else insert_size_c2.put(i_sz, insert_size_c2.get(i_sz)+1);
				} 
			}else {
				if(i==0) {
					insert_size_c1.put(int_max, insert_size_c1.get(int_max)+1);
				} else {
					insert_size_c2.put(int_max, insert_size_c2.get(int_max)+1);
				}
			}
		}
	}

	private static SAMRecord bufferRecord(final SAMRecordIterator iter, 
			final SAMRecord sam_record,
			final SAMRecord[] buff_record, 
			final float[] as,
			final boolean[] aln,
			final String record_id) {
		// TODO Auto-generated method stub

		final Set<SAMRecord> buffered_records = new HashSet<SAMRecord>();
		buffered_records.add(sam_record);
		SAMRecord tmp_record;
		while( (tmp_record=iter.hasNext()?iter.next():null)!=null && 
				tmp_record.getReadName().equals(record_id)) {
			buffered_records.add(tmp_record);
		}
		Arrays.fill(buff_record, null);
		Arrays.fill(as, -1.0f);
		Arrays.fill(aln, false);
		
		for(SAMRecord r : buffered_records) {
			if( r.getFirstOfPairFlag() && 
					!r.getSupplementaryAlignmentFlag() &&
					!r.getNotPrimaryAlignmentFlag() ) {
				buff_record[0] = r;	
				as[0] = r.getFloatAttribute("AS");
				aln[0] = !r.getReadUnmappedFlag();
			}
			if( r.getSecondOfPairFlag() && 
					!r.getSupplementaryAlignmentFlag() &&
					!r.getNotPrimaryAlignmentFlag() ) {
				buff_record[1] = r;
				as[1] = r.getFloatAttribute("AS");
				aln[1] = !r.getReadUnmappedFlag();
			}
		}
		
		return tmp_record;
	}
}
