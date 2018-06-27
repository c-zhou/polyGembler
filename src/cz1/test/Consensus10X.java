package cz1.test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Range;

import cz1.ngs.model.Sequence;
import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Executor;
import cz1.util.Utils;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Consensus10X extends Executor {

	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -s/--subject            Subject/reference sequences file in FASTA format.\n"
						+ " -a/--align              Alignment files of query sequences to the subject sequences. \n"
						+ " -#/--min-link           Minimum link to make consensus (default 3).\n"
						+ " -t/--threads            Number of threads to use.\n"
						+ " -o/--out-prefix         Prefix of the output files.\n"
						+ "\n");	
	}

	private String subject_file;
	private String[] bam_files;
	private String out_prefix;
	private int min_link = 10;
	private int num_threads = Runtime.getRuntime().availableProcessors();
	
	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-s", "--subject", true);
			myArgsEngine.add("-a", "--align", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-#", "--min-link", true);
			myArgsEngine.add("-o", "--out-prefix", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-s")) {
			this.subject_file = myArgsEngine.getString("-s");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the subject/reference file.");
		}
		
		if (myArgsEngine.getBoolean("-a")) {
			this.bam_files = myArgsEngine.getString("-a").split(":");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the BAM files.");
		}
		
		if (myArgsEngine.getBoolean("-#")) {
			this.min_link = Integer.parseInt(myArgsEngine.getString("-#"));
		}
		
		if (myArgsEngine.getBoolean("-t")) {
			int t = Integer.parseInt(myArgsEngine.getString("-t"));
			if(t<this.num_threads) this.num_threads = t;
			this.THREADS = t;
			Constants.omp_threads = this.num_threads;
			myLogger.info("OMP_THREADS = "+this.num_threads);
		}
		
		if (myArgsEngine.getBoolean("-o")) {
			this.out_prefix = myArgsEngine.getString("-o");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the prefix of output files.");
		}
	}
	
	private final static Object lock = new Object();
	private final static int min_scaff = 200;
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		final List<Sequence> sub_seqs = Sequence.parseFastaFileAsList(subject_file);
		final Map<String, List<Range<Integer>>> ranges = new HashMap<String, List<Range<Integer>>>();
		for(final Sequence sub_seq : sub_seqs) 
			ranges.put(sub_seq.seq_sn(), new ArrayList<Range<Integer>>());
		
		this.initial_thread_pool();
		for(final String bam_file : bam_files) {
			executor.submit(new Runnable() {
				private String bam_file;
				
				@Override
				public void run() {
					// TODO Auto-generated method stub
					try {
						myLogger.info("####now process "+this.bam_file);
						final BAMBarcodeIterator iter = new BAMBarcodeIterator(this.bam_file);
						while(iter.hasNext()) {
							List<SAMRecord[]> bc_records = iter.next();
							List<Molecule> mols = extractMoleculeFromList(bc_records);
							for(final Molecule mol : mols) {
								synchronized(lock) {
									ranges.get(mol.chr_id).add(Range.closed(mol.chr_start, mol.chr_end));
								}
							}
						}
						iter.close();
						myLogger.info("####"+this.bam_file+" processed");
					} catch (Exception e) {
						Thread t = Thread.currentThread();
						t.getUncaughtExceptionHandler().uncaughtException(t, e);
						e.printStackTrace();
						executor.shutdown();
						System.exit(1);
					}
				}

				public Runnable init(String bam_file) {
					// TODO Auto-generated method stub
					this.bam_file = bam_file;
					return this;
				}
				
			}.init(bam_file));
		}
		this.waitFor();
		
		myLogger.info("####calculate depth of coverage");
		final List<Sequence> scaffs = new ArrayList<Sequence>();
		for(final Sequence sub_seq : sub_seqs) {
			String sub_seqid = sub_seq.seq_sn();
			myLogger.info("####now process "+sub_seqid);
			if("Chr00".equals(sub_seqid)) continue;
			String seq_str = sub_seq.seq_str();
			int seq_ln = sub_seq.seq_ln();
			final byte[] cvg = new byte[seq_ln];
			for(Range<Integer> r : ranges.get(sub_seqid)) {
				for(int i=r.lowerEndpoint(); i<r.upperEndpoint(); i++)
					if(cvg[i]<Byte.MAX_VALUE) ++cvg[i];
			}
			final List<Integer> lowCov = new ArrayList<Integer>();
			lowCov.add(-1);
			for(int i=0; i<seq_ln; i++) {
				if(cvg[i]<min_link) lowCov.add(i);
			}
			lowCov.add(seq_ln);
			int start, end;
			for(int i=1; i<lowCov.size();i++) {
				start = lowCov.get(i-1)+1;
				end   = lowCov.get(i);
				if(end<start) continue;
				String scaff_str = seq_str.substring(start, end);
				scaff_str = scaff_str.replaceAll("^N{1,}", "");
				scaff_str = scaff_str.replaceAll("N{1,}$", "");
				if(scaff_str.length()>=min_scaff) {
					scaffs.add(new Sequence(sub_seqid+":"+start+"-"+end, scaff_str));
				}
			}
		}
		
		Collections.sort(scaffs, new Comparator<Sequence>() {

			@Override
			public int compare(Sequence arg0, Sequence arg1) {
				// TODO Auto-generated method stub
				return arg1.seq_ln()-arg0.seq_ln();
			}
		});
		
		int scaff_no = 1;
		for(final Sequence scaff : scaffs) {
			scaff.setSeqSn("scaffold"+String.format("%08d", scaff_no)+" "+scaff.seq_sn());
			++scaff_no;
		}
		
		try {
			BufferedWriter bw = Utils.getBufferedWriter(this.out_prefix+"_consensus.fa");
			for(final Sequence scaff : scaffs)
				bw.write(scaff.formatOutput());
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	private final SamReaderFactory factory =
			SamReaderFactory.makeDefault()
			.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
					SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
			.validationStringency(ValidationStringency.SILENT);
	private final static double min_qual = 1;
	private final static double max_ins  = 1000;
	private final static int max_gap = 20000;
	
	private class BAMBarcodeIterator {
		private final SamReader samReader;
		private final SAMRecordIterator iter;
		private SAMRecord samRecord = null;
		
		public BAMBarcodeIterator(String bam_file) {
			this.samReader = factory.open(new File(bam_file));
			this.iter = this.samReader.iterator();
			this.samRecord = iter.hasNext() ? iter.next() : null;
		}
		
		public boolean hasNext() {
			return samRecord != null;
		}
		
		public List<SAMRecord[]> next() {
			
			if(!this.hasNext()) throw new RuntimeException("!!!");
			
			List<SAMRecord[]> bc_records = new ArrayList<SAMRecord[]>();
			String bc = samRecord.getStringAttribute("BX");
			
			String sn;
			SAMRecord[] records = new SAMRecord[2];
			
			while( samRecord!=null && samRecord.getStringAttribute("BX").equals(bc) ) {
				sn = samRecord.getReadName();
				
				if( !samRecord.getReadUnmappedFlag() &&
						!samRecord.getNotPrimaryAlignmentFlag()&&
						!samRecord.getSupplementaryAlignmentFlag() ) {
					if(samRecord.getFirstOfPairFlag())
						records[0] = samRecord;
					else if(samRecord.getSecondOfPairFlag())
						records[1] = samRecord;
					else
						throw new RuntimeException("!!!");
				}
				
				while( (samRecord = iter.hasNext() ? iter.next() : null)!=null && samRecord.getReadName().equals(sn)) {
					if( !samRecord.getReadUnmappedFlag() &&
							!samRecord.getNotPrimaryAlignmentFlag()&&
							!samRecord.getSupplementaryAlignmentFlag() ) {
						if(samRecord.getFirstOfPairFlag())
							records[0] = samRecord;
						else if(samRecord.getSecondOfPairFlag())
							records[1] = samRecord;
						else
							throw new RuntimeException("!!!");
					}
				}

				if(records[0]!=null && records[1]!=null &&
						records[0].getReferenceIndex().intValue()==records[1].getReferenceIndex().intValue()) {
					if(records[0].getInferredInsertSize()+records[1].getInferredInsertSize()!=0)
						throw new RuntimeException("!!!");
					if(Math.abs(records[0].getInferredInsertSize())<=max_ins && 
							records[0].getMappingQuality()+records[1].getMappingQuality()>=min_qual) {
						bc_records.add(records);
					}
				}
				
				records = new SAMRecord[2];
			}
			
			return bc_records;
		}
		
		public void close() {
			try {
				this.iter.close();
				this.samReader.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	private class Molecule {
		private List<SAMRecord[]> reads_list;
		private String chr_id = null;
		private int chr_start = -1;
		private int chr_end = -1;

		public Molecule() {
			this.reads_list = new ArrayList<SAMRecord[]>();
		}

		public void add(SAMRecord[] record) {
			this.reads_list.add(record);
		}

		public void construct() {
			// TODO Auto-generated method stub
			this.chr_id = this.reads_list.get(0)[0].getReferenceName();
			int chr_start = Integer.MAX_VALUE;
			int chr_end   = Integer.MIN_VALUE;
			for(SAMRecord[] records : this.reads_list) {
				for(SAMRecord record : records) {
					if(record.getAlignmentStart()<chr_start)
						chr_start = record.getAlignmentStart();	
					if(record.getAlignmentEnd()>chr_end)
						chr_end   = record.getAlignmentEnd();
				}
			}
			this.chr_start = chr_start-1;
			this.chr_end   = chr_end;
		}
	}
	
	private double middlePoint(SAMRecord[] record) {
		// TODO Auto-generated method stub
		return (double)(Math.min(record[0].getAlignmentStart(), record[1].getAlignmentStart())+
				Math.max(record[0].getAlignmentEnd(), record[1].getAlignmentEnd()))/2;
	}
	
	private List<Molecule> extractMoleculeFromList(List<SAMRecord[]> list) {
		// TODO Auto-generated method stub
		List<Molecule> mols = new ArrayList<Molecule>();
		
		if(list.isEmpty()) return mols; 
		
		Collections.sort(list, new Comparator<SAMRecord[]>() {
			@Override
			public int compare(SAMRecord[] record0, SAMRecord[] record1) {
				// TODO Auto-generated method stub
				int f = record0[0].getReferenceIndex().intValue()-
						record1[0].getReferenceIndex().intValue();
				return f==0 ? Double.compare(middlePoint(record0), middlePoint(record1)) : f;
			}
		});
		
		Iterator<SAMRecord[]> iter = list.iterator();
		SAMRecord[] records = iter.next();
		int chr_index = records[0].getReferenceIndex();
		double chr_pos = middlePoint(records);
		Molecule mol = new Molecule();
		mol.add(records);
		while(iter.hasNext()) {
			records = iter.next();
			if(records[0].getReferenceIndex()!=chr_index ||
					middlePoint(records)-chr_pos>max_gap) {
				mol.construct();
				mols.add(mol);
				mol = new Molecule();		
			}
			chr_index = records[0].getReferenceIndex();
			chr_pos = middlePoint(records);
			mol.add(records);
		}
		mol.construct();
		mols.add(mol);
		
		return mols;
	}
	
	public static void main(String[] args) {
		Consensus10X cons = new Consensus10X();
		cons.setParameters(args);
		cons.run();
	}
}















