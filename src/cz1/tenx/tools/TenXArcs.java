package cz1.tenx.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class TenXArcs extends Executor {

	private String[] bamFiles = null;
	private String outdv = null;
	private String bcList = null;
	private int bcn = 3;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -b/--bam-file       Input BAM file.\n"
						+ " -f/--bam-fof        File of BAM file list.\n"
						+ " -w/--white-list     File of barcodes on the white list.\n"
						+ " -n/--read-number    Minimum number of read pair per barcode (defalt: 3).\n"
						+ " -t/--threads        Threads.\n"
						+ " -o/--out            Output file.\n\n");	
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-b", "--bam-file", true);
			myArgsEngine.add("-f", "--bam-fof", true);
			myArgsEngine.add("-w", "--white-list", true);
			myArgsEngine.add("-n", "--read-number", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-o", "--out", true);
			myArgsEngine.parse(args);
		}
		if (!myArgsEngine.getBoolean("-b")&&!myArgsEngine.getBoolean("-f")) {
			printUsage();
			throw new IllegalArgumentException("Please specify the BAM file or a BAM list file.");
		}
		if(myArgsEngine.getBoolean("-b")) {
			this.bamFiles = new String[]{myArgsEngine.getString("-b")};
		}
		
		if(myArgsEngine.getBoolean("-f")) {
			if(this.bamFiles!=null) 
				myLogger.warn("The program will use the BAM list file instead of the BAM file you provided.");
			this.bamFiles = myArgsEngine.getString("-f").split(",");
		}
		
		if(myArgsEngine.getBoolean("-w")) {
			this.bcList = myArgsEngine.getString("-w");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the barcode white list.");
		}
		
		if(myArgsEngine.getBoolean("-n")) {
			this.bcn = Integer.parseInt(myArgsEngine.getString("-n"));
		}
		
		this.THREADS = Runtime.getRuntime().availableProcessors();
		if(myArgsEngine.getBoolean("-t")) {
			this.THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}
		
		if(myArgsEngine.getBoolean("-o")) {
			this.outdv = myArgsEngine.getString("-o");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the output file.");
		}
	}
	

	private final SamReaderFactory factory =
			SamReaderFactory.makeDefault()
			.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
					SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
			.validationStringency(ValidationStringency.SILENT);
	
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
					bc_records.add(records);
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

	final int shift = 30;
	final int rev = 1<<shift;
	final Map<Long, Integer> arcs = new HashMap<Long, Integer>();
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		final Set<String> bc_white = new HashSet<String>();
		
		try {
			final BufferedReader br_bc = Utils.getBufferedReader(this.bcList);
			String line;
			while( (line=br_bc.readLine()) != null ) 
				bc_white.add(line.split(",")[0]);
			br_bc.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		final SamReader samReader = factory.open(new File(this.bamFiles[0]));
		final SAMSequenceDictionary seqdict = samReader.getFileHeader().getSequenceDictionary();
		try {
			samReader.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
	
		this.initial_thread_pool();
		for(final String bamFile : bamFiles) {
			this.executor.submit(new Runnable() {
				private String bamFile;
				
				@Override
				public void run() {
					// TODO Auto-generated method stub
					try {
						final BAMBarcodeIterator iter = new BAMBarcodeIterator(this.bamFile);
						List<SAMRecord[]> bc_records;
						Map<String, List<SAMRecord[]>> bc_blocks;
						while(iter.hasNext()) {
							bc_records = iter.next();
							if(bc_records.isEmpty() || 
									!bc_white.contains(bc_records.get(0)[0].getStringAttribute("BX")))
								continue;
							bc_blocks = bin(bc_records);
						}

						iter.close();
					} catch (Exception e) {
						Thread t = Thread.currentThread();
						t.getUncaughtExceptionHandler().uncaughtException(t, e);
						e.printStackTrace();
						executor.shutdown();
						System.exit(1);
					}
				}

				private Map<String, List<SAMRecord[]>> bin(List<SAMRecord[]> bc_records) {
					// TODO Auto-generated method stub
					final Map<String, List<SAMRecord[]>> bc_blocks = new HashMap<String, List<SAMRecord[]>>();
					
					return bc_blocks;
				}

				public Runnable init(String bamFile) {
					// TODO Auto-generated method stub
					this.bamFile = bamFile;
					return this;
				}
				
			}.init(bamFile));
		}
		this.waitFor();
	}

}





