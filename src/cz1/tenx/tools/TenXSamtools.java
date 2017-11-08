package cz1.tenx.tools;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.TreeMap;
import java.util.Map.Entry;

import cz1.util.ArgsEngine;
import cz1.util.Executor;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class TenXSamtools extends Executor {

	private static enum Task {sort, zzz}
	private Task task = Task.zzz;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		switch(this.task) {
		case zzz:
			myLogger.info(
					"\n\nChoose a task from listed:\n"
							+ " sort                Sort BAM file.\n\n");	
			break;
		case sort:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -i/--in-bam         Input BAM file.\n"
							+ " -b/--barcode        Sort by barcode. (BX tag).\n"
							+ " -n/--name           Sort by read name.\n"
							+ " -s/--batch-size     Batch size. (default 4000000).\n"
							+ " -t/--threads        Threads to use. (default 1).\n"
							+ " -o/--out-bam        Output files prefix.\n\n");	
			break;
		default:
			throw new RuntimeException("Undefined task!!!");
		}
	}
	
	private static enum Order {coordinate, barcode, queryname}
	private Order sort_order = Order.coordinate;
	
	private Comparator<SAMRecord> comprator = new Comparator<SAMRecord>() {
		@Override
		public int compare(SAMRecord record1, SAMRecord record2) {
			// TODO Auto-generated method stub
			if(record1==null&&record2==null) return 0;
			if(record1==null) return  1;
			if(record2==null) return -1;
			if(record1.getReadUnmappedFlag()&&record2.getReadUnmappedFlag())
				return 0;
			// unmapped record to the end
			if(record1.getReadUnmappedFlag()) return  1;
			if(record2.getReadUnmappedFlag()) return -1;
			if(record1.getAlignmentStart()-record2.getAlignmentStart()==0)
				return compareSAMRecord(record1, record2);
			return record1.getAlignmentStart()-record2.getAlignmentStart();
		}
	};
	
	private String bam_in = null;
	private String bam_out = null;
	private int batch_size = 4000000;
	
	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}

		switch(args[0].toUpperCase()) {
		case "SORT":
			this.task = Task.sort;
			break;
		default:
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");	
		}

		String[] args2 = new String[args.length-1];
		System.arraycopy(args, 1, args2, 0, args2.length);

		switch(this.task) {
		case sort:
			if (myArgsEngine == null) {
				myArgsEngine = new ArgsEngine();
				myArgsEngine.add("-i", "--in-bam", true);
				myArgsEngine.add("-b", "--barcode", false);
				myArgsEngine.add("-n", "--name", false);
				myArgsEngine.add("-s", "--batch-size", true);
				myArgsEngine.add("-t", "--threads", true);
				myArgsEngine.add("-o", "--out-bam", true);
				myArgsEngine.parse(args2);
			}
			if (myArgsEngine.getBoolean("-i")) {
				this.bam_in = myArgsEngine.getString("-i");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the input BAM file.");
			}
			if (myArgsEngine.getBoolean("-b")&&myArgsEngine.getBoolean("-n")) {
				printUsage();
				throw new IllegalArgumentException("Options -b and -n are exculsive!!!");
			}
			if (myArgsEngine.getBoolean("-b")) {
				this.sort_order = Order.barcode;
				comprator = new Comparator<SAMRecord>() {
					@Override
					public int compare(SAMRecord record1, SAMRecord record2) {
						// TODO Auto-generated method stub
						if(record1==null&&record2==null) return 0;
						if(record1==null) return  1;
						if(record2==null) return -1;
						if(record1.getStringAttribute("BX")==null&&
								record2.getStringAttribute("BX")==null)
							return 0;
						// none barcode record to the end
						if(record1.getStringAttribute("BX")==null) return  1;
						if(record2.getStringAttribute("BX")==null) return -1;
						if(record1.getStringAttribute("BX").
								compareTo(record2.getStringAttribute("BX"))==0)
							return compareSAMRecord(record1, record2);
						return record1.getStringAttribute("BX").
								compareTo(record2.getStringAttribute("BX"));
					}
				};
			}
			if (myArgsEngine.getBoolean("-n")) {
				this.sort_order = Order.queryname;
				comprator = new Comparator<SAMRecord>() {
					@Override
					public int compare(SAMRecord record1, SAMRecord record2) {
						// TODO Auto-generated method stub
						if(record1==null&&record2==null) return 0;
						if(record1==null) return  1;
						if(record2==null) return -1;
						if(record1.getReadName().compareTo(record2.getReadName())==0)
							return compareSAMRecord(record1, record2);
						return record1.getReadName().compareTo(record2.getReadName());
					}
				};
			}
			if (myArgsEngine.getBoolean("-s")) {
				this.batch_size = Integer.parseInt(myArgsEngine.getString("-s"));
			}
			if (myArgsEngine.getBoolean("-t")) {
				this.THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
			}
			if (myArgsEngine.getBoolean("-o")) {
				this.bam_out = myArgsEngine.getString("-o");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the output BAM file.");
			}
			break;
		default:
			throw new RuntimeException("!!!");	
		}
	}

	protected int compareSAMRecord(SAMRecord record1, SAMRecord record2) {
		// TODO Auto-generated method stub
		return record1.getSAMString().compareTo(record2.getSAMString());
		/**
		if(!record1.getReadName().equals(record2.getReadName()))
			return record1.getReadName().compareTo(record2.getReadName());
		if(record1.getFirstOfPairFlag()&&record2.getSecondOfPairFlag())
			return -1;
		if(record1.getSecondOfPairFlag()&&record2.getFirstOfPairFlag())
			return  1;
		if(!record1.getNotPrimaryAlignmentFlag()&&record2.getNotPrimaryAlignmentFlag())
			return -1;
		if(record1.getNotPrimaryAlignmentFlag()&&!record2.getNotPrimaryAlignmentFlag())
			return  1;
		return record1.getAlignmentStart()-record2.getAlignmentStart();
		**/
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		switch(this.task) {
		case zzz:
			myLogger.info("Task list is empty!!!");
			break;
		case sort:
			this.runSort();
			break;
		default:
			throw new RuntimeException("!!!");
		}
		return;
	}
	
	private final static Object lock = new Object();
	private static int batch = 0;
	private static long record_count = 0;
	
	private void runSort() {
		// TODO Auto-generated method stub
		final SamReaderFactory factory =
				SamReaderFactory.makeDefault()
				.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
						SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
				.validationStringency(ValidationStringency.SILENT);
		
		final SamReader inputSam = factory.open(new File(this.bam_in));
		final SAMFileHeader sort_header = inputSam.getFileHeader();
		switch(this.sort_order) {
		case coordinate:
			sort_header.setSortOrder(SortOrder.coordinate);
			break;
		case queryname:
			sort_header.setSortOrder(SortOrder.queryname);
			break;
		case barcode:
			sort_header.setSortOrder(SortOrder.unknown);
			break;
		}
		SAMRecordIterator iter=inputSam.iterator();
		long record_inCount = 0;
		
		SAMRecord[] buff = new SAMRecord[this.batch_size];
		int k = 0;
		SAMRecord temp = iter.hasNext() ? iter.next() : null;
		
		this.initial_thread_pool();
		
		while( temp!=null ) {
			buff[k++] = temp;
			record_inCount++;
			temp = iter.hasNext() ? iter.next() : null;

			if(k==this.batch_size || temp==null) {
				executor.submit(new Runnable() {
					private SAMRecord[] records;

					@Override
					public void run() {
						// TODO Auto-generated method stub
						try {
							Arrays.sort(records, comprator);
							final SAMFileWriter outputSam;
							synchronized(lock) {
								outputSam =  new SAMFileWriterFactory().
										makeSAMOrBAMWriter( sort_header,
												true, new File(bam_out+String.format("%08d", batch++)) );
							}
							int count = 0;
							for(SAMRecord record : records) {
								if(record!=null) {
									count++;
									outputSam.addAlignment(record);
								}
							}
							outputSam.close();
							synchronized(lock) {
								record_count += count;
							}
							myLogger.info("["+Thread.currentThread().getName()+"] "+record_count+" records processed.");
						} catch (Exception e) {
							Thread t = Thread.currentThread();
							t.getUncaughtExceptionHandler().uncaughtException(t, e);
							e.printStackTrace();
							executor.shutdown();
							System.exit(1);
						}
					}

					public Runnable init(SAMRecord[] buff) {
						// TODO Auto-generated method stub
						this.records = buff;
						return(this);
					}

				}.init(buff));

				k=0;
				buff = new SAMRecord[this.batch_size];
			}
		}
		iter.close();

		try {
			inputSam.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		myLogger.info(record_inCount+" records read from "+this.bam_in);
		
		this.waitFor();
		
		// merge all batches
		myLogger.info("Merge "+batch+" files.");
		final SAMFileWriter outputSam =  new SAMFileWriterFactory().
				makeSAMOrBAMWriter( sort_header, true, new File(this.bam_out) );

		final SamReader[] batchSam = new SamReader[batch];
		final SAMRecordIterator[] iterSam = new SAMRecordIterator[batch];
		final boolean[] reachFileEnd = new boolean[batch];
		final TreeMap<SAMRecord, Integer> treeMap = new TreeMap<SAMRecord, Integer>(this.comprator);
		for(int i=0; i!=batch; i++) { 
			batchSam[i] = factory.open(new File(this.bam_out+String.format("%08d", i)));
			iterSam[i] = batchSam[i].iterator();
			if(iterSam[i].hasNext()) treeMap.put(iterSam[i].next(), i);
			reachFileEnd[i] = !iterSam[i].hasNext();
		}
		
		Entry<SAMRecord, Integer> firstEntry;
		int bch, nReachFileEnd=0;
		for(boolean b : reachFileEnd) if(b) nReachFileEnd++;
		long record_outCount = 0;
		
		while( !treeMap.isEmpty() ) {
			firstEntry = treeMap.pollFirstEntry();
			outputSam.addAlignment(firstEntry.getKey());
			record_outCount++;
			bch = firstEntry.getValue();
			
			if(!reachFileEnd[bch]) {				
				treeMap.put(iterSam[bch].next(), bch);
				if(!iterSam[bch].hasNext()) {
					reachFileEnd[bch] = true;
					nReachFileEnd++;
				}
			}
			
			if(treeMap.isEmpty()&&nReachFileEnd!=batch) {
				for(int i=0; i!=batch; i++) { 
					if(!reachFileEnd[i]) {
						treeMap.put(iterSam[i].next(), i);
						if(!iterSam[i].hasNext()) {
							reachFileEnd[i] = true;
							nReachFileEnd++;
						}
					}
				}
			}
		}
		try {
			outputSam.close();
			for(int i=0; i!=batch; i++) {
				iterSam[i].close();
				batchSam[i].close();
				// new File(this.bam_out+String.format("%08d", i)).delete();
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		myLogger.info(record_outCount+" records written to "+this.bam_out);
	}
}
