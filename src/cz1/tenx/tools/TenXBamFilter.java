package cz1.tenx.tools;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader.ValidationStringency;

public class TenXBamFilter extends Executor {

	private String in_bam;
	private String out_bam;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--input-bam      Input directory containing tag sequence files.\n"
						+ " -t/--threads        Threads (default is 1).\n"
						+ " -o/--output-bam     Output directory.\n\n");
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}

		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-i", "--input-bam", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-o", "--output-bam", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-i")) {
			in_bam = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your FASTQ files.");
		}

		if (myArgsEngine.getBoolean("-t")) {
			THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}
		
		if (myArgsEngine.getBoolean("-o")) {
			out_bam = myArgsEngine.getString("-o");
		}
	}
	
	SAMFileWriter outputSam;
	private final List<Long> synchedBlock = 
			Collections.synchronizedList(new LinkedList<Long>());
	
	public void write(final long i, final SAMRecord[] records) 
			throws InterruptedException, IOException {
		synchronized (synchedBlock) {
			while (synchedBlock.get(0)!=i)
				synchedBlock.wait();
			int b = records.length;
			for(int k=0; k!=b; k++) {
				if(records[k]!=null)
					outputSam.addAlignment(records[k]);
			}
			synchedBlock.remove(0);
			synchedBlock.notifyAll();
		}
	}
	
	@Override
	public void run() {
		// TODO Auto-generated method stub

		final SAMFileReader inputSam = new SAMFileReader(new File(in_bam));
		inputSam.setValidationStringency(ValidationStringency.SILENT);
		outputSam = new SAMFileWriterFactory().
				makeSAMOrBAMWriter(inputSam.getFileHeader(), true, new File(out_bam));

		this.initial_thread_pool();

		final int block = 1000000;
		SAMRecord[] Qs = new SAMRecord[block];
		int k = 0;
		long bS = 0;
		SAMRecordIterator iter=inputSam.iterator();
		while(iter.hasNext()) {
			Qs[k++] = iter.next();

			if(k==block || !iter.hasNext()) {

				synchronized (synchedBlock) {synchedBlock.add(bS);}

				executor.submit(new Runnable() {
					private SAMRecord[] records;
					private long block_i;
					@Override
					public void run() {
						// TODO Auto-generated method stub
						try {
							final long start_i = block_i*block;
							SAMRecord record;
							Cigar cigar;
							List<CigarElement> cigar_ele;
							CigarElement e;
							int nc; 
							double l;
							int[] stats = new int[4]; // M S D I
							
							for(int i=0; i<records.length; i++) {
								if( (record=records[i])==null ) break;
								if(record.getReadUnmappedFlag() ||
										record.getMappingQuality()<30 ||
										record.getNotPrimaryAlignmentFlag() ||
										record.getSupplementaryAlignmentFlag() ||
										record.getDuplicateReadFlag() ) {
									records[i] = null;
									continue;
								}
								Arrays.fill(stats, 0);
								cigar = record.getCigar();
								l = (double) cigar.getReadLength();
								cigar_ele = cigar.getCigarElements();
								nc = cigar_ele.size();
								for(int j=0; j<nc; j++) {
									e = cigar_ele.get(j);
									switch( CigarOperator.enumToCharacter(e.getOperator()) ) {
									case 'M':
										stats[0] += e.getLength();
										break;
									case 'S':
										stats[1] += e.getLength();
										break;
									case 'D':
										stats[2] += e.getLength();
										break;
									case 'I':
										stats[3] += e.getLength();
										break;
									default:
									}
								}
								
								if(stats[0]/l<.95 || stats[1]/l>.05 || stats[2]>3 || stats[3]>3) {
									records[i] = null;
									continue;
								}
							}
							write(block_i, records);
							if(start_i%1000000==0)
								myLogger.info(start_i+" record processed");
						} catch (Exception e) {
							Thread t = Thread.currentThread();
							t.getUncaughtExceptionHandler().uncaughtException(t, e);
							e.printStackTrace();
							executor.shutdown();
							System.exit(1);
						}
					}

					public Runnable init(final SAMRecord[] Qs, final long bS) {
						this.records = Qs;
						this.block_i = bS;
						return(this);
					}
				}.init(Qs, bS));

				bS++;
				k = 0;
				Qs = new SAMRecord[block];
			}
		}
		iter.close();
		inputSam.close();
		this.waitFor();
		outputSam.close();
	}
}
