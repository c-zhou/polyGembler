package cz1.tenx.tools;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class TenXMoleculeStatsZ extends Executor {
	
	private String in_bam;
	private String out_bam;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--input-bam      Input BAM file.\n"
						+ " -o/--output-bam     Output BAM file.\n\n");
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
			myArgsEngine.add("-o", "--output-bam", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-i")) {
			in_bam = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your BAM file.");
		}
		
		if (myArgsEngine.getBoolean("-o")) {
			out_bam = myArgsEngine.getString("-o");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your BAM file.");
		}
	}
	
	private final SamReaderFactory factory =
			SamReaderFactory.makeDefault()
			.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
					SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
			.validationStringency(ValidationStringency.SILENT);
	private final static double min_qual = 20;
	private final static double max_ins  = 1000;
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		
		try {
			final SamReader inputSam = factory.open(new File(this.in_bam));

			final SAMFileWriter outputSam = new SAMFileWriterFactory().
					makeSAMOrBAMWriter(inputSam.getFileHeader(), true, new File(out_bam));

			final SAMRecordIterator iter = inputSam.iterator();
			SAMRecord record = iter.hasNext() ? iter.next() : null;
			final SAMRecord[] records = new SAMRecord[2];
			String sn;
			
			while( record!=null ) {
				sn = record.getReadName();
				
				if( !record.getReadUnmappedFlag() &&
						!record.getNotPrimaryAlignmentFlag()&&
						!record.getSupplementaryAlignmentFlag() ) {
					if(record.getFirstOfPairFlag())
						records[0] = record;
					else if(record.getSecondOfPairFlag())
						records[1] = record;
					else
						throw new RuntimeException("!!!");
				}
				
				while( (record = iter.hasNext() ? iter.next() : null)!=null && record.getReadName().equals(sn)) {
					if( !record.getReadUnmappedFlag() &&
							!record.getNotPrimaryAlignmentFlag()&&
							!record.getSupplementaryAlignmentFlag() ) {
						if(record.getFirstOfPairFlag())
							records[0] = record;
						else if(record.getSecondOfPairFlag())
							records[1] = record;
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
						outputSam.addAlignment(records[0]);
						outputSam.addAlignment(records[1]);
					}
				}
				
				Arrays.fill(records, null);
			}
			
			iter.close();
			inputSam.close();
			outputSam.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
	}
}
