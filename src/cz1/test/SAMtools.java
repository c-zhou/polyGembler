package cz1.test;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecord.SAMTagAndValue;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class SAMtools extends Executor {

	private String mySamFile;
	private String myOutput;
	private SAMFileHeader samHeader;
	
	public static void main(String[] args) {
		SAMtools samTools = new SAMtools();
		samTools.setParameters(args);
		samTools.run();
	}

	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-i", "--sam-file", true);
			myArgsEngine.add("-o", "--prefix", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-i")) {
			mySamFile = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your SAM/BAM file.");
		}
		
		if (myArgsEngine.getBoolean("-o")) {
			myOutput = myArgsEngine.getString("-o");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your output file.");
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
		final SamReader inputSam = factory.open(new File(mySamFile));
		
		samHeader = inputSam.getFileHeader();
		samHeader.setSortOrder(SortOrder.unsorted);
		SAMRecordIterator iter=inputSam.iterator();
		Set<Entry<String, String>> attr = samHeader.getAttributes();
		
		List<SAMReadGroupRecord> rgs = samHeader.getReadGroups();
		
		SAMReadGroupRecord rg = new SAMReadGroupRecord("cz1");
		rg.setSample("cz1");
		samHeader.addReadGroup(rg);
		
		
		//samHeader.setAttribute("RG", "cz1");
		final SAMFileWriter outSam = new SAMFileWriterFactory().
				makeSAMOrBAMWriter(samHeader,
						true, new File(myOutput));
		
		for(int i=0; i<100; i++) {
			SAMRecord record = iter.next();
			List<SAMTagAndValue> tags = record.getAttributes();
			record.setAttribute("RG", "cz1");
			List<SAMTagAndValue> tags2 = record.getAttributes();
			outSam.addAlignment(record);
		}
		myLogger.info("exit...");
		
		try {
			inputSam.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		outSam.close();
	}
}
