package cz1.test;

import java.io.File;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecord.SAMTagAndValue;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader.ValidationStringency;

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
		final SAMFileReader inputSam = new SAMFileReader(new File(mySamFile));
		
		
		samHeader = inputSam.getFileHeader();
		samHeader.setSortOrder(SortOrder.unsorted);
		inputSam.setValidationStringency(ValidationStringency.SILENT);
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
		
		inputSam.close();
		outSam.close();
	}
}
