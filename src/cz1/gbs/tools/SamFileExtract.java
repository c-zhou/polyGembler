package cz1.gbs.tools;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;


public class SamFileExtract extends Executor {

	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i  input bam or sam file, should be sorted.\n"
						+ " -b	bed file.\n"
						+ " -o  output directory.\n\n" );
	}

	private static String bam_in;
	private static String bed_in;
	private static String bam_out;
	
	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}

		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-i", "--input-file", true);
			myArgsEngine.add("-b", "--bed-dir", true);
			myArgsEngine.add("-o", "--output-file", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-i")) {
			bam_in = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the location of your FASTQ files.");
		}

		if (myArgsEngine.getBoolean("-b")) {
			bed_in = myArgsEngine.getString("-b");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify a barcode key file.");
		}

		if (myArgsEngine.getBoolean("-o")) {
			bam_out = myArgsEngine.getString("-o");
		} else {
			myLogger.warn("No enzyme specified.  Using enzyme listed in key file.");
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
		final SamReader inputSam = factory.open(new File(bam_in));
		final SAMFileHeader header = inputSam.getFileHeader();
		final SAMSequenceDictionary seqdic = header.getSequenceDictionary();
		final SAMFileHeader header_out = new SAMFileHeader();
		final SAMSequenceDictionary seqdic_out = new SAMSequenceDictionary();
		
		SAMRecordIterator iter=inputSam.iterator();

		File bed_file = new File(bed_in);
		
		final Set<String> extract = new HashSet<String>();
		try (BufferedReader br = new BufferedReader(
				new FileReader(bed_file))) {
			String line;
			while( (line=br.readLine()) !=null )
				extract.add(line.split("\\s+")[0]);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		header_out.setAttribute("VN", header.getAttribute("VN"));
		header_out.setAttribute("SO", header.getAttribute("SO"));
		List<SAMSequenceRecord> seqs = seqdic.getSequences();
		for(SAMSequenceRecord seq : seqs)
			if(extract.contains(seq.getSequenceName()))
				seqdic_out.addSequence(seq);
		header_out.setSequenceDictionary(seqdic_out);
		for(SAMReadGroupRecord rg : header.getReadGroups())
			header_out.addReadGroup(rg);
		for(SAMProgramRecord pg : header.getProgramRecords())
			header_out.addProgramRecord(pg);
		
		final SAMFileWriter	outputSam = new SAMFileWriterFactory().
				makeSAMOrBAMWriter(header_out,
						true, new File(bam_out));
		while(iter.hasNext()) {
			SAMRecord rec=iter.next();
			if(extract.contains(rec.getReferenceName())) 
				outputSam.addAlignment(rec);
		}
		iter.close();
		try {
			inputSam.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		outputSam.close();
		
		System.err.println(bam_in+" return true");
	}
	
}

