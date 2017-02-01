package cz1.appl;

import cz1.util.ArgsEngine;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMProgramRecord;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;


public class SamFileExtractPlugin {

	private final static String file_sep  = 
			System.getProperty("file.separator");
	private ArgsEngine myArgsEngine = null;
	private final static Logger myLogger = 
			Logger.getLogger(SamFileExtractPlugin.class);
	static {
		BasicConfigurator.configure();
	}
	
	public static void main(String[] args) {
		SamFileExtractPlugin ssp = new SamFileExtractPlugin();
		ssp.setParameters(args);
		ssp.extract();
	}

	private void printUsage() {
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i  input bam or sam file, should be sorted.\n"
						+ " -b	bed file.\n"
						+ " -o  output directory.\n\n" );
	}

	private static String bam_in;
	private static String bed_in;
	private static String bam_out;
	
	public void setParameters(String[] args) {
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

	public void extract() {
		
		final SAMFileReader inputSam = new SAMFileReader(new File(bam_in));
		inputSam.setValidationStringency(ValidationStringency.SILENT);
		final SAMFileHeader header = inputSam.getFileHeader();
		final SAMSequenceDictionary seqdic = header.getSequenceDictionary();
		final SAMFileHeader header_out = new SAMFileHeader();
		final SAMSequenceDictionary seqdic_out = new SAMSequenceDictionary();
		
		inputSam.setValidationStringency(ValidationStringency.SILENT);
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
		inputSam.close();
		outputSam.close();
		
		System.err.println(bam_in+" return true");
	}
	
}

