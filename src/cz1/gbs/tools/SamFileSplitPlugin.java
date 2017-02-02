package cz1.gbs.tools;

import cz1.util.ArgsEngine;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
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


public class SamFileSplitPlugin {

	private final static String file_sep  = 
			System.getProperty("file.separator");
	private ArgsEngine myArgsEngine = null;
	private final static Logger myLogger = 
			Logger.getLogger(SamFileSplitPlugin.class);
	static {
		BasicConfigurator.configure();
	}
	
	public static void main(String[] args) {
		SamFileSplitPlugin ssp = new SamFileSplitPlugin();
		ssp.setParameters(args);
		ssp.split();
	}

	private void printUsage() {
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i  input bam or sam file, should be sorted.\n"
						+ " -b	bed files directory (currently only support splitting by scaffolds).\n"
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
			myArgsEngine.add("-o", "--output-dir", true);
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

	public void split() {
		
		final SAMFileReader inputSam = new SAMFileReader(new File(bam_in));
		final SAMFileHeader header = inputSam.getFileHeader();
		inputSam.setValidationStringency(ValidationStringency.SILENT);
		SAMRecordIterator iter=inputSam.iterator();
		final SAMSequenceDictionary seqdic = header.getSequenceDictionary();
		
		File[] bed_files = new File(bed_in).listFiles();
		final SAMFileWriter[] outputSam = new SAMFileWriter[bed_files.length];
		Map<String, Integer> outMap = new HashMap<String, Integer>();
		for(int i=0; i<bed_files.length; i++) {
			
			Set<String> bed_seq = new HashSet<String>();
			String out = bed_files[i].getName().split("\\.")[0];
			String out_2 = new File(bam_in).getName().split("\\.")[0];
			String tmp;
			
			try (BufferedReader br = new BufferedReader(
					new FileReader(bed_files[i]))) {
				String line;
				while( (line=br.readLine()) !=null ) {
					tmp = line.split("\\s+")[0];
					bed_seq.add(tmp);
					outMap.put(tmp, i);
				}
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(1);
			}
			final SAMFileHeader header_i = new SAMFileHeader();
			final SAMSequenceDictionary seqdic_i = new SAMSequenceDictionary();
			header_i.setAttribute("VN", header.getAttribute("VN"));
			header_i.setAttribute("SO", header.getAttribute("SO"));
			List<SAMSequenceRecord> seqs = seqdic.getSequences();
			for(SAMSequenceRecord seq : seqs)
				if(bed_seq.contains(seq.getSequenceName()))
					seqdic_i.addSequence(seq);
			header_i.setSequenceDictionary(seqdic_i);
			for(SAMReadGroupRecord rg : header.getReadGroups())
				header_i.addReadGroup(rg);
			for(SAMProgramRecord pg : header.getProgramRecords())
				header_i.addProgramRecord(pg);
			
			outputSam[i] = new SAMFileWriterFactory().
					makeSAMOrBAMWriter(header_i,
							true, new File(bam_out+file_sep+
									out+file_sep+
									out_2+".sam"));
		}
		
		while(iter.hasNext()) {
			SAMRecord rec=iter.next();
			outputSam[outMap.get(rec.getReferenceName())].addAlignment(rec);
		}
		iter.close();
		inputSam.close();
		for(int i=0; i<outputSam.length; i++)
			outputSam[i].close();
		
		System.err.println(bam_in+" return true");
	}
}

