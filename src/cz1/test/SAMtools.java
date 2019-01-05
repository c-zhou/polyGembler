package cz1.test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

import cz1.ngs.model.Sequence;
import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class SAMtools extends Executor {

	private String mySamFile;
	private String myOutput;
	private String chrId;
	private int chrStt;
	private int chrEnd;
	
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
			myArgsEngine.add("-r", "--range", true);
			myArgsEngine.add("-o", "--output", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-i")) {
			mySamFile = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your SAM/BAM file.");
		}
		
		if(myArgsEngine.getBoolean("-r")) {
			String range = myArgsEngine.getString("-r");
			String[] s = range.split(":");
			chrId = s[0];
			s = s[1].split("-");
			chrStt = Integer.parseInt(s[0]);
			chrEnd = Integer.parseInt(s[1]);
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your data range.");
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
		if(!inputSam.hasIndex()) throw new RuntimeException("BAM file need to be indexed!!!");
		final SAMRecordIterator iter = inputSam.query(chrId, chrStt, chrEnd, false);
		final BufferedWriter bw = Utils.getBufferedWriter(myOutput);
		try {
			SAMRecord record;
			int readStt, readEnd;
			String readStr;
			while(iter.hasNext()) {
				record = iter.next();
				if(record.getAlignmentStart()>chrStt || 
						record.getAlignmentEnd()<chrEnd ||
						record.isSecondaryAlignment() ||
						record.getSupplementaryAlignmentFlag())
					continue;
				readStr = record.getReadString();
				readStt = record.getReadPositionAtReferencePosition(chrStt, true);
				readEnd = record.getReadPositionAtReferencePosition(chrEnd, true);
				bw.write(Sequence.formatOutput(record.getReadName(), readStr.substring(readStt-1, readEnd), 100));
			}
			iter.close();
			inputSam.close();
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
