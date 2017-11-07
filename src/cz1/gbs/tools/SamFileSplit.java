package cz1.gbs.tools;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;
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
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeUnit;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;


public class SamFileSplit extends Executor {
	
	public SamFileSplit(String bam_in, 
			String bed_in, 
			String bam_out,
			int threads) {
		this.bam_in = bam_in;
		this.bed_in = bed_in;
		this.bam_out = bam_out;
		this.THREADS = threads;
	}
	
	public SamFileSplit() {}
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--input-dir  	Input bam or sam file, should be sorted.\n"
						+ " -b/--bed-dir		Bed files directory (currently only support splitting by scaffolds).\n"
						+ " -t/--threads		Threads.\n"
						+ " -o--prefix  		Output prefix.\n\n" );
	}

	private String bam_in;
	private String bed_in;
	private String bam_out;
	
	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}

		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-i", "--input-dir", true);
			myArgsEngine.add("-b", "--bed-dir", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-o", "--prefix", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-i")) {
			bam_in = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the location of your BAM files.");
		}

		if (myArgsEngine.getBoolean("-b")) {
			bed_in = myArgsEngine.getString("-b");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the location of your BED files.");
		}
		
		if (myArgsEngine.getBoolean("-t")) {
			THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}

		if (myArgsEngine.getBoolean("-o")) {
			bam_out = myArgsEngine.getString("-o");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the location of your output BAM files.");
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		Utils.makeOutputDir(bam_out);
		final File[] beds = new File(bed_in).listFiles();
		final String[] out_prefix = new String[beds.length];
		for(int i=0; i<beds.length; i++) {
			out_prefix[i] = bam_out+"/"+beds[i].getName().replaceAll(".bed$", "");
			Utils.makeOutputDir(out_prefix[i]);
		}
		
		final File[] bams = new File(bam_in).listFiles(new FilenameFilter() {
		    @Override
		    public boolean accept(File dir, String name) {
		        return name.endsWith(".bam");
		    }
		});
		this.initial_thread_pool();
		for(File bam : bams) {
			executor.submit(new Runnable(){
				private File bam;
				@Override
				public void run() {
					// TODO Auto-generated method stub
					try {
						
						final SamReaderFactory factory =
								SamReaderFactory.makeDefault()
								.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
										SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
								.validationStringency(ValidationStringency.SILENT);
						
						final SamReader inputSam = factory.open(bam);
						final SAMFileHeader header = inputSam.getFileHeader();
						final SAMRecordIterator iter = inputSam.iterator();
						final SAMSequenceDictionary seqdic = header.getSequenceDictionary();

						final SAMFileWriter[] outputSam = new SAMFileWriter[beds.length];
						final SAMSequenceDictionary[] seqdics = new SAMSequenceDictionary[beds.length];
						final Map<String, Integer> outMap = new HashMap<String, Integer>();
						final String out = bam.getName();
						
						for(int i=0; i<beds.length; i++) {

							Set<String> bed_seq = new HashSet<String>();
							String tmp;

							BufferedReader br = new BufferedReader(
									new FileReader(beds[i]));
							String line;
							while( (line=br.readLine()) !=null ) {
								tmp = line.split("\\s+")[0];
								bed_seq.add(tmp);
								outMap.put(tmp, i);
							}
							br.close();
							
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
											true, new File(out_prefix[i]+"/"+out));
							seqdics[i] = seqdic_i;
						}

						Set<String> refs = outMap.keySet();
						String ref;
						int f;
						while(iter.hasNext()) {
							SAMRecord rec=iter.next();
							if( refs.contains(ref=rec.getReferenceName()) ) {
								f = outMap.get(ref);
								rec.setReferenceIndex(seqdics[f].getSequenceIndex(ref));
								outputSam[f].addAlignment(rec);
							}
						}
						iter.close();
						inputSam.close();
						for(int i=0; i<outputSam.length; i++)
							outputSam[i].close();

						myLogger.info(out+" return true");
					} catch (Exception e) {
						Thread t = Thread.currentThread();
						t.getUncaughtExceptionHandler().uncaughtException(t, e);
						e.printStackTrace();
						executor.shutdown();
						System.exit(1);
					}
				}
				
				public Runnable init(File bam) {
					this.bam = bam;
					return(this);
				}
			}.init(bam));
		}
		
		this.waitFor();
	}
}
