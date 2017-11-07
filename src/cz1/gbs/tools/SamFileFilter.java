package cz1.gbs.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class SamFileFilter extends Executor {
	private String input_dir;
	private String filter_file;
	private String output_dir;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--input-dir      Input sam/bam file.\n"
						+ " -x/--filter         File indicate positions to filter out.\n"
						+ " -t/--threads        Threads (default is 1).\n"
						+ " -o/--prefix         Output file directory. \n\n");
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
			myArgsEngine.add("-i", "--input-dir", true);
			myArgsEngine.add("-x", "-x/--filter", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-o", "--prefix", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-i")) {
			input_dir = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the folder of coverage files.");
		}

		if (myArgsEngine.getBoolean("-o")) {
			output_dir = myArgsEngine.getString("-o");
		}
		
		if (myArgsEngine.getBoolean("-t")) {
			THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}
		
		if (myArgsEngine.getBoolean("-x")) {
			filter_file = myArgsEngine.getString("-x");
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		
		File folder = new File(input_dir);
		File[] listOfFiles = folder.listFiles(new FilenameFilter() {
		    @Override
		    public boolean accept(File dir, String name) {
		        return name.endsWith(".bam");
		    }
		});
		int bam_n = listOfFiles.length;
		
		this.initial_thread_pool();
		
		for(int i=0; i<bam_n; i++) {
			executor.submit(new Runnable() {
				private String bam_file;
				
				@Override
				public void run() {
					// TODO Auto-generated method stub
					try{
						final SamReaderFactory factory =
								SamReaderFactory.makeDefault()
								.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
										SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
								.validationStringency(ValidationStringency.SILENT);
						
						final SamReader inputSam = factory.open(new File(input_dir+"/"+bam_file));

						final SAMFileWriter outputSam =  new SAMFileWriterFactory().
								makeSAMOrBAMWriter(inputSam.getFileHeader(),
										true, new File(output_dir+"/"+bam_file));
						final BufferedReader filter_br = Utils.getBufferedReader(filter_file);
						SAMRecordIterator iter=inputSam.iterator();
						
						String line = filter_br.readLine();
						if(line==null) {
							filter_br.close();
							iter.close();
							inputSam.close();
							outputSam.close();
							return;
						}
						String[] s = line.split("\\s+");
						String chr = s[0];
						long end_pos = Long.parseLong(s[1]);
						SAMRecord samr;
						while( iter.hasNext() ) {
							samr = iter.next();
							if( samr.getReferenceName().equals(chr) &&
									samr.getAlignmentEnd()<end_pos ) {
								outputSam.addAlignment(samr);
								continue;
							}
							if( !samr.getReferenceName().equals(chr) ||
									samr.getAlignmentStart()>end_pos) {
								line = null;
								while( (line=filter_br.readLine())!=null &&
										end_pos<samr.getAlignmentStart() ) {
									s = line.split("\\s");
									chr = s[0];
									end_pos = Long.parseLong(s[1]);
								}
								if(line!=null) {
									if( samr.getReferenceName().equals(chr) &&
											samr.getAlignmentEnd()<end_pos )
										outputSam.addAlignment(samr);
								} else {
									outputSam.addAlignment(samr);
									while(iter.hasNext())
										outputSam.addAlignment(iter.next());
									break;
								}
							}
						}
						filter_br.close();
						iter.close();
						inputSam.close();
						outputSam.close();
					} catch (Exception e) {
						Thread t = Thread.currentThread();
						t.getUncaughtExceptionHandler().uncaughtException(t, e);
						e.printStackTrace();
						executor.shutdown();
						System.exit(1);
					}
				}

				public Runnable init(String bam_file) {
					// TODO Auto-generated method stub
					this.bam_file = bam_file;
					return this;
				}
				
			}.init(listOfFiles[i].getName()));
			
			throw new RuntimeException("incorrect tools!!!!");
			
		}
		
		this.waitFor();
	}
}
