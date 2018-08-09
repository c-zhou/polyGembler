package cz1.ngs.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;

import cz1.ngs.model.Sequence;
import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Executor;
import cz1.util.Utils;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Scaffolder extends Executor {
	
	private final static boolean USE_OS_BUFFER = false;
	// 8Mb buffer size
	private final static int buffSize8Mb = 8388608;
	private final static Writer STD_OUT_BUFFER = USE_OS_BUFFER ? 
			new BufferedWriter(new OutputStreamWriter(System.out), buffSize8Mb) : new OutputStreamWriter(System.out);
			
	private String[] bamList;
	private int num_threads = Runtime.getRuntime().availableProcessors();
	private boolean debug  = false;
	private boolean ddebug = false;
	private String out_prefix = null;
	private int minQual = 20;
	private int inst;
	private int maxInst;
	private String subject_file;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -a/--align              Alignment file(s). Multiple file are separated by ':'. \n"
						+ " -aL/--align-list        Alignment file list.\n"
						+ " -s/--subject            The FASTA file contain subject/reference sequences. \n"
						+ " -q/--min-qual           Minimum alignment quality (default 20).\n"
						+ " -f/--frag-size          Insert size of the library.\n"
						+ " -t/--threads            Number of threads to use (default 16).\n"
						+ " -d/--debug              Debugging mode will have extra information printed out.\n"
						+ " -dd/--debug-debug       Debugging mode will have more information printed out than -d mode.\n"
						+ " -o/--out-prefix         Prefix of the output files.\n"
						+ "\n");	
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-a", "--align", true);
			myArgsEngine.add("-aL", "--align-list", true);
			myArgsEngine.add("-s", "--subject", true);
			myArgsEngine.add("-q", "--min-qual", true);
			myArgsEngine.add("-f", "--frag-size", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-d", "--debug", false);
			myArgsEngine.add("-dd", "--debug-debug", false);
			myArgsEngine.add("-o", "--out-prefix", true);
			myArgsEngine.parse(args);
		}

		if (!myArgsEngine.getBoolean("-a")&&!myArgsEngine.getBoolean("-aL")) {
			printUsage();
			throw new IllegalArgumentException("Please specify the alignment file(s) using -a or -aL option.");
		}
		
		if (myArgsEngine.getBoolean("-a")&&myArgsEngine.getBoolean("-aL")) {
			printUsage();
			throw new IllegalArgumentException("Options -a and -aL are exclusive.");
		}
		
		if (myArgsEngine.getBoolean("-a")) {
			this.bamList = myArgsEngine.getString("-a").trim().split(":");
		}
		
		if (myArgsEngine.getBoolean("-aL")) {
			try {
				BufferedReader br = Utils.getBufferedReader(myArgsEngine.getString("-aL").trim());
				final List<String> file_list = new ArrayList<String>();
				String line;
				while( (line = br.readLine()) != null) {
					line = line.trim();
					if(line.length()>0) file_list.add(line);
				}
				this.bamList = file_list.toArray(new String[file_list.size()]);
				br.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		}
		
		if (myArgsEngine.getBoolean("-s")) {
			this.subject_file = myArgsEngine.getString("-s");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the contig file.");
		}
		
		if (myArgsEngine.getBoolean("-q")) {
			this.minQual = Integer.parseInt(myArgsEngine.getString("-q"));
		}

		if (myArgsEngine.getBoolean("-f")) {
			this.inst = Integer.parseInt(myArgsEngine.getString("-f"));
			this.maxInst = this.inst*2;
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify insert size of the library.");
		}

		
		if (myArgsEngine.getBoolean("-t")) {
			int t = Integer.parseInt(myArgsEngine.getString("-t"));
			if(t<this.num_threads) this.num_threads = t;
			this.THREADS = t;
			Constants.omp_threads = this.num_threads;
			myLogger.info("OMP_THREADS = "+this.num_threads);
		}

		if (myArgsEngine.getBoolean("-d")) {
			this.debug = true;
		}

		if (myArgsEngine.getBoolean("-dd")) {
			this.debug  = true;
			this.ddebug = true;
		}

		if (myArgsEngine.getBoolean("-o")) {
			this.out_prefix = myArgsEngine.getString("-o");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the prefix of output files.");
		}
	}

	final SamReaderFactory factory =
			SamReaderFactory.makeDefault()
			.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
					SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
			.validationStringency(ValidationStringency.SILENT);

	private final static Object lock = new Object();

	private final static Map<Long, Integer> linkCount = new HashMap<Long, Integer>();
	private static long exceed_ins = 0;
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		myLogger.info("Reading alignments from "+this.bamList.length+" BAM file"+
				(this.bamList.length>1?"s":"")+":");
		for(String bamfile : this.bamList)
			myLogger.info(bamfile);
		myLogger.info("****");
		
		final Map<String, Sequence> sub_seqs = Sequence.parseFastaFileWithRevCmpAsMap(subject_file);
		final BidiMap<String, Integer> seq_index = new DualHashBidiMap<String, Integer>();
		final Map<String, String> symm_seqsn = new HashMap<String, String>();
		
		int index = 0;
		for(String seq : sub_seqs.keySet()) seq_index.put(seq, ++index);
		for(String seq : sub_seqs.keySet()) {
			if(!seq.endsWith("'")) {
				symm_seqsn.put(seq, seq+"'");
				symm_seqsn.put(seq+"'", seq);
			}
		}
		
		this.initial_thread_pool();
		for(String bamfile : this.bamList)
			this.executor.submit(new Runnable() {
				private String bamfile;

				@Override
				public void run() {
					// TODO Auto-generated method stub
					try {
						final SamReader in1 = factory.open(new File(bamfile));
						final SAMRecordIterator iter1 = in1.iterator();

						myLogger.info("Reading alignments from "+ this.bamfile);


						SAMRecord[] sam_records = new SAMRecord[2];
						SAMRecord sam_record;
						final Map<Long, Integer> linkCount1 = new HashMap<Long, Integer>();
						long exceed_ins1 = 0;
						final boolean[] rev = new boolean[2];
						final int[] reflen = new int[2];
						final String[] refstr = new String[2];
						int a, b;
						long refind;

						while(iter1.hasNext()) {

							sam_record = iter1.next();
							if(sam_record.getNotPrimaryAlignmentFlag() || 
									sam_record.getSupplementaryAlignmentFlag())
								continue;
							if(sam_record.getFirstOfPairFlag()) 
								sam_records[0] = sam_record;
							else sam_records[1] = sam_record;

							if(sam_records[0]!=null&&sam_records[1]!=null) {
								// so we have a confident read pair aligned to two contigs
								// check return
								if(sam_records[0].getReadUnmappedFlag() ||
										sam_records[1].getReadUnmappedFlag()) {
									Arrays.fill(sam_records, null);
									continue;
								}
								if(!sam_records[0].getReadName().
										equals(sam_records[1].getReadName()))
									throw new RuntimeException("!!!");
								if(sam_records[0].getReferenceIndex().intValue()==
										sam_records[1].getReferenceIndex().intValue()) {
									Arrays.fill(sam_records, null);
									continue;
								}
								if(sam_records[0].getMappingQuality()<minQual || 
										sam_records[1].getMappingQuality()<minQual) {
									Arrays.fill(sam_records, null);
									continue;
								}
								
								rev[0] = sam_records[0].getReadNegativeStrandFlag();
								rev[1] = sam_records[1].getReadNegativeStrandFlag();
								reflen[0]  = sub_seqs.get(sam_records[0].getReferenceName()).seq_ln();
								reflen[1]  = sub_seqs.get(sam_records[1].getReferenceName()).seq_ln();

								if(rev[0]&&rev[1]) {
									//       0                  1
									// ---------------     -------------
									//     <===                <===

									// reverse 0
									// ---------------     -------------
									//     ===>                <===
									a = reflen[0]-sam_records[0].getAlignmentStart()+1+sam_records[1].getAlignmentEnd();
									// reverse 1
									// ---------------     -------------
									//     <===                ===>
									b = sam_records[0].getAlignmentEnd()+reflen[1]-sam_records[1].getAlignmentStart()+1;

									if(ddebug)
										STD_OUT_BUFFER.write(">>>>"+Math.min(a, b)+"\n"+
												sam_records[0].getSAMString()+sam_records[1].getSAMString()+"<<<<\n");
									if(a>maxInst&&b>maxInst) {
										++exceed_ins1;
										Arrays.fill(sam_records, null);
										continue;
									}
									if(a<b) {
										// reverse 0
										refstr[0] = sam_records[0].getReferenceName()+"'";
										refstr[1] = sam_records[1].getReferenceName();
									} else {
										// reverse 1
										refstr[0] = sam_records[1].getReferenceName()+"'";
										refstr[1] = sam_records[0].getReferenceName();
									}
								} else if(rev[0]&&!rev[1]) {
									//       0                  1
									// ---------------     -------------
									//     <===                ===>
									a = sam_records[0].getAlignmentEnd()+reflen[1]-sam_records[1].getAlignmentStart()+1;

									// reverse 0 & reverse1
									// ---------------     -------------
									//     ===>                <===
									b = reflen[0]-sam_records[0].getAlignmentStart()+1+sam_records[1].getAlignmentEnd();

									if(ddebug)
										STD_OUT_BUFFER.write(">>>>"+Math.min(a, b)+"\n"+
												sam_records[0].getSAMString()+sam_records[1].getSAMString()+"<<<<\n");
									if(a>maxInst&&b>maxInst) {
										++exceed_ins1;
										Arrays.fill(sam_records, null);
										continue;
									}

									if(a<b) {
										refstr[0] = sam_records[1].getReferenceName();
										refstr[1] = sam_records[0].getReferenceName();
									} else {
										refstr[0] = sam_records[0].getReferenceName()+"'";
										refstr[1] = sam_records[1].getReferenceName()+"'";
									}

								} else if(!rev[0]&&rev[1]) {
									//       0                  1
									// ---------------     -------------
									//     ===>                <===
									a = reflen[0]-sam_records[0].getAlignmentStart()+1+sam_records[1].getAlignmentEnd();

									// reverse 0 & reverse1
									// ---------------     -------------
									//     <===                ===>
									b = sam_records[0].getAlignmentEnd()+reflen[1]-sam_records[1].getAlignmentStart()+1;

									if(ddebug)
										STD_OUT_BUFFER.write(">>>>"+Math.min(a, b)+"\n"+
												sam_records[0].getSAMString()+sam_records[1].getSAMString()+"<<<<\n");
									if(a>maxInst&&b>maxInst) {
										++exceed_ins1;
										Arrays.fill(sam_records, null);
										continue;
									}

									if(a<b) {
										refstr[0] = sam_records[0].getReferenceName();
										refstr[1] = sam_records[1].getReferenceName();
									} else {
										refstr[0] = sam_records[1].getReferenceName()+"'";
										refstr[1] = sam_records[0].getReferenceName()+"'";
									}

								} else {
									//       0                  1
									// ---------------     -------------
									//     ===>                ===>

									// reverse 0
									// ---------------     -------------
									//     <===                ===>
									a = sam_records[0].getAlignmentEnd()+reflen[1]-sam_records[1].getAlignmentStart()+1;

									// reverse 1
									// ---------------     -------------
									//     ===>                <===
									b = reflen[0]-sam_records[0].getAlignmentStart()+1+sam_records[1].getAlignmentEnd();

									if(ddebug)
										STD_OUT_BUFFER.write(">>>>"+Math.min(a, b)+"\n"+
												sam_records[0].getSAMString()+sam_records[1].getSAMString()+"<<<<\n");
									if(a>maxInst&&b>maxInst) {
										++exceed_ins1;
										Arrays.fill(sam_records, null);
										continue;
									}

									if(a<b) {
										refstr[0] = sam_records[1].getReferenceName();
										refstr[1] = sam_records[0].getReferenceName()+"'";
									} else {
										refstr[0] = sam_records[0].getReferenceName();
										refstr[1] = sam_records[1].getReferenceName()+"'";
									}
								}

								refind  = seq_index.get(refstr[0]);
								refind <<= 32;
								refind += seq_index.get(refstr[1]);
								if(linkCount1.containsKey(refind)) {
									linkCount1.put(refind, linkCount1.get(refind)+1);
								} else {
									linkCount1.put(refind, 1);
								}

								refind  = seq_index.get(symm_seqsn.get(refstr[1]));
								refind <<= 32;
								refind += seq_index.get(symm_seqsn.get(refstr[0]));
								if(linkCount1.containsKey(refind)) {
									linkCount1.put(refind, linkCount1.get(refind)+1);
								} else {
									linkCount1.put(refind, 1);
								}
								Arrays.fill(sam_records, null);
							}
						}

						iter1.close();
						in1.close();

						synchronized(lock) {
							for(long key : linkCount1.keySet()) {
								if(linkCount.containsKey(key)) {
									linkCount.put(key, linkCount.get(key)+linkCount1.get(key));
								} else {
									linkCount.put(key, linkCount1.get(key));
								}
							}
							exceed_ins += exceed_ins1;
						}

					} catch (Exception e) {
						// TODO Auto-generated catch block
						Thread t = Thread.currentThread();
						t.getUncaughtExceptionHandler().uncaughtException(t, e);
						e.printStackTrace();
						executor.shutdown();
						System.exit(1);
					}
				} 

				public Runnable init(String bamfile) {
					// TODO Auto-generated method stub
					this.bamfile = bamfile;
					return this;
				}				
			}.init(bamfile));

		this.waitFor();
		
		myLogger.info("################################################################");
		myLogger.info("#insert size exceeds "+maxInst+": "+exceed_ins);
		myLogger.info("#good links: "+linkCount.size()/2);
		myLogger.info("################################################################");
		
		try {
			BufferedWriter bw = Utils.getBufferedWriter(this.out_prefix+".links");
			String source, target;
			for(long key : linkCount.keySet()) {
				target = seq_index.getKey((int)  key     );
				source = seq_index.getKey((int) (key>>32));
				bw.write(source+"\t"+target+"\t"+linkCount.get(key)+"\n");
			}
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}




