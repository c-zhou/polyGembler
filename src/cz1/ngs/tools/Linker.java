package cz1.ngs.tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;

import cz1.ngs.model.GFA;
import cz1.ngs.model.Sequence;
import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Executor;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Linker extends Executor {

	private final static boolean USE_OS_BUFFER = true;
	// 8Mb buffer size
	private final static int buffSize8Mb = 8388608;
	private final static Writer STD_OUT_BUFFER = USE_OS_BUFFER ? 
			new BufferedWriter(new OutputStreamWriter(System.out), buffSize8Mb) : new OutputStreamWriter(System.out);

			@Override
			public void printUsage() {
				// TODO Auto-generated method stub
				myLogger.info(
						"\n\nUsage is as follows:\n"
								+ " -s/--subject            Subject/reference sequences file in FASTA format.\n"
								+ " -pe/--pe                BAM files for paired-end reads (multiple files seperated by ':').\n"
								+ " -long/--long            BAM files for long reads (multiple files seperated by ':').\n"
								+ " -g/--graph              Assembly graph (GFA) format. Currently, the program only accept \n"
								+ "                         the assembly graph format used by the assembler SPAdes (de-bruijn \n"
								+ "                         graph) or CANU (overlap). For assembly graphs in other formats: \n"
								+ "                         1. for contigs generated with de-bruijn graph, please ignore this \n"
								+ "                         option and provide k-mer size, the program is able to rebuild the \n"
								+ "                         assembly graph from contigs; and, \n"
								+ "                         2. for contigs generated with overlapping algorithms, please convert \n"
								+ "                         assembly graph to the format used by CANU.\n"
								+ "                         NOTE: it is possible to run the program without an assembly graph, \n"
								+ "                         however, the assembly might be less accurate.\n"
								+ " -t/--threads            Number of threads to use (default 16).\n"
								+ " -d/--debug              Debugging mode will have extra information printed out.\n"
								+ " -dd/--debug-debug       Debugging mode will have more information printed out than -d mode.\n"
								+ " -o/--out-prefix         Prefix of the output files.\n"
								+ "\n");	
			}

			private String subject_file = null;
			private String[] pe_files = null;
			private String[] long_files = null;
			private String graph_file = null;

			private int num_threads = Runtime.getRuntime().availableProcessors();
			private String out_prefix = null;
			private boolean debug  = false;
			private boolean ddebug = false;

			@Override
			public void setParameters(String[] args) {
				// TODO Auto-generated method stub
				if (myArgsEngine == null) {
					myArgsEngine = new ArgsEngine();
					myArgsEngine.add("-s", "--subject", true);
					myArgsEngine.add("-pe", "--pe", true);
					myArgsEngine.add("-long", "--long", true);
					myArgsEngine.add("-g","--graph", true);
					myArgsEngine.add("-t", "--threads", true);
					myArgsEngine.add("-d", "--debug", false);
					myArgsEngine.add("-dd", "--debug-debug", false);
					myArgsEngine.add("-o", "--out-prefix", true);
					myArgsEngine.parse(args);
				}

				if (myArgsEngine.getBoolean("-s")) {
					this.subject_file = myArgsEngine.getString("-s");
				} else {
					printUsage();
					throw new IllegalArgumentException("Please specify the subject/reference file.");
				}

				if (myArgsEngine.getBoolean("-pe")) {
					this.pe_files = myArgsEngine.getString("-pe").split(":");
				}

				if (myArgsEngine.getBoolean("-long")) {
					this.long_files = myArgsEngine.getString("-long").split(":");
				}

				if (myArgsEngine.getBoolean("-g")) {
					this.graph_file = myArgsEngine.getString("-g");
				}

				if (myArgsEngine.getBoolean("-o")) {
					this.out_prefix = myArgsEngine.getString("-o");
					if(new File(out_prefix+".fa").exists() || 
							new File(out_prefix+".map").exists() || 
							new File(out_prefix+"_unplaced.map").exists()) {
						throw new RuntimeException("Output files exist. Please specify a different name.");
					}
				} else {
					printUsage();
					throw new IllegalArgumentException("Please specify the prefix of output files.");
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
			}

			private Map<String, Sequence> qry_seqs;
			private Map<String, Sequence> sub_seqs;
			private GFA gfa;
			private final BidiMap<String, Integer> seq_index = new DualHashBidiMap<String, Integer>();
			private final Map<String, String> symm_seqsn = new HashMap<String, String>();

			@Override
			public void run() {
				sub_seqs = Sequence.parseFastaFileWithRevCmpAsMap(subject_file);
				int index = 0;
				for(String seq : sub_seqs.keySet()) seq_index.put(seq, ++index);
				for(String seq : sub_seqs.keySet()) {
					if(!seq.endsWith("'")) {
						symm_seqsn.put(seq, seq+"'");
						symm_seqsn.put(seq+"'", seq);
					}
				}
				gfa = new GFA(subject_file, graph_file);

				this.link_pe();
				this.link_long();
			}

			final static SamReaderFactory factory =
					SamReaderFactory.makeDefault()
					.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
							SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
					.validationStringency(ValidationStringency.SILENT);

			private final static Map<Long, Integer> linkCount = new HashMap<Long, Integer>();
			private final static int m_ins = 2000; // maximum insert size for pe read library
			private final static int m_lnk = 3;    // minimum #link to confirm a link
			private final static int m_qual = 20;  // minimum alignment quality
			private final static Object lock = new Object();

			private void link_long() {
				// TODO Auto-generated method stub
				this.initial_thread_pool();

				for(final String qf : long_files) {
					executor.submit(new Runnable() {

						private String query_file;

						@Override
						public void run() {
							// TODO Auto-generated method stub
							try {
								final SamReader in1 = factory.open(new File(query_file));
								final SAMRecordIterator iter1 = in1.iterator();
								final Map<Long, Integer> linkCount1 = new HashMap<Long, Integer>();
								final Set<String> mappedRefStr = new HashSet<String>(); 
								long refind;

								SAMRecord sam_record = iter1.hasNext()?iter1.next():null;
								String readn;
								while(sam_record!=null) {
									mappedRefStr.clear();
									if(!sam_record.getReadUnmappedFlag()&&
											!sam_record.getNotPrimaryAlignmentFlag()&&
											sam_record.getMappingQuality()>=m_qual)
										mappedRefStr.add(sam_record.getReferenceName());
									readn = sam_record.getReadName();
									while((sam_record=iter1.hasNext()?iter1.next():null)!=null && 
											readn.equals(sam_record.getReadName())) {
										if(!sam_record.getReadUnmappedFlag()&&
												!sam_record.getNotPrimaryAlignmentFlag()&&
												sam_record.getMappingQuality()>=m_qual)
											mappedRefStr.add(sam_record.getReferenceName());
									}

									if(mappedRefStr.size()<=1) continue;
									for(final String refstr1 : mappedRefStr) {
										for(final String refstr2 : mappedRefStr) {
											if(refstr1.equals(refstr2)) continue;

											refind  = seq_index.get(refstr1);
											refind <<= 32;
											refind += seq_index.get(refstr2);
											if(linkCount1.containsKey(refind)) {
												linkCount1.put(refind, linkCount1.get(refind)+1);
											} else {
												linkCount1.put(refind, 1);
											}
											
											refind  = seq_index.get(symm_seqsn.get(refstr1));
											refind <<= 32;
											refind += seq_index.get(symm_seqsn.get(refstr2));
											if(linkCount1.containsKey(refind)) {
												linkCount1.put(refind, linkCount1.get(refind)+1);
											} else {
												linkCount1.put(refind, 1);
											}

										}	
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
									myLogger.info(this.query_file+" processed.");
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

						public Runnable init(String query_file) {
							this.query_file = query_file;
							return this;
						}

					}.init(qf));
				}
				this.waitFor();
			}

			private void link_pe() {
				// TODO Auto-generated method stub
				this.initial_thread_pool();

				for(final String qf : pe_files) {
					executor.submit(new Runnable() {

						private String query_file;

						@Override
						public void run() {
							// TODO Auto-generated method stub
							try {
								final SamReader in1 = factory.open(new File(query_file));
								final SAMRecordIterator iter1 = in1.iterator();
								final SAMRecord[] sam_records = new SAMRecord[2];
								SAMRecord sam_record;
								final Map<Long, Integer> linkCount1 = new HashMap<Long, Integer>();
								final boolean[] rev = new boolean[2];
								final int[] reflen = new int[2];
								final String[] refstr = new String[2];
								int a, b;
								long refind;
								int match, refPos;
								String refStr, readStr;
								SequencePair<DNASequence, NucleotideCompound> seqPair;
								boolean ambiguous;

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
											// reads aligned to the same contig
											Arrays.fill(sam_records, null);
											continue;
										}
										if(sam_records[0].getMappingQuality()+
												sam_records[1].getMappingQuality()<m_qual) {
											Arrays.fill(sam_records, null);
											continue;
										}
										rev[0] = sam_records[0].getReadNegativeStrandFlag();
										rev[1] = sam_records[1].getReadNegativeStrandFlag();
										reflen[0]  = sub_seqs.get(sam_records[0].getReferenceName()).seq_ln();
										reflen[1]  = sub_seqs.get(sam_records[1].getReferenceName()).seq_ln();

										// to avoid distortions of the alignment caused by the redundancies
										// we realign the read pair
										ambiguous = false;
										for(int i=0; i<2; i++) {
											List<CigarElement> cigar = sam_records[i].getCigar().getCigarElements();
											match = 0;
											for(CigarElement element : cigar) {
												if(element.getOperator() == CigarOperator.M)
													match += element.getLength();
											}
											match -= sam_records[i].getIntegerAttribute("NM");


											readStr = sam_records[i].getReadString();
											if(rev[1-i]) {
												refPos = sam_records[1-i].getStart();
												refStr = sub_seqs.get(sam_records[1-i].getReferenceName()).seq_str().
														substring(Math.max(0, refPos-m_ins), refPos);
											} else {
												refPos = sam_records[1-i].getEnd();
												refStr = sub_seqs.get(sam_records[1-i].getReferenceName()).seq_str().
														substring(refPos, Math.min(refPos+m_ins, reflen[1-i]));
											}

											if(refStr.length()<0.8*match) continue;

											seqPair = Constants.localPairMatcher(refStr, readStr);

											if(seqPair.getNumIdenticals()>=0.8*match) {
												ambiguous = true;
												break;
											}
										}

										if(ambiguous) {
											Arrays.fill(sam_records, null);
											continue;
										}

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
											if(a>m_ins&&b>m_ins) {
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
											if(a>m_ins&&b>m_ins) {
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
											if(a>m_ins&&b>m_ins) {
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
											if(a>m_ins&&b>m_ins) {
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
									myLogger.info(this.query_file+" processed.");
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

						public Runnable init(String query_file) {
							this.query_file = query_file;
							return this;
						}

					}.init(qf));
				}
				this.waitFor();
			}
}


