package cz1.tenx.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;

import cz1.ngs.model.Sequence;
import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Executor;
import cz1.util.Utils;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class TenXArcs extends Executor {
	private static enum Task {link, arcs, gv, zzz}
	private Task task_list = Task.zzz;

	private String[] bamFiles = null;
	private boolean[] matePair = null; // mate-pair library
	private int[] inst = null;
	private double[] insE = null;
	private int minQual = 10;
	private String out = null;
	private String bcList = null;
	private String linkFile = null;
	private String seqFile = null;
	private int bcn = 3;

	@Override
	public void printUsage() {
		// TODO Auto-generated method stub

		switch(this.task_list) {
		case zzz:
			myLogger.info(
					"\n\nChoose a task from listed:\n"
							+ " link                  Count links.\n"
							+ " arcs                 Anchor contigs to generate scaffolds.\n"
							+ " gv                    Convert link file to assembly graph file.\n"
							+ "\n");
			break;

		case link:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -b/--bam-file          BAM file for the alignment records for read pairs.\n"
							+ " -c/--config-file        Configuration file for BAM files.\n"
							+ " -mp/--mate-pair     Mate-pair library.\n"
							+ " -q/--min-qual          Minimum alignment quality (default: 10).\n"
							+ " -i/--insert-size        Insert size of the library.\n"
							+ " -e/--ins-error          Allowed insert size error.\n" 
							+ " -t/--threads            Number of threads to use (default: all avaiable processors).\n"
							+ " -o/--out                  Output files.\n"
							+ "\n");
			break;

		case arcs:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -b/--bam-file           BAM file for the alignment records for read pairs.\n"
							+ " -f/--bam-fof            File of BAM file list.\n"
							+ " -w/--white-list         File of barcodes on the white list.\n"
							+ " -n/--read-number    Minimum number of read pair per barcode (defalt: 3).\n"
							+ " -t/--threads             Number of threads to use (default: all avaiable processors).\n"
							+ " -o/--out                  Output file.\n\n");
			break;
			
		case gv:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -l/--links                 Link file.\n"
							+ "-s/--seqs                 Sequence file. \n"
							+ " -o/--out                  Output file.\n\n");
			break;
			
		default:
			throw new RuntimeException("!!!");
		}
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}

		switch(args[0].toLowerCase()) {
		case "link":
			this.task_list = Task.link;
			break;
		case "arcs":
			this.task_list = Task.arcs;
			break;
		case "gv":
			this.task_list = Task.gv;
			break;
		default:
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");	
		}

		String[] args2 = new String[args.length-1];
		System.arraycopy(args, 1, args2, 0, args2.length);

		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-b", "--bam-file", true);
			myArgsEngine.add("-c", "--config-file", true);
			myArgsEngine.add("-mp", "--mate-pair", true);
			myArgsEngine.add("-q", "--min-qual", true);
			myArgsEngine.add("-i", "--insert-size", true);
			myArgsEngine.add("-e", "--ins-error", true);
			myArgsEngine.add("-f", "--bam-fof", true);
			myArgsEngine.add("-w", "--white-list", true);
			myArgsEngine.add("-n", "--read-number", true);
			myArgsEngine.add("-l", "--links", true);
			myArgsEngine.add("-s", "--seqs", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-o", "--out", true);
			myArgsEngine.parse(args2);
		}

		switch(this.task_list) {
		case link:
			if (!myArgsEngine.getBoolean("-b")&&!myArgsEngine.getBoolean("-c")) {
				printUsage();
				throw new IllegalArgumentException("Please specify the alignment file(s) using -b or -c option.");
			}

			if (myArgsEngine.getBoolean("-b")) {
				this.bamFiles = new String[]{myArgsEngine.getString("-b")};
				this.matePair = new boolean[1];
				if(myArgsEngine.getBoolean("-mp")) {
					this.matePair[0] = true;
					myLogger.info("a mate-pair library.");
				}
				this.inst = new int[1];
				this.insE = new double[1];
				if (myArgsEngine.getBoolean("-i")) {
					this.inst[0] = Integer.parseInt(myArgsEngine.getString("-i"));
					this.insE[0] = Double.parseDouble(myArgsEngine.getString("-e"));
				} else {
					printUsage();
					throw new IllegalArgumentException("Please specify insert size of the library (-i).");
				}
			} else {
				if(myArgsEngine.getBoolean("-mp"))
					myLogger.info("mate-pair (-mp) indicator ignored.");
				if(myArgsEngine.getBoolean("-i"))
					myLogger.info("insert size (-i) option ignored.");
			}

			if (myArgsEngine.getBoolean("-c")) {
				try {
					BufferedReader br = Utils.getBufferedReader(myArgsEngine.getString("-c").trim());
					final List<String> f = new ArrayList<String>();
					final List<Boolean> rev = new ArrayList<Boolean>();
					final List<Integer> ins = new ArrayList<Integer>();
					final List<Double>  err = new ArrayList<Double>();
					String line;
					String[] s;
					while( (line = br.readLine()) != null) {
						line = line.trim();
						s = line.split("\\s+");
						if(s.length>0) {
							f.add(s[0].trim());
							rev.add(s[1].trim().equals("0")?false:true);
							ins.add(Integer.parseInt(s[2].trim()));
							err.add(Double.parseDouble(s[3].trim()));
						}
					}
					br.close();

					final int n = f.size();
					this.bamFiles = new String[n];
					this.matePair = new boolean[n];
					this.inst = new int[n];
					this.insE = new double[n];
					for(int i=0; i<n; i++) {
						this.bamFiles[i] = f.get(i);
						this.matePair[i] = rev.get(i);
						this.inst[i] = ins.get(i);
						this.insE[i] = err.get(i);
					}
				} catch (IOException e) {
					// TODO Auto-generated catch block
					myLogger.error("configuration file with invalid format.");
					e.printStackTrace();
				}

			}

			if (myArgsEngine.getBoolean("-q")) {
				this.minQual = Integer.parseInt(myArgsEngine.getString("-q"));
			}

			break;
		case arcs:
			if (!myArgsEngine.getBoolean("-b")&&!myArgsEngine.getBoolean("-f")) {
				printUsage();
				throw new IllegalArgumentException("Please specify the BAM file or a BAM list file.");
			}
			if(myArgsEngine.getBoolean("-b")) {
				this.bamFiles = new String[]{myArgsEngine.getString("-b")};
			}

			if(myArgsEngine.getBoolean("-f")) {
				if(this.bamFiles!=null) 
					myLogger.warn("The program will use the BAM list file instead of the BAM file you provided.");
				this.bamFiles = myArgsEngine.getString("-f").split(",");
			}

			if(myArgsEngine.getBoolean("-w")) {
				this.bcList = myArgsEngine.getString("-w");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the barcode white list.");
			}

			if(myArgsEngine.getBoolean("-n")) {
				this.bcn = Integer.parseInt(myArgsEngine.getString("-n"));
			}
			break;
		case gv:
			if(myArgsEngine.getBoolean("-l")) {
				this.linkFile = myArgsEngine.getString("-l");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify link file.");
			}
			if(myArgsEngine.getBoolean("-s")) {
				this.seqFile = myArgsEngine.getString("-s");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify sequence file.");
			}
			break;
		default:
			throw new RuntimeException("!!!");
		}

		this.THREADS = Runtime.getRuntime().availableProcessors();
		if (myArgsEngine.getBoolean("-t")) {
			int t = Integer.parseInt(myArgsEngine.getString("-t"));
			if(t<this.THREADS) this.THREADS = t;
			this.THREADS = t;
			Constants.omp_threads = this.THREADS;
			myLogger.info("OMP_THREADS = "+this.THREADS);
		}

		if(myArgsEngine.getBoolean("-o")) {
			this.out = myArgsEngine.getString("-o");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the output file.");
		}
	}

	private final SamReaderFactory factory =
			SamReaderFactory.makeDefault()
			.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
					SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
			.validationStringency(ValidationStringency.SILENT);

	private class BAMBarcodeIterator {

		private final SamReader samReader;
		private final SAMRecordIterator iter;
		private SAMRecord samRecord = null;

		public BAMBarcodeIterator(String bam_file) {
			this.samReader = factory.open(new File(bam_file));
			this.iter = this.samReader.iterator();
			this.samRecord = iter.hasNext() ? iter.next() : null;
		}

		public boolean hasNext() {
			return samRecord != null;
		}

		public List<SAMRecord[]> next() {

			if(!this.hasNext()) throw new RuntimeException("!!!");

			List<SAMRecord[]> bc_records = new ArrayList<SAMRecord[]>();
			String bc = samRecord.getStringAttribute("BX");

			String sn;
			SAMRecord[] records = new SAMRecord[2];

			while( samRecord!=null && samRecord.getStringAttribute("BX").equals(bc) ) {
				sn = samRecord.getReadName();

				if( !samRecord.getReadUnmappedFlag() &&
						!samRecord.getNotPrimaryAlignmentFlag()&&
						!samRecord.getSupplementaryAlignmentFlag() ) {
					if(samRecord.getFirstOfPairFlag())
						records[0] = samRecord;
					else if(samRecord.getSecondOfPairFlag())
						records[1] = samRecord;
					else
						throw new RuntimeException("!!!");
				}

				while( (samRecord = iter.hasNext() ? iter.next() : null)!=null && samRecord.getReadName().equals(sn)) {
					if( !samRecord.getReadUnmappedFlag() &&
							!samRecord.getNotPrimaryAlignmentFlag()&&
							!samRecord.getSupplementaryAlignmentFlag() ) {
						if(samRecord.getFirstOfPairFlag())
							records[0] = samRecord;
						else if(samRecord.getSecondOfPairFlag())
							records[1] = samRecord;
						else
							throw new RuntimeException("!!!");
					}
				}

				if(records[0]!=null && records[1]!=null &&
						records[0].getReferenceIndex().intValue()==records[1].getReferenceIndex().intValue()) {
					if(records[0].getInferredInsertSize()+records[1].getInferredInsertSize()!=0)
						throw new RuntimeException("!!!");
					bc_records.add(records);
				}

				records = new SAMRecord[2];
			}

			return bc_records;
		}

		public void close() {
			try {
				this.iter.close();
				this.samReader.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	final int shift = 30;
	final int rev = 1<<shift;
	final Map<Long, Integer> arcs = new HashMap<Long, Integer>();

	public void run_arcs() {
		// TODO Auto-generated method stub
		final Set<String> bc_white = new HashSet<String>();

		try {
			final BufferedReader br_bc = Utils.getBufferedReader(this.bcList);
			String line;
			while( (line=br_bc.readLine()) != null ) 
				bc_white.add(line.split(",")[0]);
			br_bc.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		final SamReader samReader = factory.open(new File(this.bamFiles[0]));
		final SAMSequenceDictionary seqdict = samReader.getFileHeader().getSequenceDictionary();
		try {
			samReader.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		this.initial_thread_pool();
		for(final String bamFile : bamFiles) {
			this.executor.submit(new Runnable() {
				private String bamFile;

				@Override
				public void run() {
					// TODO Auto-generated method stub
					try {
						final BAMBarcodeIterator iter = new BAMBarcodeIterator(this.bamFile);
						List<SAMRecord[]> bc_records;
						Map<String, List<SAMRecord[]>> bc_blocks;
						while(iter.hasNext()) {
							bc_records = iter.next();
							if(bc_records.isEmpty() || 
									!bc_white.contains(bc_records.get(0)[0].getStringAttribute("BX")))
								continue;
							bc_blocks = bin(bc_records);
						}

						iter.close();
					} catch (Exception e) {
						Thread t = Thread.currentThread();
						t.getUncaughtExceptionHandler().uncaughtException(t, e);
						e.printStackTrace();
						executor.shutdown();
						System.exit(1);
					}
				}

				private Map<String, List<SAMRecord[]>> bin(List<SAMRecord[]> bc_records) {
					// TODO Auto-generated method stub
					final Map<String, List<SAMRecord[]>> bc_blocks = new HashMap<String, List<SAMRecord[]>>();

					return bc_blocks;
				}

				public Runnable init(String bamFile) {
					// TODO Auto-generated method stub
					this.bamFile = bamFile;
					return this;
				}

			}.init(bamFile));
		}
		this.waitFor();
	}

	final Object lock = new Object();

	private void run_link() {
		// TODO Auto-generated method stub
		final int niof = this.bamFiles.length;
		myLogger.info("Reading alignments from "+niof+" BAM file"+
				(this.bamFiles.length>1?"s":"")+":");
		for(int i=0; i<niof; i++)
			myLogger.info(this.bamFiles[i]);
		myLogger.info("****");

		
		final SamReader samReader = factory.open(new File(this.bamFiles[0]));
		final SAMSequenceDictionary seqdict = samReader.getFileHeader().getSequenceDictionary();
		try {
			samReader.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		final BidiMap<String, Integer> seq_index = new DualHashBidiMap<String, Integer>();
		final Map<String, String> symm_seqsn = new HashMap<String, String>();
		
		int index = 0;
		for(SAMSequenceRecord seq : seqdict.getSequences()) {
			seq_index.put(seq.getSequenceName(),     ++index);
			seq_index.put(seq.getSequenceName()+"'", ++index);
		}
		for(String seq : seq_index.keySet()) {
			if(!seq.endsWith("'")) {
				symm_seqsn.put(seq, seq+"'");
				symm_seqsn.put(seq+"'", seq);
			}
		}
		
		this.initial_thread_pool();
		final Map<Long, Integer> linkCount = new HashMap<Long, Integer>();

		for(int i=0; i<niof; i++) {
			myLogger.info("####Processing "+this.bamFiles[i]);

			this.executor.submit(new Runnable() {
				private String bamFile;
				private int inst;
				private double insE;
				private boolean matePair;

				@Override
				public void run() {
					// TODO Auto-generated method stub
					int readCount = 0;
					final double maxInst = inst*(1+insE);
					final double minInst = inst*(1-insE);
					
					try {
						final SamReader samReader = factory.open(new File(this.bamFile));
						final SAMRecordIterator iter = samReader.iterator();

						int inst, s1, s2, e1, e2;
						boolean rev1, rev2;
						String refstr1, refstr2;
						int reflen1, reflen2;
						SAMRecord tmp, r1, r2;
						String readName;
						tmp = iter.hasNext() ? iter.next() : null;

						while(tmp!=null) {
							synchronized(lock) {
								++readCount;
								if(readCount%1000000==0)
									myLogger.info("####"+this.bamFile+" #"+readCount+" read pairs processed.");
							}

							r1 = null;
							r2 = null;

							readName = tmp.getReadName();
							while(tmp!=null && tmp.getReadName().equals(readName)) {
								if( !tmp.getNotPrimaryAlignmentFlag() && 
										!tmp.getSupplementaryAlignmentFlag() &&
										!tmp.getReadUnmappedFlag() &&
										tmp.getMappingQuality()>=minQual) {
									if(tmp.getFirstOfPairFlag())  r1 = tmp;
									if(tmp.getSecondOfPairFlag()) r2 = tmp;
								}
								tmp = iter.hasNext()?iter.next():null;
							}

							if(r1==null||r2==null) continue;

							s1 = r1.getAlignmentStart();
							e1 = r1.getAlignmentEnd();
							rev1 = r1.getReadNegativeStrandFlag();
							s2 = r2.getAlignmentStart();
							e2 = r2.getAlignmentEnd();
							rev2 = r2.getReadNegativeStrandFlag();

							if(r1.getReferenceIndex().intValue()==r2.getReferenceIndex().intValue()) 
								continue;

							inst = Integer.MAX_VALUE;

							// mapped to different contigs
							reflen1  = seqdict.getSequence(r1.getReferenceIndex()).getSequenceLength();
							reflen2  = seqdict.getSequence(r2.getReferenceIndex()).getSequenceLength();

							if(matePair) {
								// this is a mate-pair
								if(rev1&&rev2) {
									//       0                  1
									// ---------------     -------------
									//   <===                 <===

									// reverse 0
									// ---------------     -------------
									//          ===>          <===
									// reverse 1
									// ---------------     -------------
									//   <===                    ===>

									// reverse 0 and 1 are symmetric
									refstr1 = r1.getReferenceName();
									refstr2 = r2.getReferenceName()+"'";
									inst = reflen1-s1+1+reflen2-s2+1;

								} else if(rev1&&!rev2) {
									//       0                  1
									// ---------------     -------------
									//   <===                    ===>

									// reverse 0 & reverse 1
									// ---------------     -------------
									//          ===>          <===
									refstr1 = r1.getReferenceName();
									refstr2 = r2.getReferenceName();
									inst = reflen1-s1+1+e2;

								} else if(!rev1&&rev2) {
									//       0                  1
									// ---------------     -------------
									//          ===>          <===

									// reverse 0 & reverse 1
									// ---------------     -------------
									//   <===                    ===>
									refstr1 = r2.getReferenceName();
									refstr2 = r1.getReferenceName();
									inst = e1+reflen2-s2+1;

								} else {
									//       0                  1
									// ---------------     -------------
									//         ===>               ===>

									// reverse 0
									// ---------------     -------------
									//    <===                    ===>
									// reverse 1
									// ---------------     -------------
									//         ===>          <===
									refstr1 = r1.getReferenceName()+"'";
									refstr2 = r2.getReferenceName();
									inst = e1+e2;

								}
							} else {
								// this is a paired-end
								if(rev1&&rev2) {
									//       0                  1
									// ---------------     -------------
									//   <===                 <===

									// reverse 0
									// ---------------     -------------
									//          ===>          <===
									// reverse 1
									// ---------------     -------------
									//   <===                    ===>

									// reverse 0 and 1 are symmetric
									refstr1 = r1.getReferenceName()+"'";
									refstr2 = r2.getReferenceName();
									inst = e1+e2;

								} else if(rev1&&!rev2) {
									//       0                  1
									// ---------------     -------------
									//   <===                    ===>

									// reverse 0 & reverse 1
									// ---------------     -------------
									//          ===>          <===

									refstr1 = r2.getReferenceName();
									refstr2 = r1.getReferenceName();
									inst = e1+reflen2-s2+1;

								} else if(!rev1&&rev2) {
									//       0                  1
									// ---------------     -------------
									//          ===>          <===

									// reverse 0 & reverse 1
									// ---------------     -------------
									//   <===                    ===>

									refstr1 = r1.getReferenceName();
									refstr2 = r2.getReferenceName();
									inst = reflen1-s1+1+e2;

								} else {
									//       0                  1
									// ---------------     -------------
									//         ===>               ===>

									// reverse 0
									// ---------------     -------------
									//    <===                    ===>
									// reverse 1
									// ---------------     -------------
									//         ===>          <===

									refstr1 = r1.getReferenceName();
									refstr2 = r2.getReferenceName()+"'";
									inst = reflen1-s1+1+reflen2-s2+1;
								}
							}

							if(inst<minInst || inst>maxInst) continue;
							
							long refind;
							refind  = seq_index.get(refstr1);
							refind <<= 32;
							refind += seq_index.get(refstr2);

							synchronized(lock) {
								if(linkCount.containsKey(refind)) 
									linkCount.put(refind, linkCount.get(refind)+1);
								else
									linkCount.put(refind, 1);
							}

							refind  = seq_index.get(symm_seqsn.get(refstr2));
							refind <<= 32;
							refind += seq_index.get(symm_seqsn.get(refstr1));

							synchronized(lock) {
								if(linkCount.containsKey(refind)) 
									linkCount.put(refind, linkCount.get(refind)+1);
								else
									linkCount.put(refind, 1);
							}
						}
						samReader.close();

					} catch (Exception e) {
						// TODO Auto-generated catch block
						Thread t = Thread.currentThread();
						t.getUncaughtExceptionHandler().uncaughtException(t, e);
						e.printStackTrace();
						executor.shutdown();
						System.exit(1);
					}
				}

				public Runnable init(final String bamFile, final boolean matePair, final int inst, final double insE) {
					// TODO Auto-generated method stub
					this.bamFile = bamFile;
					this.matePair = matePair;
					this.inst = inst;
					this.insE = insE;
					return this;
				}

			}.init(this.bamFiles[i], this.matePair[i], this.inst[i], this.insE[i]));
		}
		this.waitFor();
		
		try {
			BufferedWriter bw = Utils.getBufferedWriter(this.out);
			String source, target;
			for(long key : linkCount.keySet()) {
				target = seq_index.getKey((int)  key     );
				source = seq_index.getKey((int) (key>>32));
				bw.write(key+"\t"+source+"\t"+target+"\t"+linkCount.get(key)+"\n");
			}
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void run_gv() {
		// TODO Auto-generated method stub
		try {
			BufferedReader br = Utils.getBufferedReader(this.linkFile);
			String line;
			String[] s;
			final Set<String> seqs = new HashSet<String>();
			while( (line=br.readLine())!=null) {
				s = line.split("\\s+");
				seqs.add(s[1].replaceAll("'$", ""));
				seqs.add(s[2].replaceAll("'$", ""));
			}
			br.close();
			int index = 0;
			final Map<String, Integer> seq_index = new HashMap<String, Integer>();
			final Map<String, Integer> edges = new HashMap<String, Integer>();
			final BidiMap<String, String> syms = new DualHashBidiMap<String, String>();
			for(String seq : seqs) {
				syms.put(seq, seq+"'");
				syms.put(seq+"'", seq);
				seq_index.put(seq,      index);
				seq_index.put(seq+"'", index);
				++index;
			}
			final List<Sequence> org_seqs = Sequence.parseFastaFileAsList(seqFile);
			
			BufferedReader br_link = Utils.getBufferedReader(this.linkFile);
			BufferedWriter bw_origv = Utils.getBufferedWriter(this.out+"_original.gv");
			BufferedWriter bw_dstgv = Utils.getBufferedWriter(this.out+".dist.gv");
			bw_origv.write("graph G {\n");
			bw_dstgv.write("digraph arcs {\n");
			for(String seq : seq_index.keySet()) {
				if(!seq.endsWith("'"))
					bw_origv.write(seq_index.get(seq)+" [id="+seq+"];\n");
			}
			for(Sequence seq : org_seqs) {
				bw_dstgv.write("\""+seq.seq_sn()+"+\" [l="+seq.seq_ln()+"]\n");
				bw_dstgv.write("\""+seq.seq_sn()+"-\" [l="+seq.seq_ln()+"]\n");
			}
			
			String sym;
			int source, target, label = -1;
			String source_id, target_id;
			while( (line=br_link.readLine())!=null ) {
				s = line.split("\\s+");
				
				source_id = s[1].replaceAll("'$", "")+(s[1].endsWith("'")?"-":"+");
				target_id = s[2].replaceAll("'$", "")+(s[2].endsWith("'")?"-":"+");
				bw_dstgv.write("\""+source_id+"\" -> \""+target_id+"\" [d=100 e=100.0 n="+s[3]+"]\n");
				
				sym = syms.get(s[2])+"+"+syms.get(s[1]);
				if(edges.containsKey(sym)) {
					if( !(edges.get(sym).intValue()==Integer.parseInt(s[3])) )
						throw new RuntimeException("!!!");
					continue;
				}
				if(s[1].endsWith("'") &&!s[2].endsWith("'")) label = 0;
				if(s[1].endsWith("'") && s[2].endsWith("'")) label = 1;
				if(!s[1].endsWith("'")&&!s[2].endsWith("'")) label = 2;
				if(!s[1].endsWith("'")&& s[2].endsWith("'")) label = 3;
				source = seq_index.get(s[1]);
				target = seq_index.get(s[2]);
				bw_origv.write(source+"--"+target+" [label="+label+", weight="+s[3]+"];\n");
				edges.put(s[1]+"+"+s[2], Integer.parseInt(s[3]));
			}
			bw_origv.write("}\n");
			bw_dstgv.write("}\n");
			br_link.close();
			bw_origv.close();
			bw_dstgv.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		switch(this.task_list) {
		case zzz:
			myLogger.info("Task list is empty!!!");
			break;
		case link:
			this.run_link();
			break;
		case arcs:
			this.run_arcs();
			break;
		case gv:
			this.run_gv();
			break;
		default:
			throw new RuntimeException("!!!");
		}
		return;
	}

}





