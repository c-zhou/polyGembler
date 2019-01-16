package cz1.tenx.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.DirectedPseudograph;
import org.jgrapht.graph.DirectedWeightedPseudograph;

import cz1.ngs.model.OverlapEdge;
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
	private static enum Task {link, arcs, gv, dist, zzz}
	private Task task_list = Task.zzz;

	private String[] bamFiles = null;
	private boolean[] matePair = null; // mate-pair library
	private int[] inst = null;
	private double[] insE = null;
	private int minQual = 20;
	private String out = null;
	private String bcList = null;
	private String linkFile = null;
	private String seqFile = null;
	private int bcn = 3;
	private int maxGap = 30000;
	private int minLink = 5;
	private double minRatio = 0.2;
	private int radius = 30000;
	private String distFile = null;

	@Override
	public void printUsage() {
		// TODO Auto-generated method stub

		switch(this.task_list) {
		case zzz:
			myLogger.info(
					"\n\nChoose a task from listed:\n"
							+ " link                  Count links.\n"
							+ " arcs                  Anchor contigs to generate scaffolds.\n"
							+ " gv                    Convert link file to assembly graph file.\n"
							+ " dist                  Make a distance matrix for a link graph with a given radius.\n"
							+ "\n");
			break;

		case link:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -b/--bam-file         BAM file for the alignment records for read pairs.\n"
							+ " -c/--config-file      Configuration file for BAM files.\n"
							+ " -mp/--mate-pair       Mate-pair library.\n"
							+ " -q/--min-qual         Minimum alignment quality (default: 20).\n"
							+ " -i/--insert-size      Insert size of the library.\n"
							+ " -e/--ins-error        Allowed insert size error.\n" 
							+ " -t/--threads          Number of threads to use (default: all avaiable processors).\n"
							+ " -o/--out              Output files.\n\n");
			break;

		case arcs:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -b/--bam-file         BAM file for the alignment records for read pairs.\n"
							+ " -f/--bam-fof          File of BAM file list.\n"
							+ " -q/--min-qual         Minimum min average quality scores per molecule (default: 20).\n"
							+ " -w/--white-list       File of barcodes on the white list.\n"
							+ " -n/--read-number      Minimum number of read pair per barcode (defalt: 3).\n"
							+ " -g/--max-gap          Maximum gap size for two adjcent read pairs in a molecule (default: 30000).\n"
							+ " -d/--dist             Distance file.\n"
							+ " -t/--threads          Number of threads to use (default: all avaiable processors).\n"
							+ " -o/--out              Output file.\n\n");
			break;
			
		case gv:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -l/--links            Link file.\n"
							+ " -s/--seqs             Sequence file.\n"
							+ " -o/--out              Output file.\n\n");
			break;
			
		case dist:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -l/--links            Link file.\n"
							+ " -k/--min-links        Minimum #links to create a edge (default: 5).\n"
							+ " -a/--ratio            Minimum link ratio between two best contig pairs (default: 0.2).\n"
							+ "                       *higher values lead to more stringent link graph trimming*\n"
							+ "                       e.g., a node of three edges with #links=100, 20 and 10, repectively,\n"
							+ "                       -a 0.2 will keep two edges while -a 0.3 will only keep the first one.\n"
							+ " -s/--seqs             Sequence file.\n"
							+ " -r/--radius           Radius to search (default: 50000).\n"
							+ " -t/--threads          Number of threads to use (default: all avaiable processors).\n"
							+ " -o/--out              Output file.\n\n");
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
		case "dist":
			this.task_list = Task.dist;
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
			myArgsEngine.add("-g", "--max-gap", true);
			myArgsEngine.add("-d", "--dist", true);
			myArgsEngine.add("-l", "--links", true);
			myArgsEngine.add("-s", "--seqs", true);
			myArgsEngine.add("-k", "--min-links", true);
			myArgsEngine.add("-a", "--ratio", true);
			myArgsEngine.add("-r", "--radius", true);
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
			
			if(myArgsEngine.getBoolean("-g")) {
				this.maxGap = Integer.parseInt(myArgsEngine.getString("-g"));
			}
			
			if(myArgsEngine.getBoolean("-d")) {
				this.distFile = myArgsEngine.getString("-d");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the distance file.");
			}
			
			if(myArgsEngine.getBoolean("-n")) {
				this.bcn = Integer.parseInt(myArgsEngine.getString("-n"));
			}
			
			if (myArgsEngine.getBoolean("-q")) {
				this.minQual = Integer.parseInt(myArgsEngine.getString("-q"));
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
		case dist:
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
			if(myArgsEngine.getBoolean("-k")) {
				this.minLink = Integer.parseInt(myArgsEngine.getString("-k"));
			}
			if(myArgsEngine.getBoolean("-a")) {
				this.minRatio = Double.parseDouble(myArgsEngine.getString("-a"));
			}
			if(myArgsEngine.getBoolean("-r")) {
				this.radius = Integer.parseInt(myArgsEngine.getString("-r"));
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
						!samRecord.isSecondaryOrSupplementary() ) {
					if(samRecord.getFirstOfPairFlag())
						records[0] = samRecord;
					else if(samRecord.getSecondOfPairFlag())
						records[1] = samRecord;
					else
						throw new RuntimeException("!!!");
				}

				while( (samRecord = iter.hasNext() ? iter.next() : null)!=null && samRecord.getReadName().equals(sn)) {
					if( !samRecord.getReadUnmappedFlag() &&
							!samRecord.isSecondaryOrSupplementary()) {
						if(samRecord.getFirstOfPairFlag())
							records[0] = samRecord;
						else if(samRecord.getSecondOfPairFlag())
							records[1] = samRecord;
						else
							throw new RuntimeException("!!!");
					}
				}

				if(records[0]!=null && records[1]!=null &&
						records[0].getReferenceIndex().intValue()==records[1].getReferenceIndex().intValue() &&
						records[0].getReadNegativeStrandFlag()!=records[1].getReadNegativeStrandFlag() && 
						(records[0].getReadNegativeStrandFlag()&&records[0].getAlignmentStart()>records[1].getAlignmentStart() || 
								records[1].getReadNegativeStrandFlag()&&records[1].getAlignmentStart()>records[0].getAlignmentEnd())) {
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

	private class Molecule {
		private List<SAMRecord[]> reads_list;
		private String chr_id = null;
		private int chr_start = -1;
		private int chr_end = -1;
		private boolean reverse = false;
		
		public Molecule() {
			this.reads_list = new ArrayList<SAMRecord[]>();
		}

		public void add(SAMRecord[] record) {
			this.reads_list.add(record);
		}

		public void construct() {
			// TODO Auto-generated method stub
			this.chr_id = this.reads_list.get(0)[0].getReferenceName();
			int chr_start = Integer.MAX_VALUE;
			int chr_end   = Integer.MIN_VALUE;
			for(SAMRecord[] records : this.reads_list) {
				for(SAMRecord record : records) {
					if(record.getAlignmentStart()<chr_start)
						chr_start = record.getAlignmentStart();	
					if(record.getAlignmentEnd()>chr_end)
						chr_end   = record.getAlignmentEnd();
				}
			}
			this.chr_start = chr_start-1;
			this.chr_end   = chr_end;
		}

		public boolean passQualityControl() {
			// TODO Auto-generated method stub
			// 1. need at least 'bcn' read pairs
			if(this.reads_list.size()<bcn) 
				return false;
			// 2. need all read pairs aligned to the same direction
			this.reverse = this.reads_list.get(0)[0].getReadNegativeStrandFlag();
			for(SAMRecord[] rs : this.reads_list) {
				if(rs[0].getReadNegativeStrandFlag()!=this.reverse)
					return false;
			}
			// 3. need max quality score not smaller than 'minQual'
			double qual = .0, q;
			for(SAMRecord[] rs : this.reads_list) {
				q = rs[0].getMappingQuality()+rs[1].getMappingQuality();
				if(qual<q) qual = q;
			}
			if(qual/2<minQual) return false;
			return true;
		}
	}
	
	private double middlePoint(SAMRecord[] record) {
		// TODO Auto-generated method stub
		return (double)(Math.min(record[0].getAlignmentStart(), record[1].getAlignmentStart())+
				Math.max(record[0].getAlignmentEnd(), record[1].getAlignmentEnd()))/2;
	}
	
	private List<Molecule> extractMoleculeFromList(List<SAMRecord[]> list) {
		// TODO Auto-generated method stub
		List<Molecule> mols = new ArrayList<Molecule>();
		
		if(list.isEmpty()) return mols; 
		
		Collections.sort(list, new Comparator<SAMRecord[]>() {
			@Override
			public int compare(SAMRecord[] record0, SAMRecord[] record1) {
				// TODO Auto-generated method stub
				int f = record0[0].getReferenceIndex().intValue()-
						record1[0].getReferenceIndex().intValue();
				return f==0 ? Double.compare(middlePoint(record0), middlePoint(record1)) : f;
			}
		});
		
		Iterator<SAMRecord[]> iter = list.iterator();
		SAMRecord[] records = iter.next();
		int chr_index = records[0].getReferenceIndex();
		double chr_pos = middlePoint(records);
		Molecule mol = new Molecule();
		mol.add(records);
		while(iter.hasNext()) {
			records = iter.next();
			if(records[0].getReferenceIndex()!=chr_index ||
					middlePoint(records)-chr_pos>maxGap) {
				mol.construct();
				if(mol.passQualityControl()) 
					mols.add(mol);
				mol = new Molecule();		
			}
			chr_index = records[0].getReferenceIndex();
			chr_pos = middlePoint(records);
			mol.add(records);
		}
		mol.construct();
		if(mol.passQualityControl())
			mols.add(mol);
		
		return mols;
	}
	
	public void run_arcs() {
		// TODO Auto-generated method stub
		final Set<String> bc_white = new HashSet<String>();
		final Map<String, Map<String, Double>> dist = new HashMap<>();
		try {
			String line;
			final BufferedReader br_bc = Utils.getBufferedReader(this.bcList);
			while( (line=br_bc.readLine()) != null ) 
				bc_white.add(line.split(",")[0]);
			br_bc.close();
			
			final BufferedReader br_ds = Utils.getBufferedReader(this.distFile);
			String[] s, s2;
			while( (line=br_ds.readLine()) !=null ) {
				s = line.trim().split("\\s+");
				final Map<String, Double> ds = new HashMap<>();
				for(int i=1; i<s.length; i++) {
					s2 = s[i].trim().split(",");
					ds.put(s2[0], Double.parseDouble(s2[1]));
				}
				dist.put(s[0], ds);
			}
			br_ds.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		final SamReader samReader = factory.open(new File(this.bamFiles[0]));
		final SAMSequenceDictionary seqdict = samReader.getFileHeader().getSequenceDictionary();
		final int nSeq = seqdict.size();
		
		try {
			samReader.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		final Map<Long, Integer> arcs = new HashMap<Long, Integer>();
		
		this.initial_thread_pool();
		for(final String bamFile : bamFiles) {
			
            this.executor.submit(new Runnable() {
				private String bamFile;

				@Override
				public void run() {
					// TODO Auto-generated method stub
					myLogger.info("####Processing BAM file "+bamFile+".");
			        
					try {
						final BAMBarcodeIterator iter = new BAMBarcodeIterator(this.bamFile);
						List<SAMRecord[]> bc_records;
						List<Molecule> mols;
						String ref, ref1, ref2;
						long refind;

						while(iter.hasNext()) {
							bc_records = iter.next();
							if(bc_records.isEmpty() || 
									!bc_white.contains(bc_records.get(0)[0].getStringAttribute("BX")))
								continue;
							mols = extractMoleculeFromList(bc_records);
							
							final DirectedPseudograph<String, DefaultEdge> g = new DirectedPseudograph<>(DefaultEdge.class);
							for(Molecule mol : mols) {
								ref = mol.chr_id;
								if(!g.containsVertex(ref))
									g.addVertex(ref);
								if(!g.containsVertex(ref+"'"))
									g.addVertex(ref+"'");
							}
							
							for(Molecule mol1 : mols) {
								ref1 = mol1.chr_id;
								if(mol1.reverse) 
									ref1 += "'";
								for(Molecule mol2 : mols) {
									ref2 = mol2.chr_id;
									if(mol2.reverse)
										ref2 += "'";
									
								}
							}
						}
						iter.close();

                        myLogger.info("####BAM file "+bamFile+" processed.");
					} catch (Exception e) {
						Thread t = Thread.currentThread();
						t.getUncaughtExceptionHandler().uncaughtException(t, e);
						e.printStackTrace();
						executor.shutdown();
						System.exit(1);
					}
				}
				
				public Runnable init(String bamFile) {
					// TODO Auto-generated method stub
					this.bamFile = bamFile;
					return this;
				}

			}.init(bamFile));
		}
		this.waitFor();
		
		try {
			BufferedWriter bw = Utils.getBufferedWriter(this.out);
			String source, target;
			int tk;
			for(long key : arcs.keySet()) {
				tk = (int) key;
				target = seqdict.getSequence(tk<nSeq?tk:tk-nSeq).getSequenceName()+(tk<nSeq?"":"'");
				source = seqdict.getSequence((int)(key>>32)).getSequenceName();
				bw.write(key+"\t"+source+"\t"+target+"\t"+arcs.get(key)+"\n");
			}
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
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
								if( !tmp.isSecondaryOrSupplementary() &&
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
	
	private void run_dist() {
		// TODO Auto-generated method stub
		
		final DirectedWeightedPseudograph<String, DefaultWeightedEdge> link_graph = this.makeLinkGraph();
		this.trimLinkGraph(link_graph);
		
		final Map<String, Sequence> sequences = Sequence.parseFastaFileWithRevCmpAsMap(this.seqFile);
		
		final Map<String, Map<String, Double>> distMat = new HashMap<>();
		
		this.initial_thread_pool();
		for(final String source : link_graph.vertexSet()) {
			executor.submit(new Runnable() {

				private String root;
				
				@Override
				public void run() {
					// TODO Auto-generated method stub
					myLogger.info("####calulate distance array for "+root);

					final Map<String, Double> dist = new HashMap<>();
					final Map<String, Double> distPool = new HashMap<>();
					distPool.put(root, .0);
					
					final Deque<String> queue = new ArrayDeque<>();
					queue.push(root);
					
					String source, target;
					double distance, d;
					while(!queue.isEmpty()) {
						source = queue.pop();
						distance = distPool.get(source);
						for(final DefaultWeightedEdge out : link_graph.outgoingEdgesOf(source)) {
							target = link_graph.getEdgeTarget(out);
							if(!dist.containsKey(target)||dist.get(target)>distance)
								dist.put(target, distance);
							d = distance+sequences.get(target).seq_ln();
							if(d>radius) continue;
							if(!distPool.containsKey(target)||distPool.get(target)>d) {
								distPool.put(target, d);
								queue.push(target);
							}
						}
					}
					
					synchronized(lock) {
						distMat.put(root, dist);
					}
				}

				public Runnable init(String root) {
					// TODO Auto-generated method stub
					this.root = root;
					return this;
				}
				
			}.init(source));
		}
		this.waitFor();
		
		
		try {
			BufferedWriter bw = Utils.getBufferedWriter(this.out);
			StringBuilder os = new StringBuilder();
			for(Map.Entry<String, Map<String, Double>> entry : distMat.entrySet()) {
				os.setLength(0);
				os.append(entry.getKey());
				for(Map.Entry<String, Double> entry2 : entry.getValue().entrySet()) {
					os.append("\t");
					os.append(entry2.getKey());
					os.append(",");
					os.append((int) entry2.getValue().doubleValue());
				}
				os.append("\n");
				bw.write(os.toString());
			}
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
		return;
	}
	
	private void trimLinkGraph(DirectedWeightedPseudograph<String, DefaultWeightedEdge> link_graph) {
		// TODO Auto-generated method stub
		final Set<DefaultWeightedEdge> edges = new HashSet<>(link_graph.edgeSet());
		for(final String v : link_graph.vertexSet()) {
			trim(link_graph, v, true,  edges);
			trim(link_graph, v, false, edges);
		}
		link_graph.removeAllEdges(edges);
		myLogger.info("#### link graph trimming done.");
		myLogger.info("    +V "+link_graph.vertexSet().size());
		myLogger.info("    -E "+edges.size());
		myLogger.info("    +E "+link_graph.edgeSet().size());
		myLogger.info("########");
		return;
	}

	private void trim(DirectedWeightedPseudograph<String, DefaultWeightedEdge> link_graph, String v, boolean out,
			Set<DefaultWeightedEdge> edgeSet) {
		// TODO Auto-generated method stub
		DefaultWeightedEdge[] edges = out ? 
				link_graph.outgoingEdgesOf(v).toArray(new DefaultWeightedEdge[link_graph.outDegreeOf(v)]) :
					link_graph.incomingEdgesOf(v).toArray(new DefaultWeightedEdge[link_graph.inDegreeOf(v)]);
		if(edges.length==0) return;
		double[] weights = new double[edges.length];
		double maxW = 0, w;
		for(int i=0; i<weights.length; i++) {
			w = link_graph.getEdgeWeight(edges[i]);
			weights[i] = w;
			if(maxW<w) maxW = w;
		}
		if(maxW==0) throw new RuntimeException("!!!");
		for(int i=0; i<weights.length; i++) {
			if(weights[i]/maxW>=this.minRatio) 
				edgeSet.remove(edges[i]);
		}
		return;
	}

	private DirectedWeightedPseudograph<String, DefaultWeightedEdge> makeLinkGraph() {
		// TODO Auto-generated method stub
		final DirectedWeightedPseudograph<String, DefaultWeightedEdge> link_graph = new 
				DirectedWeightedPseudograph<>(DefaultWeightedEdge.class);
		try {
			String line;
			String[] s;
			BufferedReader br_seqs = Utils.getBufferedReader(seqFile);
			String seq;
			while( (line=br_seqs.readLine())!=null ) {
				if(line.startsWith(">")) {
					seq = line.split("\\s+")[0].substring(1);
					link_graph.addVertex(seq);
					link_graph.addVertex(seq+"'");
				}
			}
			br_seqs.close();
			
			BufferedReader br_link = Utils.getBufferedReader(linkFile);
			int weight;
			DefaultWeightedEdge edge;
			while( (line=br_link.readLine())!=null ) {
				s = line.split("\\s+");
				weight = Integer.parseInt(s[3]);
				if(weight<this.minLink) continue;
				edge = link_graph.addEdge(s[1], s[2]);
				link_graph.setEdgeWeight(edge, weight);
			}
			br_link.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		myLogger.info("#### link graph construction done.");
		myLogger.info("    +V "+link_graph.vertexSet().size());
		myLogger.info("    +E "+link_graph.edgeSet().size());
		myLogger.info("########");
		return link_graph;
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
		case dist:
			this.run_dist();
			break;
		default:
			throw new RuntimeException("!!!");
		}
		return;
	}

}





