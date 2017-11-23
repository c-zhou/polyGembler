package cz1.tenx.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.zip.GZIPOutputStream;

import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;

import org.renjin.eval.Context;
import org.renjin.primitives.io.serialization.RDataWriter;
import org.renjin.primitives.matrix.DoubleMatrixBuilder;
import org.renjin.sexp.ListVector;

import cz1.ngs.model.Sequence;
import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class TenXMoleculeStats extends Executor {

	private static final boolean debug = false;
	
	private String in_bam;
	private String in_vcf;
	private String out_stats;
	private boolean call_snps = false;
	// minimum minor allele frequency
	private final double min_maf = 0.1;
	private int min_depth = 30;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--input-bam      Input BAM file.\n"
						+ " -f/--input-vcf      Input VCF file.\n"
						+ " -s/--call-snps      Call SNPs.\n"
						+ " -d/--min-depth      Minimum allele depth (default 30).\n"
						+ " -t/--threads        Threads (default is 1 and will only use ).\n"
						+ " -o/--output-stats   Output file.\n\n");
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
			myArgsEngine.add("-i", "--input-bam", true);
			myArgsEngine.add("-f", "--input-vcf", true);
			myArgsEngine.add("-s", "--call-snps", false);
			myArgsEngine.add("-d", "--min-depth", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-o", "--output-stats", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-i")) {
			in_bam = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your BAM file.");
		}
		
		if (myArgsEngine.getBoolean("-s")) {
			call_snps = true;
		}
		
		if (myArgsEngine.getBoolean("-f")) {
			in_vcf = myArgsEngine.getString("-f");
		} else if(call_snps) {
			printUsage();
			throw new IllegalArgumentException("Please specify your VCF file.");
		}
		
		if (myArgsEngine.getBoolean("-d")) {
			min_depth = Integer.parseInt(myArgsEngine.getString("-d"));
		}

		if (myArgsEngine.getBoolean("-t")) {
			THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}
		
		if (myArgsEngine.getBoolean("-o")) {
			out_stats = myArgsEngine.getString("-o");
		}
	}
	
	// maximum gap size between two barcodes
	private static final int max_gap = 10000;
	// at least a barcode per
	private static final int abc_per = 10000;
	// minimum size of molecules to study
	private static final int oo_thresh = 10000;
	// maximum insert size for read pairs
	private static final int max_ins = 1000;
	private static BufferedWriter oos;
	private static BufferedWriter vcf;
	private static BufferedWriter hap;
	private static SAMFileWriter oob;
	private static final Object lock = new Object(),
			lockoos = new Object();
	private static long counter = 0, last_counter = 0;
	
	// statistics
	private long[] mappedReadsOnChromosome;
	private long     mappedReads = 0;
	private long   unmappedReads = 0;
	private long       badpaired = 0;
	private long duplicatedReads = 0;
	private final Map<String, Integer> readsOnBarcode = new HashMap<String, Integer>();
	
	private final SamReaderFactory factory =
			SamReaderFactory.makeDefault()
			.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
					SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
			.validationStringency(ValidationStringency.SILENT);
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		if(call_snps) runSNPCaller();
		else runBCSortedStats();
		this.hist();
	}
	
	private void runSNPCaller() {
		// TODO Auto-generated method stub
		final SamReader inputSam = factory.open(new File(this.in_bam));
		if(inputSam.getFileHeader().getSortOrder()!=SAMFileHeader.SortOrder.coordinate) 
			throw new RuntimeException("Sort BAM file before SNP calling!!!");
		final SAMSequenceDictionary seqDic = inputSam.getFileHeader().getSequenceDictionary();
		mappedReadsOnChromosome = new long[seqDic.size()];
		
		oos = Utils.getBufferedWriter(out_stats+".txt");
		hap = Utils.getBufferedWriter(out_stats+".haps");
		vcf = Utils.getBufferedWriter(out_stats+".snps");
		SAMFileHeader header = inputSam.getFileHeader();
		header.setSortOrder(SortOrder.unknown);
		oob = new SAMFileWriterFactory().
				makeSAMOrBAMWriter(header, true, new File(out_stats+".bam"));

		try {
			inputSam.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		this.initial_thread_pool();
		
		for(int i=0; i<seqDic.size(); i++) {
			executor.submit(new Runnable() {
				private int chr_num;
				
				@Override
				public void run() {
					// TODO Auto-generated method stub
					try {
						final SAMSequenceRecord refSeq = seqDic.getSequence(chr_num);
						final String seqnm = refSeq.getSequenceName();
						final int seqln = refSeq.getSequenceLength();
						
						final Set<Integer> snp_pos = readSNPs(in_vcf, seqnm);
						
						final SamReader inputSam = factory.open(new File(in_bam));
						final SAMRecordIterator iter = inputSam.queryOverlapping(seqnm, 0, seqln);
						
						
						// do NOT call SNPs, read it from file instead
						final Map<String, List<SAMRecord>> bc_record = new HashMap<String, List<SAMRecord>>();
						final Map<String, Integer> bc_position = new HashMap<String, Integer>();
						String bc;
						SAMRecord tmp_record;
						int pos;
						while(iter.hasNext()) {
							tmp_record = iter.next();
							if( !processRecord(tmp_record) ) 
								continue;
							bc = tmp_record.getStringAttribute("BX");
							if(bc==null) continue;
							
							pos = tmp_record.getAlignmentStart();
							
							if(mappedReads%10000==0) {
								Set<String> bc_removal = new HashSet<String>();
								for(Map.Entry<String, List<SAMRecord>> entry : bc_record.entrySet()) {
									String b = entry.getKey();
									List<SAMRecord> records = entry.getValue();
									if(bc_position.get(b)+max_gap<pos) {
										processMolecule(b, records, snp_pos);
										bc_removal.add(b);
									}
								}
								for(String b : bc_removal) {
									bc_record.remove(b);
									bc_position.remove(b);
								}
							}
							
							if(bc_record.containsKey(bc)) {
								List<SAMRecord> records = bc_record.get(bc);
								if(bc_position.get(bc)+max_gap<pos) {
									processMolecule(bc, records, snp_pos);
									records.clear();
									bc_position.remove(bc);
								}
							} else 
								bc_record.put(bc, new ArrayList<SAMRecord>());
							bc_record.get(bc).add(tmp_record);
							bc_position.put(bc, Math.max(tmp_record.getAlignmentEnd(), 
									bc_position.containsKey(bc)?bc_position.get(bc):0));
						}
						
						if(!bc_record.isEmpty()) 
							for(Map.Entry<String, List<SAMRecord>> entry : bc_record.entrySet()) 
								processMolecule(entry.getKey(), entry.getValue(), snp_pos);
						
						/***
						// do NOT call SNPs, read it from file instead
						int pos;
						SAMRecord buffer = iter.hasNext() ? iter.next() : null;
						int bufferS = buffer==null ? Integer.MAX_VALUE : buffer.getAlignmentStart();
						final Set<SAMRecord> record_pool = new HashSet<SAMRecord>();
						final Set<SAMRecord> tmp_removal = new HashSet<SAMRecord>();
						
						final int[] allele_stats = new int[5];
						
						final Map<String, TreeSet<SNP>> bc_snp = new HashMap<String, TreeSet<SNP>>();
						
						String dna_seq, bc;
						int tmp_int;
						String snp_str;
						char nucl;
						for(int p=0; p!=seqln; p++) {
							if(p%1000000==0) {
								System.err.println(seqnm+":"+p);
								// TODO
								// process bc_snp
								Set<String> bc_removal = new HashSet<String>();
								for(Map.Entry<String, TreeSet<SNP>> entry : bc_snp.entrySet()) {
									if(entry.getValue().last().position+max_gap<p) {
										processMolecule(entry.getKey(), entry.getValue());
										bc_removal.add(entry.getKey());
									}
								}
								for(String b : bc_removal) bc_snp.remove(b);
							}
							pos = p+1; // BAM format is 1-based coordination
							
							while(bufferS<=pos) {
								// TODO
								// buffer records upto this position
								while( buffer!=null && 
										buffer.getAlignmentStart()<=pos) {
									if( processRecord(buffer) ) 
										record_pool.add(buffer);
									buffer=iter.hasNext()?iter.next():null;
									bufferS = buffer==null ? Integer.MAX_VALUE : 
										buffer.getAlignmentStart();
								}
							}
							if(record_pool.isEmpty()) continue;
							
							Arrays.fill(allele_stats, 0);
							tmp_removal.clear();
							for(SAMRecord record : record_pool) {

								if(record.getAlignmentEnd()<pos) {
									tmp_removal.add(record);
									continue;
								}
								dna_seq = record.getReadString();
								tmp_int = record.getReadPositionAtReferencePosition(pos);
								nucl = tmp_int==0?'D':Character.toUpperCase(dna_seq.charAt(tmp_int-1));
								switch(nucl) {
								case 'A':
									allele_stats[0]++;
									break;
								case 'C':
									allele_stats[1]++;
									break;
								case 'G':
									allele_stats[2]++;
									break;
								case 'T':
									allele_stats[3]++;
									break;
								case 'D':
									allele_stats[4]++;
									break;
								case 'N':
									// do nothing here
									break;
								default:
									throw new RuntimeException("!!!");
								}
							}
							record_pool.removeAll(tmp_removal);
							
							snp_str = feasibleSNP(allele_stats);
							
							// not two alleles or minor allele frequency is smaller than threshold
							// don't call indels
							if(snp_str==null) continue;
							
							vcf.write(seqnm+"\t"+p+"\t"+snp_str+"\n");
							
							for(SAMRecord record : record_pool) {
								bc = record.getStringAttribute("BX");
								if(bc==null) continue;
								dna_seq = record.getReadString();
								tmp_int = record.getReadPositionAtReferencePosition(pos);
								if(tmp_int==0) continue;
								nucl = Character.toUpperCase(dna_seq.charAt(tmp_int-1));
								if(bc_snp.containsKey(bc)) { 
									TreeSet<SNP> bc_set = bc_snp.get(bc);
									if(bc_set.last().position+max_gap<p) {
										// TODO
										// gap exceeds max_gap
										// write molecule out
										processMolecule(bc, bc_set);
										bc_set.clear();
									}
								} else {
									bc_snp.put(bc, new TreeSet<SNP>(new SNP.PositionComparator()));
								}
								bc_snp.get(bc).add(new SNP(p, nucl, record));
								
							}
						}
						
						for(Map.Entry<String, TreeSet<SNP>> entry : bc_snp.entrySet()) 
							processMolecule(entry.getKey(), entry.getValue());
						**/
						
						iter.close();
						inputSam.close();
					} catch (Exception e) {
						Thread t = Thread.currentThread();
						t.getUncaughtExceptionHandler().uncaughtException(t, e);
						e.printStackTrace();
						executor.shutdown();
						System.exit(1);
					}
				}

				private void processMolecule(final String bc, 
						final List<SAMRecord> records, 
						final Set<Integer> snp_pos) {
					// TODO Auto-generated method stub
					MoleculeConstructor molCon  = new MoleculeConstructor(records);
					if(molCon.getMolecularLength()<oo_thresh) return;
					molCon.run();
				}

				/**
				private void processMolecule(final String bc, final TreeSet<SNP> snp_set) 
						throws IOException {
					// TODO Auto-generated method stub
					// write molecule haplotypes
					Set<SAMRecord> records = new HashSet<SAMRecord>();
					List<String> snp_str = new ArrayList<String>();
					SNP snp;
					int tmp_position = -1;
					while(!snp_set.isEmpty()) {
						snp = snp_set.pollFirst();
						records.add(snp.alignment);
						if(tmp_position==snp.position) {
							snp_str.remove(snp_str.size()-1);
						} else {
							snp_str.add(snp.position+","+
									snp.allele);
							tmp_position = snp.position;
						}
					}
					// write molecule statistics
					MoleculeConstructor molCon  = 
							new MoleculeConstructor(new ArrayList<SAMRecord>(records));
					if(molCon.getMolecularLength()<oo_thresh) return;
					molCon.run();
					
					if(snp_str.size()<2) return;
					StringBuilder os = new StringBuilder();
					os.append(bc);
					for(String ss : snp_str) {
						os.append(" ");
						os.append( ss);
					}
					os.append("\n");
					hap.write(os.toString());
				}
				**/

				public Runnable init(final int i) {
					// TODO Auto-generated method stub
					this.chr_num = i;
					return this;
				}

				private String feasibleSNP(int[] ints) {
					// TODO Auto-generated method stub
					if(ints[4]>0) return null;
					List<Integer> obs = new ArrayList<Integer>();
					for(int i=0; i<4; i++)
						if(ints[i]>0) obs.add(i);
					if(obs.size()!=2) return null;
					int ref = obs.get(0), alt = obs.get(1);
					int d = ints[ref]+ints[alt];
					if( d < min_depth) return null;
					double maf = (double) ints[ref]/d;
					if(maf>0.5) maf = 1-maf;
					if(maf<min_maf) return null;
					return Sequence.nucleotide[ref]+","+
							Sequence.nucleotide[alt]+"\t"+
							ints[ref]+","+ints[alt];
				}
			}.init(i));
		}

		try {
			inputSam.close();
			this.waitFor();
			oos.close();
			hap.close();
			vcf.close();
			oob.close();
		
			BufferedWriter bw_summary = Utils.getBufferedWriter(this.out_stats+".summary");
			bw_summary.write("##Reads     mapped: "+mappedReads+"\n");
			bw_summary.write("##Reads   unmapped: "+unmappedReads+"\n");
			bw_summary.write("##Reads duplicated: "+duplicatedReads+"\n");
			bw_summary.write("##Reads  bad pairs: "+badpaired+"\n");
			bw_summary.write("##Reads by pseudo-chromosomes:\n");
			for(int i=0; i<mappedReadsOnChromosome.length; i++) {
				bw_summary.write(" #"+seqDic.getSequence(i).getSequenceName()+" "+
						mappedReadsOnChromosome[i]+"\n");
			}
			bw_summary.write("##Reads by barcodes:\n");
			for(Map.Entry<String, Integer> entry : readsOnBarcode.entrySet()) {
				bw_summary.write(" #"+entry.getKey()+" "+entry.getValue()+"\n");
			}
			bw_summary.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	protected Set<Integer> readSNPs(String in_vcf2, String seqnm) {
		// TODO Auto-generated method stub
		return null;
	}

	public void runBCSortedStats() {
		// TODO Auto-generated method stub
		
		final SamReader inputSam = factory.open(new File(this.in_bam));
		final SAMSequenceDictionary seqDic = inputSam.getFileHeader().getSequenceDictionary();
		
		mappedReadsOnChromosome = new long[seqDic.size()];
		/***
		try {
			oos = Utils.getGZIPBufferedWriter(out_stats);
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		**/
		oos = Utils.getBufferedWriter(out_stats+".txt");
		oob = new SAMFileWriterFactory().
				makeSAMOrBAMWriter(inputSam.getFileHeader(), true, new File(out_stats+".bam"));
		
		SAMRecordIterator iter = inputSam.iterator();
		SAMRecord record = iter.next();
		String bc_str;		
		
		this.initial_thread_pool();
		while( record!=null ) {
			
			// if(mappedReads%1000000==0) myLogger.info(mappedReads+" mapped records processed.");

			List<SAMRecord> bc_records = new ArrayList<SAMRecord>();
			
			bc_str = record.getStringAttribute("BX");
			
			this.processRecord(bc_records, record);
			
			while(bc_str==null) {
				record = iter.hasNext() ? iter.next() : null;
				if(record==null) break;
				this.processRecord(bc_records, record);
				bc_str = record.getStringAttribute("BX");
			}
			
			if(bc_str==null) continue;
			
			while(true) {
				record = iter.hasNext() ? iter.next() : null;
				if(record==null) break;
				if(bc_str.equals(record.getStringAttribute("BX"))) {
					this.processRecord(bc_records, record);
				} else break;
			}
			
			if(bc_records.isEmpty()) continue;
			
			readsOnBarcode.put(bc_records.get(0).getStringAttribute("BX"), bc_records.size());
			
			executor.submit(new MoleculeConstructor(bc_records));
			
		}
		try {
			iter.close();
			inputSam.close();
			this.waitFor();
			oos.close();
			oob.close();
		
			BufferedWriter bw_summary = Utils.getBufferedWriter(this.out_stats+".summary");
			bw_summary.write("##Reads     mapped: "+mappedReads+"\n");
			bw_summary.write("##Reads   unmapped: "+unmappedReads+"\n");
			bw_summary.write("##Reads duplicated: "+duplicatedReads+"\n");
			bw_summary.write("##Reads  bad pairs: "+badpaired+"\n");
			bw_summary.write("##Reads by pseudo-chromosomes:\n");
			for(int i=0; i<mappedReadsOnChromosome.length; i++) {
				bw_summary.write(" #"+seqDic.getSequence(i).getSequenceName()+" "+
						mappedReadsOnChromosome[i]+"\n");
			}
			bw_summary.write("##Reads by barcodes:\n");
			for(Map.Entry<String, Integer> entry : readsOnBarcode.entrySet()) {
				bw_summary.write(" #"+entry.getKey()+" "+entry.getValue()+"\n");
			}
			bw_summary.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		myLogger.info("Completed. "+counter+" bam records processed.");
	}
	
	private boolean processRecord(SAMRecord record) {
		// TODO Auto-generated method stub
        if(record.getStringAttribute("BX")==null) return false;

		if(record.getDuplicateReadFlag()&&
				(record.getReadUnmappedFlag()||
						!record.getNotPrimaryAlignmentFlag()&&
						!record.getSupplementaryAlignmentFlag()) ) {
			synchronized (lock) {++duplicatedReads;}
			return false;
		} else if(record.getReadUnmappedFlag()) {
			// unmapped reads
			synchronized (lock) {++unmappedReads;}
			return false;
		} else if( record.getNotPrimaryAlignmentFlag() || 
				record.getSupplementaryAlignmentFlag() ) {
			// not primary or supplementary alignment
			return false;
		} else {
			synchronized (lock) {++mappedReadsOnChromosome[record.getReferenceIndex()];}
			synchronized (lock) {++mappedReads;}
			
			if(record.getMateUnmappedFlag() || 
				record.getReferenceIndex()!=
				record.getMateReferenceIndex() || 
				Math.abs(record.getInferredInsertSize())>max_ins) {
				// pair not properly aligned 
				synchronized (lock) {++badpaired;}
				return false;
			} else return true;
		}
	}
	
	private void processRecord(List<SAMRecord> record_list, SAMRecord record) {
		// TODO Auto-generated method stub
		if(processRecord(record)) record_list.add(record);
	}

	private void r() {
		// TODO Auto-generated method stub
		ScriptEngineManager manager = new ScriptEngineManager();
		ScriptEngine engine = manager.getEngineByName("Renjin"); 
		if(engine == null) { 
			throw new RuntimeException("Renjin not found!!!"); 
		}
		try {

			Process p = Runtime.getRuntime().exec("wc -l " + out_stats);
			p.waitFor();
			BufferedReader pbr = new BufferedReader(new InputStreamReader(p.getInputStream()));
			String line;
			int mol_n = 0;
			while ((line = pbr.readLine()) != null) 
				mol_n = Integer.parseInt(line.split("\\s+")[0]);
			pbr.close();
			
			DoubleMatrixBuilder mol_len = new DoubleMatrixBuilder(1, mol_n);
			final Map<String, Integer> bc_stats = new HashMap<String, Integer>();
			BufferedReader br = Utils.getBufferedReader(out_stats);
			String[] s;
			int i = 0;
			while( (line=br.readLine())!=null ) {
				s = line.split("\\s+");
				if(bc_stats.containsKey(s[0])) 
					bc_stats.put(s[0], bc_stats.get(s[0])+1);
				else bc_stats.put(s[0], 1);
				mol_len.set(0, i++, Integer.parseInt(s[2]));
			}
			br.close();
			
			int bc_n = bc_stats.size();
			DoubleMatrixBuilder bc_nstats = new DoubleMatrixBuilder(1, bc_n);
			i = 0;
			for(Map.Entry<String, Integer> entry : bc_stats.entrySet()) 
				bc_nstats.set(0, i++, entry.getValue());
			
			Context context = Context.newTopLevelContext();
			FileOutputStream fos = new FileOutputStream(out_stats+".RData");
			GZIPOutputStream zos = new GZIPOutputStream(fos);
			RDataWriter writer = new RDataWriter(context, zos);
			
			ListVector.NamedBuilder Rdat = new ListVector.NamedBuilder();
			Rdat.add("mol_stats", mol_len.build());
			Rdat.add("bc_stats", bc_nstats.build());
			writer.save(Rdat.build());
			writer.close();
			
		} catch (IOException | InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private void hist() {
		// TODO Auto-generated method stub
		
		try {
			if(!new File(out_stats+".txt").exists()) return;
			BufferedReader br = Utils.getBufferedReader(out_stats+".txt");
			String line;
			String[] s;
			final Map<Integer, Integer> mol_hist = new HashMap<Integer, Integer>();
			final Map<Integer, Long> mol_L = new HashMap<Integer, Long>();
			final Map<String, Integer> bc_hist = new HashMap<String, Integer>();
			mol_hist.put(100,0);
			mol_hist.put(200,0);
			mol_hist.put(500,0);
			mol_hist.put(1000,0);
			mol_hist.put(2000,0);
			mol_hist.put(5000,0);
			mol_hist.put(10000,0);
			mol_hist.put(20000,0);
			mol_hist.put(50000,0);
			mol_hist.put(100000,0);
			mol_hist.put(Integer.MAX_VALUE,0);
			mol_L.put(100,0L);
			mol_L.put(200,0L);
			mol_L.put(500,0L);
			mol_L.put(1000,0L);
			mol_L.put(2000,0L);
			mol_L.put(5000,0L);
			mol_L.put(10000,0L);
			mol_L.put(20000,0L);
			mol_L.put(50000,0L);
			mol_L.put(100000,0L);
			mol_L.put(Integer.MAX_VALUE,0L);
			int mol_len;
			counter = 0;
			while( (line=br.readLine())!=null ) {
				if(++counter%1000000==0) myLogger.info(counter+" items processed.");;
				s = line.split("\\s+");
				if(bc_hist.containsKey(s[0])) 
					bc_hist.put(s[0], bc_hist.get(s[0])+1);
				else bc_hist.put(s[0], 1);
				if(s[1].startsWith("Chr00")) continue;
				mol_len = Integer.parseInt(s[2]);
				if(mol_len<100) {
					mol_hist.put(100, mol_hist.get(100)+1);
					mol_L.put(100, mol_L.get(100)+mol_len);
				} else if(mol_len<200) {
					mol_hist.put(200, mol_hist.get(200)+1);
					mol_L.put(200, mol_L.get(200)+mol_len);
				} else if(mol_len<500) {
					mol_hist.put(500, mol_hist.get(500)+1);
					mol_L.put(500, mol_L.get(500)+mol_len);
				} else if(mol_len<1000) {
					mol_hist.put(1000, mol_hist.get(1000)+1);
					mol_L.put(1000, mol_L.get(1000)+mol_len);
				} else if(mol_len<2000) {
					mol_hist.put(2000, mol_hist.get(2000)+1);
					mol_L.put(2000, mol_L.get(2000)+mol_len);
				} else if(mol_len<5000) {
					mol_hist.put(5000, mol_hist.get(5000)+1);
					mol_L.put(5000, mol_L.get(5000)+mol_len);
				} else if(mol_len<10000) {
					mol_hist.put(10000, mol_hist.get(10000)+1);
					mol_L.put(10000, mol_L.get(10000)+mol_len);
				} else if(mol_len<20000) {
					mol_hist.put(20000, mol_hist.get(20000)+1);
					mol_L.put(20000, mol_L.get(20000)+mol_len);
				} else if(mol_len<50000) {
					mol_hist.put(50000, mol_hist.get(50000)+1);
					mol_L.put(50000, mol_L.get(50000)+mol_len);
				} else if(mol_len<100000) {
					mol_hist.put(100000, mol_hist.get(100000)+1);
					mol_L.put(100000, mol_L.get(100000)+mol_len);
				} else {
					mol_hist.put(Integer.MAX_VALUE, mol_hist.get(Integer.MAX_VALUE)+1);
					mol_L.put(Integer.MAX_VALUE, mol_L.get(Integer.MAX_VALUE)+mol_len);
				}
			}
			br.close();
			
			BufferedWriter mol_bw = Utils.getBufferedWriter(out_stats+".molL");
			mol_bw.write("<100bp\t"+mol_hist.get(100)+"\t"+mol_L.get(100)+"bp\n");
			mol_bw.write("<200bp\t"+mol_hist.get(200)+"\t"+mol_L.get(200)+"bp\n");
			mol_bw.write("<500bp\t"+mol_hist.get(500)+"\t"+mol_L.get(500)+"bp\n");
			mol_bw.write("<1Kbp\t"+mol_hist.get(1000)+"\t"+mol_L.get(1000)+"bp\n");
			mol_bw.write("<2Kbp\t"+mol_hist.get(2000)+"\t"+mol_L.get(2000)+"bp\n");
			mol_bw.write("<5Kbp\t"+mol_hist.get(5000)+"\t"+mol_L.get(5000)+"bp\n");
			mol_bw.write("<10Kbp\t"+mol_hist.get(10000)+"\t"+mol_L.get(10000)+"bp\n");
			mol_bw.write("<20Kbp\t"+mol_hist.get(20000)+"\t"+mol_L.get(20000)+"bp\n");
			mol_bw.write("<50Kbp\t"+mol_hist.get(50000)+"\t"+mol_L.get(50000)+"bp\n");
			mol_bw.write("<100Kbp\t"+mol_hist.get(100000)+"\t"+mol_L.get(100000)+"bp\n");
			mol_bw.write(">100Kbp\t"+mol_hist.get(Integer.MAX_VALUE)+"\t"+mol_L.get(Integer.MAX_VALUE)+"bp\n");
			mol_bw.close();
			
			BufferedWriter bc_bw = Utils.getBufferedWriter(out_stats+".bc");
			for(Map.Entry<String, Integer> entry : bc_hist.entrySet()) 
				bc_bw.write(entry.getKey()+"\t"+entry.getValue()+"\n");
			bc_bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	

	private class MoleculeConstructor implements Runnable {

		private final List<SAMRecord> bc_records;
		private final int bc_size;
		
		private final StringBuilder os = new StringBuilder();
		private String barcode;
		private String chr;
		private int p0, p1, read_count;
		private final List<Integer> gap = new ArrayList<Integer>();
		private final List<SAMRecord> os_record = new ArrayList<SAMRecord>();

		public MoleculeConstructor(final List<SAMRecord> bc_records) {
			this.bc_records = bc_records;
			Collections.sort(bc_records, new Comparator<SAMRecord>() {
				@Override
				public int compare(SAMRecord r, SAMRecord r2) {
					// TODO Auto-generated method stub
					int f = r.getReferenceIndex() - r2.getReferenceIndex();
					return f==0 ? 
							(r.getAlignmentStart()-r2.getAlignmentStart()) : f;
				}
			});
			this.bc_size = bc_records.size();
		}
		
		public int getMolecularLength() {
				int a = Integer.MAX_VALUE, b = Integer.MIN_VALUE;
				for(SAMRecord r : bc_records) {
					if(r.getAlignmentStart()<a) a=r.getAlignmentStart();
					if(r.getAlignmentEnd()>b) b=r.getAlignmentEnd();
				}
				return b-a+1;
		}
		
		@Override
		public void run() {
			// TODO Auto-generated method stub
			try {
				
				synchronized (lock) {
					counter += bc_records.size();
					if(counter-last_counter>1000000) {
						myLogger.info("["+Thread.currentThread().getName()+"] "+counter+" records processed.");
						last_counter = counter;
					}
				}

				barcode = bc_records.get(0).getStringAttribute("BX");

				SAMRecord record;
				if(debug) {
					myLogger.info("["+Thread.currentThread().getName()+"] "+barcode+", "+bc_records.size());
					for(int i=0; i!=bc_size; i++) {
						record = bc_records.get(i);
						myLogger.info("    "+record.getReferenceName()+": "+record.getAlignmentStart());
					}
				}

				Iterator<SAMRecord> it = bc_records.iterator();
				record = it.next();
				this.initNewMolecule(record);
				int g;
				int molen;

				while( it.hasNext() ) {
					record = it.next();
					if( record.getAlignmentStart()-p1>max_gap || 
							!chr.equals(record.getReferenceName()) ) {
						if( (molen=p1-p0+1)>=oo_thresh && os_record.size()*abc_per>=molen ) {
							this.addNewMoleculeRecord();
							synchronized (lockoos) {
								for(SAMRecord rc : os_record) oob.addAlignment(rc);
							}
						}
						this.initNewMolecule(record);
					} else {
						g = record.getAlignmentStart()-p1;
						gap.add(g);
						p1 = Math.max(record.getAlignmentEnd(), p1);
						os_record.add(record);
						++read_count;
					}
				}

				if( (molen=p1-p0+1)>=oo_thresh && os_record.size()*abc_per>=molen ) {
					this.addNewMoleculeRecord();
					synchronized (lockoos) {
						for(SAMRecord rc : os_record) oob.addAlignment(rc);
					}
				}

				oos.write(os.toString());
			} catch (Exception e) {
				Thread t = Thread.currentThread();
				t.getUncaughtExceptionHandler().uncaughtException(t, e);
				e.printStackTrace();
				executor.shutdown();
				System.exit(1);
			}
		}

		private void initNewMolecule(SAMRecord record) {
			// TODO Auto-generated method stub
			chr = record.getReferenceName();
			p0 = record.getAlignmentStart();
			p1 = record.getAlignmentEnd();
			read_count = 1;
			gap.clear();
			os_record.clear();
			os_record.add(record);
		}

		private void addNewMoleculeRecord() throws IOException {
			// TODO Auto-generated method stub
			int maxG, l, gS;
			double cvg;
			l = p1-p0+1;
			if(gap.isEmpty()) {
				maxG = 0;
				gS = 0;
				cvg = 1;
			} else {
				maxG = gap.get(0);
				gS = 0;
				for(Integer i : gap) {
					if(i>maxG) maxG=i;
					gS += i;
				}
				cvg = (l-gS)/(double) l;
			}
			//os.setLength(0);
			os.append(barcode);
			os.append("\t");
			os.append(chr);
			os.append(":");
			os.append(p0);
			os.append("-");
			os.append(p1);
			os.append("\t");
			os.append(l);
			os.append("\t");
			os.append(read_count);
			os.append("\t");
			os.append(gS);
			os.append("\t");
			os.append(cvg);
			os.append("\t");
			os.append(maxG);
			os.append("\t");
			if(gap.isEmpty()) os.append("-");
			else {
				os.append(gap.get(0));
				int n = gap.size();
				for(int i=1; i!=n; i++) {
					os.append(",");
					os.append(gap.get(i));
				}
			}
			os.append("\n");
			//oos.write(os.toString());
		}
	}

	private static final class SNP {
		private final int position;
		private final char allele;
		
		public SNP(final int position,
				final char allele) {
			this.position = position;
			this.allele = Character.toUpperCase(allele);
		}
		
		@Override
		public boolean equals(Object obj) {
			if ((obj instanceof SNP) && 
					(((SNP) obj).position == this.position) &&
					(((SNP) obj).allele == this.allele) ) {
				return true;
			} else {
				return false;
			}
		}
		
		@Override
	    public int hashCode() {
	        int result = Integer.hashCode(this.position);
	        result <<= 2;
	        switch(this.allele) {
	        case 'A':
	        	// result += 0;
	        	break;
	        case 'C':
	        	result += 1;
	        	break;
	        case 'G':
	        	result += 2;
	        	break;
	        case 'T':
	        	result += 3;
	        	break;
	        default:
	        	throw new RuntimeException("!!!");
	        }
	        return result;
	    }
		
		private static class PositionComparator implements Comparator<SNP> {

			@Override
			public int compare(SNP snp0, SNP snp1) {
				// TODO Auto-generated method stub
				return snp0.position==snp1.position ? 
						snp0.allele-snp1.allele : snp0.position-snp1.position;
			}
		}
	}
	
    public static void main(String[] args) {
        TenXMoleculeStats stats = new TenXMoleculeStats();
        stats.setParameters(new String[]{
                        "-i", "C:\\Users\\chenxi.zhou\\Desktop\\10x_igv\\440166_tanzania_306_2M.bam",
                        "-t", "3",
                        "-s",
                        "-f", "null",
                        "-o", "C:\\Users\\chenxi.zhou\\Desktop\\10x_igv\\out4"
        });
        stats.run();
        //stats.r();
        stats.hist();
    }
}
