package cz1.tenx.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NavigableMap;
import java.util.Set;
import java.util.TreeMap;

import cz1.ngs.model.Sequence;
import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class GenomeComparisonZZ extends Executor {

	private String in_bam1;
	private String in_bam2;
	private String in_ref1;
	private String in_ref2;
	private String in_vcf;
	private String out_prefix;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i1/--input-bam1    Input BAM file 1.\n"
						+ " -i2/--input-bam2    Input BAM file 2.\n"
						+ " -r1/--input-ref1    Input reference genome 1.\n"
						+ " -r2/--input-ref2    Input reference genome 2.\n"
						+ " -x/--variants       Input VCF file.\n"
						+ " -o/--out-prefix     Output file prefix.\n\n");
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
			myArgsEngine.add("-i1", "--input-bam1", true);
			myArgsEngine.add("-i2", "--input-bam2", true);
			myArgsEngine.add("-r1", "--input-ref1", true);
			myArgsEngine.add("-r2", "--input-ref2", true);
			myArgsEngine.add("-x", "--variants", true);
			myArgsEngine.add("-o", "--out-prefix", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-i1")) {
			in_bam1 = myArgsEngine.getString("-i1");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your BAM file.");
		}

		if (myArgsEngine.getBoolean("-i2")) {
			in_bam2 = myArgsEngine.getString("-i2");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your BAM file.");
		}

		if (myArgsEngine.getBoolean("-r1")) {
			in_ref1 = myArgsEngine.getString("-r1");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your reference file.");
		}
		
		if (myArgsEngine.getBoolean("-r2")) {
			in_ref2 = myArgsEngine.getString("-r2");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your reference file.");
		}
		
		if (myArgsEngine.getBoolean("-x")) {
			in_vcf = myArgsEngine.getString("-x");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your VCF file.");
		}
		
		if (myArgsEngine.getBoolean("-o")) {
			out_prefix = myArgsEngine.getString("-o");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your output file prefix.");
		}
	}
	
	private final SamReaderFactory factory =
			SamReaderFactory.makeDefault()
			.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
					SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
			.validationStringency(ValidationStringency.SILENT);
	private final static double min_qual = 20;
	private final static double max_ins  = 1000;
	private final static int max_gap = 10000;
	private static final int min_mol = 10000;
	private static final int abc_per = 10000;
	
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
					if(Math.abs(records[0].getInferredInsertSize())<=max_ins && 
							records[0].getMappingQuality()+records[1].getMappingQuality()>=min_qual) {
						bc_records.add(records);
					}
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
		private List<String> reads_set;
		private List<SAMRecord[]> reads_list;
		private String chr_id = null;
		private int chr_start = -1;
		private int chr_end = -1;
		private int mol_sz = -1;

		public Molecule() {
			this.reads_set = new ArrayList<String>();
			this.reads_list = new ArrayList<SAMRecord[]>();
		}

		public void add(SAMRecord[] record) {
			this.reads_list.add(record);
			this.reads_set.add(record[0].getReadName());
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
			this.mol_sz    = chr_end-chr_start;
		}
	}
	
	private final static class Variant {
		private final String refA;
		private final String altA;
		private final boolean rev;
		private final int altPos;
		
		public Variant(final String refA, final String altA, final int altPos, final boolean rev) {
			this.refA   = refA;
			this.altA   = altA;
			this.altPos = altPos;
			this.rev    = rev;
		}
	}
	
	private final static Map<String, TreeMap<Integer, Variant>> variants = new HashMap<String, TreeMap<Integer, Variant>>();

	private Map<String, Sequence> refSequence1;
	private Map<String, Sequence> refSequence2;
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		
		try {
			refSequence1 = Sequence.parseFastaFileAsMap(in_ref1);
			refSequence2 = Sequence.parseFastaFileAsMap(in_ref2);
			
			final BufferedReader br_vcf = Utils.getBufferedReader(in_vcf);
			String line;
			String[] s;
			long pcount = 0;
			while( (line=br_vcf.readLine())!=null ) {
				if(line.startsWith("#")) continue;
				s = line.split("\\s+");
				if(!variants.containsKey(s[0]))
					variants.put(s[0], new TreeMap<Integer, Variant>());
				variants.get(s[0]).put(Integer.parseInt(s[1]), new Variant(s[3], s[4], Integer.parseInt(s[2]), s[5].equals("1")?false:true));
				++pcount;
				if(pcount%1000000==0) myLogger.info(pcount+" variants loaded into memory.");
			}
			br_vcf.close();
			myLogger.info("Variants loaded from "+in_vcf);
			for(String key : variants.keySet()) 
				myLogger.info(key+": "+variants.get(key).size());
			
			final BAMBarcodeIterator iter1 = new BAMBarcodeIterator(this.in_bam1);
			final BAMBarcodeIterator iter2 = new BAMBarcodeIterator(this.in_bam2);
			BufferedWriter bw_con = new BufferedWriter(new FileWriter(out_prefix+".txt"));
			
			while(iter1.hasNext()&&iter2.hasNext()) {
				List<SAMRecord[]> bc_records1 = iter1.next();
				List<SAMRecord[]> bc_records2 = iter2.next();

				List<Molecule> mols1 = extractMoleculeFromList(bc_records1);
				List<Molecule> mols2 = extractMoleculeFromList(bc_records2);
				
				int h;
				int n1 = mols1.size(), n2 = mols2.size();
				// this is for homologous molecules
				for(int i=0; i!=n1; i++) {
					for(int j=0; j!=n2; j++) {
						if( (h=homologous(mols1.get(i), mols2.get(j)))!=0 ) {
							bw_con.write("@MOLECULE | "+mols1.get(i).chr_id+":"+mols1.get(i).chr_start+"-"+mols1.get(i).chr_end+"\t"+
									mols2.get(j).chr_id+":"+mols2.get(j).chr_start+"-"+mols2.get(j).chr_end+"\t1.0\t"+
									compare(mols1.get(i), mols2.get(j), h)+"\n");
						}
					}
				}
			}
			
			iter1.close();
			iter2.close();
			bw_con.close();
			
			myLogger.info("#forward molecule, "+fwd_homo);
			myLogger.info("#reverse molecule, "+rev_homo);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private String compare(Molecule mol1, Molecule mol2, int dirc) {
		// TODO Auto-generated method stub
		final List<SAMRecord[]> r1 = mol1.reads_list;
		final List<SAMRecord[]> r2 = mol2.reads_list;
		int n1 = r1.size(), n2 = r2.size();
		if(n1!=n2) throw new RuntimeException("!!!");
	
		int startv1, endv1, startv2, endv2, refposv, altposv, keyv, seql, a;
		String refv, a0, a1, dnaseq1, dnaseq2, prefix;
		SAMRecord[] sams1, sams2;
		SAMRecord sam1, sam2;
		NavigableMap<Integer, Variant> mapv;
		Variant var;
		boolean fwdaln, sim1, sim2;
		final int[] A = new int[2], count = new int[3], alnpos = new int[2];
		final String[] dnaseq = new String[2];
		
		int c1 = 0, c2 = 0, c = 0;
		
		final String chrSeq1 = refSequence1.get(r1.get(0)[0].getReferenceName()).seq_str();
		final String chrSeq2 = refSequence2.get(r2.get(0)[0].getReferenceName()).seq_str();
		
		final String mol = mol1.chr_id+":"+mol1.chr_start+"-"+mol1.chr_end+" | "+
				mol2.chr_id+":"+mol2.chr_start+"-"+mol2.chr_end;
		
		final StringBuilder diff = new StringBuilder();
		diff.append("==========>\n");
		
		for(int i=0; i<n1; i++) {
			sams1 = r1.get(i);
			sams2 = r2.get(dirc==1?i:(n1-1-i));
			
			Arrays.fill(count, 0);
			
			for(int j=0; j<2; j++) {
				sam1 = sams1[j];
				sam2 = sams2[j];
				
				if(!sam1.getReadName().equals(sam2.getReadName()))
					throw new RuntimeException("!!!");
				
				prefix = "@FASTQ | "+sam1.getReadName()+" | "+j+" | "+(dirc==1?"FWD":"REV")+" | "+mol;
				
				diff.append(sam1.getSAMString());
				diff.append(chrSeq1.substring(sam1.getAlignmentStart()-1, sam1.getAlignmentEnd())+"\n");
				diff.append(sam2.getSAMString());
				diff.append(chrSeq2.substring(sam2.getAlignmentStart()-1, sam2.getAlignmentEnd())+"\n");
				
				refv    = sam1.getReferenceName();
				dnaseq1 = sam1.getReadString();
				dnaseq2 = sam2.getReadString();
				seql    = dnaseq1.length();
				
				startv1 = sam1.getAlignmentStart();
				endv1   = sam1.getAlignmentEnd();
				startv2 = sam2.getAlignmentStart();
				endv2   = sam2.getAlignmentEnd();
				
				fwdaln  = sam1.getReadNegativeStrandFlag()==sam2.getReadNegativeStrandFlag();
				mapv   = variants.get(refv).subMap(startv1, true, endv1, true);
				
				for(final Map.Entry<Integer, Variant> mv : mapv.entrySet()) {
					keyv = mv.getKey();
					var  = mv.getValue();
					
					if(var.rev==fwdaln) continue;
					
					altposv = var.altPos;
					
					if(altposv<startv2||altposv>endv2) 
						continue;
					
					refposv = sam1.getReadPositionAtReferencePosition(keyv)-1;
					altposv = sam2.getReadPositionAtReferencePosition(altposv)-1;
					
					diff.append(prefix+" | "+keyv+" "+var.altPos+" | "+var.refA+" "+var.altA+" | vv1,"+refposv+"; vv2,"+altposv+" | ");
					
					if(var.refA.equals(".")) {
						if(altposv<0) {
							++count[2];
							diff.append(2);
						} else {
							if(sim(var.altA, dnaseq2.substring(altposv, Math.min(altposv+var.altA.length(), seql)))) {
								count[1] += var.altA.length();
								diff.append(1);
							} else {
								++count[2];
								diff.append(2);
							}
							diff.append(" | ");

							diff.append(".,.");
							diff.append("; ");
							diff.append(var.altA+","+dnaseq2.substring(altposv, Math.min(altposv+var.altA.length(), seql)));
						}
					} else if(var.altA.equals(".")) {
						if(refposv<0) {
							++count[2];
							diff.append(2);
						} else {
							if(sim(var.refA, dnaseq1.substring(refposv, Math.min(refposv+var.refA.length(), seql)))) { 
								count[0] += var.refA.length();
								diff.append(0);
							} else {
								++count[2];
								diff.append(2);
							}
							diff.append(" | ");

							diff.append(var.refA+","+dnaseq1.substring(refposv, Math.min(refposv+var.refA.length(), seql)));
							diff.append("; ");
							diff.append(".,.");
						}
					} else if(altposv>=0 && refposv>=0) {
						if(var.refA.length()<var.altA.length()) {
							a0 = var.altA;
							a1 = var.refA;
							A[0] = 1;
							A[1] = 0;
							alnpos[0] = altposv;
							alnpos[1] = refposv;
							dnaseq[0] = dnaseq2;
							dnaseq[1] = dnaseq1;
						} else {
							a0 = var.refA;
							a1 = var.altA;
							A[0] = 0;
							A[1] = 1;
							alnpos[0] = refposv;
							alnpos[1] = altposv;
							dnaseq[0] = dnaseq1;
							dnaseq[1] = dnaseq2;
						}
						
						sim1 = sim(a0, dnaseq[0].substring(alnpos[0], Math.min(alnpos[0]+a0.length(), seql)));
						sim2 = sim(a1, dnaseq[1].substring(alnpos[1], Math.min(alnpos[1]+a1.length(), seql)));
						a = 2;
						if(sim1&&!sim2) {
							a = A[0];
							count[a] += a0.length();
						} else if(!sim1&&sim2) {
							a = A[1];
							count[a] += a1.length();
						} else {
							a = 2; 
							++count[2];
						}
						diff.append(a);
						diff.append(" | ");
						
						diff.append(var.refA+","+dnaseq1.substring(refposv, Math.min(refposv+var.refA.length(), seql)));
						diff.append("; ");
						diff.append(var.altA+","+dnaseq2.substring(altposv, Math.min(altposv+var.altA.length(), seql)));
					} else {
						++count[2];
						diff.append("2");
					}
					
					diff.append("\n");
				}
			}
			
			if(count[0]>count[1]) {
				++c1;
			} else if(count[0]<count[1]) {
				++c2;
			} else {
				++c;
			}
			
			diff.append(count[0]);
			diff.append("/");
			diff.append(count[1]);
			diff.append("/");
			diff.append(count[2]);
			// diff.append(",");
			diff.append("\n");
		}
		
		// return n1+"\t"+n2+"\t"+c1+"\t"+c2+"\t"+c+"\t"+diff.toString().replaceAll(",$", "");
		
		return n1+"\t"+n2+"\t"+c1+"\t"+c2+"\t"+c+"\n"+diff.toString();
	}

	private final static double min_simularity = 0.9;
	
	private boolean sim(String target, String query) {
		// TODO Auto-generated method stub
		int n = Math.min(target.length(), query.length());
		if(n==0) return false;
		int eq = 0;
		for(int i=0; i<n; i++) {
			if(target.charAt(i)==query.charAt(i)) 
				++eq;
		}
		return eq>=min_simularity*n;
	}

	private static long fwd_homo = 0;
	private static long rev_homo = 0;
	
	private int homologous(Molecule mol1, Molecule mol2) {
		// TODO Auto-generated method stub
		// return: 0, non-homologous
		//         1, forward
		//         2, reverse
		List<String> r1 = mol1.reads_set;
		List<String> r2 = mol2.reads_set;
		int n = r1.size();
		if(n!=r2.size()) return 0;
		boolean h1 = true, h2 = true;
		for(int i=0; i<n; i++) {
			if(!r1.get(i).equals(r2.get(i))) {
				h1 = false;
				break;
			}
		}
		for(int i=0; i<n; i++) {
			if(!r1.get(i).equals(r2.get(n-1-i))) {
				h2 = false;
				break;
			}
		}
		if(h1) ++fwd_homo;
		if(h2) ++rev_homo;
		if(h1) return  1;
		if(h2) return -1;
		return 0;
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
					middlePoint(records)-chr_pos>max_gap) {
				mol.construct();
				if(mol.mol_sz>=min_mol&&mol.reads_set.size()*abc_per>=mol.mol_sz&&!mol.chr_id.equals("Chr00")) mols.add(mol);
				mol = new Molecule();		
			}
			chr_index = records[0].getReferenceIndex();
			chr_pos = middlePoint(records);
			mol.add(records);
		}
		mol.construct();
		if(mol.mol_sz>=min_mol&&mol.reads_set.size()*abc_per>=mol.mol_sz&&!mol.chr_id.equals("Chr00")) mols.add(mol);
		
		return mols;
	}
}
