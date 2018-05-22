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
		private final int altPos;
		
		public Variant(final String refA, final String altA, final int altPos) {
			this.refA   = refA;
			this.altA   = altA;
			this.altPos = altPos;
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
				variants.get(s[0]).put(Integer.parseInt(s[1]), new Variant(s[3], s[4], Integer.parseInt(s[2])));
				++pcount;
				if(pcount%1000000==0) myLogger.info(pcount+" variants loaded into memory.");
			}
			br_vcf.close();
			myLogger.info("Variants loaded from "+in_vcf);
			for(String key : variants.keySet()) 
				myLogger.info(key+"-"+variants.get(key).size());
			
			final BAMBarcodeIterator iter1 = new BAMBarcodeIterator(this.in_bam1);
			final BAMBarcodeIterator iter2 = new BAMBarcodeIterator(this.in_bam2);
			BufferedWriter bw_con = new BufferedWriter(new FileWriter(out_prefix+".txt"));
			
			while(iter1.hasNext()&&iter2.hasNext()) {
				List<SAMRecord[]> bc_records1 = iter1.next();
				List<SAMRecord[]> bc_records2 = iter2.next();

				List<Molecule> mols1 = extractMoleculeFromList(bc_records1);
				List<Molecule> mols2 = extractMoleculeFromList(bc_records2);
				
				int n1 = mols1.size(), n2 = mols2.size();
				// this is for homologous molecules
				for(int i=0; i!=n1; i++) {
					for(int j=0; j!=n2; j++) {
						if( homologous(mols1.get(i), mols2.get(j)) ) {
							bw_con.write(mols1.get(i).chr_id+":"+mols1.get(i).chr_start+"-"+mols1.get(i).chr_end+"\t"+
									mols2.get(j).chr_id+":"+mols2.get(j).chr_start+"-"+mols2.get(j).chr_end+"\t1.0\t"+
									compare(mols1.get(i), mols2.get(j))+"\n");
						}
					}
				}
			}
			
			iter1.close();
			iter2.close();
			bw_con.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private String compare(Molecule mol1, Molecule mol2) {
		// TODO Auto-generated method stub
		final List<SAMRecord[]> r1 = mol1.reads_list;
		final List<SAMRecord[]> r2 = mol2.reads_list;
		int n1 = r1.size(), n2 = r2.size();
		if(n1!=n2) throw new RuntimeException("!!!");
	
		int startv1, endv1, startv2, endv2, refposv, altposv, keyv, seql, a;
		String refv, a0, a1, dnaseq;
		SAMRecord[] sams1, sams2;
		SAMRecord sam1, sam2;
		NavigableMap<Integer, Variant> mapv;
		Variant var;
		final int[] A = new int[2], count = new int[3], alnpos = new int[2];
		
		int c1 = 0, c2 = 0, c = 0;
		
		final String chrSeq1 = refSequence1.get(r1.get(0)[0].getReferenceName()).seq_str();
		final String chrSeq2 = refSequence2.get(r2.get(0)[0].getReferenceName()).seq_str();
		
		final StringBuilder diff = new StringBuilder();
		diff.append("==========>\n");
		
		for(int i=0; i<n1; i++) {
			sams1 = r1.get(i);
			sams2 = r2.get(i);
			
			Arrays.fill(count, 0);
			
			for(int j=0; j<2; j++) {
				sam1 = sams1[j];
				sam2 = sams2[j];
				
				diff.append(sam1.getSAMString());
				diff.append(chrSeq1.substring(sam1.getAlignmentStart()-1, sam1.getAlignmentEnd()));
				diff.append(sam2.getSAMString());
				diff.append(chrSeq2.substring(sam2.getAlignmentStart()-1, sam2.getAlignmentEnd()));
				
				refv   = sam1.getReferenceName();
				dnaseq = sam1.getReadString();
				seql   = dnaseq.length();
				
				startv1 = sam1.getAlignmentStart();
				endv1   = sam1.getAlignmentEnd();
				startv2 = sam2.getAlignmentStart();
				endv2   = sam2.getAlignmentEnd();
				
				mapv   = variants.get(refv).subMap(startv1, true, endv1, true);
				
				for(final Map.Entry<Integer, Variant> mv : mapv.entrySet()) {
					keyv = mv.getKey();
					var     = mv.getValue();
					altposv = var.altPos;
					
					if(altposv<startv2||altposv>endv2) 
						continue;
					
					refposv = sam1.getReadPositionAtReferencePosition(keyv)-1;
					altposv = sam2.getReadPositionAtReferencePosition(altposv)-1;
					
					diff.append(keyv+" "+var.refA+" "+var.altA+" "+var.altPos+"; vv1,"+refposv+"; vv2,"+altposv+" | ");
					
					if(refposv<0) {
						if(var.refA.equals(".")) ++count[0];
						
						diff.append("0");
						
					} else if(altposv<0) {
						if(var.altA.equals(".")) ++count[1];
						
						diff.append("1");
						
					} else if(!var.refA.equals(".")&&!var.altA.equals(".")) {
						if(var.refA.length()<var.altA.length()) {
							a0 = var.altA;
							a1 = var.refA;
							A[0] = 1;
							A[1] = 0;
							alnpos[0] = altposv;
							alnpos[1] = refposv;
						} else {
							a0 = var.refA;
							a1 = var.altA;
							A[0] = 0;
							A[1] = 1;
							alnpos[0] = refposv;
							alnpos[1] = altposv;
						}
						
						a = 2;
						if(dnaseq.substring(alnpos[0], Math.min(alnpos[0]+a0.length(), seql)).equals(a0)) {
							a = A[0];
						} else if(dnaseq.substring(alnpos[1], Math.min(alnpos[1]+a1.length(), seql)).equals(a1)) {
							a = A[1];
						}
						++count[a];
						
						diff.append(a);
						diff.append(" | ");
						
						diff.append(a0+","+dnaseq.substring(alnpos[0], Math.min(alnpos[0]+a0.length(), seql)));
						diff.append("; ");
						diff.append(a1+","+dnaseq.substring(alnpos[1], Math.min(alnpos[1]+a1.length(), seql)));
						
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

	private boolean homologous(Molecule mol1, Molecule mol2) {
		// TODO Auto-generated method stub
		List<String> r1 = mol1.reads_set;
		List<String> r2 = mol2.reads_set;
		if(r1.size()!=r2.size()) return false;
		for(int i=0; i<r1.size(); i++)
			if(!r1.get(i).equals(r2.get(i)))
				return false;
		return true;
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
				if(mol.mol_sz>=min_mol&&mol.reads_set.size()*abc_per>=mol.mol_sz) mols.add(mol);
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
