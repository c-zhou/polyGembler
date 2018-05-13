package cz1.tenx.tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class GenomeComparisonZ extends Executor {
	
	private String in_bam1;
	private String in_bam2;
	private String out_prefix;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i1/--input-bam1    Input BAM file 1.\n"
						+ " -i2/--input-bam2    Input BAM file 2.\n"
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
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		try {
			final BAMBarcodeIterator iter1 = new BAMBarcodeIterator(this.in_bam1);
			final BAMBarcodeIterator iter2 = new BAMBarcodeIterator(this.in_bam2);

			BufferedWriter bw_con = new BufferedWriter(new FileWriter(out_prefix+".txt"));
			BufferedWriter bw_1 = new BufferedWriter(new FileWriter(out_prefix+".1txt"));
			BufferedWriter bw_2 = new BufferedWriter(new FileWriter(out_prefix+".2txt"));

			while(iter1.hasNext()&&iter2.hasNext()) {
				List<SAMRecord[]> bc_records1 = iter1.next();
				List<SAMRecord[]> bc_records2 = iter2.next();

				List<Molecule> mols1 = extractMoleculeFromList(bc_records1);
				List<Molecule> mols2 = extractMoleculeFromList(bc_records2);

				int n1 = mols1.size(), n2 = mols2.size();
				Set<Integer> m1 = new HashSet<Integer>();
				Set<Integer> m2 = new HashSet<Integer>();
				for(int i=0; i!=n1; i++) {
					for(int j=0; j!=n2; j++) {
						if( homologous(mols1.get(i), mols2.get(j)) ) {
							bw_con.write(mols1.get(i).chr_id+":"+mols1.get(i).chr_start+"-"+mols1.get(i).chr_end+"\t"+
									mols2.get(j).chr_id+":"+mols2.get(j).chr_start+"-"+mols2.get(j).chr_end+"\t1.0\t"+
									compare(mols1.get(i), mols2.get(j))+"\n");
							m1.add(i);
							m2.add(j);
						}
					}
				}
				for(int i=0; i!=n1; i++) {
					if(!m1.contains(i)) 
						bw_1.write(mols1.get(i).chr_id+":"+mols1.get(i).chr_start+"-"+mols1.get(i).chr_end+"\n");
				}
				for(int i=0; i!=n2; i++) {
					if(!m2.contains(i)) 
						bw_2.write(mols2.get(i).chr_id+":"+mols2.get(i).chr_start+"-"+mols2.get(i).chr_end+"\n");
				}
			}

			iter1.close();
			iter2.close();

			bw_1.close();
			bw_2.close();
			bw_con.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
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

	private String compare(Molecule mol1, Molecule mol2) {
		// TODO Auto-generated method stub
		
		int a12 = 0, a21 = 0, aeq = 0;
		float as1, as2;
		SAMRecord[] s1, s2;
		
		final List<SAMRecord[]> r1 = mol1.reads_list;
		final List<SAMRecord[]> r2 = mol2.reads_list;
		int n1 = r1.size(), n2 = r2.size();
		if(n1!=n2) throw new RuntimeException("!!!");
		
		for(int i=0; i<n1; i++) {
			s1 = r1.get(i);
			as1 = s1[0].getFloatAttribute("AS")+s1[1].getFloatAttribute("AS");
			s2 = r2.get(i);
			as2 = s2[0].getFloatAttribute("AS")+s2[1].getFloatAttribute("AS");
			if(as1>as2) {
				++a12;
			} else if(as1<as2) {
				++a21;
			} else {
				++aeq;
			}
		}
		StringBuilder diff = new StringBuilder();
		for(int i=0; i<n1-1; i++) {
			s1 = r1.get(i);
			diff.append(s1[0].getFloatAttribute("AS")+s1[1].getFloatAttribute("AS"));
			diff.append("/");
			s2 = r2.get(i);
			diff.append(s2[0].getFloatAttribute("AS")+s2[1].getFloatAttribute("AS"));
			diff.append(",");
		}
		s1 = r1.get(n1-1);
		diff.append(s1[0].getFloatAttribute("AS")+s1[1].getFloatAttribute("AS"));
		diff.append("/");
		s2 = r2.get(n2-1);
		diff.append(s2[0].getFloatAttribute("AS")+s2[1].getFloatAttribute("AS"));
	
		return n1+"\t"+n2+"\t"+a12+"\t"+a21+"\t"+aeq+"\t"+diff.toString();
	}
	
	private double middlePoint(SAMRecord[] record) {
		// TODO Auto-generated method stub
		return (double)(Math.min(record[0].getAlignmentStart(), record[1].getAlignmentStart())+
				Math.max(record[0].getAlignmentEnd(), record[1].getAlignmentEnd()))/2;
	}
	
	private List<Molecule> extractMoleculeFromList(List<SAMRecord[]> list) {
		// TODO Auto-generated method stub
		Collections.sort(list, new Comparator<SAMRecord[]>() {
			@Override
			public int compare(SAMRecord[] record0, SAMRecord[] record1) {
				// TODO Auto-generated method stub
				int f = record0[0].getReferenceIndex().intValue()-
						record1[0].getReferenceIndex().intValue();
				return f==0 ? Double.compare(middlePoint(record0), middlePoint(record1)) : f;
			}
		});
		
		List<Molecule> mols = new ArrayList<Molecule>();
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
		if(mol.mol_sz>=min_mol&&mol.reads_set.size()*abc_per>=mol.mol_sz) mols.add(mol);
		
		return mols;
	}
}
