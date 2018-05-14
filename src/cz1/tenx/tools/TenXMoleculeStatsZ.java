package cz1.tenx.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
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

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class TenXMoleculeStatsZ extends Executor {
	
	private String in_bam;
	private String in_vcf;
	private String out_prefix;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--input-bam      Input BAM file.\n"
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
			myArgsEngine.add("-i", "--input-bam", true);
			myArgsEngine.add("-x", "--variants", true);
			myArgsEngine.add("-o", "--out-prefix", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-i")) {
			in_bam = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your BAM file.");
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
		private List<SAMRecord[]> reads_list;
		private String chr_id = null;
		private int chr_start = -1;
		private int chr_end = -1;

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
	}
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		
		try {
			final Map<String, TreeMap<Integer, String[]>> variants = new HashMap<String, TreeMap<Integer, String[]>>();
			final BufferedReader br_vcf = Utils.getBufferedReader(in_vcf);
			String line;
			String[] s;
			while( (line=br_vcf.readLine())!=null ) {
				if(line.startsWith("#")) continue;
				s = line.split("\\s+");
				if(!variants.containsKey(s[0]))
					variants.put(s[0], new TreeMap<Integer, String[]>());
				variants.get(s[0]).put(Integer.parseInt(s[1]), new String[]{s[3],s[4]});
			}
			br_vcf.close();
			myLogger.info("Variants loaded from "+in_vcf);
			for(String key : variants.keySet()) 
				myLogger.info(key+"-"+variants.get(key).size());
			
			final BAMBarcodeIterator iter = new BAMBarcodeIterator(this.in_bam);
			BufferedWriter bw_mol = new BufferedWriter(new FileWriter(out_prefix+".mol"));

			int startv, endv, posv, keyv, seql, count, a;
			String refv, a0, a1, dnaseq;
			String[] alleles;
			NavigableMap<Integer, String[]> mapv;
			final Map<Integer, Set<Integer>> molv = new HashMap<Integer, Set<Integer>>(); 
			final StringBuilder os = new StringBuilder();
			while(iter.hasNext()) {
				List<SAMRecord[]> bc_records = iter.next();
				List<Molecule> mols = extractMoleculeFromList(bc_records);
				
				for(final Molecule molecule : mols) {
					final List<SAMRecord[]> samList = molecule.reads_list;
					molv.clear();
					
					for(final SAMRecord[] sams : samList) {
						for(final SAMRecord sam : sams) {
							refv   = sam.getReferenceName();
							startv = sam.getAlignmentStart();
							endv   = sam.getAlignmentEnd();
							dnaseq = sam.getReadString();
							seql   = dnaseq.length();
							mapv   = variants.get(refv).subMap(startv, true, endv, true);
							for(final Map.Entry<Integer, String[]> mv : mapv.entrySet()) {
								keyv = mv.getKey();
								posv = sam.getReadPositionAtReferencePosition(keyv)-1;
								if(posv>0) {
									alleles = mv.getValue();
									if(alleles[0].length()<alleles[1].length()) {
										a0 = alleles[1];
										a1 = alleles[0];
									} else {
										a0 = alleles[0];
										a1 = alleles[1];
									}
									a = -1;
									if(dnaseq.substring(posv, Math.min(posv+a0.length(), seql)).equals(a0)) {
										a = 0;
									} else if(dnaseq.substring(posv, Math.min(posv+a1.length(), seql)).equals(a1)) {
										a = 1;
									}
									if(a!=-1) {
										if(!molv.containsKey(keyv)) molv.put(keyv, new HashSet<Integer>());
										molv.get(keyv).add(a);
									}
								}
							}
						}
					}
					
					os.setLength(0);	
					os.append(molecule.chr_id);
					os.append("\t");
					os.append(molecule.chr_start);
					os.append("\t");
					os.append(molecule.chr_end);
					os.append("\t");
					
					count = 0;
					final List<Integer> keys = new ArrayList<Integer>(molv.keySet());
					Collections.sort(keys);
					for(Integer key : keys) {
						Set<Integer> ale = molv.get(key);
						if(ale.size()>1) continue;
						for(Integer e : ale) {
							os.append(key);
							os.append(",");
							os.append(e);
							os.append(";");
							++count;
						}
					}
					
					if(count>1) {
						bw_mol.write(os.toString().replaceAll(";$", ""));
						bw_mol.write("\n");
					}
				}
			}

			iter.close();
			bw_mol.close();
		} catch (IOException e) {
			e.printStackTrace();
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
					middlePoint(records)-chr_pos>max_gap) {
				mol.construct();
				mols.add(mol);
				mol = new Molecule();		
			}
			chr_index = records[0].getReferenceIndex();
			chr_pos = middlePoint(records);
			mol.add(records);
		}
		mol.construct();
		mols.add(mol);
		
		return mols;
	}
}
