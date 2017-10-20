package cz1.tenx.tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.renjin.repackaged.guava.collect.Sets;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;

public class GenomeComparison extends Executor {

	private String in_bam1 = null;
	private String in_bam2 = null;
	private String out_f = null;

	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i1/--input-bam-1   Input BAM file 1.\n"
						+ " -i2/--input-bam-2   Input BAM file 2.\n"
						+ " -t/--threads        Threads (default is 1).\n"
						+ " -o/--prefix         Output files prefix.\n\n");	
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
			myArgsEngine.add("-i1", "--input-bam-1", true);
			myArgsEngine.add("-i2", "--input-bam-2", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-o", "--prefix", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-i1")) {
			in_bam1 = myArgsEngine.getString("-i1");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your FASTQ files.");
		}

		if (myArgsEngine.getBoolean("-i2")) {
			in_bam2 = myArgsEngine.getString("-i2");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your FASTQ files.");
		}

		if (myArgsEngine.getBoolean("-t")) {
			THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}

		if (myArgsEngine.getBoolean("-o")) {
			out_f = myArgsEngine.getString("-o");
		}	
	}

	private class Molecule {
		private Set<String> reads_set;
		private List<SAMRecord> reads_list;
		private String chr_id = null;
		private int chr_start = -1;
		private int chr_end = -1;

		public Molecule() {
			this.reads_set = new HashSet<String>();
			this.reads_list = new ArrayList<SAMRecord>();
		}

		public void add(SAMRecord record) {
			this.reads_list.add(record);
			this.reads_set.add(record.getReadName().substring(25));
		}

		public void construct() {
			// TODO Auto-generated method stub
			Collections.sort(reads_list, new Comparator<SAMRecord>() {
				@Override
				public int compare(SAMRecord r, SAMRecord r2) {
					// TODO Auto-generated method stub
					return r.getAlignmentStart()-r2.getAlignmentStart();
				}
			});
			this.chr_id = this.reads_list.get(0).getReferenceName();
			this.chr_start = this.reads_list.get(0).getAlignmentStart();
			//TODO: chr_end is not right, need to iterate through the whole SAMRecord list
			this.chr_end = this.reads_list.get(this.reads_list.size()-1).getAlignmentEnd();
		}
	}

	BufferedWriter bw_con, bw_1, bw_2;

	@Override
	public void run() {
		// TODO Auto-generated method stub
		final SAMFileReader inputSam1 = new SAMFileReader(new File(in_bam1));
		inputSam1.setValidationStringency(ValidationStringency.SILENT);
		SAMRecordIterator iter1 = inputSam1.iterator();

		final SAMFileReader inputSam2 = new SAMFileReader(new File(in_bam2));
		inputSam2.setValidationStringency(ValidationStringency.SILENT);
		SAMRecordIterator iter2 = inputSam2.iterator();

		this.initial_thread_pool();

		try {
			bw_con = new BufferedWriter(new FileWriter(out_f+".txt"));
			bw_1 = new BufferedWriter(new FileWriter(out_f+".1txt"));
			bw_2 = new BufferedWriter(new FileWriter(out_f+".2txt"));

			SAMRecord record1 = iter1.next(), 
					record2 = iter2.next();
			List<SAMRecord> barcoded_records1 = new ArrayList<SAMRecord>(),
					barcoded_records2 = new ArrayList<SAMRecord>();
			barcoded_records1.add(record1);
			barcoded_records2.add(record2);
			String barcode1 = record1.getReadName().split(":")[0],
					barcode2 = record2.getReadName().split(":")[0];

			jobs:
			while(true) {
				if(barcode1.compareTo(barcode2)<0) {
					//output barcode 1 up to barcode 1
					while(iter1.hasNext()) {
						while( (record1=iter1.next()).getReadName().split(":")[0].equals(barcode1) ) { 
							barcoded_records1.add(record1);
							if(!iter1.hasNext()) break;
						}

						//TODO barcoded records 1 output
						Molecule[] mols = extractMoleculeFromList(barcoded_records1);
						for(int i=0; i!=mols.length; i++) 
							bw_1.write(mols[i].chr_id+":"+mols[i].chr_start+"-"+mols[i].chr_end+"\n");

						barcoded_records1.clear();

						if(iter1.hasNext()) {
							barcoded_records1.add(record1);
							barcode1 = record1.getReadName().split(":")[0];
						} else break jobs;

						if(barcode1.compareTo(barcode2)>=0) break;
					}
					if(!barcode1.equals(barcode2)) continue;

				} else if(barcode1.compareTo(barcode2)>0) {
					//output barcode 2 up to barcode 2
					while(iter2.hasNext()) {

						while( (record2=iter2.next()).getReadName().split(":")[0].equals(barcode2) ) { 
							barcoded_records2.add(record2);
							if(!iter2.hasNext()) break;	
						}

						//TODO barcoded records 2 output
						Molecule[] mols = extractMoleculeFromList(barcoded_records2);
						for(int i=0; i!=mols.length; i++) 
							bw_2.write(mols[i].chr_id+":"+mols[i].chr_start+"-"+mols[i].chr_end+"\n");

						barcoded_records2.clear();
						if(iter2.hasNext()) {
							barcoded_records2.add(record2);
							barcode2 = record2.getReadName().split(":")[0];
						} else break jobs;
						if(barcode2.compareTo(barcode1)>=0) break;
					}
					if(!barcode1.equals(barcode2)) continue;
				} 
				
				while( (record1=iter1.next()).getReadName().split(":")[0].equals(barcode1) ) { 
					barcoded_records1.add(record1);
					if(!iter1.hasNext()) break;
				}

				while( (record2=iter2.next()).getReadName().split(":")[0].equals(barcode2) ) {
					barcoded_records2.add(record2);
					if(!iter2.hasNext()) break;
				}

				executor.submit(new Runnable() {
					List<SAMRecord> list1;
					List<SAMRecord> list2;

					@Override
					public void run() {
						// TODO Auto-generated method stub
						try {
							Molecule[] mols1 = extractMoleculeFromList(list1);
							Molecule[] mols2 = extractMoleculeFromList(list2);
							int n1 = mols1.length, n2 = mols2.length;
							double insec;
							for(int i=0; i!=n1; i++) {
								for(int j=0; j!=n2; j++) {
									if( (insec=intersect(mols1[i], mols2[j]))>=0.5 )
										bw_con.write(mols1[i].chr_id+":"+mols1[i].chr_start+"-"+mols1[i].chr_end+"\t"+
												mols2[j].chr_id+":"+mols2[j].chr_start+"-"+mols2[j].chr_end+"\t"+insec+"\n");
								}
							}
						} catch (Exception e) {
							Thread t = Thread.currentThread();
							t.getUncaughtExceptionHandler().uncaughtException(t, e);
							e.printStackTrace();
							executor.shutdown();
							System.exit(1);
						}
					}

					public Runnable init(List<SAMRecord> list1, List<SAMRecord> list2) {
						// TODO Auto-generated method stub
						this.list1 = list1;
						this.list2 = list2;
						return this;
					}

				}.init(new ArrayList<SAMRecord>(barcoded_records1), 
						new ArrayList<SAMRecord>(barcoded_records2)));

				if(!iter1.hasNext() || !iter2.hasNext()) break jobs;
				
				barcoded_records1.clear();
				barcoded_records2.clear();
				barcoded_records1.add(record1);
				barcoded_records2.add(record2);
				barcode1 = record1.getReadName().split(":")[0];
				barcode2 = record2.getReadName().split(":")[0];			
			}

			
			if(!iter1.hasNext()) {
				//TODO output remaining iter2
				while(iter2.hasNext()) {
					while( (record2=iter2.next()).getReadName().split(":")[0].equals(barcode2) ) { 
						barcoded_records2.add(record2);
						if(!iter2.hasNext()) break;
					}
					
					Molecule[] mols = extractMoleculeFromList(barcoded_records2);
					for(int i=0; i!=mols.length; i++) 
						bw_2.write(mols[i].chr_id+":"+mols[i].chr_start+"-"+mols[i].chr_end+"\n");
				
					barcoded_records2.clear();
					if(iter2.hasNext()) {
						barcoded_records2.add(record2);
						barcode2 = record2.getReadName().split(":")[0];
					}
				}
			}

			if(!iter2.hasNext()) {
				//TODO output remaining iter1
				while(iter1.hasNext()) {
					while( (record1=iter1.next()).getReadName().split(":")[0].equals(barcode1) ) { 
						barcoded_records1.add(record1);
						if(!iter1.hasNext()) break;
					}
					
					Molecule[] mols = extractMoleculeFromList(barcoded_records1);
					for(int i=0; i!=mols.length; i++) 
						bw_1.write(mols[i].chr_id+":"+mols[i].chr_start+"-"+mols[i].chr_end+"\n");
				
					barcoded_records1.clear();
					if(iter1.hasNext()) {
						barcoded_records1.add(record1);
						barcode1 = record1.getReadName().split(":")[0];
					}
				}
			}
			
			this.waitFor();
			inputSam1.close();
			inputSam2.close();
			bw_con.close();
			bw_1.close();
			bw_2.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private double intersect(Molecule mol1, Molecule mol2) {
		// TODO Auto-generated method stub
		int a = Sets.union(mol1.reads_set, mol2.reads_set).size();
		int b = Sets.intersection(mol1.reads_set, mol2.reads_set).size();
		return b/(double) a;
	}

	private Molecule[] extractMoleculeFromList(List<SAMRecord> list) {
		// TODO Auto-generated method stub
		int n = Integer.parseInt(list.get(list.size()-1).getReadName().split(":")[1])+1;
		Molecule[] mols = new Molecule[n];
		for(int i=0; i!=n; i++) 
			mols[i] = new Molecule();
		int k;
		for(SAMRecord record : list) {
			k = Integer.parseInt(record.getReadName().substring(19, 24));
			mols[k].add(record);
		}
		for(int i=0; i!=n; i++) {
			mols[i].construct();
		}
		return mols;
	}

	public static void main(String[] args) {
		GenomeComparison gc = new GenomeComparison();
		gc.setParameters(args);
		gc.run();
	}
}
