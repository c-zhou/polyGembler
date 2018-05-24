package cz1.tenx.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.NavigableMap;
import java.util.TreeMap;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class GenomeComparisonZZ2 extends Executor {

	private String in_bam1;
	private String in_bam2;
	private String in_vcf;
	private String out_prefix;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i1/--input-bam1    Input BAM file 1.\n"
						+ " -i2/--input-bam2    Input BAM file 2.\n"
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
	
	private final static class Variant {
		private final String refA;
		private final String altA;
		private final boolean rev;
		private final int altPos;
		private final String keyStr;
		
		public Variant(final String refA, final String altA, final int altPos, final boolean rev, final String keyStr) {
			this.refA   = refA;
			this.altA   = altA;
			this.altPos = altPos;
			this.rev    = rev;
			this.keyStr = keyStr;
		}
	}

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
		
		public SAMRecord[] next() {
			
			if(!this.hasNext()) throw new RuntimeException("!!!");
			
			final SAMRecord[] records = new SAMRecord[2];
			String sn = samRecord.getReadName();

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
					return records;
				}
			}

			return null;
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

	
	private final static Map<String, TreeMap<Integer, Variant>> variants = new HashMap<String, TreeMap<Integer, Variant>>();
	private final static Map<String, int[]> depthStats = new HashMap<String, int[]>();
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		
		try {
			final BufferedReader br_vcf = Utils.getBufferedReader(in_vcf);
			String line;
			String[] s;
			long pcount = 0;
			while( (line=br_vcf.readLine())!=null ) {
				if(line.startsWith("#")) continue;
				s = line.split("\\s+");
				if(!variants.containsKey(s[0]))
					variants.put(s[0], new TreeMap<Integer, Variant>());
				variants.get(s[0]).put(Integer.parseInt(s[1]), new Variant(s[3], s[4], Integer.parseInt(s[2]), s[5].equals("1")?false:true, line));
				depthStats.put(line, new int[2]);
				++pcount;
				if(pcount%1000000==0) myLogger.info(pcount+" variants loaded into memory.");
			}
			br_vcf.close();
			myLogger.info("Variants loaded from "+in_vcf);
			for(String key : variants.keySet()) 
				myLogger.info(key+": "+variants.get(key).size());
			
			final BAMBarcodeIterator iter1 = new BAMBarcodeIterator(this.in_bam1);
			final BAMBarcodeIterator iter2 = new BAMBarcodeIterator(this.in_bam2);
			
			while(iter1.hasNext()&&iter2.hasNext()) {
				SAMRecord[] records1 = iter1.next();
				SAMRecord[] records2 = iter2.next();	
				if(records1==null||records2==null||
						!records1[0].getReferenceName().equals(records2[0].getReferenceName())||
						records1[0].getReferenceName().equals("Chr00")) 
					continue;
				if(!records1[0].getReadName().equals(records2[0].getReadName()))
					throw new RuntimeException("!!!");
				compare(records1, records2);
			}
			
			iter1.close();
			iter2.close();
			BufferedWriter bw_con = new BufferedWriter(new FileWriter(out_prefix+".dp1"));
			int[] dp;
			for(String keyStr : depthStats.keySet()) {
				dp = depthStats.get(keyStr);
				bw_con.write(keyStr+"\t"+dp[0]+"\t"+dp[1]+"\t");
			}
			bw_con.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void compare(SAMRecord[] sams1, SAMRecord[] sams2) {
		// TODO Auto-generated method stub
		
		int startv1, endv1, startv2, endv2, refposv, altposv, keyv, seql;
		String refv, a0, a1, dnaseq1, dnaseq2;
		SAMRecord sam1, sam2;
		NavigableMap<Integer, Variant> mapv;
		Variant var;
		boolean fwdaln;
		boolean sim1, sim2;
		final int[] A = new int[2], alnpos = new int[2];
		final String[] dnaseq = new String[2];
		
		for(int j=0; j<2; j++) {
			sam1 = sams1[j];
			sam2 = sams2[j];

			if(!sam1.getReadName().equals(sam2.getReadName()))
				throw new RuntimeException("!!!");
			
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

				if(var.refA.equals(".")) {
					if(altposv>=0&&sim(var.altA, dnaseq2.substring(altposv, Math.min(altposv+var.altA.length(), seql)))) 
						++depthStats.get(var.keyStr)[1];
				} else if(var.altA.equals(".")) {
					if(refposv>=0&&sim(var.refA, dnaseq1.substring(refposv, Math.min(refposv+var.refA.length(), seql))))  
						++depthStats.get(var.keyStr)[0];
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
					if(sim1&&!sim2) {
						++depthStats.get(var.keyStr)[A[0]];
					} else if(!sim1&&sim2) {
						++depthStats.get(var.keyStr)[A[1]];
					}
				}
			}
		}
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
}
