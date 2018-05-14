package cz1.tenx.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.inference.TestUtils;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

public class TenXVcftools extends Executor {

	private static enum Task {filter, zzz}
	private Task task = Task.zzz;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		switch(this.task) {
		case zzz:
			myLogger.info(
					"\n\nChoose a task from listed:\n"
							+ " filter               Filter VCF file.\n\n");	
			break;
		case filter:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -i/--in-vcf         Input VCF file.\n"
							+ " -p/--ploidy         Ploidy. \n"
							+ " -o/--out-vcf        Output VCF file.\n\n");	
			break;
		default:
			throw new RuntimeException("Undefined task!!!");
		}
	}
	
	private String vcf_in;
	private String vcf_out;
	private int ploidy = 6;
	private double[][] config = null;
	
	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}

		switch(args[0].toUpperCase()) {
		case "FILTER":
			this.task = Task.filter;
			break;
		default:
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");	
		}

		String[] args2 = new String[args.length-1];
		System.arraycopy(args, 1, args2, 0, args2.length);

		switch(this.task) {
		case filter:
			if (myArgsEngine == null) {
				myArgsEngine = new ArgsEngine();
				myArgsEngine.add("-i", "--in-vcf", true);
				myArgsEngine.add("-p", "--ploidy", true);
				myArgsEngine.add("-o", "--out-vcf", true);
				myArgsEngine.parse(args2);
			}
			if (myArgsEngine.getBoolean("-i")) {
				this.vcf_in = myArgsEngine.getString("-i");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the input VCF file.");
			}
			if (myArgsEngine.getBoolean("-p")) {
				this.ploidy = Integer.parseInt(myArgsEngine.getString("-p"));
			}
			configure();
			if (myArgsEngine.getBoolean("-o")) {
				this.vcf_out = myArgsEngine.getString("-o");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the output VCF file.");
			}
			
			break;
		default:
			throw new RuntimeException("!!!");	
		}
	}

	private void configure() {
		// TODO Auto-generated method stub
		this.config = new double[this.ploidy-1][2];
		for(int i=1; i<this.ploidy; i++) {
			config[i-1][0] = (double)i/this.ploidy;
			config[i-1][1] = 1.0-config[i-1][0];
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		switch(this.task) {
		case zzz:
			myLogger.info("Task list is empty!!!");
			break;
		case filter:
			this.filter();
			break;
		default:
			throw new RuntimeException("!!!");
		}
		return;
	}
	
	private void filter() {
		// TODO Auto-generated method stub
		try {
			BufferedReader br_vcf = Utils.getBufferedReader(vcf_in);
			BufferedWriter bw_vcf = Utils.getBufferedWriter(vcf_out);
			String line;
			String[] s;
			String chrom, pos, ref, alt, qual;
			
			final long[] allele_depth = new long[2];
			final double[] chisq_p = new double[this.ploidy-1];
			double pval;
			long ada;
			
			while( (line=br_vcf.readLine())!=null ) {
				if(line.startsWith("#")) {
					bw_vcf.write(line+"\n");
					continue;
				}
				s = line.split("\\s+");
				chrom = s[0];
				pos   = s[1];
				ref   = s[3];
				alt   = s[4];
				qual  = s[5];
				if(alt.contains(",")||Double.parseDouble(qual)<20) continue;
				s = s[9].split(":");
				allele_depth[0] = Long.parseLong(s[4]);
				allele_depth[1] = Long.parseLong(s[6]);
				
				ada = allele_depth[0]+allele_depth[1];
				if(allele_depth[0]<5||allele_depth[1]<5||ada<30||ada>500) continue;

				for(int z=0; z<this.ploidy-1; z++) 
					chisq_p[z] = TestUtils.chiSquareTest(config[z], allele_depth);
				
				pval = StatUtils.max(chisq_p);
				if(pval>=0.05)
					bw_vcf.write(chrom+"\t"+pos+"\t.\t"+ref+"\t"+alt+"\t"+qual+"\t.\tPV="+String.format("%.3f", pval)+"\tAD\t"+s[4]+","+s[6]+"\n");
			}
			bw_vcf.close();
			br_vcf.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
	}
}
