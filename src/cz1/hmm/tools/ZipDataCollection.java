package cz1.hmm.tools;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

public class ZipDataCollection extends Executor {

	String vcf_in = null;
	String id = null;
	String out_file = null;
	
	public ZipDataCollection() {}
	
	public ZipDataCollection(String vcf_in, 
			String id, 
			String out_file) {
		this.vcf_in = vcf_in;
		this.id = id;
		this.out_file = out_file;
	}
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--vcf            Input VCF file.\n"
						+ " -s/--id             Unique id of this run (default: input VCF file name prefix).\n"
						+ " -o/--prefix         Prefix for output (default: input VCF file folder).\n\n");
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
			myArgsEngine.add("-i", "--vcf", true);
			myArgsEngine.add("-s", "--id", true);
			myArgsEngine.add("-o", "--prefix", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-i")) {
			vcf_in = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your VCF file.");
		}
		
		if (myArgsEngine.getBoolean("-s")) {
			id = myArgsEngine.getString("-s");
		} else {
			id = new File(vcf_in).getName().replaceAll(".vcf.gz$", "").
					replaceAll(".vcf$", "");
		}
		
		if (myArgsEngine.getBoolean("-o")) {
			out_file = myArgsEngine.getString("-o");
		} else {
			out_file = new File(vcf_in).getParent();
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		String zipFilePath = out_file+"/"+id+".zip";
		try {
			
			ZipOutputStream out = new ZipOutputStream(new BufferedOutputStream(new 
					FileOutputStream(zipFilePath), 65536));
			
			BufferedReader br = Utils.getBufferedReader(vcf_in);
			String line;
			while( (line=br.readLine()).startsWith("##") ) {}
			String[] s = line.split("\\s+");
			List<String> sList = Arrays.asList(s);
			int _idx_start = sList.indexOf("FORMAT")+1,
					_idx_pos = sList.indexOf("POS"),
					_idx_ref = sList.indexOf("REF"),
					_idx_alt = sList.indexOf("ALT"),
					_idx_format = sList.indexOf("FORMAT"),
					_idx_gt = -1,
					_idx_ad = -1,
					_idx_pl = -1;
					
			out.putNextEntry(new ZipEntry("samples"));
			Set<String> allSamples = new HashSet<>();
			String sample;
			for(int i=_idx_start; i<s.length; i++) {
				sample = s[i];
				if(sample.contains(":")) {
					sample = sample.split(":")[0];
					myLogger.warn("sample renamed: "+s[i]+" -> "+sample);
				}
				if(allSamples.contains(sample)) {
					br.close();
					out.close();
					throw new RuntimeException("Duplicate samples "+sample+"!!!");
				}
				allSamples.add(sample);
				out.write((sample+"\n").getBytes());
			}
			List<Contig> contigs = new ArrayList<Contig>();
			line = br.readLine();
			if(line == null) {
				br.close();
				out.close();
				throw new RuntimeException("Empty VCF file!!!");
			}
			s = line.split("\\s+");
			final Set<String> fields = new HashSet<String>();
			String[] fs = s[_idx_format].split(":");
			for(int i=0; i<fs.length; i++) {
				fields.add(fs[i]);
				if(fs[i].equals("GT")) {
					_idx_gt = i;
				} else if(fs[i].equals("AD")) {
					_idx_ad = i;
				} else if(fs[i].equals("PL")) {
					_idx_pl = i;
				}
			}
			while(line != null) {
				int i=1;
				String contig = line.split("\\s+")[0];
				List<String[]> snps = new ArrayList<String[]>();
				snps.add(line.split("\\s+"));
				while( (line=br.readLine())!=null && 
						line.split("\\s+")[0].equals(contig)) {
					i++;
					snps.add(line.split("\\s+"));
				}
				contigs.add(new Contig(contig,i));
				
				out.putNextEntry(new ZipEntry(contig+"/"));
				
				out.putNextEntry(new ZipEntry(contig+"/position"));
				for(int j=0; j<snps.size(); j++)
					out.write((snps.get(j)[_idx_pos]+"\n").
							getBytes());
				
				out.putNextEntry(new ZipEntry(contig+"/allele"));
				for(int j=0; j<snps.size(); j++)
					out.write((snps.get(j)[_idx_ref]+
							"\t"+snps.get(j)[_idx_alt]+
							"\n").
							getBytes());
				
				if(fields.contains("GT")) {
					out.putNextEntry(new ZipEntry(contig+"/GT"));
					for(int j=0; j<snps.size(); j++) {
						for(int k=_idx_start; k<snps.get(j).length; k++)
							out.write((snps.get(j)[k].split(":")[_idx_gt]+"\t")
									.getBytes());
						out.write("\n".getBytes());
					}
				}
				if(fields.contains("AD")) {
					out.putNextEntry(new ZipEntry(contig+"/AD"));
					for(int j=0; j<snps.size(); j++) {
						for(int k=_idx_start; k<snps.get(j).length; k++)
							try {
								out.write((snps.get(j)[k].split(":")[_idx_ad]+"\t")
										.getBytes());
							} catch (ArrayIndexOutOfBoundsException e) {
								out.write(".\t".getBytes());
							}
						out.write("\n".getBytes());
					}
				}
				if(fields.contains("PL")) {
					out.putNextEntry(new ZipEntry(contig+"/PL"));
					for(int j=0; j<snps.size(); j++) {
						for(int k=_idx_start; k<snps.get(j).length; k++)
							try {
								out.write((snps.get(j)[k].split(":")[_idx_pl]+"\t")
										.getBytes());
							} catch (ArrayIndexOutOfBoundsException e) {
								out.write(".\t".getBytes());
							}
						out.write("\n".getBytes());
					}
				}
			}
			br.close();
			
			Collections.sort(contigs);
			Collections.reverse(contigs);
			
			out.putNextEntry(new ZipEntry("contig"));
			for(Contig contig : contigs)
				out.write((contig.id+"\t"+contig.markers+"\n").getBytes());
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private class Contig implements Comparable<Contig> {
		private String id;
		private int markers;
		
		public Contig(String id,
				int markers) {
			this.id = id;
			this.markers = markers;
		}
		
		@Override
		public int compareTo(Contig contig) {
			// TODO Auto-generated method stub
			return this.markers-contig.markers;
		}
	}
}
