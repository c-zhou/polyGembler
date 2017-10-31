package cz1.ngs.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.DiscreteDomain;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

import cz1.ngs.model.Blast6Record;
import cz1.ngs.model.Sequence;
import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

public class Redundas extends Executor {

	private String query_file = null;
	private String blast_out = null;
	private double min_ident = 0.95;
	private int min_match = 100;
	private int max_overhang = Integer.MAX_VALUE;
	private double min_frac = 0.95;
	private String out_prefix = null;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -q/--query              Query sequences to remove redundancy.\n"
						+ " -b/--blast              Self-to-self BLAST output (-outfmt 6) of query sequences.\n"
						+ " -i/--min-identity       Minimum identity between the query and subject sequences \n"
						+ "                         to mark a redundancy (default 0.95).\n"
						+ " -m/--min-match          Minimum match length between the query and subject sequences \n"
						+ "                         for an alignment to consider (default 100).\n"
						+ " -x/--max-overhang       Maximum distance to the end of the subject sequence for an \n"
						+ "                         alignment record to consider (default no restriction).\n"
						+ " -f/--min-fraction       Minimum alignment fraction of the query sequence for an \n"
						+ "                         alignment record to consider (default 0.95).\n"
						+ " -o/--out-prefix         Prefix of the output FASTQ file.\n"
						+ "\n");	
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-q", "--query", true);
			myArgsEngine.add("-b", "--blast", true);
			myArgsEngine.add("-i", "--min-identity", true);
			myArgsEngine.add("-m", "--min-match", true);
			myArgsEngine.add("-x", "--max-overhang", true);
			myArgsEngine.add("-f", "--min-fraction", true);
			myArgsEngine.add("-o", "--out-prefix", true);
			myArgsEngine.parse(args);
		}
		
		if (myArgsEngine.getBoolean("-q")) {
			this.query_file = myArgsEngine.getString("-q");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the query file.");
		}
		
		if (myArgsEngine.getBoolean("-b")) {
			this.blast_out = myArgsEngine.getString("-b");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the BLAST file.");
		}
		
		if (myArgsEngine.getBoolean("-i")) {
			this.min_ident = Double.parseDouble(myArgsEngine.getString("-i"));
		}
		
		if (myArgsEngine.getBoolean("-m")) {
			this.min_match = Integer.parseInt(myArgsEngine.getString("-m"));
		}
		
		if (myArgsEngine.getBoolean("-x")) {
			this.max_overhang = Integer.parseInt(myArgsEngine.getString("-x"));
		}
		
		if (myArgsEngine.getBoolean("-f")) {
			this.min_frac = Double.parseDouble(myArgsEngine.getString("-f"));
		}
		
		if (myArgsEngine.getBoolean("-o")) {
			this.out_prefix = myArgsEngine.getString("-o");
			if(new File(out_prefix+".fa").exists() || 
					new File(out_prefix+"_redundas.fa").exists()) {
				throw new RuntimeException("Output files exist. Please specify a different name.");
			}
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the prefix of output files.");
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		Map<String, Sequence> sequence_map = Sequence.parseFastaFileAsMap(query_file);
		
		String line;
		String[] s;
		Set<String> seq_rm = new HashSet<String>();
		try {
			BufferedReader br_blast6 = new BufferedReader(new FileReader(blast_out));
			Set<Blast6Record> buffer_b6 = new HashSet<Blast6Record>();
			line = br_blast6.readLine();
			String qseqid;
			while( line!=null ) {
				s = line.split("\\s+");
				qseqid = s[0];
				
				buffer_b6.clear();
				
				buffer_b6.add(Blast6Record.blast6Record(line));
				while( (line=br_blast6.readLine())!=null && 
						line.startsWith(qseqid) ) 
					buffer_b6.add(Blast6Record.blast6Record(line));		
				
				int sz = sequence_map.get(qseqid).seq_ln();
			
				RangeSet<Integer> range_covered = TreeRangeSet.create();
				for(Blast6Record record : buffer_b6) {
					if( !seq_rm.contains(record.sseqid()) &&
							!record.qseqid().equals(record.sseqid()) &&
							record.pident()>=min_ident && 
							record.length()>=min_match &&
							(record.qstart()<=max_overhang ||
							record.qend()>sz-max_overhang) )
					range_covered.add(Range.closedOpen(record.qstart(), record.qend()).canonical(DiscreteDomain.integers()));
				}

				int unique_cvg = 0;
				for(Range<Integer> range : range_covered.asRanges()) 
					unique_cvg += range.upperEndpoint()-range.lowerEndpoint();

				if( (double)unique_cvg/sz>=min_frac ) {
					seq_rm.add(qseqid);
					myLogger.info("Redundant sequence "+qseqid+"\t"+sz+"\t"+unique_cvg+"\t"+(double)unique_cvg/sz);
				}
			}
			br_blast6.close();
		
			BufferedWriter bw_unique = Utils.getBufferedWriter(this.out_prefix+".fa");
			BufferedWriter bw_redundas = Utils.getBufferedWriter(this.out_prefix+"_redundas.fa");
			for(Map.Entry<String, Sequence> entry : sequence_map.entrySet()) {
				if(seq_rm.contains(entry.getKey())) 
					bw_redundas.write(entry.getValue().formatOutput());
				else bw_unique.write(entry.getValue().formatOutput());
			}
			bw_unique.close();
			bw_redundas.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
