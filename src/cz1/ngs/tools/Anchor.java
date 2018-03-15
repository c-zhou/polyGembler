package cz1.ngs.tools;

import java.io.File;

import cz1.util.ArgsEngine;
import cz1.util.Executor;

public class Anchor extends Executor {
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -s/--subject            Subject/reference sequences file in FASTA format.\n"
						+ " -q/--query              Query sequences to anchor in FASTA format.\n"
						+ " -a/--align              Alignment file of query sequences to the subject sequences. \n"
						+ "                         IMPORTANT: The alignment file need to be grouped by the query sequences,\n"
						+ "                         i.e., the program assume all the alignment records for a query sequence \n"
						+ "                         can be continuously read from the alignment file.\n"
						+ " -blast/--blast          Alignment file is BLAST output format 6 (default).\n"
						+ " -bam/--bam              Alignment file is BAM format.\n"
						+ " -mummer/--mummer        Alignment file is MUMmer output format (not implemented yet).\n"
						+ " -g/--graph              Assembly graph (GFA) format. Currently, the program only accept \n"
						+ "                         the assembly graph format used by the assembler SPAdes (de-bruijn \n"
						+ "                         graph) or CANU (overlap). For assembly graphs in other formats: \n"
						+ "                         1. for contigs generated with de-bruijn graph, please ignore this \n"
						+ "                         option and provide k-mer size, the program is able to rebuild the \n"
						+ "                         assembly graph from contigs; and, \n"
						+ "                         2. for contigs generated with overlapping algorithms, please convert \n"
						+ "                         assembly graph to the format used by CANU.\n"
						+ "                         NOTE: it is possible to run the program without an assembly graph, \n"
						+ "                         however, the assembly might be less accurate.\n"
						+ " -k/--kmer-size          Kmer size used for assembly graph construction. Which is supposed \n"
						+ "                         to be the same to the k-mer size used for contig construction. \n"
						+ " -i/--min-identity       Minimum identity between the query and subject sequences \n"
						+ "                         for an alignment record to consider (default 0.90).\n"
						+ " -f/--min-fraction       Minimum alignment fraction of the query sequence (default 0.5).\n"
						+ " -di/--diff-identity     Threshold of the difference of the identity between the primary and secondary \n"
						+ "                         alignments. If the difference is smaller than this value, the query \n"
						+ "                         sequence will be regarded as duplications. Otherwise, the secondary \n"
						+ "                         alignments will be discared (default 0.01).\n"
						+ " -df/--diff-fraction     Threshold of the difference of the alignment fraction between the primary and \n"
						+ "                         secondary alignments. If the difference is smaller than this value, the query \n"
						+ "                         sequence will be regarded as duplications. Otherwise, the secondary \n"
						+ "                         alignments will be discared (default 0.05).\n"
						+ " -o/--out-prefix         Prefix of the output files.\n"
						+ "\n");	
	}

	private String subject_file = null;
	private String query_file = null;
	private String blast_out = null;
	private String asm_graph = null;
	private double min_ident = 90;
	private double min_frac = 0.5;
	private double diff_ident = 0.01;
	private double diff_frac = 0.05;
	private int min_overlap = 10;
	// maximum shift distance for two collinear alignment segments
	// 50% of the smaller segment size
	private double collinear_shift = 0.5;
	private String out_prefix = null;
	private int kmer_size = -1;
	
	private enum ALN_type {blast, bam, mummer};
	
	ALN_type aln_type = null;
	
	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-s", "--subject", true);
			myArgsEngine.add("-q", "--query", true);
			myArgsEngine.add("-a", "--align", true);
			myArgsEngine.add("-blast","--blast", false);
			myArgsEngine.add("-bam","--bam", false);
			myArgsEngine.add("-mummer","--mummer", false);
			myArgsEngine.add("-g","--graph", false);
			myArgsEngine.add("-k", "--kmer-size", true);
			myArgsEngine.add("-i", "--min-identity", true);
			myArgsEngine.add("-f", "--min-fraction", true);
			myArgsEngine.add("-di", "--diff-identity", true);
			myArgsEngine.add("-df", "--diff-fraction", true);
			myArgsEngine.add("-o", "--out-prefix", true);
			myArgsEngine.parse(args);
		}
		
		if (myArgsEngine.getBoolean("-s")) {
			this.subject_file = myArgsEngine.getString("-s");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the subject/reference file.");
		}

		if (myArgsEngine.getBoolean("-q")) {
			this.query_file = myArgsEngine.getString("-q");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the query file.");
		}
		
		if (myArgsEngine.getBoolean("-a")) {
			this.blast_out = myArgsEngine.getString("-a");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the alignment file of query sequences to subject sequences.");
		}
		
		aln_type = ALN_type.blast;
		
		int k_type = (myArgsEngine.getBoolean("-blast") ? 1 : 0) +
				(myArgsEngine.getBoolean("-bam") ? 1 : 0) +
				(myArgsEngine.getBoolean("-mummer") ? 1 : 0) ;
		if(k_type>1) throw new IllegalArgumentException("-blast, -bam and -mummer options are exclusive.");
		
		if(myArgsEngine.getBoolean("-bam")) {
			this.aln_type = ALN_type.bam;
		}
		
		if(myArgsEngine.getBoolean("-mummer")) {
			this.aln_type = ALN_type.mummer;
		}
		
		if (myArgsEngine.getBoolean("-g")) {
			this.asm_graph = myArgsEngine.getString("-g");
		}
		
		if (myArgsEngine.getBoolean("-k")) {
			this.kmer_size = Integer.parseInt(myArgsEngine.getString("-k"));
		}
		
		if (myArgsEngine.getBoolean("-o")) {
			this.out_prefix = myArgsEngine.getString("-o");
			if(new File(out_prefix+".fa").exists() || 
					new File(out_prefix+".map").exists() || 
					new File(out_prefix+"_unplaced.map").exists()) {
				throw new RuntimeException("Output files exist. Please specify a different name.");
			}
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the prefix of output files.");
		}
		
		if (myArgsEngine.getBoolean("-i")) {
			this.min_ident = 100*Math.min(1.0, Double.parseDouble(myArgsEngine.getString("-i")));
		}
		
		if (myArgsEngine.getBoolean("-f")) {
			this.min_frac = Double.parseDouble(myArgsEngine.getString("-f"));
		}
		
		if (myArgsEngine.getBoolean("-di")) {
			this.diff_ident = Double.parseDouble(myArgsEngine.getString("-di"));
		}
		
		if (myArgsEngine.getBoolean("-df")) {
			this.diff_frac = Double.parseDouble(myArgsEngine.getString("-df"));
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		
		
	}
	
}
