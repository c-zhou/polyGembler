package cz1.ngs.tools;

import java.awt.Dimension;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;
import org.jgrapht.DirectedGraph;
import org.jgrapht.ListenableGraph;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.ext.JGraphXAdapter;
import org.jgrapht.graph.AsSubgraph;
import org.jgrapht.graph.DefaultListenableGraph;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.DirectedWeightedPseudograph;

import com.google.common.collect.DiscreteDomain;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;
import com.mxgraph.layout.mxIGraphLayout;
import com.mxgraph.layout.mxOrganicLayout;
import com.mxgraph.swing.mxGraphComponent;

import cz1.ngs.model.AlignmentSegment;
import cz1.ngs.model.Blast6Segment;
import cz1.ngs.model.DirectedWeightedOverlapPseudograph;
import cz1.ngs.model.GFA;
import cz1.ngs.model.OverlapEdge;
import cz1.ngs.model.Sequence;
import cz1.ngs.model.TraceableAlignmentSegment;
import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Executor;
import cz1.util.Utils;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Scaffolder extends Executor {

	private static enum Task {all, link, graphmap, parse, anchor, zzz}
	private Task task_list = Task.zzz;
	
	private final static boolean USE_OS_BUFFER = false;
	// 8Mb buffer size
	private final static int buffSize8Mb = 8388608;
	private final static Writer STD_OUT_BUFFER = USE_OS_BUFFER ? 
			new BufferedWriter(new OutputStreamWriter(System.out), buffSize8Mb) : new OutputStreamWriter(System.out);
			
	private String[][] bamList;
	private boolean[] matePair; // mate-pair library
	private String asm_graph; // assembly graph (GFA) format
	private int num_threads = Runtime.getRuntime().availableProcessors();
	private boolean debug  = false;
	private boolean ddebug = false;
	private String out_prefix = null;
	private int minQual = 0;
	private int jump = 0;
	private int[] inst;
	private int[] maxInst;
	private String subject_file;
	private String link_pref;
	private String graph_pref;
	private String blast_file;
	private String[] fastq_file;
	private double min_ident = 90;
	private double min_cov = 0.5;
	private double diff_ident = 5;
	private double diff_cov = 0.1;
	
	private final static String link_lib = "-l([0-9]+)";
	private final static Pattern link_pat = Pattern.compile(link_lib);
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		switch(this.task_list) {
		case zzz:
			myLogger.info(
					"\n\nChoose a task from listed:\n"
							+ " link                    Count links.\n"
							+ " graphmap                Count links for long reads. \n"
							+ " parse                   Parse scaffolds. \n"
							+ " anchor                  Anchor contigs to generate scaffolds. \n"
							+ " all                     Count links and then parse scaffolds.\n"
							+ "\n");
			break;

			
		case link:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -1                      BAM file for the alignment records of the 1st reads in read pairs.\n"
							+ " -2                      BAM file for the alignment records of the 2nd reads in read pairs.\n"
							+ " -12                     Configuration file for BAM files.\n"
							+ " -mp/--mate-pair         Mate-pair library.\n"
							+ " -s/--subject            The FASTA file contain subject/reference sequences.\n"
							+ " -g/--graph              Assembly graph file (GFA format).\n"
							+ " -j/--jump               Jump (Default: 0).\n"
							+ " -q/--min-qual           Minimum alignment quality (default 0).\n"
							+ " -f/--frag-size          Insert size of the library.\n"
							+ " -t/--threads            Number of threads to use (default 16).\n"
							+ " -d/--debug              Debugging mode will have extra information printed out.\n"
							+ " -dd/--debug-debug       Debugging mode will have more information printed out than -d mode.\n"
							+ " -o/--out-prefix         Prefix of the output files.\n"
							+ "\n");
			break;
			
		case graphmap:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -b/--blast-file         Alignment file generated with BLAST.\n"
							+ " -fq/--fastq             FASTQ file(s). Multiple files seperated with commas.\n"
							+ " -s/--subject            The FASTA file contain subject/reference sequences.\n"
							+ " -g/--graph              Assembly graph file (GFA format).\n"
							+ " -t/--threads            Number of threads to use (default 16).\n"
							+ " -d/--debug              Debugging mode will have extra information printed out.\n"
							+ " -dd/--debug-debug       Debugging mode will have more information printed out than -d mode.\n"
							+ " -o/--out-prefix         Prefix of the output files.\n"
							+ "\n");
			break;
		case parse:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -s/--subject            The FASTA file contain subject/reference sequences. \n"
							+ " -g/--graph              Assembly graph file (GFA format).\n"
							+ " -a                      The file prefix for graphs.\n"
							+ " -l                      The file prefix for links.\n"
							+ " -d/--debug              Debugging mode will have extra information printed out.\n"
							+ " -dd/--debug-debug       Debugging mode will have more information printed out than -d mode.\n"
							+ " -o/--out-prefix         Prefix of the output files.\n"
							+ "\n");
			break;
		case anchor:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -s/--subject            The FASTA file contain subject/reference sequences. \n"
							+ " -g/--graph              Assembly graph file (GFA format).\n"
							+ " -b/--blast-file         Alignment file generated with BLAST.\n"
							+ " -i/--min-identity       Minimum identity between the query and subject sequences \n"
							+ "                         for an alignment record to consider (default 0.90).\n"
							+ " -c/--min-coverage       Minimum alignment coverage of the query sequence (default 0.5).\n"
							+ " -di/--diff-identity     Threshold of the difference of the identity between the primary and secondary \n"
							+ "                         alignments. If the difference is smaller than this value, the query \n"
							+ "                         sequence will be regarded as duplications. Otherwise, the secondary \n"
							+ "                         alignments will be discared (default 0.05).\n"
							+ " -dc/--diff-coverage     Threshold of the difference of the alignment coverage between the primary and \n"
							+ "                         secondary alignments. If the difference is smaller than this value, the query \n"
							+ "                         sequence will be regarded as duplications. Otherwise, the secondary \n"
							+ "                         alignments will be discared (default 0.1).\n"
							+ " -d/--debug              Debugging mode will have extra information printed out.\n"
							+ " -dd/--debug-debug       Debugging mode will have more information printed out than -d mode.\n"
							+ " -o/--out-prefix         Prefix of the output files.\n"
							+ "\n");
			break;
		default:
			throw new RuntimeException("!!!");
		}
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}
		
		switch(args[0].toLowerCase()) {
		case "link":
			this.task_list = Task.link;
			break;
		case "graphmap":
			this.task_list = Task.graphmap;
			break;
		case "parse":
			this.task_list = Task.parse;
			break;
		case "all":
			this.task_list = Task.all;
			break;
		default:
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");	
		}
		
		String[] args2 = new String[args.length-1];
		System.arraycopy(args, 1, args2, 0, args2.length);
		
		
		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-1", null, true);
			myArgsEngine.add("-2", null, true);
			myArgsEngine.add("-12", null, true);
			myArgsEngine.add("-b", "--blast-file", true);
			myArgsEngine.add("-mp", "--mate-pair", false);
			myArgsEngine.add("-s", "--subject", true);
			myArgsEngine.add("-fq", "--fastq", true);
			myArgsEngine.add("-g", "--graph", true);
			myArgsEngine.add("-a", null, true);
			myArgsEngine.add("-l", null, true);
			myArgsEngine.add("-j", "--jump", true);
			myArgsEngine.add("-q", "--min-qual", true);
			myArgsEngine.add("-f", "--frag-size", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-i", "--min-identity", true);
			myArgsEngine.add("-c", "--min-coverage", true);
			myArgsEngine.add("-di", "--diff-identity", true);
			myArgsEngine.add("-dc", "--diff-coverage", true);
			myArgsEngine.add("-d", "--debug", false);
			myArgsEngine.add("-dd", "--debug-debug", false);
			myArgsEngine.add("-o", "--out-prefix", true);
			myArgsEngine.parse(args2);
		}

		switch(this.task_list) {
		case link:
			if (!(myArgsEngine.getBoolean("-1")&&myArgsEngine.getBoolean("-2"))&&!myArgsEngine.getBoolean("-12")) {
				printUsage();
				throw new IllegalArgumentException("Please specify the alignment file(s) using (-1, -2) or -12 option.");
			}
			
			if ((myArgsEngine.getBoolean("-1")||myArgsEngine.getBoolean("-2"))&&myArgsEngine.getBoolean("-12")) {
				printUsage();
				throw new IllegalArgumentException("Options (-1, -2) and -12 are exclusive.");
			}
			
			if (myArgsEngine.getBoolean("-1")!=myArgsEngine.getBoolean("-2")) {
				printUsage();
				throw new IllegalArgumentException("Options (-1, -2) need to be used in pair.");
			}
			
			if (myArgsEngine.getBoolean("-1")&&myArgsEngine.getBoolean("-2")) {
				this.bamList = new String[1][2];
				this.bamList[0][0] = myArgsEngine.getString("-1");
				this.bamList[0][1] = myArgsEngine.getString("-2");
				this.matePair = new boolean[1];
				if(myArgsEngine.getBoolean("-mp")) {
					this.matePair[0] = true;
					myLogger.info("a mate-pair library.");
				}
				this.inst = new int[1];
				this.maxInst = new int[1];
				if (myArgsEngine.getBoolean("-f")) {
					this.inst[0] = Integer.parseInt(myArgsEngine.getString("-f"));
					this.maxInst[0] = 2*this.inst[0];
				} else {
					printUsage();
					throw new IllegalArgumentException("Please specify insert size of the library (-f).");
				}
			} else {
				if(myArgsEngine.getBoolean("-mp"))
					myLogger.info("mate-pair (-mp) indicator ignored.");
				if(myArgsEngine.getBoolean("-f"))
					myLogger.info("insert size (-f) option ignored.");
			}
			
			if (myArgsEngine.getBoolean("-12")) {
				try {
					BufferedReader br = Utils.getBufferedReader(myArgsEngine.getString("-12").trim());
					final List<String> f1 = new ArrayList<String>();
					final List<String> f2 = new ArrayList<String>();
					final List<Boolean> rev = new ArrayList<Boolean>();
					final List<Integer> ins = new ArrayList<Integer>();
					String line;
					String[] s;
					while( (line = br.readLine()) != null) {
						line = line.trim();
						s = line.split("\\s+");
						if(s.length>0) {
							f1.add(s[0].trim());
							f2.add(s[1].trim());
							rev.add(s[2].trim().equals("0")?false:true);
							ins.add(Integer.parseInt(s[3].trim()));
						}
					}
					br.close();
					
					final int n = f1.size();
					this.bamList = new String[n][2];
					this.matePair = new boolean[n];
					this.inst = new int[n];
					this.maxInst = new int[n];
					for(int i=0; i<n; i++) {
						this.bamList[i][0] = f1.get(i);
						this.bamList[i][1] = f2.get(i);
						this.matePair[i] = rev.get(i);
						this.inst[i] = ins.get(i);
						this.maxInst[i] = 2*this.inst[i];
					}
				} catch (IOException e) {
					// TODO Auto-generated catch block
					myLogger.error("configuration file with invalid format.");
					e.printStackTrace();
				}
				
			}
			
			if (myArgsEngine.getBoolean("-q")) {
				this.minQual = Integer.parseInt(myArgsEngine.getString("-q"));
			}
			
			if (myArgsEngine.getBoolean("-j")) {
				this.jump = Integer.parseInt(myArgsEngine.getString("-j"));
			}

			if (myArgsEngine.getBoolean("-t")) {
				int t = Integer.parseInt(myArgsEngine.getString("-t"));
				if(t<this.num_threads) this.num_threads = t;
				this.THREADS = t;
				Constants.omp_threads = this.num_threads;
				myLogger.info("OMP_THREADS = "+this.num_threads);
			}
			break;
		case graphmap:
			if (myArgsEngine.getBoolean("-b")) {
				this.blast_file = myArgsEngine.getString("-b");
			} else {
				throw new IllegalArgumentException("Please specify the BLAST alignment file using -b option.");
			}
			if (myArgsEngine.getBoolean("-fq")) {
				this.fastq_file = myArgsEngine.getString("-fq").split(",");
			} else {
				throw new IllegalArgumentException("Please specify the FASTQ file(s) using -fq option.");
			}
			break;
		case parse:
			if (myArgsEngine.getBoolean("-a")) {
				this.graph_pref = myArgsEngine.getString("-a");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the graph file prefix (-a).");
			}
			
			if (myArgsEngine.getBoolean("-l")) {
				this.link_pref = myArgsEngine.getString("-l");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the link file prefix (-l).");
			}
			break;
		case anchor:
			if (myArgsEngine.getBoolean("-i")) {
				this.min_ident = 100*Double.parseDouble(myArgsEngine.getString("-i"));
			}
			if (myArgsEngine.getBoolean("-di")) {
				this.diff_ident = 100*Double.parseDouble(myArgsEngine.getString("-di"));
			}
			if (myArgsEngine.getBoolean("-c")) {
				this.min_cov = Double.parseDouble(myArgsEngine.getString("-c"));
			}
			if (myArgsEngine.getBoolean("-di")) {
				this.diff_cov = Double.parseDouble(myArgsEngine.getString("-dc"));
			}
			break;
		case all:
			break;
		default:
			throw new RuntimeException("!!!");
		}
		
		if (myArgsEngine.getBoolean("-s")) {
			this.subject_file = myArgsEngine.getString("-s");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the contig file.");
		}

		if (myArgsEngine.getBoolean("-g")) {
			this.asm_graph = myArgsEngine.getString("-g");
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify the graph file.");
		}
		
		if (myArgsEngine.getBoolean("-d")) {
			this.debug = true;
		}

		if (myArgsEngine.getBoolean("-dd")) {
			this.debug  = true;
			this.ddebug = true;
		}

		if (myArgsEngine.getBoolean("-o")) {
			this.out_prefix = myArgsEngine.getString("-o");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the prefix of output files.");
		}
	}

	/***
	private void parseLinkLibrary(String[] args) {
		// TODO Auto-generated method stub
		final TreeMap<Integer, String> lib = new TreeMap<Integer, String>();
		for(int i=0; i<args.length; i++) {
			String arg = args[i];
			if(arg.matches(link_lib)) {
				Matcher m = link_pat.matcher(arg);
				m.find();
				int k = Integer.parseInt(m.group(1));
				if(lib.containsKey(k)) {
					this.printUsage();
					throw new RuntimeException("multiple link library share the same id.");
				}
				lib.put(k, args[++i]);
			}
		}
		final int n = lib.size();
		this.link_files = new String[n];
		int i = 0;
		for(Map.Entry<Integer, String> ent : lib.entrySet()) link_files[i++] = ent.getValue();
	}
	***/
	
	final SamReaderFactory factory =
			SamReaderFactory.makeDefault()
			.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
					SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
			.validationStringency(ValidationStringency.SILENT);

	private final static Object lock = new Object();

	@Override
	public void run() {
		switch(this.task_list) {
		case zzz:
			myLogger.info("Task list is empty!!!");
			break;
		case link:
			this.run_link();
			break;
		case graphmap:
			this.run_graphmap();
			break;
		case parse:
			this.run_parse();
			break;
		case anchor:
			this.run_anchor();
			break;
		case all:
			this.run_all();
			break;
		default:
			throw new RuntimeException("!!!");
		}
		return;
			
	}
	
	private void run_anchor() {
		// TODO Auto-generated method stub
		sub_seqs = Sequence.parseFastaFileWithRevCmpAsMap(subject_file);
		seq_index = new DualHashBidiMap<String, Integer>();
		symm_seqsn = new HashMap<String, String>();
		
		int index = 0;
		for(String seq : sub_seqs.keySet()) seq_index.put(seq, ++index);
		for(String seq : sub_seqs.keySet()) {
			if(!seq.endsWith("'")) {
				symm_seqsn.put(seq, seq+"'");
				symm_seqsn.put(seq+"'", seq);
			}
		}
		
		final GFA gfa = new GFA(this.subject_file, this.asm_graph);
		for(OverlapEdge edge : gfa.edgeSet()) gfa.setEdgeWeight(edge, .0);
		this.maxCC(gfa.gfa());
		myLogger.info("Loading assembly graph done.");
		
		try {
			BufferedReader br = Utils.getBufferedReader(this.blast_file);
			
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void run_all() {
		// TODO Auto-generated method stub
		sub_seqs = Sequence.parseFastaFileWithRevCmpAsMap(subject_file);
		seq_index = new DualHashBidiMap<String, Integer>();
		symm_seqsn = new HashMap<String, String>();
		
		int index = 0;
		for(String seq : sub_seqs.keySet()) seq_index.put(seq, ++index);
		for(String seq : sub_seqs.keySet()) {
			if(!seq.endsWith("'")) {
				symm_seqsn.put(seq, seq+"'");
				symm_seqsn.put(seq+"'", seq);
			}
		}
		
		final GFA gfa = new GFA(this.subject_file, this.asm_graph);
		for(OverlapEdge edge : gfa.edgeSet()) gfa.setEdgeWeight(edge, .0);
		this.maxCC(gfa.gfa());
		myLogger.info("Loading assembly graph done.");
		
		try {
			BufferedWriter bw = Utils.getBufferedWriter(this.out_prefix+".fa");
			String scaff;
			scaff = this.getScaffold(new String[]{"tig18432716","tig18432717","tig18432718","tig64182101","tig18305290","tig18305291","tig01088712'","tig18236803"}, gfa);
			bw.write(Sequence.formatOutput("Seq1 len="+scaff.length(), scaff));
			//scaff = this.getScaffold(new String[]{"tig01889442","tig18236802","tig18236803"}, gfa);
			//bw.write(Sequence.formatOutput("Seq2 len="+scaff.length(), scaff));
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		final Set<String> vset = new HashSet<String>();
		final String source = "tig64138778"; //"tig18270178";
		final Queue<String> vertices = new LinkedList<String>();
		vertices.add(source);
		vertices.add(symm_seqsn.get(source));
		String vertice;
		while(!vertices.isEmpty()) {
			vertice = vertices.poll();
			if(vset.contains(vertice))
				continue;
			vset.add(vertice);
			for(final OverlapEdge edge : gfa.outgoingEdgesOf(vertice)) {
				vertice = gfa.getEdgeTarget(edge);
				vertices.add(vertice);
				vertices.add(symm_seqsn.get(vertice));
			}
			for(final OverlapEdge edge : gfa.incomingEdgesOf(vertice)) {
				vertice = gfa.getEdgeSource(edge);
				vertices.add(vertice);
				vertices.add(symm_seqsn.get(vertice));
			}
		}
		
		final ListenableGraph<String, OverlapEdge> subgraph = new DefaultListenableGraph<>(new AsSubgraph<>(gfa.gfa(), vset));
		myLogger.info("####Subgraph #V "+vset.size()+", #E "+subgraph.edgeSet().size());
		
		SwingUtilities.invokeLater(new Runnable() {
            public void run(){
            	
                JGraphXAdapter<String, OverlapEdge> graphAdapter = 
                        new JGraphXAdapter<String, OverlapEdge>(subgraph);
                mxOrganicLayout layout = new mxOrganicLayout(graphAdapter);
                
                layout.setOptimizeBorderLine(true);
                layout.setOptimizeEdgeCrossing(true);
                layout.setOptimizeEdgeDistance(true);
                layout.setOptimizeEdgeLength(true);
                layout.setOptimizeNodeDistribution(true);
                layout.setFineTuning(true);
                layout.setMaxIterations(10000);
                layout.execute(graphAdapter.getDefaultParent());

            	JFrame frame = new JFrame("Graphmap");
                frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
                frame.add(new mxGraphComponent(graphAdapter));
                frame.pack();
                frame.setSize(800, 800);
                frame.setVisible(true);
                
            }
        });
	}

	private String getScaffold(String[] scaff, GFA gfa) {
		// TODO Auto-generated method stub
		int olap;
		StringBuilder a = new StringBuilder();
		a.setLength(0);
		a.append(sub_seqs.get(scaff[0]).seq_str());
		for(int i=1; i<scaff.length; i++) {
			olap = (int) gfa.getEdge(scaff[i-1],scaff[i]).olapR();
			a.append(sub_seqs.get(scaff[i]).seq_str().substring(olap));
		}
		return a.toString();
	}

	private final int max_out  = 3;
	private final int min_link = 1;
	DirectedWeightedPseudograph<String, DefaultWeightedEdge> g0, g1;
	
	private void run_parse() {
		// TODO Auto-generated method stub
		sub_seqs = Sequence.parseFastaFileWithRevCmpAsMap(subject_file);
		seq_index = new DualHashBidiMap<String, Integer>();
		symm_seqsn = new HashMap<String, String>();
		
		int index = 0;
		for(String seq : sub_seqs.keySet()) seq_index.put(seq, ++index);
		for(String seq : sub_seqs.keySet()) {
			if(!seq.endsWith("'")) {
				symm_seqsn.put(seq, seq+"'");
				symm_seqsn.put(seq+"'", seq);
			}
		}
		
		// final GFA gfa = new GFA(this.subject_file, this.asm_graph);
		// for(OverlapEdge edge : gfa.edgeSet()) gfa.setEdgeWeight(edge, .0);
		// maxCC(gfa.gfa());
		// myLogger.info("Loading assembly graph done.");
		
		g0 = this.readLinkGraph(this.graph_pref+".0.gfa", this.link_pref+".0.links");
		this.stats(g0);
		this.maxCC(g0);
		
		g1 = this.readLinkGraph(this.graph_pref+".1.gfa", this.link_pref+".1.links");
		this.stats(g1);
		this.maxCC(g1);
		
		this.trim(g0);
		this.trim(g1);
		
		this.maxCC(g0);
		this.maxCC(g1);

		final Set<String> outs = new HashSet<String>();
		int round = 0;
		final Set<String> processed = new HashSet<String>();
		final List<List<String>> scaffs = new ArrayList<List<String>>();
		while(!g0.vertexSet().isEmpty()) {
			++round;
			
			outs.clear();
			for(final String v : g0.vertexSet()) {
				if(g0.incomingEdgesOf(v).isEmpty()) 
					outs.add(v);
			}

			final int N = outs.size();
			if(N==0) break;
			
			myLogger.info("==> Parse #round "+round+", "+"#V "+g0.vertexSet().size()+", #outs "+N+".");
			int n = 0;
			
			processed.clear();
			for(final String out : outs) {
				final List<String> scaff = this.findScaff(out);
				++n;
				scaffs.add(scaff);
				for(String s : scaff) processed.add(s);
				if(n%10000==0) myLogger.info("#round "+round+", #outs "+n+"/"+N+" completed.");
			}
			g0.removeAllVertices(processed);
			
			myLogger.info("####parse #round "+round+" completed. Processed #V "+processed.size()+".");
		}
		
		myLogger.info("####Parsing scaffolds completed #scaffolds "+scaffs.size()+
				", unprocessed #contigs "+g0.vertexSet().size());
		
		myLogger.info("######################################");
		try {
			BufferedWriter bw = Utils.getBufferedWriter(this.out_prefix+".scaffs");
			for(final List<String> scaff : scaffs) {
				int L = 0;
				for(final String s : scaff) { 
					bw.write(s+"->");
					L += this.sub_seqs.get(s).seq_ln();
				}	
				bw.write(L+"\n");
			}
			for(final String s : g0.vertexSet()) {
				bw.write(s+"->"+this.sub_seqs.get(s).seq_ln()+"\n");	
			}
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
		myLogger.info("######################################");
		
		/***
		final int nL = this.link_files.length;
		for(int i=0; i<nL; i++) {
			myLogger.info("reading link file "+this.link_files[i]+"...");
			try {
				BufferedReader br = Utils.getBufferedReader(this.link_files[i]);
				String[] s;
				String line;
				int link;
				OverlapEdge edge;
				while( (line=br.readLine())!=null ) {
					s = line.split("\\s+");
					link = Integer.parseInt(s[2]);
					if(gfa.containsEdge(s[0], s[1])) {
						edge = gfa.getEdge(s[0], s[1]);
						gfa.setEdgeWeight(edge, gfa.getEdgeWeight(edge)+link);
					}
				}
				br.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			myLogger.info("done.");
		}

		myLogger.info("######################################");
		final Set<OverlapEdge> edges = new HashSet<OverlapEdge>();
		for(OverlapEdge edge : gfa.edgeSet()) 
			if(gfa.getEdgeWeight(edge)<min_link) 
				edges.add(edge);
		gfa.removeAllEdges(edges);
		myLogger.info("link graph edge trimming (#link>"+(min_link-1)+") completed.");
		maxCC(gfa.gfa(), 10);
		
		
		myLogger.info("######################################");
		final ConnectivityInspector<String,OverlapEdge> connInsp = new ConnectivityInspector<String,OverlapEdge>(gfa.gfa());
		final List<Set<String>> conns = connInsp.connectedSets();
		Collections.sort(conns, new Comparator<Set<String>>() { // sort conns by size: decreasing
			@Override
			public int compare(Set<String> arg0, Set<String> arg1) {
				// TODO Auto-generated method stub
				return arg1.size()-arg0.size();
			}	
		});

		AsSubgraph<String, OverlapEdge> maxCC = new AsSubgraph<String, OverlapEdge>(gfa.gfa(), conns.get(0));
		myLogger.info("maximum connected component,");
		myLogger.info("#V "+maxCC.vertexSet().size());
		myLogger.info("#E "+maxCC.edgeSet().size());
		
		int countv = 0;
		for(String vertex : maxCC.vertexSet()) {
		//for(String vertex : new String[]{"tig13758488"}) {
			if(maxCC.incomingEdgesOf(vertex).isEmpty()) {
				final Set<OverlapEdge> outs = maxCC.outgoingEdgesOf(vertex);
				if(outs.isEmpty()) continue;
				final Set<String> visitedv = new HashSet<String>();
				visitedv.add(vertex);
				//if(ddebug) myLogger.info("#root "+vertex);
				
				final Set<String> nextv = new HashSet<String>();
				String v;
				for(OverlapEdge out : outs) { 
					v = maxCC.getEdgeTarget(out);
					if(!visitedv.contains(v)) {
						nextv.add(v);
						//if(ddebug) myLogger.info(dash(1)+vertex+": "+v);
					}
				}
				
				int counte = 0;
				//while(!nextv.isEmpty()) {
				while(nextv.size()==1) {
					++counte;
					final Set<String> childs = new HashSet<String>();
					for(String next : nextv) {
						visitedv.add(next);
						final Set<OverlapEdge> child = maxCC.outgoingEdgesOf(next);
						for(OverlapEdge c : child) { 
							v = maxCC.getEdgeTarget(c);
							if(!visitedv.contains(v)) {
								childs.add(v);
								//if(ddebug) myLogger.info(dash(counte+1)+next+": "+v);
							}
						}		
					}
					nextv.clear();
					nextv.addAll(childs);
				}
				
				++countv;
				myLogger.info("#"+countv+" vertex "+vertex+": "+counte);
			}
			if(countv==1000) break;
		}
		
		myLogger.info("######################################");
		try {
			BufferedWriter bw = Utils.getBufferedWriter(this.out_prefix+".g");
			for(OverlapEdge edge : gfa.edgeSet()) {
				bw.write(gfa.getEdgeSource(edge));
				bw.write("\t");
				bw.write(gfa.getEdgeTarget(edge));
				bw.write("\t");
				bw.write(scoring.get(edge)+"\n");
			}
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
		
		int[] mstats = new int[10];

		myLogger.info("######################################");
		stats(gfa.gfa(), mstats);
		
		while(mstats[1]!=0||mstats[2]!=0) {
			myLogger.info("######################################");
			final Set<String> vertices = new HashSet<String>();
			for(String vertex : gfa.vertexSet()) {
				if(gfa.outgoingEdgesOf(vertex).isEmpty()||
						gfa.incomingEdgesOf(vertex).isEmpty())
					vertices.add(vertex);
			}
			gfa.removeAllVertices(vertices);
			myLogger.info("remove #"+vertices.size()+" vertices.");
			maxCC(gfa.gfa(), 10);
			stats(gfa.gfa(), mstats);
		}
		**/
	}


	private List<String> findScaff(String first) {
		// TODO Auto-generated method stub
		
		final List<String> scaff = new ArrayList<String>();
		scaff.add(first);
		final List<String> nexts = new ArrayList<String>();
		final List<String> prevs = new ArrayList<String>();
		final Set<String> skips = new HashSet<String>();
		String next = null, prev = null, out;
		
		findNext:
			while(true) {
				out = scaff.get(scaff.size()-1);
				nexts.clear();
				for(final DefaultWeightedEdge edge : g0.outgoingEdgesOf(out)) 
					nexts.add(g0.getEdgeTarget(edge));
				
				next = null;
				if(nexts.size()==0) {
					// no next, break
					
				} else if(nexts.size()==1) {
					// only one next, continue
					next = nexts.get(0);
				} else {
					// ambiguous next, decide if next can be decided
					skips.clear();
					for(final DefaultWeightedEdge edge : g1.outgoingEdgesOf(out))
						skips.add(g1.getEdgeTarget(edge));
					// no skip found
					if(skips.isEmpty()) break;
					
					for(String s : nexts) {
						for(final DefaultWeightedEdge edge : g0.outgoingEdgesOf(s)) {
							if(skips.contains(g0.getEdgeTarget(edge))) {
								if(next==null) next = s;
								else break findNext;
							}
						}
					}
				}
				if(next==null) break;
				
				// now for next we check if no ambiguous prev
				prevs.clear();
				for(final DefaultWeightedEdge edge : g0.incomingEdgesOf(next))
					prevs.add(g0.getEdgeSource(edge));
				prev = null;
				if(prevs.isEmpty()) {
					// not right
					// 'out' will be definitely included here
					throw new RuntimeException("!!!");
				} else if(prevs.size()==1) {
					prev = prevs.get(0);
				} else{
					// found more than one incoming edges
					for(final String s : prevs) {
						skips.clear();
						for(final DefaultWeightedEdge edge : g1.outgoingEdgesOf(s)) 
							skips.add(g1.getEdgeTarget(edge));
						for(final String n : skips) {
							if(g0.containsEdge(next, n)) {
								if(prev==null||prev.equals(s))
									prev = s;
								else break findNext;
							}
						}
					}
				}
				// prev is not out
				if(prev==null || !prev.equals(out)) break;
				
				scaff.add(next);
			}
		
		return scaff;
	}

	private final double max_diff = 0.7;
	private final double max_link = 10d;
	
	private void trim(DirectedWeightedPseudograph<String, DefaultWeightedEdge> graph) {
		// TODO Auto-generated method stub
		/***
		myLogger.info("######################################");
		final Map<DefaultWeightedEdge, Double> scoring = new HashMap<DefaultWeightedEdge, Double>();
		for(final DefaultWeightedEdge edge : graph.edgeSet()) scoring.put(edge, .0);
		for(String vertex : graph.vertexSet()) {
			for(int d=0; d<2; d++) {
				final List<DefaultWeightedEdge> conn = new ArrayList<DefaultWeightedEdge>(
						d==0?graph.outgoingEdgesOf(vertex):graph.incomingEdgesOf(vertex));
				if(conn.isEmpty()) continue; 
				// we score these edges
				final int n = conn.size();
				final double[] score = new double[n];
				for(int i=0; i<n; i++) score[i] = graph.getEdgeWeight(conn.get(i));
				// softmax rescaling
				final double m = StatUtils.max(score);
				for(int i=0; i<n; i++) score[i] = Math.exp(score[i]-m);
				final double s = StatUtils.sum(score);
				for(int i=0; i<n; i++) score[i] /= s;
				for(int i=0; i<n; i++) {
					final DefaultWeightedEdge edge = conn.get(i);
					scoring.put(edge, score[i]+scoring.get(edge));		
				}
			}
		}
		myLogger.info("link graph edge scoring completed.");
		
		myLogger.info("######################################");
		final Set<DefaultWeightedEdge> edges = new HashSet<DefaultWeightedEdge>();
		for(String vertex : graph.vertexSet()) {
			for(int d=0; d<2; d++) {
				final List<DefaultWeightedEdge> conn = new ArrayList<DefaultWeightedEdge>(
						d==0?graph.outgoingEdgesOf(vertex):graph.incomingEdgesOf(vertex));
				if(conn.isEmpty()) continue; 
				double m = Double.NEGATIVE_INFINITY;
				for(DefaultWeightedEdge edge : conn) 
					if(m<scoring.get(edge)) m = scoring.get(edge);
				m *= max_diff;
				for(DefaultWeightedEdge edge : conn)
					if(m<=scoring.get(edge)) edges.add(edge);
			}
		}
		final Set<DefaultWeightedEdge> edge2 = new HashSet<DefaultWeightedEdge>();
		for(DefaultWeightedEdge edge : graph.edgeSet())
			if(!edges.contains(edge)) edge2.add(edge);
		graph.removeAllEdges(edge2);
		
		myLogger.info("link graph edge trimming completed.");
		myLogger.info(edge2.size()+" edges removed from link graph.");
		***/
		
		myLogger.info("######################################");
		final Set<DefaultWeightedEdge> edges = new HashSet<DefaultWeightedEdge>();
		for(String vertex : graph.vertexSet()) {
			for(int d=0; d<2; d++) {
				final List<DefaultWeightedEdge> conn = new ArrayList<DefaultWeightedEdge>(
						d==0?graph.outgoingEdgesOf(vertex):graph.incomingEdgesOf(vertex));
				if(conn.isEmpty()) continue; 
				double m = Double.NEGATIVE_INFINITY;
				for(DefaultWeightedEdge edge : conn) 
					if(m<graph.getEdgeWeight(edge)) m = graph.getEdgeWeight(edge);
				m = Math.min(max_link, m*max_diff);
				for(DefaultWeightedEdge edge : conn)
					if(m<=graph.getEdgeWeight(edge)) edges.add(edge);
			}
		}
		final Set<DefaultWeightedEdge> edge2 = new HashSet<DefaultWeightedEdge>();
		for(DefaultWeightedEdge edge : graph.edgeSet())
			if(!edges.contains(edge)) edge2.add(edge);
		graph.removeAllEdges(edge2);
		
		myLogger.info("link graph edge trimming completed.");
		myLogger.info(edge2.size()+" edges removed from link graph.");
	}

	private DirectedWeightedPseudograph<String, DefaultWeightedEdge> readLinkGraph(final String graph, final String link) {
		// TODO Auto-generated method stub
		final DirectedWeightedPseudograph<String, DefaultWeightedEdge> g = 
				new DirectedWeightedPseudograph<String, DefaultWeightedEdge>(DefaultWeightedEdge.class);
		
		try {
			BufferedReader br_graph = Utils.getBufferedReader(graph);
			String line;
			String[] s;
			while( (line=br_graph.readLine())!=null ) {
				s = line.trim().split("\\s+");
				if(s[0].equals("V"))
					g.addVertex(s[1]);
				else if(s[0].equals("E")||!g.containsEdge(s[1], s[2]))
					g.addEdge(s[1], s[2]);
				else
					throw new RuntimeException("!!!");
			}
			br_graph.close();
			
			BufferedReader br_link  = Utils.getBufferedReader(link);
			while( (line=br_link.readLine())!=null ) {
				s = line.trim().split("\\s+");
				g.setEdgeWeight(g.getEdge(s[0], s[1]), Double.parseDouble(s[2]));
			}
			br_link.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return g;
	}

	private String dash(int n) {
		// TODO Auto-generated method stub
		StringBuilder buff = new StringBuilder();
		for(int i=0; i<n; i++)
			buff.append("-");
		return buff.toString().trim();
	}

	private <V, E, G extends DirectedWeightedPseudograph<V, E>> void stats(G graph, int[] mstats) {
		// TODO Auto-generated method stub
		Arrays.fill(mstats, 0);
		int in1, out;
		for(V vertex : graph.vertexSet()) {
			in1 = graph.outgoingEdgesOf(vertex).size();
			out = graph.incomingEdgesOf(vertex).size();
			if(in1==1&&out==1) ++mstats[0];
			if(out==0) ++mstats[1];
			if(in1==0) ++mstats[2];
		}
		myLogger.info("#vertex   : "+graph.vertexSet().size());
		myLogger.info("#edge     : "+graph.edgeSet().size());
		myLogger.info("#in & out : "+mstats[0]);
		myLogger.info("#in  only : "+mstats[1]);
		myLogger.info("#out only : "+mstats[2]);
	}
	
	private <V, E, G extends DirectedWeightedPseudograph<V, E>> void stats(G graph) {
		// TODO Auto-generated method stub
		stats(graph, new int[10]);
	}

	private <V, E, G extends DirectedWeightedPseudograph<V, E>> void maxCC(G graph, int C) {
		// TODO Auto-generated method stub
		final ConnectivityInspector<V, E> connInsp = new ConnectivityInspector<V, E>(graph);
		final List<Set<V>> conns = connInsp.connectedSets();
		Collections.sort(conns, new Comparator<Set<V>>() { // sort conns by size: decreasing
			@Override
			public int compare(Set<V> arg0, Set<V> arg1) {
				// TODO Auto-generated method stub
				return arg1.size()-arg0.size();
			}	
		});
		myLogger.info("#V: "+graph.vertexSet().size());
		myLogger.info("#E: "+graph.edgeSet().size());
		myLogger.info("#Connected Components: "+conns.size());
		for(int j=0; j<C; j++) {
			long s = 0;
			int n = 0;
			final Set<V> conn = conns.get(j);
			for(V v : conn) {
				s += sub_seqs.get(v).seq_ln();
				if(graph.inDegreeOf(v)==0) ++n;
			}
			myLogger.info("Maximum CC ["+j+"] #V: "+conns.get(j).size()+" "+n+" "+s+" "+conns.get(j).iterator().next());
		}
	}

	private <V, E, G extends DirectedWeightedPseudograph<V, E>> void maxCC(G graph) {
		maxCC(graph, 10);
	}
	
	private void offer(LinkedList<DefaultWeightedEdge> maxE, 
			DefaultWeightedEdge E, 
			DirectedWeightedPseudograph<String, DefaultWeightedEdge> G) {
		// TODO Auto-generated method stub
		int i = 0;
		for(DefaultWeightedEdge e : maxE) {
			if(G.getEdgeWeight(E)<G.getEdgeWeight(e))
				break;
			++i;
		}
		maxE.add(i, E);
		return;
	}

	private class BAMBarcodeIterator {
		private final SamReader[] samReader = new SamReader[2];
		private final SAMRecordIterator[] iter = new SAMRecordIterator[2];
		private final SAMRecord[] samRecord = new SAMRecord[2];

		public BAMBarcodeIterator(final String f1, final String f2) {
			this.samReader[0] = factory.open(new File(f1));
			this.samReader[1] = factory.open(new File(f2));
			this.iter[0] = this.samReader[0].iterator();
			this.iter[1] = this.samReader[1].iterator();
			this.samRecord[0] = iter[0].hasNext() ? iter[0].next() : null;
			this.samRecord[1] = iter[1].hasNext() ? iter[1].next() : null;
		}

		public boolean hasNext() {
			return samRecord[0]!=null && samRecord[1]!=null;
		}

		public List<Set<SAMRecord>> next() {

			if(!this.hasNext()) throw new RuntimeException("!!!");

			List<Set<SAMRecord>> records = new ArrayList<Set<SAMRecord>>();
			String rn = samRecord[0].getReadName();
			if(!rn.equals(samRecord[1].getReadName()))
				throw new RuntimeException("!!!");

			for(int i=0; i<2; i++) {
				final Set<SAMRecord> r = new HashSet<SAMRecord>();
				while( samRecord[i]!=null && rn.equals(samRecord[i].getReadName()) ) {
					r.add(samRecord[i]);
					samRecord[i] = iter[i].next();
				}
				records.add(r);
			}
			
			return records;
		}

		public void close() {
			try {
				for(int i=0; i<2; i++) {
					this.iter[i].close();
					this.samReader[i].close();
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
	
	private class SAMPair {
		private final SAMRecord  first;
		private final SAMRecord second;
		private final String  firstr;
		private final String secondr;
		private final int ins;
		private final int J;
		
		public SAMPair(final SAMRecord first,
				final SAMRecord second,
				final String  firstr,
				final String secondr,
				final int ins,
				final int J) {
			this.first   =   first;
			this.second  =  second;
			this.firstr  =  firstr;
			this.secondr = secondr;
			this.ins = ins;
			this.J = J;
		}
	}
	
	private final static double m_clip = 0.2d; // max clip size (%) to treat an alignment end-to-end
	private final Map<String, Integer> qry_seqs = new HashMap<String, Integer>();
	private final static boolean plot = true;
	
	private void run_graphmap() {
		// TODO Auto-generated method stub
		final GFA gfa = new GFA(this.subject_file, this.asm_graph);
		myLogger.info("Loading assembly graph done.");
		
		sub_seqs = Sequence.parseFastaFileWithRevCmpAsMap(subject_file);
		seq_index = new DualHashBidiMap<String, Integer>();
		symm_seqsn = new HashMap<String, String>();
		
		int index = 0;
		for(String seq : sub_seqs.keySet()) seq_index.put(seq, ++index);
		for(String seq : sub_seqs.keySet()) {
			if(!seq.endsWith("'")) {
				symm_seqsn.put(seq, seq+"'");
				symm_seqsn.put(seq+"'", seq);
			}
		}
		try {
			String line, fq_id;
			double tL = 0;
			for(String fq : this.fastq_file) {
				BufferedReader br_fq = Utils.getBufferedReader(fq);
				while( (line=br_fq.readLine())!=null ) {
					fq_id = line.split("\\s+")[0].substring(1);
					if(qry_seqs.containsKey(fq_id)) throw new RuntimeException("!!!");
					line = br_fq.readLine();
					qry_seqs.put(fq_id, line.length());
					tL += line.length();
					br_fq.readLine();
					br_fq.readLine();
				}
				br_fq.close();
			}
			myLogger.info("####Query sequence #"+qry_seqs.size()+", avg length "+(int)(tL/qry_seqs.size()));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		try {
			BufferedReader br_blast = Utils.getBufferedReader(blast_file);
			Blast6Segment tmp_record = Blast6Segment.blast6Segment(br_blast.readLine());
			Blast6Segment primary_record, secondary_record;
			String qry_sn, sub_sn;
			double aln_frac, clip;
			int qry_ln, sub_ln, qstart, qend, sstart, send, tmp_int;
			final List<Blast6Segment> buff = new ArrayList<Blast6Segment>();
			final List<TraceableAlignmentSegment> alignments = new ArrayList<TraceableAlignmentSegment>();
			TraceableAlignmentSegment alignment;
			
			while(tmp_record!=null) {
				buff.clear();
				buff.add(tmp_record);
				qry_sn = tmp_record.qseqid();
				while( (tmp_record=Blast6Segment.blast6Segment(br_blast.readLine()))!=null
						&&
						tmp_record.qseqid().equals(qry_sn) ) {
					buff.add(tmp_record);
				}
								
				Collections.sort(buff, new AlignmentSegment.QueryCoordinationComparator());
				
				qry_ln = qry_seqs.get(qry_sn);
				alignments.clear();
				for(final Blast6Segment seg : buff) {
					qstart = seg.qstart();
					qend   = seg.qend();
					if(seg.sstart()<seg.send()) {
						sstart = seg.sstart();
						send   = seg.send();
						sub_sn = seg.sseqid();
					} else {
						sstart = seg.send();
						send   = seg.sstart();
						sub_sn = this.symm_seqsn.get(seg.sseqid());
					}
					alignment = new TraceableAlignmentSegment(qry_sn, sub_sn, qstart, qend, sstart, send);
					clip = Math.max(100, Math.max((send-sstart+1)*m_clip, (qend-qstart+1)*m_clip));
					sub_ln = sub_seqs.get(sub_sn).seq_ln();
					alignment.setEndToEnd(qstart<=clip, qend+clip>=qry_ln, sstart<=clip, send+clip>=sub_ln);
					alignment.setClip(qstart-1, qry_ln-qend, sstart-1, sub_ln-send);
					alignment.calcScore();
					
					if(alignment.getEndToEnd()) alignments.add(alignment);
				}
				
				myLogger.info(qry_sn+" "+alignments.size());
				//myLogger.info("####debug break point");
				if(alignments.isEmpty()) {
					myLogger.info(qry_sn+" "+0+"%");
					continue;
				}
				
				if(plot) {
					final DirectedWeightedPseudograph<String, DefaultWeightedEdge> graph = 
							new DirectedWeightedPseudograph<String, DefaultWeightedEdge>(DefaultWeightedEdge.class);
					RangeSet<Integer> cov = TreeRangeSet.create();
					for(final TraceableAlignmentSegment a : alignments)
						cov.add(Range.closed(a.qstart(), a.qend()).canonical(DiscreteDomain.integers()));
					double tcov = 0;
					for(final Range<Integer> r : cov.asRanges()) 
						tcov += r.upperEndpoint()-r.lowerEndpoint();
					myLogger.info(tcov/qry_ln);
					
					
				} else {
					SwingUtilities.invokeLater(new Runnable() {
			            public void run(){
			            	
			                final ListenableGraph<String, DefaultWeightedEdge> graph = 
									new DefaultListenableGraph<>(new DirectedWeightedPseudograph<>(DefaultWeightedEdge.class));
			                
			                final List<String> vs = new ArrayList<String>();
			                String sub_sn;
							for(final TraceableAlignmentSegment a : alignments) {
								sub_sn = a.sseqid();
								vs.add(sub_sn);
								if(!graph.containsVertex(sub_sn)) graph.addVertex(sub_sn);	
							}

							String source, target;
							for(int i=0; i<alignments.size(); i++) {
								source = vs.get(i);
								for(int j=i+1; j<alignments.size(); j++) {
									target = vs.get(j);
									if(gfa.containsEdge(source, target))
										graph.addEdge(source, target);
								}
							}

			                JGraphXAdapter<String, DefaultWeightedEdge> graphAdapter = 
			                        new JGraphXAdapter<String, DefaultWeightedEdge>(graph);
			                mxOrganicLayout layout = new mxOrganicLayout(graphAdapter);
			                
			                layout.setOptimizeBorderLine(true);
			                layout.setOptimizeEdgeCrossing(true);
			                layout.setOptimizeEdgeDistance(true);
			                layout.setOptimizeEdgeLength(true);
			                layout.setOptimizeNodeDistribution(true);
			                layout.setFineTuning(true);
			                layout.setMaxIterations(10000);
			                layout.execute(graphAdapter.getDefaultParent());

			            	JFrame frame = new JFrame("Graphmap");
			                frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			                frame.add(new mxGraphComponent(graphAdapter));
			                frame.pack();
			                frame.setSize(800, 800);
			                frame.setVisible(true);
			                
			            }
			        });
				}
				
				myLogger.info("####debug break point");
			}
			br_blast.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	private final static List<Map<Long, Double>> linkCount = new ArrayList<Map<Long, Double>>(); // link count
	private static long exceed_ins = 0;
	private final static double olap_ext = 0.3;
	private static long readCount = 0;
	
	private Map<String, Sequence> sub_seqs;
	private BidiMap<String, Integer> seq_index;
	private Map<String, String> symm_seqsn;
	
	private void run_link() {
		// TODO Auto-generated method stub
		final int niof = this.bamList.length;
		myLogger.info("Reading alignments from "+niof+" BAM file"+
				(this.bamList.length>1?"s":"")+":");
		for(int i=0; i<niof; i++)
			myLogger.info(this.bamList[i][0]+" paired with "+this.bamList[i][1]+" "+
					(this.matePair[i]?"mate-pair":"paired-end")+" library.");
		myLogger.info("****");
		
		final GFA gfa = new GFA(this.subject_file, this.asm_graph);
		myLogger.info("Loading assembly graph done.");
		
		sub_seqs = Sequence.parseFastaFileWithRevCmpAsMap(subject_file);
		seq_index = new DualHashBidiMap<String, Integer>();
		symm_seqsn = new HashMap<String, String>();
		
		int index = 0;
		for(String seq : sub_seqs.keySet()) seq_index.put(seq, ++index);
		for(String seq : sub_seqs.keySet()) {
			if(!seq.endsWith("'")) {
				symm_seqsn.put(seq, seq+"'");
				symm_seqsn.put(seq+"'", seq);
			}
		}
		
		final List<DirectedWeightedOverlapPseudograph<String>> gfas = new ArrayList<DirectedWeightedOverlapPseudograph<String>>();
		final List<Map<OverlapEdge, Set<String>>> routes = new ArrayList<Map<OverlapEdge, Set<String>>>();
		
		makeJumpGraph(gfa.gfa(), this.jump, gfas, routes);
		myLogger.info("Construction of jump assembly graph done.");
		for(int j=0; j<=jump; j++) {
			myLogger.info("@jump library #"+j+":");
			myLogger.info("#V: "+gfas.get(j).vertexSet().size());
			myLogger.info("#E: "+gfas.get(j).edgeSet().size());
		}
		
		try {
			for(int j=0; j<=jump; j++) {
				BufferedWriter bw = Utils.getBufferedWriter(this.out_prefix+"."+j+".gfa");
				final DirectedWeightedOverlapPseudograph<String> g = gfas.get(j);
				final Map<OverlapEdge, Set<String>> route = routes.get(j);
				for(final String v : g.vertexSet()) 
					bw.write("V\t"+v+"\n");
				for(final OverlapEdge e : g.edgeSet()) {
					bw.write("E\t"+g.getEdgeSource(e)+"\t"+g.getEdgeTarget(e)+"\t"+g.getEdgeOverlap(e)+"\t");
					bw.write(";");
					if(route!=null) 
						for(final String v : route.get(e)) 
							bw.write(v+";");
					bw.write("\n");
				}
				bw.close();
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		for(int j=0; j<=jump; j++) 
			linkCount.add(new HashMap<Long, Double>());
		
		for(int i=0; i<niof; i++) {
			myLogger.info("####Processing "+this.bamList[i][0]+" paired with "+this.bamList[i][1]);
			this.initial_thread_pool();
			final BAMBarcodeIterator iter = new BAMBarcodeIterator(this.bamList[i][0], this.bamList[i][1]);
			while(iter.hasNext()) {
				this.executor.submit(new Runnable() {
					private Set<SAMRecord> sr1;
					private Set<SAMRecord> sr2;
					private boolean matePair;
					private int maxInst;
					
					@Override
					public void run() {
						// TODO Auto-generated method stub
						if(sr1.isEmpty()||sr2.isEmpty()) return;
						
						final List<SAMPair> sampairs = new ArrayList<SAMPair>();
						SAMRecord tmp;
						int inst, s1, s2, e1, e2, olap, j;
						final boolean[] rev = new boolean[2];
						final int[] reflen = new int[2];
						final String[] ref = new String[2];
						final String[] refstr = new String[2];
						
                        synchronized(lock) {
                            ++readCount;
                            if(readCount%1000000==0)
                                myLogger.info("# "+readCount+" reads processed.");
                        }

						try {
							if(ddebug) myLogger.info(sr1.iterator().next().getReadName()+": "+sr1.size()+"x"+sr2.size());
							for(SAMRecord r1 : sr1) {
								// check return
								if(r1.getReadUnmappedFlag() ||
										r1.getMappingQuality()<minQual) 
									continue;
								s1 = r1.getAlignmentStart();
								e1 = r1.getAlignmentEnd();
								rev[0] = r1.getReadNegativeStrandFlag();
								ref[0] = r1.getReferenceName();
								
								for(SAMRecord r2 : sr2) {
									// check return
									if(r2.getReadUnmappedFlag() ||
											r2.getMappingQuality()<minQual) 
										continue;
									s2 = r2.getAlignmentStart();
									e2 = r2.getAlignmentEnd();
									rev[1] = r2.getReadNegativeStrandFlag();
									ref[1] = r2.getReferenceName();
									
									if(!r1.getReadName().
											equals(r2.getReadName()))
										throw new RuntimeException("!!!");
									
									inst = Integer.MAX_VALUE;
									j = -1;
									if(r1.getReferenceIndex().intValue()==r2.getReferenceIndex().intValue()) {
										// mapped to the same contig
										if(rev[0]!=rev[1]) {
											// two reads need to be mapped in reverse direction
											if(rev[1]) {
												tmp = r1;
												r1 = r2;
												r2 = tmp;
											}
											if(s1<s2&&matePair) {
												inst = e2-s1;
											} else if(s1>s2&&!matePair) {
												inst = e1-s2;
											}
											if(inst<=maxInst) return;
										}
									} else {
										// mapped to different contigs
										reflen[0]  = sub_seqs.get(r1.getReferenceName()).seq_ln();
										reflen[1]  = sub_seqs.get(r2.getReferenceName()).seq_ln();

										if(matePair) {
											// this is a mate-pair
											if(rev[0]&&rev[1]) {
												//       0                  1
												// ---------------     -------------
												//   <===                 <===

												// reverse 0
												// ---------------     -------------
												//          ===>          <===
												// reverse 1
												// ---------------     -------------
												//   <===                    ===>

												// reverse 0 and 1 are symmetric
												refstr[0] = r1.getReferenceName();
												refstr[1] = r2.getReferenceName()+"'";
												if((j=getJump(refstr[0], refstr[1]))==-1) continue;
												olap = (int) gfas.get(j).getEdge(refstr[0], refstr[1]).olap();
												if( reflen[0]-s1+1-olap<olap_ext*(r1.getAlignmentEnd()-r1.getAlignmentStart()) ||
														reflen[1]-s2+1-olap<olap_ext*(r2.getAlignmentEnd()-r2.getAlignmentStart()) )
													continue;
												inst = reflen[0]-s1+1+reflen[1]-s2+1-olap;
												
											} else if(rev[0]&&!rev[1]) {
												//       0                  1
												// ---------------     -------------
												//   <===                    ===>

												// reverse 0 & reverse 1
												// ---------------     -------------
												//          ===>          <===
												refstr[0] = r1.getReferenceName();
												refstr[1] = r2.getReferenceName();
												if((j=getJump(refstr[0], refstr[1]))==-1) continue;
												olap = (int) gfas.get(j).getEdge(refstr[0], refstr[1]).olap();
												if( reflen[0]-s1+1-olap<olap_ext*(r1.getAlignmentEnd()-r1.getAlignmentStart()) ||
														e2-olap<olap_ext*(r2.getAlignmentEnd()-r2.getAlignmentStart()) )
													continue;
												inst = reflen[0]-s1+1+e2-olap;
												
											} else if(!rev[0]&&rev[1]) {
												//       0                  1
												// ---------------     -------------
												//          ===>          <===

												// reverse 0 & reverse 1
												// ---------------     -------------
												//   <===                    ===>
												refstr[0] = r2.getReferenceName();
												refstr[1] = r1.getReferenceName();
												if((j=getJump(refstr[0], refstr[1]))==-1) continue;
												olap = (int) gfas.get(j).getEdge(refstr[0], refstr[1]).olap();
												if( e1-olap<olap_ext*(r1.getAlignmentEnd()-r1.getAlignmentStart()) ||
														reflen[1]-s2+1-olap<olap_ext*(r2.getAlignmentEnd()-r2.getAlignmentStart()) )
													continue;
												inst = e1+reflen[1]-s2+1-olap;

											} else {
												//       0                  1
												// ---------------     -------------
												//         ===>               ===>

												// reverse 0
												// ---------------     -------------
												//    <===                    ===>
												// reverse 1
												// ---------------     -------------
												//         ===>          <===
												refstr[0] = r1.getReferenceName()+"'";
												refstr[1] = r2.getReferenceName();
												if((j=getJump(refstr[0], refstr[1]))==-1) continue;
												olap = (int) gfas.get(j).getEdge(refstr[0], refstr[1]).olap();
												if( e1-olap<olap_ext*(r1.getAlignmentEnd()-r1.getAlignmentStart()) ||
														e2-olap<olap_ext*(r2.getAlignmentEnd()-r2.getAlignmentStart()) )
													continue;
												inst = e1+e2-olap;

											}
										} else {
											// this is a paired-end
											if(rev[0]&&rev[1]) {
												//       0                  1
												// ---------------     -------------
												//   <===                 <===

												// reverse 0
												// ---------------     -------------
												//          ===>          <===
												// reverse 1
												// ---------------     -------------
												//   <===                    ===>
												
												// reverse 0 and 1 are symmetric
												refstr[0] = r1.getReferenceName()+"'";
												refstr[1] = r2.getReferenceName();
												if((j=getJump(refstr[0], refstr[1]))==-1) continue;
												olap = (int) gfas.get(j).getEdge(refstr[0], refstr[1]).olap();
												if( e1-olap<olap_ext*(r1.getAlignmentEnd()-r1.getAlignmentStart()) ||
														e2-olap<olap_ext*(r2.getAlignmentEnd()-r2.getAlignmentStart()) )
													continue;
												inst = e1+e2-olap;

											} else if(rev[0]&&!rev[1]) {
												//       0                  1
												// ---------------     -------------
												//   <===                    ===>

												// reverse 0 & reverse 1
												// ---------------     -------------
												//          ===>          <===

												refstr[0] = r2.getReferenceName();
												refstr[1] = r1.getReferenceName();
												if((j=getJump(refstr[0], refstr[1]))==-1) continue;
												olap = (int) gfas.get(j).getEdge(refstr[0], refstr[1]).olap();
												if( e1-olap<olap_ext*(r1.getAlignmentEnd()-r1.getAlignmentStart()) ||
														reflen[1]-s2+1-olap<olap_ext*(r2.getAlignmentEnd()-r2.getAlignmentStart()) )
													continue;
												inst = e1+reflen[1]-s2+1-olap;
												
											} else if(!rev[0]&&rev[1]) {
												//       0                  1
												// ---------------     -------------
												//          ===>          <===

												// reverse 0 & reverse 1
												// ---------------     -------------
												//   <===                    ===>

												refstr[0] = r1.getReferenceName();
												refstr[1] = r2.getReferenceName();
												if((j=getJump(refstr[0], refstr[1]))==-1) continue;
												olap = (int) gfas.get(j).getEdge(refstr[0], refstr[1]).olap();
												if( reflen[0]-s1+1-olap<olap_ext*(r1.getAlignmentEnd()-r1.getAlignmentStart()) ||
														e2-olap<olap_ext*(r2.getAlignmentEnd()-r2.getAlignmentStart()) )
													continue;
												inst = reflen[0]-s1+1+e2-olap;

											} else {
												//       0                  1
												// ---------------     -------------
												//         ===>               ===>

												// reverse 0
												// ---------------     -------------
												//    <===                    ===>
												// reverse 1
												// ---------------     -------------
												//         ===>          <===

												refstr[0] = r1.getReferenceName();
												refstr[1] = r2.getReferenceName()+"'";
												if((j=getJump(refstr[0], refstr[1]))==-1) continue;
												olap = (int) gfas.get(j).getEdge(refstr[0], refstr[1]).olap();
												if( reflen[0]-s1+1-olap<olap_ext*(r1.getAlignmentEnd()-r1.getAlignmentStart()) ||
														reflen[1]-s2+1-olap<olap_ext*(r2.getAlignmentEnd()-r2.getAlignmentStart()) )
													continue;
												inst = reflen[0]-s1+1+reflen[1]-s2+1-olap;
											}
										}

										if(ddebug)
											STD_OUT_BUFFER.write(">>>>"+inst+"\n"+
													r1.getSAMString()+r2.getSAMString()+"<<<<\n");
										
										if(inst>maxInst) continue;
										sampairs.add(new SAMPair(r1, r2, refstr[0], refstr[1], inst, j));
									}
								}
							}
							
							if(sampairs.isEmpty()) return;
							
							double w = 1.0/sampairs.size();
							long refind;
							for(SAMPair sampair : sampairs) {
								
								final Map<Long, Double> linkCount1 = linkCount.get(sampair.J);
								
								refind  = seq_index.get(sampair.firstr);
								refind <<= 32;
								refind += seq_index.get(sampair.secondr);

								synchronized(lock) {
									if(linkCount1.containsKey(refind)) 
										linkCount1.put(refind, linkCount1.get(refind)+w);
									else
										linkCount1.put(refind, w);
								}

								refind  = seq_index.get(symm_seqsn.get(sampair.secondr));
								refind <<= 32;
								refind += seq_index.get(symm_seqsn.get(sampair.firstr));

								synchronized(lock) {
									if(linkCount1.containsKey(refind)) 
										linkCount1.put(refind, linkCount1.get(refind)+w);
									else
										linkCount1.put(refind, w);
								}
							}
						} catch (Exception e) {
							// TODO Auto-generated catch block
							Thread t = Thread.currentThread();
							t.getUncaughtExceptionHandler().uncaughtException(t, e);
							e.printStackTrace();
							executor.shutdown();
							System.exit(1);
						}
					}

					private int getJump(String source, String target) {
						// TODO Auto-generated method stub
						int j = -1;
						for(int i=0; i<=jump; i++) {
							if(gfas.get(i).containsEdge(source, target)) {
								if(j!=-1) 
									throw new RuntimeException("!!!");
								j = i;
							}
						}
						return j;
					}

					public Runnable init(final List<Set<SAMRecord>> records, final boolean matePair, final int maxInst) {
						// TODO Auto-generated method stub
						this.sr1 = records.get(0);
						this.sr2 = records.get(1);
						this.matePair = matePair;
						this.maxInst = maxInst;
						return this;
					}
					
				}.init(iter.next(), this.matePair[i], this.maxInst[i]));
			}
			iter.close();
			this.waitFor();
		}
		
		myLogger.info("################################################################");
		for(int j=0; j<=jump; j++) 
			myLogger.info("#good links jump library #"+j+": "+linkCount.get(j).size()/2);
		myLogger.info("################################################################");
		
		try {
			for(int j=0; j<=jump; j++) {
				BufferedWriter bw = Utils.getBufferedWriter(this.out_prefix+"."+j+".links");
				String source, target;
				final Map<Long, Double> linkCount1 = linkCount.get(j);
				for(long key : linkCount1.keySet()) {
					target = seq_index.getKey((int)  key     );
					source = seq_index.getKey((int) (key>>32));
					bw.write(source+"\t"+target+"\t"+linkCount1.get(key)+"\n");
				}
				bw.close();
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void makeJumpGraph(DirectedWeightedOverlapPseudograph<String> gfa, 
			int jump, 
			List<DirectedWeightedOverlapPseudograph<String>> gfas,
			List<Map<OverlapEdge, Set<String>>> routes) {
		// TODO Auto-generated method stub
		gfas.add(gfa);
		routes.add(null);
		for(int j=0; j<jump; j++) {
			myLogger.info("####"+j);
			final DirectedWeightedOverlapPseudograph<String> g = new DirectedWeightedOverlapPseudograph<String>(OverlapEdge.class);
			final Map<OverlapEdge, Set<String>> route = new HashMap<OverlapEdge, Set<String>>();
			
			final DirectedWeightedOverlapPseudograph<String> prevg = gfas.get(j);
			for(final String v : gfa.vertexSet()) g.addVertex(v);
			int countv = 0;
			OverlapEdge nexte;
			double olap;
			for(final String v : gfa.vertexSet()) {
				++countv;
				if(countv%1000==0) myLogger.info("@graph "+j+": V,"+countv+"; "+"E,"+g.edgeSet().size());
				for(OverlapEdge e : prevg.outgoingEdgesOf(v)) {
					String prevv = prevg.getEdgeTarget(e);
					outer:
						for(OverlapEdge ne : gfa.outgoingEdgesOf(prevv)) {
							String nv = gfa.getEdgeTarget(ne);
							for(int k=0; k<=j; k++) {
								if(gfas.get(k).containsEdge(v, nv))
									continue outer;
							}

							olap = prevg.getEdgeOverlap(e)+gfa.getEdgeOverlap(ne)-sub_seqs.get(prevv).seq_ln();
							if(g.containsEdge(v, nv)) {
								route.get(g.getEdge(v, nv)).add(prevv);
								if(olap<g.getEdgeOverlap(v, nv)) g.setEdgeOverlap(v, nv, olap);
							} else {
								nexte = g.addEdge(v, nv);
								nexte.setOlap(olap);
								final Set<String> r  = new HashSet<String>();
								r.add(prevv);
								route.put(nexte, r);
							}
						}
				}
			}
			gfas.add(g);
			routes.add(route);
		}
		return;
	}
}




