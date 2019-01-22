package cz1.ngs.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
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
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

import org.jgrapht.GraphPath;
import org.jgrapht.alg.shortestpath.BellmanFordShortestPath;
import org.jgrapht.graph.DefaultDirectedWeightedGraph;
import org.jgrapht.graph.DefaultWeightedEdge;

import com.google.common.collect.DiscreteDomain;
import com.google.common.collect.ImmutableRangeSet;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

import cz1.ngs.model.AlignmentSegment;
import cz1.ngs.model.CompoundAlignmentSegment;
import cz1.ngs.model.Sequence;
import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Executor;
import cz1.util.Utils;
import edu.umd.marbl.mhap.impl.MatchResult;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Stitcher extends Executor {

	private static enum Task {all, parse, anchor, stitch, gapclose, zzz}
	private Task task_list = Task.zzz;

	private String align_file;
	private boolean debug  = false;
	private boolean ddebug = false;
	private String out_prefix = null;
	private String coordinate_file;
	private String subject_file;
	private String query_file;
	private String overlap_file = null;
	private int minQual = 20;

	private int num_threads = Runtime.getRuntime().availableProcessors();

	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		switch(this.task_list) {
		case zzz:
			myLogger.info(
					"\n\nChoose a task from listed:\n"
							+ " parse                   Count links.\n"
							+ " stitch                  Anchor contigs to generate scaffolds. \n"
							+ " all                     Run parse followed by stitch.\n"
							+ "\n");
			break;

		case parse:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -s/--subject            The FASTA file contain subject/reference sequences. \n"
							+ " -a/--align              The alignment file. \n"
							+ " -t/--threads            Number of threads to use (default: all available cores).\n"
							+ " -d/--debug              Debugging mode will have extra information printed out.\n"
							+ " -dd/--debug-debug       Debugging mode will have more information printed out than -d mode.\n"
							+ " -o/--out-prefix         Prefix of the output files.\n"
							+ "\n");
			break;
		case anchor:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -s/--subject            The FASTA file contain subject/reference sequences. \n"
							+ " -q/--query              The FASTA file contain query/scaffold/contig sequences. \n"
							+ " -a/--align              The alignment file with query sequences aligned to subject sequences.\n"
							+ " -t/--threads            Number of threads to use (default: all available cores).\n"
							+ " -d/--debug              Debugging mode will have extra information printed out.\n"
							+ " -dd/--debug-debug       Debugging mode will have more information printed out than -d mode.\n"
							+ " -o/--out-prefix         Prefix of the output files.\n"
							+ "\n");
			break;
		case stitch:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -q/--query              The FASTA file contain query/scaffold/contig sequences. \n"
							+ " -c/--coordinate         The coordinate file for pseudomolecules (int AGP format).\n"
							+ " -t/--threads            Number of threads to use (default: all available cores).\n"
							+ " -d/--debug              Debugging mode will have extra information printed out.\n"
							+ " -dd/--debug-debug       Debugging mode will have more information printed out than -d mode.\n"
							+ " -o/--out-prefix         Prefix of the output files.\n"
							+ "\n");
			break;
		case gapclose:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -q/--query              The FASTA file contain query/scaffold/contig sequences. \n"
							+ " -c/--coordinate         The coordinate file for pseudomolecules (int AGP format).\n"
							+ " -t/--threads            Number of threads to use (default: all available cores).\n"
							+ " -ovl/--overlap          Overlap file generated with MHAP.\n"
							+ " -d/--debug              Debugging mode will have extra information printed out.\n"
							+ " -dd/--debug-debug       Debugging mode will have more information printed out than -d mode.\n"
							+ " -o/--out-prefix         Prefix of the output files.\n"
							+ "\n");
			break;
		case all:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -s/--subject            The FASTA file contain subject/reference sequences. \n"
							+ " -q/--query              The FASTA file contain query/scaffold/contig sequences. \n"
							+ " -a/--align              The alignment file with query sequences aligned to subject sequences.\n"
							+ " -ovl/--overlap          Overlap file generated with MHAP.\n"
							+ " -t/--threads            Number of threads to use (default: all available cores).\n"
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
		case "parse":
			this.task_list = Task.parse;
			break;
		case "anchor":
			this.task_list = Task.anchor;
			break;
		case "stitch":
			this.task_list = Task.stitch;
			break;
		case "gapclose":
			this.task_list = Task.gapclose;
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
			myArgsEngine.add("-a", "--align", true);
			myArgsEngine.add("-s", "--subject", true);
			myArgsEngine.add("-q", "--query", true);
			myArgsEngine.add("-c", "--coordinate", true);
			myArgsEngine.add("-ovl", "--overlap", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-d", "--debug", false);
			myArgsEngine.add("-dd", "--debug-debug", false);
			myArgsEngine.add("-o", "--out-prefix", true);
			myArgsEngine.parse(args2);
		}

		switch(this.task_list) {
		case parse:
			if (myArgsEngine.getBoolean("-s")) {
				this.subject_file = myArgsEngine.getString("-s");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the subject file.");
			}
			
			if(myArgsEngine.getBoolean("-a")) {
				this.align_file = myArgsEngine.getString("-a");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the alignment BAM file.");
			}
			break;
		case anchor:
		case all:
			if (myArgsEngine.getBoolean("-s")) {
				this.subject_file = myArgsEngine.getString("-s");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the subject file.");
			}
			
			if(myArgsEngine.getBoolean("-a")) {
				this.align_file = myArgsEngine.getString("-a");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the alignment BAM file.");
			}
			if (myArgsEngine.getBoolean("-q")) {
				this.query_file = myArgsEngine.getString("-q");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the query file.");
			}
			if (myArgsEngine.getBoolean("-ovl")) {
				this.overlap_file = myArgsEngine.getString("-ovl");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the overlap file.");
			}
			break;
		case stitch:
			if (myArgsEngine.getBoolean("-q")) {
				this.query_file = myArgsEngine.getString("-q");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the query file.");
			}
			if (myArgsEngine.getBoolean("-c")) {
				this.coordinate_file = myArgsEngine.getString("-c");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the coordinate file.");
			}
			break;
		case gapclose:
			if (myArgsEngine.getBoolean("-q")) {
				this.query_file = myArgsEngine.getString("-q");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the query file.");
			}
			if (myArgsEngine.getBoolean("-c")) {
				this.coordinate_file = myArgsEngine.getString("-c");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the coordinate file.");
			}
			if (myArgsEngine.getBoolean("-ovl")) {
				this.overlap_file = myArgsEngine.getString("-ovl");
			}
			break;
		default:
			throw new RuntimeException("!!!");
		}
		
		if (myArgsEngine.getBoolean("-t")) {
			int t = Integer.parseInt(myArgsEngine.getString("-t"));
			if(t<this.num_threads) {
				this.num_threads = t;
			}
		}
		this.THREADS = this.num_threads;
		Constants.omp_threads = this.num_threads;
		
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

	@Override
	public void run() {
		// TODO Auto-generated method stub
		switch(this.task_list) {
		case parse:
			this.parse();
			break;
		case anchor:
			this.anchor();
			break;
		case stitch:
			this.stitch();
			break;
		case gapclose:
			this.gapclose();
			break;
		case all:
			if(!new File(this.out_prefix+"_parsed.bam").exists()) {
				this.parse();
			} else {
				myLogger.info("####parse BAM file previously done.");
			}
			this.align_file = this.out_prefix+"_parsed.bam";
			if(!new File(this.out_prefix+"_anchored.agp").exists()) {
				this.anchor();	
			} else {
				myLogger.info("####anchor scaffolds previously done.");
			}
			this.coordinate_file = this.out_prefix+"_anchored.agp";
			if(!new File(this.out_prefix+"_stitch.agp").exists()) {
				this.stitch();
			} else {
				myLogger.info("####stitch scaffolds previously done.");
			}
			this.coordinate_file = this.out_prefix+"_stitch.agp";
			if(!new File(this.out_prefix+"_gapclosed.agp").exists()||
					!new File(this.out_prefix+"_gapclosed.fasta").exists()) {
				this.gapclose();
			} else {
				myLogger.info("####gap closing previously done.");
			}
			break;
		default:
			printUsage();
			throw new IllegalArgumentException("\n\nUnkonwn task list.\n\n");	
		}
	}

	
	private final static double maxOvlErr = 0.10;
	private final static int minOvlSize = 100;
	private final static double scfReusePenalty = 6.0;
	private final static double mxPathSize = 100000;
	private static int bellmanFord = 0;
	
	private void gapclose() {
		// TODO Auto-generated method stub
		final Map<String, Sequence> qry_seqs = Sequence.parseFastaFileWithRevCmpAsMap(query_file);
		final Map<String, int[]> ovlMap = new HashMap<>();
		
		/**
		try {
			BufferedReader br = Utils.getBufferedReader(this.overlap_file);
			String line;
			String[] s;
			while( (line=br.readLine())!=null ) {
				s = line.split("\\s+");
				ovlMap.put(s[0], new int[]{Integer.parseInt(s[1]),
						Integer.parseInt(s[2]),
						Integer.parseInt(s[3]),
						Integer.parseInt(s[4])});
			}
			myLogger.info("#overlaps "+ovlMap.size());
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		**/
		
		try {
			BufferedReader br = Utils.getBufferedReader(this.overlap_file);
			String match, seq1, seq2; 
			String[] parseMatch;
			int s1, s2, e1, e2, l1, l2, d1, d2, tmp;
			double olap, ovlErr;
			int[] concat = new int[2];
			String key, rkey;
			
			while( (match=br.readLine())!=null ){
				parseMatch = match.split("\\s+");
				ovlErr = Double.parseDouble(parseMatch[2]);
				if(ovlErr>maxOvlErr) continue;
				
				seq1 = parseMatch[0];
				l1 = Integer.parseInt(parseMatch[7]);
				d1 = Integer.parseInt(parseMatch[4]);
				s1 = Integer.parseInt(parseMatch[5]);
				e1 = Math.min(l1, Integer.parseInt(parseMatch[6])+11); //MHap default --ordered-kmer-size=12
				
				if(d1==1) {
					tmp = s1;
					s1 = l1-e1;
					e1 = l1-tmp;
				}
				
				seq2 = parseMatch[1];
				l2 = Integer.parseInt(parseMatch[11]);
				d2 = Integer.parseInt(parseMatch[8]);
				s2 = Integer.parseInt(parseMatch[9]);
				e2 = Math.min(l2, Integer.parseInt(parseMatch[10])+11); //MHap default --ordered-kmer-size=12
				
				if(d2==1) {
					tmp = s2;
					s2 = l2-e2;
					e2 = l2-tmp;
				}
				
				olap = (e1-s1+e2-s2)/2.0;
				if(olap<minOvlSize) continue;
				
				Arrays.fill(concat, 0);
				if( s1/olap<=mclip&&(l2-e2)/olap<=mclip ) concat[0]=1;
				if( (l1-e1)/olap<=mclip&&s2/olap<=mclip ) concat[1]=1;
				
				if(concat[0]==concat[1]) continue;
				if(concat[0]==1) {
					// seq1     seq2
					// ======   =======
					// ->            ->
					key  = seq2+(d2==1?"'":"")+"->"+seq1+(d1==1?"'":"");
					rkey = seq1+(d1==1?"":"'")+"->"+seq2+(d2==1?"":"'");
					if(ovlMap.containsKey(key)||ovlMap.containsKey(rkey))
						throw new RuntimeException("!!!");
					ovlMap.put(key,  new int[]{s2, e2, s1, e1});
					ovlMap.put(rkey, new int[]{l1-e1, l1-s1, l2-e2, l2-s2});
				} else if(concat[1]==1) {
					// seq1     seq2
					// ======   =======
					//     ->   ->
					key  = seq1+(d1==1?"'":"")+"->"+seq2+(d2==1?"'":"");
					rkey = seq2+(d2==1?"":"'")+"->"+seq1+(d1==1?"":"'");
					if(ovlMap.containsKey(key)||ovlMap.containsKey(rkey))
						throw new RuntimeException("!!!");
					ovlMap.put(key,  new int[]{s1, e1, s2, e2});
					ovlMap.put(rkey, new int[]{l2-e2, l2-s2, l1-e1, l1-s1});
				}
			}
			br.close();
			
			myLogger.info("#overlaps "+ovlMap.size());
		
			BufferedWriter bw = Utils.getBufferedWriter(this.out_prefix+"_gapclosed.link");
			for(Map.Entry<String, int[]> entry : ovlMap.entrySet()) 
				bw.write(entry.getKey()+"\t"+Utils.cat(entry.getValue()," ")+"\n");
			bw.close();
		
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		final Set<String> unplaced_scf = new HashSet<String>();
		unplaced_scf.addAll(qry_seqs.keySet());
		
		final Set<String> gapclose_jobs = new HashSet<>();
		try {
			BufferedReader br = Utils.getBufferedReader(this.coordinate_file);
			String prev, next;
			String[] ps, ns;
			while((prev=br.readLine())!=null&&(prev.startsWith("#")||prev.contains("GAP"))) {}
			ps = prev.split("\\s+");
			unplaced_scf.remove(ps[5]    );
			unplaced_scf.remove(ps[5]+"'");
			while((next=br.readLine())!=null) {
				while(next!=null&&(next.startsWith("#")||next.contains("GAP"))) 
					next = br.readLine();

				ns = next.split("\\s+");
				if(ps[0].equals(ns[0])) 
					gapclose_jobs.add(ps[5]+(ps[8].equals("+")?"":"'")+"____"+
							ns[5]+(ns[8].equals("+")?"":"'"));
				prev = next;
				ps = ns;
				unplaced_scf.remove(ps[5]    );
				unplaced_scf.remove(ps[5]+"'");
			}
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		myLogger.info("#unplaced scaffolds "+unplaced_scf.size()/2);
		myLogger.info("#gap close jobs "+gapclose_jobs.size());
		
		final DefaultDirectedWeightedGraph<String, DefaultWeightedEdge> ovlGraph = 
				new DefaultDirectedWeightedGraph<>(DefaultWeightedEdge.class);
		DefaultWeightedEdge edge = null;
		for(String key : qry_seqs.keySet()) {
			ovlGraph.addVertex(key+"@5");
			ovlGraph.addVertex(key+"@3");
			edge = ovlGraph.addEdge(key+"@5", key+"@3");
			if(unplaced_scf.contains(key)) {
				ovlGraph.setEdgeWeight(edge, qry_seqs.get(key).seq_ln());
			} else {
				ovlGraph.setEdgeWeight(edge, qry_seqs.get(key).seq_ln()*scfReusePenalty);
			}
		}
		
		String[] vs;
		int[] ovl;
		for(String key : ovlMap.keySet()) {
			vs = key.split("->");
			ovl = ovlMap.get(key);
			edge = ovlGraph.addEdge(vs[0]+"@3", vs[1]+"@5");
			ovlGraph.setEdgeWeight(edge, (ovl[0]-ovl[1]+ovl[2]-ovl[3])/2.0);
		}
		myLogger.info("#OVL graph #V "+ovlGraph.vertexSet().size());
		myLogger.info("#OVL graph #E "+ovlGraph.edgeSet().size()  );
		
		final Map<String, GraphPath<String, DefaultWeightedEdge>> results_collector = new HashMap<>();
		this.initial_thread_pool();
		for(final String gapclose_job : gapclose_jobs) {
			this.executor.submit(new Runnable() {
				private String gapclose_job;
				
				@Override
				public void run() {
					// TODO Auto-generated method stub
					try {
						synchronized(lock) {
							++bellmanFord;
							if(bellmanFord%1000==0)
								myLogger.info("#Bellman-Ford shortest path completed "+bellmanFord);
						}
						
						String[] parse = this.gapclose_job.split("____");
						final String source = parse[0]+"@3";
						final String target = parse[1]+"@5";
						if( !ovlGraph.containsVertex(source) ||
								!ovlGraph.containsVertex(target) ||
								ovlGraph.containsEdge(source, target)) 
							return;
						GraphPath<String, DefaultWeightedEdge> path = 
								BellmanFordShortestPath.findPathBetween(ovlGraph, source, target);
						if(path!=null) {
							synchronized(lock) {
								results_collector.put(gapclose_job, path);
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

				public Runnable init(String gapclose_job) {
					// TODO Auto-generated method stub
					this.gapclose_job = gapclose_job;
					return this;
				}
			}.init(gapclose_job));
		}
		this.waitFor();
		myLogger.info("#results in collector "+results_collector.size());
	
		try {
			BufferedReader br = Utils.getBufferedReader(this.coordinate_file);
			BufferedWriter bw_agp = Utils.getBufferedWriter(this.out_prefix+"_gapclosed.agp");
			BufferedWriter bw_scf = Utils.getBufferedWriter(this.out_prefix+"_gapclosed.fasta");
			
			final Map<String, List<String[]>> blocks = new HashMap<>();
			String line;
			String[] s, s1;
			while( (line=br.readLine())!=null ) {
				if(line.startsWith("#")) {
					bw_agp.write(line+"\n");
					continue;
				}
				s = line.split("\\s+");
				if(!blocks.containsKey(s[0]))
					blocks.put(s[0], new ArrayList<>());
				blocks.get(s[0]).add(s);
			}
			
			final List<String> chrs = new ArrayList<>();
			chrs.addAll(blocks.keySet());
			Collections.sort(chrs);
			
			final StringBuilder seq = new StringBuilder();
			String prev, next, ovl_key, source, target;
			int ext, qst, len, chrPos, n, bs;
            double mxgap;
			List<String[]> block;
			final Set<String> placed = new HashSet<>();
		    GraphPath<String, DefaultWeightedEdge> bfPath;
            List<DefaultWeightedEdge> edgeList;
            int baseClosed = 0;
            int gapClosed = 0;
            boolean fill;

			for(final String chr : chrs) {
				seq.setLength(0);
				block = blocks.get(chr);
				n = block.size();
				s = block.get(0);
				next = s[5]+(s[8].equals("+")?"":"'");
				len = qry_seqs.get(next).seq_ln();
				chrPos = 1;
				ext = len;
				bs = 1;
				bw_agp.write(chr+"\t"+chrPos+"\t"+(chrPos+ext-1)+"\t"+bs+
						"\tW\t"+s[5]+"\t1\t"+len+"\t"+s[8]+"\n");
				seq.append(qry_seqs.get(next).seq_str());
				chrPos += ext;
				bs += 1;
				placed.add(next.replaceAll("'$", ""));
				
				for(int i=2; i<n; i+=2) {
					s1 = block.get(i-1);
					s = block.get(i);
					prev = next;
					next = s[5]+(s[8].equals("+")?"":"'");
					ovl_key = prev+"____"+next;
			        fill = false;
                    if(results_collector.containsKey(ovl_key)) {
                        bfPath = results_collector.get(ovl_key);
                        mxgap = Math.max(Integer.parseInt(s1[6])*6.0, mxPathSize);
                        
                        if(bfPath.getWeight()<=mxgap) fill = true;
                    }

                    if(fill) {
                        target = prev;
                        edgeList = results_collector.get(ovl_key).getEdgeList();
					
                        for(int j=0; j<edgeList.size(); j+=2) {
							source = target;
							target = ovlGraph.getEdgeTarget(edgeList.get(j)).replaceAll("@5$", "");
							
							len = qry_seqs.get(target).seq_ln();
							qst = ovlMap.get(source+"->"+target)[3];
							if(len==qst) continue;
							
							ext = mgap;
							bw_agp.write(chr+"\t"+chrPos+"\t"+(chrPos+ext-1)+"\t"+bs+
									"\tN\tGAP\t"+ext+"\tyes\t+\n");
							seq.append(Sequence.polyN(ext));
							chrPos += ext;
							bs += 1;
							
							ext = len-qst;
							bw_agp.write(chr+"\t"+chrPos+"\t"+(chrPos+ext-1)+"\t"+bs+
									"\tW\t"+target.replaceAll("'$", "")+"\t"+(qst+1)+"\t"+
									len+"\t"+(target.endsWith("'")?"-":"+")+"\n");
							seq.append(qry_seqs.get(target).seq_str().substring(qst));
							chrPos += ext;
							bs += 1;
							
							placed.add(target.replaceAll("'$", ""));
						}

                        gapClosed += 1;
                        baseClosed += Integer.parseInt(s1[6]);
					} else {
						ext = Integer.parseInt(s1[6]);
						bw_agp.write(chr+"\t"+chrPos+"\t"+(chrPos+ext-1)+"\t"+bs+
								"\tN\tGAP\t"+ext+"\tyes\t+\n");
						seq.append(Sequence.polyN(ext));
						chrPos += ext;
						bs += 1;
						
						len = qry_seqs.get(next).seq_ln();
						ext = len;
						bw_agp.write(chr+"\t"+chrPos+"\t"+(chrPos+ext-1)+"\t"+bs+
								"\tW\t"+s[5]+"\t1\t"+len+"\t"+s[8]+"\n");
						seq.append(qry_seqs.get(next).seq_str());
						chrPos += ext;
						bs += 1;
					}
					
					placed.add(next.replaceAll("'$", ""));
				}
				bw_scf.write(Sequence.formatOutput(chr, seq.toString()));
			}
			
            myLogger.info("#gap closed "+gapClosed);
            myLogger.info("#base closed "+baseClosed+"bp");

			for(String seqid : qry_seqs.keySet()) {
				if(!seqid.endsWith("'")&&!placed.contains(seqid))
					bw_scf.write(qry_seqs.get(seqid).formatOutput());
			}
			
			br.close();
			bw_agp.close();
			bw_scf.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

	private void stitch() {
		// TODO Auto-generated method stub
		final Map<String, Sequence> qry_seqs = Sequence.parseFastaFileWithRevCmpAsMap(query_file);
		
		final Set<String> overlap_jobs = new HashSet<>();
		try {
			BufferedReader br = Utils.getBufferedReader(this.coordinate_file);
			String prev, next;
			String[] ps, ns;
			while((prev=br.readLine())!=null&&(prev.startsWith("#")||prev.contains("GAP"))) {}
			ps = prev.split("\\s+");
			while((next=br.readLine())!=null) {
				while(next!=null&&(next.startsWith("#")||next.contains("GAP"))) 
					next = br.readLine();

				ns = next.split("\\s+");
				if(ps[0].equals(ns[0])) 
					overlap_jobs.add(ps[5]+(ps[8].equals("+")?"":"'")+"____"+
							ns[5]+(ns[8].equals("+")?"":"'"));
				prev = next;
				ps = ns;
			}
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		final MHapWrapper mhap = new MHapWrapper(maxOvlErr);
		myLogger.info("####overlap jobs "+overlap_jobs.size());
		final Map<String, String> overlaps = new HashMap<>();

		this.initial_thread_pool();
		for(final String overlap_job : overlap_jobs) {

			this.executor.submit(new Runnable() {
				private String overlap_job;

				@Override
				public void run() {
					// TODO Auto-generated method stub

					try {
						String[] seqPair = overlap_job.split("____");

						final Sequence seq1 = qry_seqs.get(seqPair[0]);
						final Sequence seq2 = qry_seqs.get(seqPair[1]);

						List<MatchResult> matches = mhap.search(seq1, seq2);
						if(matches.isEmpty()) return;

						if(matches.size()>1) Collections.sort(matches);
						// parse overlap results two find overlaps
						MatchResult match = matches.get(0);
						if(!seq1.seq_sn().equals(match.getFromId().getHeader())||
								!seq2.seq_sn().equals(match.getToId().getHeader()))
							throw new RuntimeException("!!!");
						String[] parseMatch = match.toString().split("\\s+");
						int a1 = Integer.parseInt(parseMatch[5]);
						int a2 = Math.min(Integer.parseInt(parseMatch[7]), 
								Integer.parseInt(parseMatch[6]) +11); //MHap default --ordered-kmer-size=12
						int b1 = Integer.parseInt(parseMatch[9]);
						int b2 = Math.min(Integer.parseInt(parseMatch[11]),
								Integer.parseInt(parseMatch[10])+11); //MHap default --ordered-kmer-size=12
						int clipAt3End = seq1.seq_ln()-a2;
						int clipAt5End = b1;
						double olap = (a2-a1+b2-b1)/2.0;
						if(olap<minOvlSize||clipAt3End/olap>mclip||clipAt5End/olap>mclip)
							return;
						synchronized(lock) {
							overlaps.put(this.overlap_job, match.toString());
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

				public Runnable init(final String overlap_job) {
					// TODO Auto-generated method stub
					this.overlap_job = overlap_job;
					return this;
				}

			}.init(overlap_job));
		}
		this.waitFor();

		String overlap_file = this.out_prefix+"_stitch.ovl";
		try {
			BufferedWriter bw = Utils.getBufferedWriter(overlap_file);
			for(String key : overlaps.keySet())
				bw.write(overlaps.get(key)+"\n");
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		final Map<String, Integer> ovlMap = new HashMap<>();
		try {
			BufferedReader br = Utils.getBufferedReader(overlap_file);
			String line;
			String[] s;
			while((line=br.readLine())!=null) {
				s = line.split("\\s+");
				ovlMap.put(s[0]+"____"+s[1], Math.min(qry_seqs.get(s[1]).seq_ln(), 
						Integer.parseInt(s[10])+11)); //MHap default --ordered-kmer-size=12
			}
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		try {
			BufferedReader br = Utils.getBufferedReader(this.coordinate_file);
			BufferedWriter bw_agp = Utils.getBufferedWriter(this.out_prefix+"_stitch.agp");
			BufferedWriter bw_scf = Utils.getBufferedWriter(this.out_prefix+"_stitch.fasta");
			
			final Map<String, List<String[]>> blocks = new HashMap<>();
			String line;
			String[] s, s1;
			while( (line=br.readLine())!=null ) {
				if(line.startsWith("#")) {
					bw_agp.write(line+"\n");
					continue;
				}
				s = line.split("\\s+");
				if(!blocks.containsKey(s[0]))
					blocks.put(s[0], new ArrayList<>());
				blocks.get(s[0]).add(s);
			}
			
			final List<String> chrs = new ArrayList<>();
			chrs.addAll(blocks.keySet());
			Collections.sort(chrs);
			
			final StringBuilder seq = new StringBuilder();
			String prev, next, ovl_key;
			int ext, qst, len, chrPos, n;
			List<String[]> block;
			final Set<String> placed = new HashSet<>();
			
			for(final String chr : chrs) {
				seq.setLength(0);
				block = blocks.get(chr);
				n = block.size();
				s = block.get(0);
				next = s[5]+(s[8].equals("+")?"":"'");
				len = qry_seqs.get(next).seq_ln();
				chrPos = 1;
				ext = len;
				bw_agp.write(s[0]+"\t"+chrPos+"\t"+(chrPos+ext-1)+"\t"+s[3]+
						"\t"+s[4]+"\t"+s[5]+"\t1\t"+len+"\t"+s[8]+"\n");
				seq.append(qry_seqs.get(next).seq_str());
				chrPos += ext;
				placed.add(next.replaceAll("'$", ""));
				
				for(int i=2; i<n; i+=2) {
					s1 = block.get(i-1);
					s = block.get(i);
					prev = next;
					next = s[5]+(s[8].equals("+")?"":"'");
					ovl_key = prev+"____"+next;
					len = qry_seqs.get(next).seq_ln();
					if(ovlMap.containsKey(ovl_key)) {
						qst = ovlMap.get(ovl_key);
						if(len==qst) continue;
						
						ext = mgap;
						bw_agp.write(s1[0]+"\t"+chrPos+"\t"+(chrPos+ext-1)+"\t"+s1[3]+
								"\t"+s1[4]+"\t"+s1[5]+"\t"+ext+"\tyes\t+\n");
						seq.append(Sequence.polyN(ext));
						chrPos += ext;
						
						ext = len-qst;
						bw_agp.write(s[0]+"\t"+chrPos+"\t"+(chrPos+ext-1)+"\t"+s[3]+
								"\t"+s[4]+"\t"+s[5]+"\t"+(qst+1)+"\t"+len+"\t"+s[8]+"\n");
						seq.append(qry_seqs.get(next).seq_str().substring(qst));
						chrPos += ext;
					} else {
						ext = Integer.parseInt(s1[6]);
						bw_agp.write(s1[0]+"\t"+chrPos+"\t"+(chrPos+ext-1)+"\t"+s1[3]+
								"\t"+s1[4]+"\t"+s1[5]+"\t"+ext+"\tyes\t+\n");
						seq.append(Sequence.polyN(ext));
						chrPos += ext;
						
						ext = len;
						bw_agp.write(s[0]+"\t"+chrPos+"\t"+(chrPos+ext-1)+"\t"+s[3]+
								"\t"+s[4]+"\t"+s[5]+"\t1\t"+len+"\t"+s[8]+"\n");
						seq.append(qry_seqs.get(next).seq_str());
						chrPos += ext;
					}
					
					placed.add(next.replaceAll("'$", ""));
				}
				bw_scf.write(Sequence.formatOutput(chr, seq.toString()));
			}
			
			for(String seqid : qry_seqs.keySet()) {
				if(!seqid.endsWith("'")&&!placed.contains(seqid))
					bw_scf.write(qry_seqs.get(seqid).formatOutput());
			}
			
			br.close();
			bw_agp.close();
			bw_scf.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private final SamReaderFactory factory =
			SamReaderFactory.makeDefault()
			.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
					SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
			.validationStringency(ValidationStringency.SILENT);

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

		public List<SAMRecord> next() {

			if(!this.hasNext()) throw new RuntimeException("!!!");

			List<SAMRecord> records = new ArrayList<SAMRecord>();
			String query = samRecord.getReadName();

			while( samRecord!=null && samRecord.getReadName().equals(query) ) {
				if(!samRecord.getReadUnmappedFlag() && 
						!samRecord.isSecondaryAlignment())
					records.add(samRecord);
				samRecord = iter.hasNext()?iter.next():null;
			}

			return records;
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

	private final static double mxgr = 3.0;
	private final static double diff = 0.05;
	private final static int molap = 300;
	private final static double mclip = 0.05;
	private final static int mgap = 11;
	private final static int[] bs = new int[]{10000, 5000, 1000, 500, 100};
	
	private void anchor() {
		// TODO Auto-generated method stub
		final Map<String, Sequence> qry_seqs = Sequence.parseFastaFileWithRevCmpAsMap(query_file);
		final Map<String, Sequence>  sub_seqs = Sequence.parseFastaFileAsMap(subject_file);
		
		final Map<String, ImmutableRangeSet<Integer>> sub_gaps = new HashMap<>();

		for(Map.Entry<String, Sequence> entry : sub_seqs.entrySet()) {
			String seq_sn = entry.getKey();
			String seq_str = entry.getValue().seq_str();

			final RangeSet<Integer> range_set = TreeRangeSet.create();
			char c;
			for(int j=0; j<seq_str.length(); j++) {
				c = seq_str.charAt(j);
				if(c=='N'||c=='n')
					range_set.add(Range.closed(j, j).
							canonical(DiscreteDomain.integers()));
			}
			sub_gaps.put(seq_sn, ImmutableRangeSet.copyOf(range_set));
		}
		
		final Map<String, List<CompoundAlignmentSegment>> initPlace = new HashMap<>();

		// read alignment file and place the query sequences
		final BAMBarcodeIterator iter = new BAMBarcodeIterator(this.align_file);

		List<SAMRecord> query_recs;
		String qry_id, sub_id;
		final Map<String, List<SAMRecord>> blocks = new HashMap<>();
		List<SAMRecord> block;
		RangeSet<Integer> qryCov, subCov;
		int qstart, qend, sstart, send;
		Range<Integer> span;
		String bx;
		
		while(iter.hasNext()) {
			query_recs = iter.next();
			if(query_recs.isEmpty()) continue;
			
			blocks.clear();
			for(SAMRecord rec : query_recs) {
				bx = rec.getStringAttribute("BX");
				if(!blocks.containsKey(bx)) 
					blocks.put(bx, new ArrayList<>());
				blocks.get(bx).add(rec);
			}
			
			for(String key : blocks.keySet()) {
				block = blocks.get(key);
				qryCov = getQryCovs(block);
				subCov = getSubCovs(block);
				qry_id = block.get(0).getReadName();
				if(reverse(block)) qry_id+="'";
				sub_id = block.get(0).getReferenceName();

				if(!initPlace.containsKey(sub_id))
					initPlace.put(sub_id, new ArrayList<>());
				span = qryCov.span();
				qstart = span.lowerEndpoint();
				qend   = span.upperEndpoint()-1;
				span = subCov.span();
				sstart = span.lowerEndpoint();
				send   = span.upperEndpoint()-1;
				initPlace.get(sub_id).add(new CompoundAlignmentSegment(
						key,      // a universal key for the compount alignment segment
						qry_id,   // query (e.g., gene) sequence id
						sub_id,   // subject (e.g., reference genome) sequence id
						qstart,   // start of alignment in query
						qend,     // end of alignment in query
						sstart,   // start of alignment in subject
						send,     // end of alignment in subject
						subCov,
						qryCov,
						getCovFromCovs(subCov),
						getCovFromCovs(qryCov),
						block));
			}

		}
		iter.close();

		// now process each subject
		
		final BufferedWriter bw_out = Utils.getBufferedWriter(this.out_prefix+"_anchored.agp");
		try {
			bw_out.write("##agp-version  2.0\n");
			bw_out.write("#ASSEMBLY NAME: batatas\n");
			bw_out.write("#DESCRIPTION: Pseudochromosome assembly\n");
			bw_out.write("#PROGRAM: PolyGembler stitcher\n");
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		this.initial_thread_pool();
		
		final List<String> sub_ids = new ArrayList<>();
		sub_ids.addAll(initPlace.keySet());
		Collections.sort(sub_ids);
		for(final String sub : sub_ids) {
		//for(final String sub : new String[]{"chr01_itf"}) {
			if(!sub.startsWith("chr")) continue;
	 		this.executor.submit(new Runnable() {
				private String sub_id;

				@Override
				public void run() {
					// TODO Auto-generated method stub
					try {
						final List<CompoundAlignmentSegment> init = initPlace.get(sub_id);
						final TreeMap<Range<Integer>, ScaffBlock> asmb = new TreeMap<>(
								new Comparator<Range<Integer>>() {
									@Override
									public int compare(Range<Integer> r1, Range<Integer> r2) {
										// TODO Auto-generated method stub
										int ld = Integer.compare(r1.lowerEndpoint(), r2.lowerEndpoint());
										return ld==0?Integer.compare(r1.upperEndpoint(),r2.upperEndpoint()):ld;
									}
								});
						final RangeSet<Integer> asmSubCov = TreeRangeSet.create();
						final int subL = sub_seqs.get(this.sub_id).seq_ln();
						final int[] subDepth = new int[subL];
						RangeSet<Integer> cov;
						int lb, ub;
						for(CompoundAlignmentSegment as : init) {
							cov = as.getSubCov();
							for(Range<Integer> r : cov.asRanges()) {
								lb = r.lowerEndpoint()-1;
								ub = r.upperEndpoint()-1;
								for(int i=lb; i<ub; i++)
									++subDepth[i];
							}
						}

						Collections.sort(init, new AlignmentSegment.SubjectCoordinationComparator());
						// remove all contained segments
						CompoundAlignmentSegment prev, next;
						int pstart, pend, nstart;
						ImmutableRangeSet<Integer> pSubCov, nSubCov;
						outerloop:
							for(int i=0; i<init.size(); i++) {
								prev = init.get(i);
								if(prev==null) continue;
								pstart = prev.sstart();
								pend   = prev.send();
								pSubCov = ImmutableRangeSet.copyOf(prev.getSubCov());
								for(int j=i+1; j<init.size(); j++) {
									next = init.get(j);
									if(next==null) continue;
									nstart = next.sstart();
									nSubCov = ImmutableRangeSet.copyOf(next.getSubCov());
									if(nstart>pend) 
										continue outerloop;
									if(pstart==nstart) {
										// it's possible that prev is contained
										if(contained(nSubCov, pSubCov)) {
											init.set(i, null);
											continue outerloop;
										}
									}
									if(contained(pSubCov, nSubCov)) 
										init.set(j, null);
								}
							}
						
						
						for(int iter = 0; iter<bs.length; iter++) {
							// block size from large to small
							// extract segments with size not smaller than b
							final int b = bs[iter];
							myLogger.info("Iteration #"+iter+", block size "+b);
							final List<CompoundAlignmentSegment> cands = new ArrayList<>();
							for(int i=0; i<init.size(); i++) {
								CompoundAlignmentSegment cand = init.get(i);
								if(cand!=null&&cand.getQryLen()>=b) {
									cands.add(cand);
									init.set(i, null);
								}
							}
							myLogger.info("  #candidates "+cands.size());
							
							// calculate score for each candidate 
							final ImmutableRangeSet<Integer> complAsmSubCov = 
									ImmutableRangeSet.copyOf(asmSubCov.complement());
							double score;
							int upper, lower;
							for(CompoundAlignmentSegment cand : cands) {
								score = .0;
								for(Range<Integer> r : complAsmSubCov.intersection(cand.getSubCov()).asRanges()) {
									lower = r.lowerEndpoint()-1;
									upper = r.upperEndpoint()-1;
									for(int i=lower; i<upper; i++) {
										score += 1.0/subDepth[i];
									}
								}
								if(Double.isInfinite(score)) throw new RuntimeException("infinite score!!!");
								cand.setScore(score);
							}
							
							Collections.sort(cands, new Comparator<CompoundAlignmentSegment>() {

								@Override
								public int compare(CompoundAlignmentSegment o1, CompoundAlignmentSegment o2) {
									// TODO Auto-generated method stub
									return Double.compare(o2.getScore(), o1.getScore());
								}
								
							});
							
							for(CompoundAlignmentSegment cand : cands) {
								myLogger.info("  "+cand.qseqid()+"\t"+String.format("%.2f", cand.getScore())+"\t"+
										cand.getSubCov().toString());
							}
							
							// try to fill gaps on the molecule with these segments
							ImmutableRangeSet<Integer> subCov, subExt;
							int sstart, send, qstart, qend, subExtLen, olap, clip3, clip5, qlen, pos;
							Range<Integer> subRange, floorRange, ceilingRange, range, floorInter, ceilingInter;
							CompoundAlignmentSegment cand, bseg;
							SortedMap<Range<Integer>, ScaffBlock> subMap;
							List<Range<Integer>> interactions = new ArrayList<>();
							ScaffBlock scaffBlock, nextBlock;
							for(int i=0; i<cands.size(); i++) {
								cand = cands.get(i);
								if(cand.getKey().equals("b00079575"))
									myLogger.info("");
								myLogger.info("Prcessing candidate "+cand.getKey()+": "+cand.qseqid()+"\t"+
										String.format("%.2f", cand.getScore())+"\t"+cand.getSubCov().toString());
								
								qlen   = qry_seqs.get(cand.qseqid()).seq_ln();
								qstart = cand.qstart();
								qend   = cand.qend();
								sstart = cand.sstart();
								send   = cand.send();
								subCov = ImmutableRangeSet.copyOf(cand.getSubCov());
								subExt = subCov.intersection(asmSubCov.complement());
								subExtLen = getCovFromCovs(subExt);
								
								if(subExtLen==0) continue;
								
								subRange = Range.closed(sstart, send).canonical(DiscreteDomain.integers());
								// now find keys in asmb interacts subRange
								floorRange = asmb.floorKey(subRange);
								subMap = floorRange==null?asmb.tailMap(subRange):asmb.tailMap(floorRange);
								interactions.clear();
								for(Range<Integer> key : subMap.keySet()) {
									if(key.lowerEndpoint()>send) break;
									if(subRange.isConnected(key)&&!subRange.intersection(key).isEmpty())
										interactions.add(key);
								}
								boolean assembled = false;
								if(interactions.isEmpty()) {
									// no interactions
									scaffBlock = new ScaffBlock();
									scaffBlock.segments.add(cand);
									scaffBlock.ranges.add(Range.closed(1, qlen).canonical(DiscreteDomain.integers()));
									scaffBlock.clipAt5End = qstart-1;
									scaffBlock.clipAt3End = qlen-qend;
									asmb.put(subRange, scaffBlock);
									assembled = true;
								} else if(interactions.size()==1) {
									
									if(sstart<=interactions.get(0).upperEndpoint() && send>=interactions.get(0).upperEndpoint()) {
										// 5'-end overlap
										floorRange = interactions.get(0);
										floorInter = subRange.intersection(floorRange);
										scaffBlock = asmb.get(floorRange);
										
										olap  = floorInter.upperEndpoint().intValue()-floorInter.lowerEndpoint().intValue();
										if(olap<=molap) {
											// very small overlap
											// do not clip scaffs
											scaffBlock.ranges.add(Range.closed(1, qlen).canonical(DiscreteDomain.integers()));
											scaffBlock.segments.add(cand);
											scaffBlock.clipAt3End = qlen-qend;
											asmb.remove(floorRange);
											asmb.put(Range.closed(floorRange.lowerEndpoint().intValue(), send).canonical(DiscreteDomain.integers()),
													scaffBlock);
											assembled = true;
										} else {
											clip3 = scaffBlock.clipAt3End;
											clip5 = qstart-1;
											if((double)clip3/qlen<=mclip&&(double)clip5/qlen<=mclip) {
												// 5'-end and 3'-end clip very small
												// clip scaffs and join
												int z = scaffBlock.ranges.size()-1;
												range = scaffBlock.ranges.get(z);
												bseg  = scaffBlock.segments.get(z);
												pos   = getReadPositionAtReferencePosition(bseg, floorInter.upperEndpoint().intValue()-1, true);
												if(range.lowerEndpoint().intValue()>pos)
													pos = getReadPositionAtReferencePosition(bseg, floorInter.upperEndpoint().intValue()-1, false);
												scaffBlock.ranges.set( z, 
														Range.closed(range.lowerEndpoint().intValue(), pos).canonical(DiscreteDomain.integers()) );
												pos   = getReadPositionAtReferencePosition(cand, floorInter.upperEndpoint().intValue()-1, false);
												if(pos>qlen)
													pos = getReadPositionAtReferencePosition(cand, floorInter.upperEndpoint().intValue()-1, true);
												scaffBlock.ranges.add( Range.closed(pos, qlen).canonical(DiscreteDomain.integers()) );
												scaffBlock.segments.add(cand);
												scaffBlock.clipAt3End = qlen-qend;
												asmb.remove(floorRange);
												asmb.put(Range.closed(floorRange.lowerEndpoint().intValue(), send).canonical(DiscreteDomain.integers()),
														scaffBlock);
												assembled = true;
											}
										}
									} else if(sstart<=interactions.get(0).lowerEndpoint() && send>=interactions.get(0).lowerEndpoint()) {
										// 3'-end overlap
										ceilingRange = interactions.get(0);
										ceilingInter = subRange.intersection(ceilingRange);
										scaffBlock = asmb.get(ceilingRange);
										
										olap  = ceilingInter.upperEndpoint().intValue()-ceilingInter.lowerEndpoint().intValue();
										if(olap<=molap) {
											// very small overlap
											// do not clip scaffs
											scaffBlock.ranges.add(0, Range.closed(1, qlen).canonical(DiscreteDomain.integers()));
											scaffBlock.segments.add(0, cand);
											scaffBlock.clipAt5End = qstart-1;
											asmb.remove(ceilingRange);
											asmb.put(Range.closed(sstart, ceilingRange.upperEndpoint().intValue()).canonical(DiscreteDomain.integers()),
													scaffBlock);
											assembled = true;
										} else {
											clip3 = qlen-qend;
											clip5 = scaffBlock.clipAt5End;
											if((double)clip3/qlen<=mclip&&(double)clip5/qlen<=mclip) {
												// 5'-end and 3'-end clip very small
												// clip scaffs and join
												int z = 0;
												range = scaffBlock.ranges.get(z);
												bseg  = scaffBlock.segments.get(z);
												pos   = getReadPositionAtReferencePosition(bseg, ceilingInter.lowerEndpoint().intValue(), false);
												if(pos>range.upperEndpoint().intValue())
													pos = getReadPositionAtReferencePosition(bseg, ceilingInter.lowerEndpoint().intValue(), true);
												scaffBlock.ranges.set( z, 
														Range.closed(pos, range.upperEndpoint().intValue()).canonical(DiscreteDomain.integers()) );
												pos   = getReadPositionAtReferencePosition(cand, ceilingInter.lowerEndpoint().intValue(), true);
												if(1>pos)
													pos   = getReadPositionAtReferencePosition(cand, ceilingInter.lowerEndpoint().intValue(), false);
												scaffBlock.ranges.add( 0, Range.closed(1, pos).canonical(DiscreteDomain.integers()) );
												scaffBlock.clipAt5End = qstart-1;
												scaffBlock.segments.add(0, cand);
												asmb.remove(ceilingRange);
												asmb.put(Range.closed(sstart, ceilingRange.upperEndpoint().intValue()).canonical(DiscreteDomain.integers()), 
														scaffBlock);
												assembled = true;
											}
										}
									}
								} else if(interactions.size()>1) {
									// will fill the gap if assembled
									floorRange   = interactions.get(0);
									ceilingRange = interactions.get(interactions.size()-1);
									floorInter   = subRange.intersection(floorRange);
									ceilingInter = subRange.intersection(ceilingRange);
									
									if(floorInter.isEmpty()||ceilingInter.isEmpty())
										throw new RuntimeException("!!!");
									
									final String[] actions = new String[2];
									
									scaffBlock = asmb.get(floorRange);
									olap  = floorInter.upperEndpoint().intValue()-floorInter.lowerEndpoint().intValue();
									if(olap<=molap) {
										actions[0] = "join";
									} else {
										clip3 = scaffBlock.clipAt3End;
										clip5 = qstart-1;
										if((double)clip3/qlen<=mclip&&(double)clip5/qlen<=mclip) {
											actions[0] = "clip";
										}
									}
									
									scaffBlock = asmb.get(ceilingRange);
									olap  = ceilingInter.upperEndpoint().intValue()-ceilingInter.lowerEndpoint().intValue();
									if(olap<=molap) {
										actions[1] = "join";
									} else {
										clip3 = qlen-qend;
										clip5 = scaffBlock.clipAt5End;
										if((double)clip3/qlen<=mclip&&(double)clip5/qlen<=mclip) {
											actions[1] = "clip";
										}
									}
									
									if(actions[0]!=null&&actions[1]!=null) {
										// will fill the gap
										scaffBlock = asmb.get(floorRange);
										if(actions[0].equals("join")) {
											scaffBlock.ranges.add(Range.closed(1, qlen).canonical(DiscreteDomain.integers()));
											scaffBlock.segments.add(cand);
											scaffBlock.clipAt3End = qlen-qend;
										} else if(actions[0].equals("clip")) {
											int z = scaffBlock.ranges.size()-1;
											range = scaffBlock.ranges.get(z);
											bseg  = scaffBlock.segments.get(z);
											pos   = getReadPositionAtReferencePosition(bseg, floorInter.upperEndpoint().intValue()-1, true);
											if(range.lowerEndpoint().intValue()>pos)
												pos = getReadPositionAtReferencePosition(bseg, floorInter.upperEndpoint().intValue()-1, false);
											scaffBlock.ranges.set( z, 
													Range.closed(range.lowerEndpoint().intValue(), pos).canonical(DiscreteDomain.integers()) );
											pos   = getReadPositionAtReferencePosition(cand, floorInter.upperEndpoint().intValue()-1, false);
											if(pos>qlen)
												getReadPositionAtReferencePosition(cand, floorInter.upperEndpoint().intValue()-1, true);
											scaffBlock.ranges.add(Range.closed(pos, qlen).canonical(DiscreteDomain.integers()));
											scaffBlock.segments.add(cand);
											scaffBlock.clipAt3End = qlen-qend;
										}
										
										nextBlock = asmb.get(ceilingRange);
										if(actions[1].equals("join")) {
											scaffBlock.ranges.addAll(nextBlock.ranges);
											scaffBlock.segments.addAll(nextBlock.segments);
											scaffBlock.clipAt3End = nextBlock.clipAt3End;
										} else if(actions[1].equals("clip")) {
											int z = scaffBlock.ranges.size()-1;
											range = scaffBlock.ranges.get(z);
											bseg  = scaffBlock.segments.get(z);
											pos   = getReadPositionAtReferencePosition(bseg, ceilingInter.lowerEndpoint().intValue(), true); 
											if(range.lowerEndpoint().intValue()>pos)
												pos = getReadPositionAtReferencePosition(bseg, ceilingInter.lowerEndpoint().intValue(), false);
											scaffBlock.ranges.set( z, 
													Range.closed(range.lowerEndpoint().intValue(), pos).canonical(DiscreteDomain.integers()) );
											z = 0;
											range = nextBlock.ranges.get(z);
											bseg  = nextBlock.segments.get(z);
											pos   = getReadPositionAtReferencePosition(bseg, ceilingInter.lowerEndpoint().intValue(), false);
											if(pos>range.upperEndpoint().intValue())
												pos = getReadPositionAtReferencePosition(bseg, ceilingInter.lowerEndpoint().intValue(), true);
											nextBlock.ranges.set( z, 
													Range.closed(pos, range.upperEndpoint().intValue()).canonical(DiscreteDomain.integers()) );
											scaffBlock.ranges.addAll(nextBlock.ranges);
											scaffBlock.segments.addAll(nextBlock.segments);
											scaffBlock.clipAt3End = nextBlock.clipAt3End;
										}
										
										for(Range<Integer> key : interactions) asmb.remove(key); 
										asmb.put(Range.closed(floorRange.lowerEndpoint().intValue(), ceilingRange.upperEndpoint().intValue()).
												canonical(DiscreteDomain.integers()), scaffBlock);
										assembled = true;
									}
								} else {
									throw new RuntimeException("!!!");
								}
								
								
								if(assembled) {
									asmSubCov.add(Range.closed(sstart, send).canonical(DiscreteDomain.integers()));
									myLogger.info("candidate "+cand.getKey()+" assembled.");
								} else {
									myLogger.info("candidate "+cand.getKey()+" discarded.");
								}		
							}
							
							long lenCov = 0, lenScf = 0;
							for(Range<Integer> r : asmSubCov.asRanges()) 
								lenCov += r.upperEndpoint()-r.lowerEndpoint();
							for(Map.Entry<Range<Integer>, ScaffBlock> entry : asmb.entrySet()) {
								for(Range<Integer> r : entry.getValue().ranges)
								lenScf += r.upperEndpoint().intValue()-r.lowerEndpoint().intValue();
							}
							myLogger.info("#STATS "+this.sub_id+" "+iter+"\n"
									+ "Iteration #"+iter+", block size "+b+" done.\n"
									+ "  #assembled segments "+asmb.size()+"\n"
									+ "  #assembled length   "+lenScf+" ("+(double)lenScf/sub_seqs.get(this.sub_id).seq_ln()+")\n"
									+ "  #assembled coverage "+lenCov+" ("+(double)lenCov/sub_seqs.get(this.sub_id).seq_ln()+")\n");
						}
						
						synchronized(lock) {
							int[] pseudo = new int[]{1, 1};
							int gap;
							Range<Integer> prevKey, nextKey;
							ScaffBlock prevBlock, nextBlock;
							final Iterator<Map.Entry<Range<Integer>, ScaffBlock>> iter =  asmb.entrySet().iterator();
							Map.Entry<Range<Integer>, ScaffBlock> entry = iter.hasNext() ? iter.next() : null;
							if(entry!=null) {
								nextKey   = entry.getKey();
								nextBlock = entry.getValue();
								while(iter.hasNext()) {
									prevKey   = nextKey;
									prevBlock = nextBlock;
									
									entry = iter.next();
									nextKey   = entry.getKey();
									nextBlock = entry.getValue();
									
									this.writeScaffBlock(prevBlock, bw_out, pseudo);
									gap = Math.max(mgap, nextKey.lowerEndpoint().intValue()-prevKey.upperEndpoint()+1);
									bw_out.write(this.sub_id+"\t"+pseudo[1]+"\t"+(pseudo[1]+gap-1)+"\t"+pseudo[0]+"\tN\tGAP\t"+gap+"\tyes\t+\n");
									pseudo[0] += 1;
									pseudo[1] += gap;
								}
								
								this.writeScaffBlock(nextBlock, bw_out, pseudo);
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

				private void writeScaffBlock(ScaffBlock scaffBlock, BufferedWriter bw_out, int[] pseudo) {
					// TODO Auto-generated method stub
					try {
						final List<CompoundAlignmentSegment> segments = scaffBlock.segments;
						final List<Range<Integer>> ranges = scaffBlock.ranges;
						CompoundAlignmentSegment segment;
						Range<Integer> range;
						int qstart, qend, qlen;
						String qseqid, orien;
						for(int i=0; i<segments.size(); i++) {
							segment = segments.get(i);
							range = ranges.get(i);
							qstart = range.lowerEndpoint().intValue();
							qend   = range.upperEndpoint().intValue()-1;
							qlen   = qend-qstart+1;
							qseqid = segment.qseqid();
							orien  = "+";
							if(qseqid.endsWith("'")) {
								qseqid = qseqid.replaceAll("'$", "");
								orien  = "-";
							}
							bw_out.write(this.sub_id+"\t"+pseudo[1]+"\t"+(pseudo[1]+qlen-1)+"\t"+pseudo[0]+"\tW\t"+
									qseqid+"\t"+qstart+"\t"+qend+"\t"+orien+"\n");
							pseudo[0] += 1;
							pseudo[1] += qlen;
							
							if(i<segments.size()-1){
								bw_out.write(this.sub_id+"\t"+pseudo[1]+"\t"+(pseudo[1]+mgap-1)+"\t"+pseudo[0]+"\tN\tGAP\t"+mgap+"\tyes\t+\n");
								pseudo[0] += 1;
								pseudo[1] += mgap;
							}
						}
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					return;
				}

				public Runnable init(String sub_id) {
					// TODO Auto-generated method stub
					this.sub_id = sub_id;
					return this;
				}

			}.init(sub));
		}
		this.waitFor();
		
		try {
			bw_out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	protected int getReadPositionAtReferencePosition(CompoundAlignmentSegment as, final int pos, boolean toLeft) {
		// TODO Auto-generated method stub
		final RangeSet<Integer> subCov = TreeRangeSet.create();
		for(SAMRecord record : as.getSamRecords()) {
			subCov.add(Range.closed(record.getAlignmentStart(), 
					record.getAlignmentEnd()).canonical(DiscreteDomain.integers()));
		}
		int pos1 = -1;
		if(!subCov.contains(pos)) {
			int tmp;
			if(toLeft) {
				for(Range<Integer> r : subCov.asRanges()) 
					if( (tmp=r.upperEndpoint().intValue()-1)<pos) 
						pos1 = tmp;
				if(pos1>pos) throw new RuntimeException("!!!");
			} else {
				for(Range<Integer> r : subCov.asDescendingSetOfRanges()) 
					if( (tmp=r.lowerEndpoint().intValue())>pos) 
						pos1 = tmp;
				if(pos1<pos) throw new RuntimeException("!!!");
			}
		} else {
			pos1 = pos;
		}
		
		for(SAMRecord record : as.getSamRecords()) {
			if(record.getAlignmentStart()>pos1||record.getAlignmentEnd()<pos1)
				continue;
			int offset = getHardClipAt5End(record);
			return record.getReadPositionAtReferencePosition(pos1,true)+offset;
		}
		
		throw new RuntimeException("!!!");
	}

	protected boolean contained(ImmutableRangeSet<Integer> cov1, ImmutableRangeSet<Integer> cov2) {
		// TODO Auto-generated method stub
		return cov1.complement().intersection(cov2).isEmpty();
	}
	
	private final Object lock = new Object();
	private long molCount = 0;
	private final static int mext = 30;

	private void parse() {
		// TODO Auto-generated method stub
				final Map<String, Sequence>  sub_seqs = Sequence.parseFastaFileAsMap(subject_file);
		final Map<String, ImmutableRangeSet<Integer>> sub_gaps = new HashMap<>();

		for(Map.Entry<String, Sequence> entry : sub_seqs.entrySet()) {
			String seq_sn = entry.getKey();
			String seq_str = entry.getValue().seq_str();

			final RangeSet<Integer> range_set = TreeRangeSet.create();
			char c;
			for(int j=0; j<seq_str.length(); j++) {
				c = seq_str.charAt(j);
				if(c=='N'||c=='n')
					range_set.add(Range.closed(j, j).
							canonical(DiscreteDomain.integers()));
			}
			sub_gaps.put(seq_sn, ImmutableRangeSet.copyOf(range_set));
		}

		final BAMBarcodeIterator iter = new BAMBarcodeIterator(this.align_file);

		final SAMFileWriter outputSam =  new SAMFileWriterFactory().
				makeSAMOrBAMWriter(iter.samReader.getFileHeader(), true, 
						new File(this.out_prefix+"_parsed.bam") );

		List<SAMRecord> query_recs;

		this.initial_thread_pool();
		while(iter.hasNext()) {
			query_recs = iter.next();
			if(query_recs.isEmpty()) continue;

			this.executor.submit(new Runnable() {
				private List<SAMRecord> query_recs;

				@Override
				public void run() {
					// TODO Auto-generated method stub

					try {
						final List<List<SAMRecord>> blocks = new ArrayList<>();
						List<SAMRecord> block;
						SAMRecord prev, next;
						RangeSet<Integer> qryCov, subCov, subGap, qryExt, subExt;
						int cov, gap, offset, qryExtLen, subExtLen;
						Range<Integer> qrySpan, subSpan;
						
						sortBySubPos(query_recs);
						
						blocks.clear();
						next = query_recs.get(0);
						blocks.add(new ArrayList<>());
						blocks.get(blocks.size()-1).add(next);
						subCov = TreeRangeSet.create();
						qryCov = TreeRangeSet.create();
						subExt = TreeRangeSet.create();
						qryExt = TreeRangeSet.create();
						subGap = TreeRangeSet.create();

						for(int i=1; i<query_recs.size(); i++) {
							prev = next;
							subCov.add(Range.closed(prev.getAlignmentStart(), prev.getAlignmentEnd()).canonical(DiscreteDomain.integers()));
							offset = getHardClipAt5End(prev);
							qryCov.add(Range.closed(prev.getReadPositionAtReferencePosition(prev.getAlignmentStart(),true)+offset, 
									prev.getReadPositionAtReferencePosition(prev.getAlignmentEnd(),true)+offset).canonical(DiscreteDomain.integers()));
							next =  query_recs.get(i);

							if(prev.getReferenceIndex().intValue()==next.getReferenceIndex().intValue()) {
								offset = getHardClipAt5End(next);
								qrySpan = Range.closed(next.getReadPositionAtReferencePosition(next.getAlignmentStart(),true)+offset, 
										next.getReadPositionAtReferencePosition(next.getAlignmentEnd(),true)+offset).canonical(DiscreteDomain.integers());
								
								qryExt.clear();
								qryExt.add(qrySpan);
								qryExtLen = getCovFromCovs(ImmutableRangeSet.copyOf(qryExt).complement().intersection(qryCov));

								subSpan = Range.closed(next.getAlignmentStart(), next.getAlignmentEnd()).canonical(DiscreteDomain.integers());
								subExt.clear();
								subExt.add(subSpan);
								subExtLen = getCovFromCovs(ImmutableRangeSet.copyOf(subExt).complement().intersection(subCov));

								if(qryExtLen>=mext&&subExtLen>=mext&&extendIsInRightOrder(qryCov, qrySpan)) {
									gap = next.getAlignmentStart()-prev.getAlignmentEnd();

									if(gap>0) {
										cov = Math.min(getCovFromCovs(subCov), next.getAlignmentEnd()-next.getAlignmentStart());
										subGap.clear();
										subGap.add(Range.closed(prev.getAlignmentEnd(), next.getAlignmentStart()).canonical(DiscreteDomain.integers()));
										gap -= getCovFromCovs(sub_gaps.get(prev.getReferenceName()).intersection(subGap));
										if(gap>cov*mxgr) {
											mergeBlocks(blocks, sub_gaps);
											blocks.add(new ArrayList<>());
											subCov.clear();
											qryCov.clear();
										}
									}
								} else {
									mergeBlocks(blocks, sub_gaps);
									blocks.add(new ArrayList<>());
									subCov.clear();
									qryCov.clear();
								}
							} else {
								mergeBlocks(blocks, sub_gaps);
								blocks.add(new ArrayList<>());
								subCov.clear();
								qryCov.clear();
							}
							blocks.get(blocks.size()-1).add(next);
						}

						final double[] covs = new double[blocks.size()];
						final List<RangeSet<Integer>> qryCovs = new ArrayList<>();
						double ext = Double.MIN_VALUE;
						for(int i=0; i<blocks.size(); i++) {
							block = blocks.get(i);
							qryCov = getQryCovs(block);
							cov = getCovFromCovs(qryCov);
							covs[i] = cov;
							qryCovs.add(qryCov);
							if(cov>ext) {
								ext = cov;
							}
						}

						ext = ext*(1-diff);
						synchronized(lock) {
							for(int i=0; i<blocks.size(); i++) {
								if(covs[i]<ext) continue;
								block = blocks.get(i);
								++molCount;
								for(SAMRecord rec : block) {
									rec.setAttribute("BX", "b"+String.format("%08d", molCount));
									outputSam.addAlignment(rec);
								}
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

				public Runnable init(List<SAMRecord> query_recs) {
					// TODO Auto-generated method stub
					this.query_recs = query_recs;
					return this;
				}

			}.init(query_recs));
		}

		this.waitFor();
		iter.close();
		outputSam.close();
	}

	protected boolean extendIsInRightOrder(RangeSet<Integer> cov, Range<Integer> ext) {
		// TODO Auto-generated method stub
		// need to be a right side extension
		Range<Integer> span = cov.span();
		return span.lowerEndpoint().intValue()<=ext.lowerEndpoint().intValue()&&
				span.upperEndpoint().intValue()<=ext.upperEndpoint().intValue();
	}

	private void mergeBlocks(List<List<SAMRecord>> blocks, Map<String, ImmutableRangeSet<Integer>> sub_gaps) {
		// TODO Auto-generated method stub
		List<SAMRecord> prev, next;
		final int n = blocks.size();
		RangeSet<Integer> prevSubCov, nextSubCov, prevQryCov, nextQryCov, subGap = TreeRangeSet.create();
		Range<Integer> prevSubSpan, nextSubSpan, nextQrySpan;
		int cov, gap, qryExtLen, subExtLen;

		outerloop:
		for(int i=n-1; i>0; ) {
			next = blocks.get(i);
			nextSubCov = getSubCovs(next);
			nextSubSpan = nextSubCov.span();
			nextQryCov = getQryCovs(next);
			nextQrySpan = nextQryCov.span();
			
			for(int j=i-1; j>=0; j--) {
				prev = blocks.get(j);
				if(prev.get(0).getReferenceIndex().intValue()!=
						next.get(0).getReferenceIndex().intValue())
					return;
				prevSubCov = getSubCovs(prev);
				prevSubSpan = prevSubCov.span();
				prevQryCov = getQryCovs(prev);

				subExtLen = getCovFromCovs(ImmutableRangeSet.copyOf(prevSubCov).complement().intersection(nextSubCov));
				qryExtLen = getCovFromCovs(ImmutableRangeSet.copyOf(prevQryCov).complement().intersection(nextQryCov));
				
				if(qryExtLen>=mext&&subExtLen>=mext&&extendIsInRightOrder(prevQryCov, nextQrySpan)) {
					
					gap = nextSubSpan.lowerEndpoint()-prevSubSpan.upperEndpoint();
					if(gap<0) {
						for(int k=i; k>j; k--) {
							prev.addAll(blocks.get(k));
							blocks.remove(k);
						}
						sortBySubPos(prev);
						i = j;
						continue outerloop;
					} else {
						cov = Math.min(getCovFromCovs(prevSubCov), getCovFromCovs(nextSubCov));
						subGap.clear();
						subGap.add(Range.closed(prevSubSpan.upperEndpoint(), nextSubSpan.lowerEndpoint()).canonical(DiscreteDomain.integers()));
						gap -= getCovFromCovs(sub_gaps.get(prev.get(0).getReferenceName()).intersection(subGap));
						if(gap<=cov*mxgr) {
							for(int k=i; k>j; k--) {
								prev.addAll(blocks.get(k));
								blocks.remove(k);
							}
							sortBySubPos(prev);
							i = j;
							continue outerloop;
						}
					}
				}
			}
			return;
		}
	}
	
	protected void sortBySubPos(List<SAMRecord> recs) {
		// TODO Auto-generated method stub
		Collections.sort(recs, new Comparator<SAMRecord>() {

			@Override
			public int compare(SAMRecord o1, SAMRecord o2) {
				// TODO Auto-generated method stub
				int r1 = o1.getReferenceIndex().intValue(),
						r2 = o2.getReferenceIndex().intValue();
				return r1==r2 ? o1.getAlignmentStart()-o2.getAlignmentStart() : r1-r2;
			}

		});

	}

	private boolean reverse(List<SAMRecord> sams) {
		// TODO Auto-generated method stub
		RangeSet<Integer> forward = TreeRangeSet.create();
		RangeSet<Integer> reverse = TreeRangeSet.create();
		for(SAMRecord sam : sams) {
			int offset = getHardClipAt5End(sam);
			int lower = sam.getReadPositionAtReferencePosition(sam.getAlignmentStart(), true);
			int upper = sam.getReadPositionAtReferencePosition(sam.getAlignmentEnd(),   true);
			if(sam.getReadNegativeStrandFlag()) {
				reverse.add(Range.closed(lower+offset, upper+offset).canonical(DiscreteDomain.integers()));
			} else {
				forward.add(Range.closed(lower+offset, upper+offset).canonical(DiscreteDomain.integers()));
			}
		}

		return getCovFromCovs(forward)<getCovFromCovs(reverse);
	}

	private int getHardClipAt5End(SAMRecord sam) {
		// TODO Auto-generated method stub
		CigarElement cigar = sam.getCigar().getCigarElement(0);
		int clip = 0;
		if(cigar.getOperator()==CigarOperator.HARD_CLIP)
			clip += cigar.getLength();
		return clip;
	}

	private int getCovFromCovs(RangeSet<Integer> rs) {
		// TODO Auto-generated method stub
		int cov = 0;
		for(Range<Integer> r : rs.asRanges())
			cov += r.upperEndpoint()-r.lowerEndpoint();
		return cov;
	}

	private int getSpanFromCovs(RangeSet<Integer> rs) {
		// TODO Auto-generated method stub
		Range<Integer> span = rs.span();
		return span.upperEndpoint()-span.lowerEndpoint();
	}

	private RangeSet<Integer> getQryCovs(List<SAMRecord> recs) {
		// TODO Auto-generated method stub
		RangeSet<Integer> covs = TreeRangeSet.create();
		for(SAMRecord rec : recs) {
			int offset = getHardClipAt5End(rec);
			int lower = rec.getReadPositionAtReferencePosition(rec.getAlignmentStart(), true);
			int upper = rec.getReadPositionAtReferencePosition(rec.getAlignmentEnd(),   true);
			covs.add(Range.closed(lower+offset, upper+offset).canonical(DiscreteDomain.integers()));
		}
		return covs;
	}

	private RangeSet<Integer> getSubCovs(List<SAMRecord> recs) {
		// TODO Auto-generated method stub
		RangeSet<Integer> covs = TreeRangeSet.create();
		for(SAMRecord rec : recs) 
			covs.add(Range.closed(rec.getAlignmentStart(), rec.getAlignmentEnd())
					.canonical(DiscreteDomain.integers()));
		return covs;
	}

	private int distance(RangeSet<Integer> span, Range<Integer> ext) {
		// TODO Auto-generated method stub

		Range<Integer> s = span.span();
		int ls = s.lowerEndpoint(), us = s.upperEndpoint(),
				le = ext.lowerEndpoint(), ue = ext.upperEndpoint();
		if(ue<=ls) 
			return ls-ue;
		else if(us<=le)
			return le-us;
		return 0;
	}
	
	private class ScaffBlock {
		// this field is to maintain the segments
		final List<CompoundAlignmentSegment> segments = new ArrayList<>();
		// these field are to maintain 
		// the ranges on the query sequences
		final List<Range<Integer>> ranges = new ArrayList<>();
		// clips of left end of the left-most segment
		int clipAt5End;
		// clips of the right end of the right-most segment
		int clipAt3End;
		
		public void setClipAt5End(int clip) {
			this.clipAt5End = clip;
		}
		
		public void setClipAt3End(int clip) {
			this.clipAt3End = clip;
		}
	}
}














