package cz1.hmm.tools;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleDirectedWeightedGraph;

import cz1.hmm.model.ModelReader;
import cz1.util.ArgsEngine;
import cz1.util.Utils;

public class SinglePointAnalysis extends RFUtils {

	private final static Logger myLogger = LogManager.getLogger(SinglePointAnalysis.class);
	
	private String out_prefix;
	private int wbp = 30000;
	private int wnm = 30;
	private String map_func = "kosambi";
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--hap-file               Directory with input haplotype files.\n"
						+ " -o/--prefix                 Output file prefix.\n"
						+ " -wbp/-windows-bp            Window size (#basepairs) (default 30000).\n"
						+ " -wnm/-windows-nm            Window size (#markers) (default 30).\n"
						+ " -ex/--experiment-id         Common prefix of haplotype files for this experiment.\n"
						+ " -nb/--best                  The most likely nb haplotypes will be used (default 10).\n"
						+ " -phi/--skew-phi             For a haplotype inference, the frequencies of parental \n"
						+ "                             haplotypes need to be in the interval [1/phi, phi], \n"
						+ "                             otherwise will be discared (default 2).\n"
						+ " -nd/--drop                  At least nd haplotype inferences are required for \n"
						+ "                             a contig/scaffold to be analysed (default 1).\n"
						+ " -t/--threads                #threads (default 1).\n"	
				);
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
			myArgsEngine.add( "-i", "--hap-file", true);
			myArgsEngine.add( "-o", "--prefix", true);
			myArgsEngine.add( "-wbp", "--windows-bp", true);
			myArgsEngine.add( "-wnm", "--windows-nm", true);
			myArgsEngine.add( "-ex", "--experiment-id", true);
			myArgsEngine.add( "-nb", "--best", true);
			myArgsEngine.add( "-phi", "--skew-phi", true);
			myArgsEngine.add( "-nd", "--drop", true);
			myArgsEngine.add( "-t", "--threads", true);
			myArgsEngine.parse(args);
		}
		
		if(myArgsEngine.getBoolean("-i")) {
			in_haps = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your input zip file.");
		}

		if(myArgsEngine.getBoolean("-o")) {
			out_prefix = myArgsEngine.getString("-o");
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify your output file prefix.");
		}
		
		if(myArgsEngine.getBoolean("-ex")) {
			expr_id = myArgsEngine.getString("-ex");
		}  else {
			expr_id = guessExperimentId();
			myLogger.warn("No experiment prefix provided, I guess it's "+expr_id+". Please\n"
					+ "specify it with -ex/--experiment-id option if it's incorrect.");
		}
		
		best_n = 10;
		if(myArgsEngine.getBoolean("-nb")) {
			best_n = Integer.parseInt(myArgsEngine.getString("-nb"));
		}
		
		if(myArgsEngine.getBoolean("-wbp")) {
			wbp = Integer.parseInt(myArgsEngine.getString("-wbp"));
		}
		
		if(myArgsEngine.getBoolean("-wnm")) {
			wnm = Integer.parseInt(myArgsEngine.getString("-wnm"));
		}
		
		if(myArgsEngine.getBoolean("-t")) {
			THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}
		
		if(myArgsEngine.getBoolean("-phi")) {
			skew_phi = Integer.parseInt(myArgsEngine.getString("-phi"));
		}
		
		if(myArgsEngine.getBoolean("-nd")) {
			drop_thres = Integer.parseInt(myArgsEngine.getString("-nd"));
		}
	}
	
	public SinglePointAnalysis (String in_haps, 
			String out_prefix,
			String expr_id, 
			int threads,
			double skew_phi,
			int drop_thres,
			int best_n) { 
		this.in_haps = in_haps;
		this.out_prefix = out_prefix;
		this.expr_id = expr_id;
		THREADS = threads;
		this.skew_phi = skew_phi;
		this.drop_thres = drop_thres;
		this.best_n = best_n;
	}
	
	public SinglePointAnalysis() {
		// TODO Auto-generated constructor stub
		super();
	}

	BufferedWriter rfWriter;
	@Override
	public void run() {
		// TODO Auto-generated method stub
		super.initialise();
		
		rfWriter = Utils.getBufferedWriter(this.out_prefix+".map");
		this.initial_thread_pool();
		for(String scaff : fileObj.keySet()) 
			executor.submit(new RfCalculator(scaff));
		this.waitFor();
		try {
			rfWriter.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		//render();
		
		myLogger.info("["+Utils.getSystemTime()+"] DONE.");
	}
	
	private class RfCalculator implements Runnable {
		private final String scaff;
		
		public RfCalculator(String scaff) {
			// TODO Auto-generated constructor stub
			this.scaff = scaff;
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub
			try {
				List<FileObject> objs = fileObj.get(scaff);
				final int n = objs.size();
				ModelReader modelReader = new ModelReader(objs.get(0).file);
				int[] pos = modelReader.getSnpPosition(scaff);
				modelReader.close();
				int m = pos.length;
				int[] dists = new int[m-1];
				for(int i=0; i<dists.length; i++) dists[i] = pos[i+1]-pos[i];
				
				final List<Map<String, char[][]>> haplotypes = new ArrayList<>();
				final int[] hapn = new int[n];
				for(int i=0; i<n; i++) {
					FileObject obj = objs.get(i);
					modelReader = new ModelReader(obj.file);
					Map<String, char[][]> haps = modelReader.getHaplotypeByPositionRange(obj.position, ploidy);
					modelReader.close();
					for(String f : parents) if(f!=null) haps.remove(f);
					for(char[][] hap : haps.values()) 
						if(hap[0][0]!='*') ++hapn[i];
					hapn[i] *= ploidy;
					haplotypes.add(haps);
				}
				
				final SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> jumpGraph = new SimpleDirectedWeightedGraph<>(DefaultWeightedEdge.class);
				for(int j=0; j<m; j++) jumpGraph.addVertex(j);
				int dist, kb, extn;
				double[] jump;
				double[][] jumps;
				char[] h;
				Map<String, char[][]> haps;
				DefaultWeightedEdge edge;
				for(int j=0; j<m-1; j++) {
					kb = Math.min(m-1, j+wnm);
					dist = 0;
					for(int k=j; k<kb&&dist<wbp; k++)
						dist += dists[k];
					while(dist<wbp&&kb<m-1) {
						dist += dists[kb];
						++kb;
					}
					
					extn=kb-j;
					jumps = new double[n][];
					for(int i=0; i<n; i++) {
						jump = new double[extn];
						haps = haplotypes.get(i);
						for(char[][] hap : haps.values()) {
							if(hap[0][0]=='*') continue;
							
							for(int p=0; p<ploidy; p++) {
								h = hap[p];
								for(int k=j+1; k<=kb; k++) {
									if(h[j]!=h[k]) 
										++jump[k-j-1];
								}
							}
						}
						for(int k=0; k<extn; k++) jump[k]/=hapn[i];
						jumps[i] = jump;
					}
					double[] means = new double[extn];
					for(int i=0; i<n; i++)
						sum(means, jumps[i], means);
					for(int i=0; i<extn; i++) means[i] /= n;
					
					for(int k=j+1; k<=kb; k++) {
						edge = jumpGraph.addEdge(j, k);
						jumpGraph.setEdgeWeight(edge, RFUtils.geneticDistance(means[k-j-1], map_func));
					}
				}

				/***
				myLogger.info("Jump graph (before pruning): ");
				myLogger.info("  #V, "+jumpGraph.vertexSet().size());
				myLogger.info("  #E, "+jumpGraph.edgeSet().size());
				double rf_thresh = 0.1;
				Set<DefaultWeightedEdge> edges2remove = new HashSet<>();
				Set<DefaultWeightedEdge> edges2keep = new HashSet<>();
				Set<DefaultWeightedEdge> edgeSet;
				int z;
				double r, w;
				for(int vertex : jumpGraph.vertexSet()) {
					for(int k=0; k<2; k++) {
						edgeSet = k==0 ? jumpGraph.incomingEdgesOf(vertex) : jumpGraph.outgoingEdgesOf(vertex);
						z = 0;
						edge = null;
						w = Double.MAX_VALUE;
						for(DefaultWeightedEdge e : edgeSet) {
							if( (r=jumpGraph.getEdgeWeight(e))>rf_thresh) {
								edges2remove.add(e);
								++z;
								if(r<w) {
									w = r;
									edge = e;
								}
							}
						}
						if(z==edgeSet.size()) {
							// means no edge left
							// keep at least one edge
							edges2keep.add(edge);
						}
					}
				}
				edges2remove.removeAll(edges2keep);
				jumpGraph.removeAllEdges(edges2remove);
				
				myLogger.info("Jump graph (after pruning): ");
				myLogger.info("  #V, "+jumpGraph.vertexSet().size());
				myLogger.info("  #E, "+jumpGraph.edgeSet().size());
				***/
				
				List<DefaultWeightedEdge> path = findPathBetween(jumpGraph, 0, m-1);
				final int p = path.size();
				double[] rs = new double[p];
				int[] ds = new int[p];
				int source, target;
				for(int i=0; i<p; i++) {
					edge = path.get(i);
					source = jumpGraph.getEdgeSource(edge);
					target = jumpGraph.getEdgeTarget(edge);
					int d = 0;
					for(int j=source; j<target; j++) d += dists[j];
					ds[i] = d;
					rs[i] = RFUtils.inverseGeneticDistance(jumpGraph.getEdgeWeight(edge), map_func);
				}
				
				if(sum(ds)!=sum(dists)) throw new RuntimeException("!!!");
				
				synchronized(lock) {
					rfWriter.write("C "+scaff+"\n");
					rfWriter.write("D "+Utils.paste(ds, ",")+"\n");
					rfWriter.write(Utils.paste(rs, ",")+"\n");
				}
				
			} catch (Exception e) {
				Thread t = Thread.currentThread();
				t.getUncaughtExceptionHandler().uncaughtException(t, e);
				e.printStackTrace();
				executor.shutdown();
				System.exit(1);
			}
		}

		private List<DefaultWeightedEdge> findPathBetween(SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> graph, int source, int target) {
			// TODO Auto-generated method stub
			// graph is topological sorted directed acyclic graph
			final int m = graph.vertexSet().size();
			double[] dist = new double[m];
			int[] predecessor = new int[m];
			for(int v=0; v<m; v++) {
				dist[v] = Double.POSITIVE_INFINITY;
				predecessor[v] = -1;
			}
			dist[source] = 0;
			int u;
			double d;
			for(int v=0; v<m; v++) {
				for(DefaultWeightedEdge e : graph.outgoingEdgesOf(v)) {
					u = graph.getEdgeTarget(e);
					if((d=dist[v]+graph.getEdgeWeight(e))<dist[u]) {
						dist[u] = d;
						predecessor[u] = v;
					}
				}
			}
			List<DefaultWeightedEdge> path = new ArrayList<>();
			u = target;
			while(u!=source) {
				if(predecessor[u]==-1) return null;
				path.add(graph.getEdge(predecessor[u], u));
				u = predecessor[u];
			}
			
			Collections.reverse(path);
			return path;
		}
	}
}
