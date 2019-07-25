package cz1.hmm.tools;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.SimpleDirectedWeightedGraph;

import cz1.hmm.model.ModelReader;
import cz1.util.ArgsEngine;
import cz1.util.Utils;

public class AssemblyError extends RFUtils {

	private final static Logger myLogger = Logger.getLogger(AssemblyError.class);
	
	private String out_prefix;
	private double rf_thresh = 0.1;
	private int wbp = 30000;
	private int wnm = 30;
	private String map_func = "haldane";
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--hap-file               Directory with input haplotype files.\n"
						+ " -o/--prefix                 Output file prefix.\n"
						+ " -r/--rf-thresh              Recombination frequency threshold for assembly error detection (default 0.1).\n"
						+ " -wbp                        Window size (#basepairs) for type II assembly error dectection (default 30000).\n"
						+ " -wnm                        Window size (#markers) for type II assembly error dectection (default 30).\n"
						+ " -ex/--experiment-id         Common prefix of haplotype files for this experiment.\n"
						+ " -nb/--best                  The most likely nb haplotypes will be used (default 10).\n"
						+ " -phi/--skew-phi             For a haplotype inference, the frequencies of parental \n"
						+ "                             haplotypes need to be in the interval [1/phi, phi], \n"
						+ "                             otherwise will be discared (default 2).\n"
						+ " -nd/--drop                  At least nd haplotype inferences are required for \n"
						+ "                             a contig/scaffold to be analysed (default 1).\n"
						+ " -t/--threads                Threads (default 1).\n"	
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
			myArgsEngine.add( "-ex", "--experiment-id", true);
			myArgsEngine.add( "-i", "--hap-file", true);
			myArgsEngine.add( "-o", "--prefix", true);
			myArgsEngine.add( "-nb", "--best", true);
			myArgsEngine.add( "-t", "--threads", true);
			myArgsEngine.add( "-phi", "--skew-phi", true);
			myArgsEngine.add( "-nd", "--drop", true);
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
		
		best_n = 30; // use as much as possible up to 30
		if(myArgsEngine.getBoolean("-nb")) {
			best_n = Integer.parseInt(myArgsEngine.getString("-nb"));
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
	
	public AssemblyError (String in_haps, 
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
	
	public AssemblyError() {
		// TODO Auto-generated constructor stub
		super();
	}

	BufferedWriter rfWriter;
	@Override
	public void run() {
		// TODO Auto-generated method stub
		super.initialise();
		
		rfWriter = Utils.getBufferedWriter(this.out_prefix+".err");
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
				int[] distance = modelReader.getDistance();
				modelReader.close();
				int lb = objs.get(0).position[0], 
						ub = objs.get(0).position[1];
				final int m = ub-lb+1;
				final int deno = n*nF1*ploidy;
				int[] position = new int[m];
				for(int i=0; i<m; i++) position[i] = lb+i;
				int[] dists = new int[m-1];
				System.arraycopy(distance, lb, dists, 0, m-1);
				
				final List<Map<String, char[][]>> haplotypes = new ArrayList<>();
				for(int i=0; i<n; i++) {
					FileObject obj = objs.get(i);
					modelReader = new ModelReader(obj.file);
					Map<String, char[][]> haps = modelReader.getHaplotypeByPosition(position, ploidy);
					modelReader.close();
					for(String f : parents) haps.remove(f);
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
							for(int p=0; p<ploidy; p++) {
								h = hap[p];
								for(int k=j+1; k<=kb; k++) {
									if(h[j]!=h[k]) 
										++jump[k-j-1];
								}
							}
						}
						jumps[i] = jump;
					}
					double[] means = new double[extn];
					for(int i=0; i<n; i++)
						sum(means, jumps[i], means);
					for(int i=0; i<extn; i++) means[i] /= deno;
					
					for(int k=j+1; k<=kb; k++) {
						edge = jumpGraph.addEdge(j, k);
						jumpGraph.setEdgeWeight(edge, RFUtils.geneticDistance(means[k-j-1], map_func));
					}
				}
				
				myLogger.info("Jump graph (before pruning): ");
				myLogger.info("  #V, "+jumpGraph.vertexSet().size());
				myLogger.info("  #E, "+jumpGraph.edgeSet().size());
				
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
				
				List<DefaultWeightedEdge> path = findPathBetween(jumpGraph, 0, m-1);
				
				myLogger.info(scaff+" ("+0+"->"+(m-1)+"): "+(path.size()+1)+", "+getWeight(jumpGraph, path)*100+"cm");
				
				final int p = path.size();
				double[] rs = new double[p];
				for(int i=0; i<p; i++) {
					rs[i] = RFUtils.inverseGeneticDistance(jumpGraph.getEdgeWeight(path.get(i)), map_func);
					myLogger.info(path.get(i).toString()+" "+rs[i]);
				}
				
				// now detect assembly errors
				for(int i=0; i<p; i++) {
					if(rs[i]>rf_thresh) {
						// assembly error
						synchronized(lock) {
							edge = path.get(i);
							rfWriter.write(scaff+" "+edge.toString()+" "+rs[i]+"\n");
						}
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

		private double getWeight(SimpleDirectedWeightedGraph<Integer, DefaultWeightedEdge> graph,
				List<DefaultWeightedEdge> path) {
			// TODO Auto-generated method stub
			double w = 0;
			for(DefaultWeightedEdge e : path) w += graph.getEdgeWeight(e);
			return w;
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

	public Set<String> split(String string, String string2) {
		// TODO Auto-generated method stub
		return null;
	}

	public Map<String, int[][]> errs() {
		// TODO Auto-generated method stub
		return null;
	}
	
	/***
	private void render() {
		// TODO Auto-generated method stub
		for(int i=0; i<dc.length; i++) {
			if(dc[i][0]==null) continue;
			String scaff = dc[i][0].markers[0].replaceAll("_[0-9]{1,}$", "");
			double[][] rfs = mapCalc.get(scaff);
			double[] rf = new double[rfs[0].length];
			int bound = 0;
			for(bound=0; bound<rfs.length; bound++) 
				if(rfs[bound]==null) break;
			for(int j=0; j<rf.length; j++) 
				for(int k=0; k<bound; k++)
					rf[j] += rfs[k][j];
			String[] markers = dc[i][0].markers;
			final List<int[]> breakage_pos = new ArrayList<int[]>();
			for(int j=0; j<rf.length; j++) {
 				rf[j] /= bound;
 				if(rf[j]>=breakage_thres) {
 					String[] s = markers[j].split("_");
 					int x = Integer.parseInt(s[s.length-1]);
 					s = markers[j+1].split("_");
 					int x2 = Integer.parseInt(s[s.length-1]);
 					breakage_pos.add(new int[]{x, x2}); 
 				}
 			}
			if(breakage_pos.size()==0) continue;
			int[][] ls = new int[breakage_pos.size()][2];
			for(int j=0; j<ls.length; j++) ls[j] = breakage_pos.get(j);
			errs.put(scaff, ls);
		}
		return;
	}
	
	public Set<String> split(String in_vcf, String out_vcf) {
		// TODO Auto-generated method stub
		final Set<String> scaff_breakge = new HashSet<String>();
		try {
			BufferedReader br = Utils.getBufferedReader(in_vcf);
			BufferedWriter bw = Utils.getBufferedWriter(out_vcf);
			 
			String[] s;
			String line = br.readLine();
			while( line!=null ) {
				if(line.startsWith("#")) {
					bw.write(line+"\n");
					line = br.readLine();
					continue;
				}
				s = line.split("\\s+");
				if(!this.errs.containsKey(s[0])) {
					bw.write(line+"\n");
					line = br.readLine();
				} else {
					String scaff = s[0];
					int[][] tmp = this.errs.get(scaff);
					double[] breakage_pos = new double[tmp.length+1];
					for(int i=0; i<tmp.length; i++) 
						breakage_pos[i] = (tmp[i][0]+tmp[i][1])/2.0;
					breakage_pos[tmp.length] = Double.POSITIVE_INFINITY;
					int sub = 1;
					scaff_breakge.add(scaff+"_"+1);
					while( line!=null ) {
						s = line.split("\\s+");
						if( !scaff.equals(s[0]) ) break;
						if( Double.parseDouble(s[1])>breakage_pos[sub-1] )
							scaff_breakge.add(scaff+"_"+(++sub));
						bw.write(line.replaceAll("^"+scaff, scaff+"_"+sub)+"\n");
						line = br.readLine();
					}
				}
			}
			br.close();
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return scaff_breakge;
	}
	
	**/
}
