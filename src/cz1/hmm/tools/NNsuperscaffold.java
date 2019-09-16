package cz1.hmm.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.log4j.Logger;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

public class NNsuperscaffold extends Executor {

	private final static Logger myLogger = Logger.getLogger(NNsuperscaffold.class);
	
	private String rf_file = null;
	private String RLibPath = null;
	private double lod_thres = 3;
	private int hs;
	private int nn = 2;
	private String out_prefix = null;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--rf                     Recombination frequency file.\n"
						+ " -h/--hap-size               #haplotypes (popsize*ploidy).\n"
						+ " -l/--lod                    LOD score threshold (default: 3).\n"
						+ " -n/--neighbour              #nearest neighbours (default: 2).\n"
						+ " -rlib/--R-external-libs     R external library path.\n"
						+ " -o/--prefix                 Output file prefix.\n\n"
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
			myArgsEngine.add("-i", "--rf", true);
			myArgsEngine.add("-h", "--hap-size", true);
			myArgsEngine.add("-l", "--lod", true);
			myArgsEngine.add("-n", "--neighbour", true);
			myArgsEngine.add("-o", "--prefix", true);
			myArgsEngine.add("-rlib", "--R-external-libs", true);
			myArgsEngine.parse(args);
		}
		
		if(myArgsEngine.getBoolean("-i")) {
			rf_file = myArgsEngine.getString("-i");
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify your recombinatio frequency file.");
		}
		
		if(myArgsEngine.getBoolean("-h")) {
			hs = Integer.parseInt(myArgsEngine.getString("-h"));
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the number of haplotypes (popsize*ploidy).");
		}

		if(myArgsEngine.getBoolean("-l")) {
			lod_thres = Double.parseDouble(myArgsEngine.getString("-l"));
		}
		
		if(myArgsEngine.getBoolean("-n")) {
			nn = Integer.parseInt(myArgsEngine.getString("-n"));
		}
		
		if(myArgsEngine.getBoolean("-o")) {
			out_prefix = myArgsEngine.getString("-o");
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify your output file prefix.");
		}

		if (myArgsEngine.getBoolean("-rlib")) {
			RLibPath = myArgsEngine.getString("-rlib");
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		final String temfile_prefix = Utils.makeTempDir();
		final boolean tmpdirCreated = Utils.makeOutputDir(new File(temfile_prefix));
		final String concorde_path =
				RFUtils.makeExecutable("cz1/hmm/executable/concorde", temfile_prefix);
		new File(concorde_path).setExecutable(true, true);
		final String nnss_path =
				RFUtils.makeExecutable("cz1/hmm/scripts/make_nnsuperscaffold.R", temfile_prefix);
		RFUtils.makeExecutable("cz1/hmm/scripts/include.R", temfile_prefix);
		RFUtils.makeRMatrix(rf_file, out_prefix+".RData", hs);
		double max_r = Math.min(RFUtils.calcRfFromLOD(lod_thres, hs), RFUtils.inverseGeneticDistance(0.5, "kosambi"));
		myLogger.info("Using recombination frequency threshold: "+max_r+".");
		final String command =
				"Rscript "+nnss_path+" "
						+ "-i "+out_prefix+".RData "
						+ "-r "+max_r+" "
						+ "-o "+out_prefix+".nns "
						+ "-n "+nn+" "
						+ "--concorde "+new File(concorde_path).getParent()+" "
						+ (RLibPath==null ? "" : "--include "+RLibPath+" ")
						+ "--tmpdir "+new File(temfile_prefix).getAbsolutePath();
		this.consume(this.bash(command));

		if(tmpdirCreated) Utils.deleteDirectory(new File(temfile_prefix));
	}
	
	public void nj(double max_r) {
		// TODO Auto-generated method stub

		final BidiMap<String, Integer> scaffs = new DualHashBidiMap<>();
		final BidiMap<Cluster, Integer> clusts = new DualHashBidiMap<>();
		final TreeMap<Double, Set<ClustPair>> minRfs = new TreeMap<>();
		try {
			BufferedReader br1 = Utils.getBufferedReader(rf_file);
			int c = 0;
			String line, scaff;
			String[] s;
			double minf;
			double[] allf;
			Cluster c1, c2;
			int i1, i2;
			ClustPair pair;
			long line_count = 0;
			while((line=br1.readLine())!=null) {
				if(line.startsWith("#")) {
					scaff = line.trim().replaceAll("^##", "");
					scaffs.put(scaff, c);
					c1 = new Cluster(c);
					clusts.put(c1, c);
					++c;
				} else {
					++line_count;
					if(line_count%1000000==0) myLogger.info("#lines loaded: "+line_count);
					s = line.trim().split("\\s+");
					minf = Double.parseDouble(s[0]);
					i1 = scaffs.get(s[5]);
					i2 = scaffs.get(s[6]);
					allf = new double[4];
					allf[0] = Double.parseDouble(s[1]);
					allf[3] = Double.parseDouble(s[4]);
					if(i1>i2) {
						c1 = clusts.getKey(i1);
						c2 = clusts.getKey(i2);
						allf[1] = Double.parseDouble(s[2]);
						allf[2] = Double.parseDouble(s[3]);
					} else {
						c1 = clusts.getKey(i2);
						c2 = clusts.getKey(i1);
						allf[1] = Double.parseDouble(s[3]);
						allf[2] = Double.parseDouble(s[2]);
					}
					c1.minf.put(c2, minf);
					c1.allf.put(c2, allf);
					pair = new ClustPair(c1, c2);
					if(!minRfs.containsKey(minf))
						minRfs.put(minf, new HashSet<>());
					minRfs.get(minf).add(pair);
				}
			}
			myLogger.info("####total lines loaded: "+line_count);
			br1.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		int nclust = clusts.size();
		Map.Entry<Double, Set<ClustPair>> entry;
		double f, extf;
		Set<ClustPair> pairs;
		ClustPair pair;
		Cluster c1, c2, c, clust;
		boolean j1, j2;
		Map<Cluster, Double> minf;
		Map<Cluster, double[]> allf;
		int i1, i2, x1, x2;
		double[] fs, fz;
		Set<Double> fs2cl = new HashSet<>();

		while(!minRfs.isEmpty()) {

			entry = minRfs.firstEntry();
			f = entry.getKey();
			pairs = entry.getValue();
			
			if(f>max_r) break;
			myLogger.info("####"+nclust+" "+f+" "+clusts.size());

			// all pairs have same recombination frequencies
			// need to decide which pair to join
			// before that need to decide the orientation of join
			pair = null;
			extf = Double.MAX_VALUE;
			for(ClustPair p : pairs) {
				p.join();
				if(p.extf<extf) 
					pair = p;
			}
			if(pair==null) pair = pairs.iterator().next();

			// now join the pair
			c1 = pair.c1;
			c2 = pair.c2;
			j1 = pair.j1;
			j2 = pair.j2;
			i1 = clusts.get(c1);
			i2 = clusts.get(c2);

			c = new Cluster(i1, i2, c1, c2, j1, j2, f);
			allf = c.allf;
			minf = c.minf;

			// now update clusts
			clusts.removeValue(i1);
			clusts.removeValue(i2);

			// now update allf
			x1 = j1?2:0;
			x2 = j2?2:0;
			for(int i : clusts.values()) {
				clust = clusts.getKey(i);
				fs = new double[4];
				fz = i<i1 ? c1.allf.get(clust) : clust.allf.get(c1);
				System.arraycopy(fz, x1, fs, 0, 2);
				fz = i<i2 ? c2.allf.get(clust) : clust.allf.get(c2);
				System.arraycopy(fz, x2, fs, 2, 2);
				allf.put(clust, fs);
			}

			// now update minf
			for(Map.Entry<Cluster, double[]> ent : allf.entrySet())
				minf.put(ent.getKey(), StatUtils.min(ent.getValue()));

			// now update minRfs
			pairs.remove(pair);

			// remove pairs with c1 c2
			for(Map.Entry<Cluster, Integer> ent : clusts.entrySet()) {
				clust = ent.getKey();
				int i = ent.getValue();
				if(i<i1) {
					minRfs.get(c1.minf.get(clust)).remove(new ClustPair(c1, clust));
				} else {
					minRfs.get(clust.minf.get(c1)).remove(new ClustPair(clust, c1));
				}
				if(i<i2) {
					minRfs.get(c2.minf.get(clust)).remove(new ClustPair(c2, clust));
				} else {
					minRfs.get(clust.minf.get(c2)).remove(new ClustPair(clust, c2));
				}
			}

			// add pairs with c
			for(Map.Entry<Cluster, Double> ent : minf.entrySet()) {
				f = ent.getValue();
				if(!minRfs.containsKey(f)) minRfs.put(f, new HashSet<>());
				minRfs.get(f).add(new ClustPair(c, ent.getKey()));
			}

			// clear empty entries
			fs2cl.clear();
			for(Map.Entry<Double, Set<ClustPair>> ent : minRfs.entrySet())
				if(ent.getValue().isEmpty()) fs2cl.add(ent.getKey());
			for(double f2cl : fs2cl) minRfs.remove(f2cl);

			// add new cluster c and update nclust
			clusts.put(c, nclust);
			++nclust;
		}
		
		myLogger.info("####clusters: "+clusts.size());
		StringBuilder out = new StringBuilder();
		List<Integer> ids;
		List<Double> dists;
		List<Boolean> joins;
		
		for(Map.Entry<Cluster, Integer> ent : clusts.entrySet()) {
			c = ent.getKey();
			out.setLength(0);
			ids = c.ids;
			dists = c.dists;
			joins = c.joins;
			
			out.append("-c ");
			out.append(scaffs.getKey(ids.get(0)));
			for(int i=1; i<ids.size(); i++) {
				out.append(":");
				out.append(scaffs.getKey(ids.get(i)));	
			}
			
			if(dists.size()>0) {
				out.append(" -s ");
				out.append(dists.get(0));
				for(int i=1; i<dists.size(); i++) {
					out.append(":");
					out.append(dists.get(i));	
				}		
			}
			
			if(joins.size()>1) {
				out.append(" -r ");
				out.append(joins.get(0));
				for(int i=1; i<joins.size(); i++) {
					out.append(":");
					out.append(joins.get(i));	
				}		
			}
			
			myLogger.info("#"+ent.getValue()+"\t"+out.toString());
		}
		
	}

	private final class ClustPair { 
		private final Cluster c1;
		private final Cluster c2;
		private boolean j1;
		private boolean j2;
		private double extf = Double.NaN;

		public ClustPair(Cluster c1, Cluster c2) {
			this.c1 = c1;
			this.c2 = c2;
		}

		public void join() {
			// TODO Auto-generated method stub
			if(!Double.isNaN(extf)) return;
			double minf = c1.minf.get(c2);
			double ef, ef1, ef2;
			double[] allf = c1.allf.get(c2);
			double[] fs;
			int x1 = -1, x2 = -1;
			extf = Double.MAX_VALUE;

			for(int i=0; i<4; i++) {
				if(allf[i]==minf) {
					// free ends
					switch(i) {
					case 0:
						x1 = 2;
						x2 = 2;
						break;
					case 1:
						x1 = 2;
						x2 = 0;
						break;
					case 2:
						x1 = 0;
						x2 = 2;
						break;
					case 3:
						x1 = 0;
						x2 = 0;
						break;
					}
					ef1 = Double.MAX_VALUE;
					for(Map.Entry<Cluster, double[]> ent : c1.allf.entrySet()) {
						if(ent.getKey().equals(c2)) continue;
						fs = ent.getValue();
						for(int j=x1; j<x1+2; j++)
							if(ef1>fs[j]) ef1 = fs[j];
					}

					ef2 = Double.MAX_VALUE;
					for(Map.Entry<Cluster, double[]> ent : c2.allf.entrySet()) {
						if(ent.getKey().equals(c1)) continue;
						fs = ent.getValue();
						for(int j=x2; j<x2+2; j++)
							if(ef2>fs[j]) ef2 = fs[j];
					}

					if( (ef=ef1+ef2)<extf ) {
						extf = ef;
						j1 = x1!=0;
						j2 = x2!=0;
					}
				}
			}
			return;
		}

		@Override
		public int hashCode() {
			int hash = 17;
			hash = hash*31+c1.first;
			hash = hash*31+c1.last;
			hash = hash*31+c2.first;
			hash = hash*31+c2.last;
			return hash;   
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			ClustPair other = (ClustPair) obj;
			return this.c1.equals(other.c1) && 
					this.c2.equals(other.c2);
		}
	}

	private final class Cluster {
		private final int first;
		private final int last;
		private final Map<Cluster, Double> minf = new HashMap<>();
		private final Map<Cluster, double[]> allf = new HashMap<>();
		private final List<Integer> ids = new ArrayList<>();
		private final List<Double> dists = new ArrayList<>();
		private final List<Boolean> joins = new ArrayList<>();

		public Cluster(int first, int last, Cluster c1, Cluster c2, boolean j1, boolean j2, double j) {
			this.first = first;
			this.last  = last;
			
			if(!j1) {
				Collections.reverse(c1.ids);
				Collections.reverse(c1.dists);
				Collections.reverse(c1.joins);
				for(int i=0; i<c1.joins.size(); i++)
					c1.joins.set(i, !c1.joins.get(i));
			}
			if(!j2) {
				Collections.reverse(c2.ids);
				Collections.reverse(c2.dists);
				Collections.reverse(c2.joins);
				for(int i=0; i<c2.joins.size(); i++)
					c2.joins.set(i, !c2.joins.get(i));
			}
			ids.addAll(c1.ids);
			ids.addAll(c2.ids);
			dists.addAll(c1.dists);
			dists.add(j);
			dists.addAll(c2.dists);
			joins.addAll(c1.joins);
			joins.addAll(c2.joins);
		}

		public Cluster(int first) {
			// TODO Auto-generated constructor stub
			this.first = first;
			this.last  = first;
			ids.add(first);
			joins.add(true);
		}
		
		@Override
		public int hashCode() {
			int hash = 17;
			hash = hash*31+first;
			hash = hash*31+last;
			return hash;   
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Cluster other = (Cluster) obj;
			return this.first==other.first && 
					this.last==other.last;
		}
	}
}
