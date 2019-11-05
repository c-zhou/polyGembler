package cz1.hmm.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
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
	private double rf_thres = RFUtils.inverseGeneticDistance(0.5, "kosambi");
	private String out_prefix = null;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--input                  Recombination frequency file.\n"
						+ " -r/--rf                     Recombination frequency threshold (default: 0.38).\n"
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
			myArgsEngine.add("-i", "--input", true);
			myArgsEngine.add("-h", "--hap-size", true);
			myArgsEngine.add("-l", "--lod", true);
			myArgsEngine.add("-r", "--rf", true);
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
		
		if(myArgsEngine.getBoolean("-r")) {
			rf_thres = Double.parseDouble(myArgsEngine.getString("-r"));
		}
		
		if(myArgsEngine.getBoolean("-o")) {
			out_prefix = myArgsEngine.getString("-o");
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify your output file prefix.");
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		myLogger.info("Using recombination frequency threshold: "+rf_thres+".");
		nj(rf_thres);
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
					c1 = clusts.getKey(i1);
					c2 = clusts.getKey(i2);
					
					// recomb freqs arranged in i1<i2 order
					allf = new double[4];
					allf[0] = Double.parseDouble(s[1]);
					allf[3] = Double.parseDouble(s[4]);
					if(i1<i2) {
						allf[1] = Double.parseDouble(s[2]);
						allf[2] = Double.parseDouble(s[3]);
					} else {
						allf[1] = Double.parseDouble(s[3]);
						allf[2] = Double.parseDouble(s[2]);	
					}
					c1.minf.put(c2, minf);
					c1.allf.put(c2, allf);
					c2.minf.put(c1, minf);
					c2.allf.put(c1, allf);
					
					// cluster pair arranged in i1<i2 order
					pair = i1<i2?new ClustPair(c1, c2):new ClustPair(c2, c1);
					minRfs.putIfAbsent(minf, new HashSet<>());
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
		Cluster c1, c2, c, cc;
		boolean j1, j2;
		Map<Cluster, Double> minf;
		Map<Cluster, double[]> allf;
		int i1, i2;
		double[] fs, fs1, fs2;
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
			i1 = c1.clust;
			i2 = c2.clust;

			c = new Cluster(nclust, c1, c2, j1, j2, f);
			allf = c.allf;
			minf = c.minf;

			// now update clusts
			clusts.remove(c1);
			clusts.remove(c2);

			// now update allf and minf
			for(int i : clusts.values()) {
				cc = clusts.getKey(i);
				fs = new double[4];
				//** for c1
				fs1 = c1.allf.get(cc);
				if(i<i1) {
					if(j1) {
						fs[0] = fs1[1];
						fs[2] = fs1[3];
					} else {
						fs[0] = fs1[0];
						fs[2] = fs1[2];
					}
				} else {
					if(j1) {
						fs[0] = fs1[2];
						fs[2] = fs1[3];
					} else {
						fs[0] = fs1[0];
						fs[2] = fs1[1];
					}
				}
				
				//** for c2
				fs2 = c2.allf.get(cc);
				if(i<i2) {
					if(j2) {
						fs[1] = fs2[0];
						fs[3] = fs2[2];
					} else {
						fs[1] = fs2[1];
						fs[3] = fs2[3];
					}
				} else {
					if(j2) {
						fs[1] = fs2[0];
						fs[3] = fs2[1];
					} else {
						fs[1] = fs2[2];
						fs[3] = fs2[3];
					}
				}
				
				f = StatUtils.min(fs);
				allf.put(cc, fs);
				minf.put(cc, f);
				cc.allf.put(c, fs);
				cc.minf.put(c, f);
			}

			// now update minf
			for(Map.Entry<Cluster, double[]> ent : allf.entrySet())
				minf.put(ent.getKey(), StatUtils.min(ent.getValue()));

			// now update minRfs
			pairs.remove(pair);

			// remove pairs with c1 c2
			for(Map.Entry<Cluster, Integer> ent : clusts.entrySet()) {
				cc = ent.getKey();
				int i = cc.clust;
				if(i1<i) {
					minRfs.get(c1.minf.get(cc)).remove(new ClustPair(c1, cc));
				} else {
					minRfs.get(cc.minf.get(c1)).remove(new ClustPair(cc, c1));
				}
				if(i2<i) {
					minRfs.get(c2.minf.get(cc)).remove(new ClustPair(c2, cc));
				} else {
					minRfs.get(cc.minf.get(c2)).remove(new ClustPair(cc, c2));
				}
			}

			// add pairs with c
			for(Map.Entry<Cluster, Double> ent : minf.entrySet()) {
				f = ent.getValue();
				minRfs.putIfAbsent(f, new HashSet<>());
				minRfs.get(f).add(new ClustPair(ent.getKey(), c));
			}

			// clear empty entries in minRfs
			fs2cl.clear();
			for(Map.Entry<Double, Set<ClustPair>> ent : minRfs.entrySet())
				if(ent.getValue().isEmpty()) fs2cl.add(ent.getKey());
			for(double f2cl : fs2cl) minRfs.remove(f2cl);
			
			// clear c1 c2 from minf and allf
			for(Cluster c0 : clusts.keySet()) {
				c0.minf.remove(c1);
				c0.minf.remove(c2);
				c0.allf.remove(c1);
				c0.allf.remove(c2);
			}
			
			// add new cluster c and update nclust
			clusts.put(c, nclust);
			++nclust;
		}
		
		myLogger.info("####clusters: "+clusts.size());
		StringBuilder out = new StringBuilder();
		List<Integer> ids;
		List<Double> dists;
		List<Boolean> joins;
		
		try {
			BufferedWriter bw = Utils.getBufferedWriter(this.out_prefix+".nns");
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
				for(int i=0; i<ids.size()-1;i++)
					myLogger.info(scaffs.getKey(ids.get(i))+"\t"+dists.get(i)+"\t"+joins.get(i));
				myLogger.info(scaffs.getKey(ids.get(ids.size()-1))+"\t\t\t\t"+joins.get(joins.size()-1));
				
				out.append("\n");
				bw.write(out.toString());
			}
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
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
			boolean j1, j2;
			
			Cluster cc;
			extf = Double.MAX_VALUE;

			for(int i=0; i<4; i++) {
				if(allf[i]==minf) {
					ef1 = Double.MAX_VALUE;
					ef2 = Double.MAX_VALUE;
					
					switch(i) {
					case 0:
						j1 = true;
						j2 = false;
						break;
					case 1:
						j1 = true;
						j2 = true;
						break;
					case 2:
						j1 = false;
						j2 = false;
						break;
					case 3:
						j1 = false;
						j2 = true;
						break;
					default:
						throw new RuntimeException("!!!");
					}
					
					if(j1) {
						for(Map.Entry<Cluster, double[]> ent : c1.allf.entrySet()) {
							cc = ent.getKey();
							if(cc.equals(c2)) continue;
							fs = ent.getValue();
							if(c1.clust<cc.clust) {
								if(ef1>fs[2]) ef1 = fs[2];
							} else {
								if(ef1>fs[1]) ef1 = fs[1];
							}
							if(ef1>fs[3]) ef1 = fs[3];
						}
					} else {
						for(Map.Entry<Cluster, double[]> ent : c1.allf.entrySet()) {
							cc = ent.getKey();
							if(cc.equals(c2)) continue;
							fs = ent.getValue();
							if(c1.clust<cc.clust) {
								if(ef1>fs[1]) ef1 = fs[1];
							} else {
								if(ef1>fs[2]) ef1 = fs[2];
							}
							if(ef1>fs[0]) ef1 = fs[0];
						}
					}
					
					if(j2) {
						for(Map.Entry<Cluster, double[]> ent : c2.allf.entrySet()) {
							cc = ent.getKey();
							if(cc.equals(c1)) continue;
							fs = ent.getValue();
							if(c2.clust<cc.clust) {
								if(ef2>fs[1]) ef2 = fs[1];
							} else {
								if(ef2>fs[2]) ef2 = fs[2];
							}
							if(ef2>fs[0]) ef2 = fs[0];
						}
					} else {
						for(Map.Entry<Cluster, double[]> ent : c2.allf.entrySet()) {
							cc = ent.getKey();
							if(cc.equals(c1)) continue;
							fs = ent.getValue();
							if(c2.clust<cc.clust) {
								if(ef2>fs[2]) ef2 = fs[2];
							} else {
								if(ef2>fs[1]) ef2 = fs[1];
							}
							if(ef2>fs[3]) ef2 = fs[3];
						}
					}
					
					if( (ef=ef1+ef2)<extf ) {
						this.extf = ef;
						this.j1 = j1;
						this.j2 = j2;
					}
				}
			}
			return;
		}

		@Override
		public int hashCode() {
			int hash = 17;
			hash = hash*31+c1.clust;
			hash = hash*31+c2.clust;
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
		private final int clust;
		private final Map<Cluster, Double> minf = new HashMap<>();
		private final Map<Cluster, double[]> allf = new HashMap<>();
		private final List<Integer> ids = new ArrayList<>();
		private final List<Double> dists = new ArrayList<>();
		private final List<Boolean> joins = new ArrayList<>();

		public Cluster(int clust, Cluster c1, Cluster c2, boolean j1, boolean j2, double j) {
			this.clust = clust;
			if(j1) {
				Collections.reverse(c1.ids);
				Collections.reverse(c1.dists);
				Collections.reverse(c1.joins);
				for(int i=0; i<c1.joins.size(); i++)
					c1.joins.set(i, !c1.joins.get(i));
			}
			if(j2) {
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

		public Cluster(int clust) {
			// TODO Auto-generated constructor stub
			this.clust = clust;
			ids.add(clust);
			joins.add(false);
		}
		
		@Override
		public int hashCode() {
			return clust;
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
			return this.clust==other.clust;
		}
	}
}
