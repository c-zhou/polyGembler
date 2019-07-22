package cz1.hmm.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPOutputStream;

import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.TreeBidiMap;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.apache.commons.math3.stat.inference.GTest;
import org.apache.log4j.Logger;
import org.renjin.eval.Context;
import org.renjin.primitives.io.serialization.RDataWriter;
import org.renjin.primitives.matrix.DoubleMatrixBuilder;
import org.renjin.sexp.ListVector;
import org.renjin.sexp.StringArrayVector;
import org.renjin.sexp.StringVector;

import cz1.hmm.model.ModelReader;
import cz1.util.Executor;
import cz1.util.Utils;

public abstract class RFUtils extends Executor {
	private final static Logger myLogger = Logger.getLogger(RFUtils.class);
	
	protected String in_haps;
	protected String expr_id = null;
	protected String[] parents;
	protected String[] progeny;
	protected int ploidy;
	protected int drop_thres = 1;
	protected double skew_phi = 2.0;
	protected int best_n = 1;
	protected final String goodness_of_fit = "fraction";

	protected double[] probs_uniform = new double[]{.5,.5,.5,.5};

	protected NumberFormat formatter = new DecimalFormat("#0.000");
	protected int nF1;
	// protected final Set<String> conjPair = new HashSet<>();
	protected final Map<String, List<FileObject>> fileObj = new HashMap<>();
	
	protected final static Object lock = new Object();

	protected class FileObject {
		protected final String file;
		protected final String[] markers; 
		protected final int[] position;

		public FileObject(String file, 
				String[] markers,
				int[] position) {
			this.file = file;
			this.markers = markers;
			this.position = position;
		}
	}

	protected class FileLoader implements Runnable {
		private final String[] files;

		public FileLoader(String[] files) {
			this.files = files;
		}
		
		@Override
		public void run() {
			// TODO Auto-generated method stub
			try {
				String id = new File(files[0]).getName().replace(expr_id,"experiment").split("\\.")[1];
				
				double[] ll = new double[files.length];
				for(int i=0; i<ll.length; i++) {
					ModelReader modelReader = new ModelReader(files[i]);
					ll[i] = modelReader.getLoglik();
					modelReader.close();
				}
				
				int[] maxN = maxN(ll);
				StringBuilder oos = new StringBuilder(id+"\n");
				boolean[] drop = new boolean[maxN.length];
				int dropped = 0;
				for(int k=0; k<maxN.length; k++) {
					ModelReader modelReader = new ModelReader(files[maxN[k]]);
					int[] haps_observed = modelReader.getHapCounts();
					modelReader.close();
					
					long[] observed;
					double p;
					switch(goodness_of_fit) {
					case "fraction":
						double[] phases = new double[ploidy*2];
						for(int z=0; z<phases.length; z++) 
							phases[z] = (double) haps_observed[z];
						double expected = StatUtils.sum(phases)/ploidy/2;
						double maf = StatUtils.max(phases)/expected, 
								mif = StatUtils.min(phases)/expected;
						if( maf>skew_phi || mif<1/skew_phi) {
							myLogger.info(files[maxN[k]]+
									" was dropped due to large haploptype frequency variance. (" +
									cat(phases, ",") +")");
							drop[k] = true;
						}
						if(drop[k]) dropped++;
						oos.append("["+(drop[k]?"drop](maf,":"keep](maf,")+maf+";mif,"+mif+") "+
								cat(haps_observed,",")+"\t"+
								files[maxN[k]]+"\n");
						break;
					case "chisq":
						observed = new long[ploidy*2];
						for(int z=0; z<observed.length; z++) 
							observed[z] = (long) haps_observed[z];
						p = new ChiSquareTest().chiSquareTest(probs_uniform, observed);
						if(p<skew_phi) drop[k] = true;
						if(drop[k]) dropped++;
						oos.append("["+(drop[k]?"drop](p,":"keep](p,")+formatter.format(p)+") "+
								cat(haps_observed,",")+"\t"+
								files[maxN[k]]+"\n");
						break;
					case "gtest":
						observed = new long[ploidy*2];
						for(int z=0; z<observed.length; z++) 
							observed[z] = (long) haps_observed[z];
						p = new GTest().gTest(probs_uniform, observed);
						if(p<skew_phi) drop[k] = true;
						if(drop[k]) dropped++;
						oos.append("["+(drop[k]?"drop](p,":"keep](p,")+formatter.format(p)+") "+
								cat(haps_observed,",")+"\t"+
								files[maxN[k]]+"\n");
						break;
					default:
						throw new RuntimeException("Goodness-of-fit test should be fraction, chisq or gTest.");
					}
				}
				myLogger.info(oos.toString());
				myLogger.info(id+" - dropped "+dropped);
				if( drop.length-dropped<drop_thres ) {
					myLogger.info("Scaffold "+id+" dropped.");
					return;
				}
				
				List<String> selected = new ArrayList<>(); 
				int z=0;
				for(int k=0; k<drop.length; k++) {
					if(!drop[k]) {
						selected.add(files[maxN[k]]);
						z++;
					}
					if(z>=best_n) break;
				}
				
				List<String> scaff_all = new ArrayList<String>();
				String[][] markers = null;
				int[][] positions = null;
				int scaff_n = 0;
				
				ModelReader modelReader = new ModelReader(selected.get(0));
				List<String> snpId = modelReader.getSnpId();
				modelReader.close();
				List<List<String>> markers_all = new ArrayList<List<String>>();

				String marker = snpId.get(0);
				String scaff_prev = marker.replaceAll("_[0-9]{1,}$", ""), scaff;
				scaff_all.add(scaff_prev);
				markers_all.add(new ArrayList<String>());
				int n = 0;
				markers_all.get(n).add(marker);
				for(int i=1; i<snpId.size(); i++) {
					marker = snpId.get(i);
					scaff = marker.replaceAll("_[0-9]{1,}$", "");
					if(scaff.equals(scaff_prev))
						markers_all.get(n).add(marker);
					else {
						markers_all.add(new ArrayList<String>());
						n++;
						markers_all.get(n).add(marker);
						scaff_prev = scaff;
						scaff_all.add(scaff_prev);
					}
				}
				
				int cuv = 0;
				scaff_n = scaff_all.size();
				markers = new String[scaff_n][];
				positions = new int[scaff_n][2];
				for(int i=0; i<scaff_n; i++) {
					markers[i] = new String[markers_all.get(i).size()];
					markers_all.get(i).toArray(markers[i]);
					int s = Integer.parseInt(markers[i][0].
							replaceAll(".*[^\\d](\\d+).*", "$1")),
							e = Integer.parseInt(markers[i][1].
									replaceAll(".*[^\\d](\\d+).*", "$1"));
					if(s<=e) {
						positions[i][0] = cuv;
						cuv += markers[i].length;
						positions[i][1] = cuv-1;
					} else {
						positions[i][1] = cuv;
						cuv += markers[i].length;
						positions[i][0] = cuv-1;
					}
				}
				
				// conjunctive pairs
				/***
				for(int i=0; i<scaff_n; i++) {
					String scaff_i = scaff_all.get(i);
					for(int j=i+1; j<scaff_n; j++) {
						String scaff_j = scaff_all.get(j);
						conjPair.add(scaff_i+Constants.collapsed_str+scaff_j);
						conjPair.add(scaff_j+Constants.collapsed_str+scaff_i);
					}
				}
				**/

				for(String file : selected) {
					for(int j=0; j<scaff_n; j++) {
						scaff = scaff_all.get(j); 
						synchronized(lock) {
							if(!fileObj.containsKey(scaff))
								fileObj.put(scaff, new ArrayList<>());
							fileObj.get(scaff).add(new FileObject(
									file,
									markers[j],
									positions[j]));
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
		
		protected int[] maxN(double[] ll) {
			double[] ll0 = Arrays.copyOf(ll, ll.length);
			int n = ll.length;
			int[] maxN = new int[n];
			Arrays.fill(maxN, -1);
			for(int k=0; k<n; k++) {
				if(k>=ll0.length) return maxN;
				int p = 0;
				double e = Double.NEGATIVE_INFINITY;
				for(int s=0; s<ll0.length; s++)
					if(ll0[s]>e) {
						e = ll0[s];
						p = s;
					}
				maxN[k] = p;
				ll0[p] = Double.NEGATIVE_INFINITY;
			}
			return maxN;
		}

		protected int[] maxN(double[] ll, int N) {
			double[] ll0 = Arrays.copyOf(ll, ll.length);
			int[] maxN = new int[N];
			Arrays.fill(maxN, -1);
			for(int k=0; k<N; k++) {
				if(k>=ll0.length) return maxN;
				int p = 0;
				double e = Double.NEGATIVE_INFINITY;
				for(int s=0; s<ll0.length; s++)
					if(ll0[s]>e) {
						e = ll0[s];
						p = s;
					}
				maxN[k] = p;
				ll0[p] = Double.NEGATIVE_INFINITY;
			}
			return maxN;
		}
	}
	
	protected static String cat(double[] array, String sep) {
		String s = ""+array[0];
		for(int i=1; i<array.length; i++)
			s += sep+array[i];
		return s;
	}
	
	protected String cat(long[] array, String sep) {
		// TODO Auto-generated method stub
		String s = ""+array[0];
		for(int i=1; i<array.length; i++)
			s += sep+array[i];
		return s;
	}

	protected static String cat(int[] array, String sep) {
		String s = ""+array[0];
		for(int i=1; i<array.length; i++)
			s += sep+array[i];
		return s;
	}
	
	protected String guessExperimentId() {
		// TODO Auto-generated method stub
		File in_dir = new File(in_haps);
		File[] haps = in_dir.listFiles();
		Map<String, Integer> stats = new HashMap<String, Integer>();
		for(File hap : haps) {
			String h = hap.getName().split("\\.")[0];
			if(!stats.keySet().contains(h))
				stats.put(h, 0);
			stats.put(h, stats.get(h)+1);
		}
		String expr_id = null;
		int count = 0;
		for(String i : stats.keySet()) {
			if(stats.get(i)>count) {
				expr_id = i;
				count = stats.get(i);
			}
		}
		return expr_id;
	}
	
	protected void initialise() {
		// TODO Auto-generated catch block
		
		File folder = new File(in_haps);
		File[] listFiles = folder.listFiles();
		ModelReader modelReader = new ModelReader(listFiles[0].getAbsolutePath());
		nF1 = modelReader.getSampleNo()-2;
		ploidy = modelReader.getPloidy();
		parents = modelReader.getParents();
		progeny = modelReader.getProgeny();
		modelReader.close();
		probs_uniform = new double[ploidy*2];
		Arrays.fill(probs_uniform, .5/ploidy);
		myLogger.info(nF1+" F1 samples in the experiment.");

		Map<String, List<String>> map = new HashMap<String, List<String>>();
		List<String> list;
		String[] s;
		for(File file:listFiles) {
			String name = file.getName();
			if( name.startsWith(expr_id) ) {
				name = name.replace(expr_id,"experiment");
				s = name.split("\\.");

				if(map.get(s[1])==null) {
					list = new ArrayList<String>();
					list.add(file.getAbsolutePath());
					map.put(s[1], list);
				} else{
					map.get(s[1]).add(file.getAbsolutePath());
				}
			}
		}
		
		String[] keys = new String[map.keySet().size()];
		map.keySet().toArray(keys);

		this.initial_thread_pool();
		for(int i=0; i<keys.length; i++) {
			List<String> files = map.get(keys[i]);
			executor.submit(new FileLoader(
					files.toArray(new String[files.size()])));
		}
		this.waitFor();
		
		myLogger.info("["+Utils.getSystemTime()+"] LOADING FILES DONE.");
		myLogger.info("["+Utils.getSystemTime()+"] READING LOG LIKELIHOOD DONE.");
		myLogger.info(map.keySet().size());
	}
	
	/** 
	 * make RData object from Renjin 
	 * renjin-script-engine-*-with dependencies.jar required
	 * https://nexus.bedatadriven.com/content/groups/public/org/renjin/renjin-script-engine/
	**/
	protected static void makeRMatrix(String in_rf, String out_Rmat) {
		// TODO Auto-generated method stub
		ScriptEngineManager manager = new ScriptEngineManager();
		ScriptEngine engine = manager.getEngineByName("Renjin"); 
		if(engine == null) { 
			throw new RuntimeException("Renjin not found!!!"); 
		}
		try {
			BufferedReader br = Utils.getBufferedReader(in_rf);
			final BidiMap<Integer, String> scaffs = 
					new TreeBidiMap<Integer, String>();
			String line;
			String s[];
			int w = 0;
			while( (line=br.readLine())!=null &&
					line.startsWith("##")) {
				scaffs.put(w++, line.replaceAll("^##", ""));
			}
			
			int n = scaffs.size();
			int A = w*(w-1)/2;
			DoubleMatrixBuilder dMat = new DoubleMatrixBuilder(n,n);
			DoubleMatrixBuilder iMat = new DoubleMatrixBuilder(n,n);
			DoubleMatrixBuilder dAllMat = new DoubleMatrixBuilder(A*2,4);
			
			w = 0;
			while( line!=null ) {
				s = line.split("\\s+");
				int i=scaffs.getKey(s[5]),
						j=scaffs.getKey(s[6]);
				double d = Double.parseDouble(s[0]);
				dMat.set(i,j,d);
				dMat.set(j,i,d);
				iMat.set(i,j,w+1);
				iMat.set(j,i,w+1+A);
				for(int k=0; k<4; k++) {
					d = Double.parseDouble(s[k+1]);
					dAllMat.set(w, k, d);
					dAllMat.set(w+A, (k==0||k==3)?k:(3-k), d);
				}
				w++;
				line = br.readLine();
			}
			br.close();
			StringVector scf = new StringArrayVector(scaffs.values());
			dMat.setRowNames(scf);
			dMat.setColNames(scf);
			iMat.setRowNames(scf);
			iMat.setColNames(scf);
			
			Context context = Context.newTopLevelContext();
			FileOutputStream fos = new FileOutputStream(out_Rmat);
			GZIPOutputStream zos = new GZIPOutputStream(fos);
			RDataWriter writer = new RDataWriter(context, zos);
			
			ListVector.NamedBuilder Rdat = new ListVector.NamedBuilder();
			Rdat.add("scaffs", scf);
			Rdat.add("n", n);
			Rdat.add("A", A);
			Rdat.add("distanceAll", dAllMat.build());
			Rdat.add("indexMat", iMat.build());
			Rdat.add("distanceMat", dMat.build());
			writer.save(Rdat.build());
			writer.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static double calcRfFromLOD(double lod_thres, int n) {
		// TODO Auto-generated method stub
		double lb = 0, ub = 0.5;
		double rf;
		double lod;
		while(true) {
			rf = (lb+ub)/2;
			lod = calcLODFromRf(rf, n);
			if(lod<lod_thres) {
				ub = rf;
			} else if(lod-lod_thres<1e-6) {
				return rf;
			} else {
				lb = rf;
			}
		}
	}
	
	public static double calcLODFromRf(double theta, int n) {
		// TODO Auto-generated method stub
		double r = n*theta;
		return (n-r)*Math.log10(1-theta)+r*Math.log10(theta)-n*Math.log10(0.5);
	}
	
	public static double geneticDistance(double r, String mapFunc) {
		// TODO Auto-generated method stub
		switch(mapFunc.toLowerCase()) {
		case "kosambi":
			return .25*Math.log((1+2*r)/(1-2*r));
		case "haldane":
			return -.5*Math.log(1-2*r);	
		default:
			throw new RuntimeException("Undefined genetic mapping function.");
		}
	}
	
	public static double inverseGeneticDistance(double d, String mapFunc) {
		// TODO Auto-generated method stub
		switch(mapFunc.toLowerCase()) {
		case "kosambi":
			return .5*(Math.exp(4*d)-1)/(Math.exp(4*d)+1);
		case "haldane":
			return .5*(1-Math.exp(-2*d));	
		default:
			throw new RuntimeException("Undefined genetic mapping function.");
		}
	}
}
