package cz1.hmm.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
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
import cz1.util.Constants;
import cz1.util.Executor;
import cz1.util.Utils;

public abstract class RFUtils extends Executor {
	private final static Logger myLogger = Logger.getLogger(RFUtils.class);
	
	protected final static double RF_INF = 0.4999999;
	
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
	protected final Set<String> conjPair = new HashSet<>();
	protected final Map<String, List<FileObject>> fileObj = new HashMap<>();
	
	protected final static Object lock = new Object();

	protected class FileObject {
		protected final String file;
		protected final int[] position;
		protected final double loglik;
		
		public FileObject(String file, 
				int[] position,
				double loglik) {
			this.file = file;
			this.position = position;
			this.loglik = loglik;
		}
	}

	protected class FileLoader implements Runnable {
		private final String file;

		public FileLoader(String file) {
			this.file = file;
		}
		
		@Override
		public void run() {
			// TODO Auto-generated method stub
			try {
				ModelReader modelReader = new ModelReader(file);
				int[] haps_observed = modelReader.getHapCounts();
				modelReader.close();
					
				long[] observed;
				double p;
				boolean drop = false;
				switch(goodness_of_fit) {
				case "fraction":
					double[] phases = new double[ploidy*2];
					for(int z=0; z<phases.length; z++) 
						phases[z] = (double) haps_observed[z];
					double expected = StatUtils.sum(phases)/ploidy/2;
					double maf = StatUtils.max(phases)/expected, 
							mif = StatUtils.min(phases)/expected;
					if( maf>skew_phi || mif<1/skew_phi) drop = true;
					myLogger.info("["+(drop?"drop](maf,":"keep](maf,")+maf+";mif,"+mif+") "+
							Utils.paste(haps_observed,",")+"\t"+file);
					break;
				case "chisq":
					observed = new long[ploidy*2];
					for(int z=0; z<observed.length; z++) 
						observed[z] = (long) haps_observed[z];
					p = new ChiSquareTest().chiSquareTest(probs_uniform, observed);
					if(p<skew_phi) drop = true;
					myLogger.info("["+(drop?"drop](p,":"keep](p,")+formatter.format(p)+") "+
							Utils.paste(haps_observed,",")+"\t"+file);
					break;
				case "gtest":
					observed = new long[ploidy*2];
					for(int z=0; z<observed.length; z++) 
						observed[z] = (long) haps_observed[z];
					p = new GTest().gTest(probs_uniform, observed);
					if(p<skew_phi) drop = true;
					myLogger.info("["+(drop?"drop](p,":"keep](p,")+formatter.format(p)+") "+
							Utils.paste(haps_observed,",")+"\t"+file);
					break;
				default:
					throw new RuntimeException("Goodness-of-fit test should be fraction, chisq or gTest.");
				}
				
				if(drop) return;

				modelReader = new ModelReader(file);
				double[] ll = modelReader.getModelLoglik();
				String[] chrs = modelReader.getChrs();
				boolean[] chrs_rev = modelReader.getChrsRev();
				int[] model_length = modelReader.getModelLength();
				modelReader.close();
				
				int chrs_n = chrs.length;
				int[][] positions = new int[chrs_n][2];
				
				int cuv = 0;
				for(int i=0; i<chrs_n; i++) {
					if(chrs_rev[i]) {
						positions[i][1] = cuv;
						cuv += model_length[i];
						positions[i][0] = cuv-1;
					} else {
						positions[i][0] = cuv;
						cuv += model_length[i];
						positions[i][1] = cuv-1;
					}
				}
				
				// conjunctive pairs
				for(int i=0; i<chrs_n; i++) {
					String scaff_i = chrs[i];
					for(int j=i+1; j<chrs_n; j++) {
						String scaff_j = chrs[j];
						conjPair.add(scaff_i+Constants.collapsed_str+scaff_j);
						conjPair.add(scaff_j+Constants.collapsed_str+scaff_i);
					}
				}

				String chr;
				for(int j=0; j<chrs_n; j++) {
					chr = chrs[j]; 
					synchronized(lock) {
						if(!fileObj.containsKey(chr))
							fileObj.put(chr, new ArrayList<>());
						fileObj.get(chr).add(new FileObject(
								file,
								positions[j],
								ll[j]));
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
	
	protected int[] max(int[] arr1, int[] arr2) {
		// TODO Auto-generated method stub
		int m1 = max(arr1), m2 = max(arr2);
		if(m1>m2) return arr1;
		else if(m1<m2) return arr2;
		else {
			int s1 = sum(arr1), s2 = sum(arr2);
			if(s1>=s2) return arr1;
			else return arr2;
		}
	}

	protected int max(int[] arr) {
		// TODO Auto-generated method stub
		int m = Integer.MIN_VALUE;
		for(int a : arr) 
			if(a>m) m = a;
		return m;
	}
	
	protected void sum(int[] arr1, int[] arr2, int[] arr) {
		// TODO Auto-generated method stub
		for(int i=0; i<arr.length; i++)
			arr[i] = arr1[i]+arr2[i];
	}
	
	protected void sum(double[] arr1, double[] arr2, double[] arr) {
		// TODO Auto-generated method stub
		for(int i=0; i<arr.length; i++)
			arr[i] = arr1[i]+arr2[i];
	}

	protected int sum(int[] arr) {
		// TODO Auto-generated method stub
		int s = 0;
		for(int a : arr) s+=a;
		return s;
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

		this.initial_thread_pool();
		for(File f : listFiles) 
			executor.submit(new FileLoader(f.getAbsolutePath()));
		this.waitFor();
		
		Set<String> drops = new HashSet<>();
		int size;
		for(String chr : fileObj.keySet()) {
			List<FileObject> obj = fileObj.get(chr);
			size = obj.size();
			if(size<drop_thres) {
				myLogger.info("[drop] "+chr+" ("+size+")");
				drops.add(chr);	
			} else {
				myLogger.info("[keep] "+chr+" ("+size+")");
				Collections.sort(obj, new Comparator<FileObject>() {

					@Override
					public int compare(FileObject obj0, FileObject obj1) {
						// TODO Auto-generated method stub
						return Double.compare(obj1.loglik, obj0.loglik);
					}
					
				});
				obj.subList(Math.min(best_n, size), size).clear();
			}
		}
		fileObj.keySet().removeAll(drops);
	}
	
	/** 
	 * make RData object from Renjin 
	 * renjin-script-engine-*-with dependencies.jar required
	 * https://nexus.bedatadriven.com/content/groups/public/org/renjin/renjin-script-engine/
	**/
	protected static void makeRMatrix(String in_rf, String out_Rmat, int hs) {
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
			double d, l;
			while( (line=br.readLine())!=null &&
					line.startsWith("##")) {
				scaffs.put(w++, line.replaceAll("^##", ""));
			}
			
			int n = scaffs.size();
			int A = w*(w-1)/2;
			DoubleMatrixBuilder dMat = new DoubleMatrixBuilder(n,n);
			DoubleMatrixBuilder lMat = new DoubleMatrixBuilder(n,n);
			DoubleMatrixBuilder iMat = new DoubleMatrixBuilder(n,n);
			DoubleMatrixBuilder dAllMat = new DoubleMatrixBuilder(A*2,4);
			DoubleMatrixBuilder lAllMat = new DoubleMatrixBuilder(A*2,4);
			
			w = 0;
			while( line!=null ) {
				s = line.split("\\s+");
				int i=scaffs.getKey(s[5]),
						j=scaffs.getKey(s[6]);
				d = Math.min(RF_INF, Double.parseDouble(s[0]));
				l = calcLODFromRf(d, hs);
				dMat.set(i,j,d);
				dMat.set(j,i,d);
				iMat.set(i,j,w+1);
				iMat.set(j,i,w+1+A);
				lMat.set(i,j,l);
				lMat.set(j,i,l);
				for(int k=0; k<4; k++) {
					d = Math.min(RF_INF, Double.parseDouble(s[k+1]));
					l = calcLODFromRf(d, hs);
					dAllMat.set(w, k, d);
					dAllMat.set(w+A, (k==0||k==3)?k:(3-k), d);
					lAllMat.set(w, k, l);
					lAllMat.set(w+A, (k==0||k==3)?k:(3-k), l);
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
			lMat.setRowNames(scf);
			lMat.setColNames(scf);
			
			Context context = Context.newTopLevelContext();
			FileOutputStream fos = new FileOutputStream(out_Rmat);
			GZIPOutputStream zos = new GZIPOutputStream(fos);
			RDataWriter writer = new RDataWriter(context, zos);
			
			ListVector.NamedBuilder Rdat = new ListVector.NamedBuilder();
			Rdat.add("scaffs", scf);
			Rdat.add("n", n);
			Rdat.add("A", A);
			Rdat.add("distanceMat", dMat.build());
			Rdat.add("distanceAll", dAllMat.build());
			Rdat.add("lodMat", lMat.build());
			Rdat.add("loadAll", lAllMat.build());
			Rdat.add("indexMat", iMat.build());
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
		double lod = Double.NEGATIVE_INFINITY, lod1;
		while(true) {
			rf = (lb+ub)/2;
			lod1 = lod;
			lod = calcLODFromRf(rf, n);
			if(Math.abs(lod-lod_thres)<1e-12 || 
					Math.abs(lod-lod1)<1e-12 ||
					rf>=RF_INF) {
				return Math.min(RF_INF, rf);
			} else if(lod<lod_thres) {
				ub = rf;
			} else {
				lb = rf;
			}
		}
	}
	
	public static double calcLODFromRf(double theta, int n) {
		// TODO Auto-generated method stub
		theta = Math.min(theta, RF_INF);
		double r = n*theta;
		return (n-r)*Math.log10(1-theta)+r*Math.log10(theta)-n*Math.log10(0.5);
	}
	
	public static double geneticDistance(double r, String mapFunc) {
		// TODO Auto-generated method stub
		// if(r>=0.5) return Double.POSITIVE_INFINITY;
		r = Math.min(r, RF_INF);
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
			return Math.min(RF_INF, .5*(Math.exp(4*d)-1)/(Math.exp(4*d)+1));
		case "haldane":
			return Math.min(RF_INF, .5*(1-Math.exp(-2*d)));	
		default:
			throw new RuntimeException("Undefined genetic mapping function.");
		}
	}
}
