package cz1.hmm.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipFile;

import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.TreeBidiMap;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.apache.commons.math3.stat.inference.GTest;
import org.renjin.eval.Context;
import org.renjin.primitives.io.serialization.RDataWriter;
import org.renjin.primitives.matrix.DoubleMatrixBuilder;
import org.renjin.sexp.ListVector;
import org.renjin.sexp.StringArrayVector;
import org.renjin.sexp.StringVector;

import cz1.util.Constants;
import cz1.util.Executor;
import cz1.util.Utils;

public abstract class RFUtils extends Executor {

	protected static final double kb = 1024;
	protected static final double mb = 1024*1024;
	protected static final double gb = 1024*1024*1024;

	protected String in_haps;
	protected String out_prefix;
	protected String[] founder_haps;
	protected String expr_id = null;
	protected int drop_thres = 1;
	protected double skew_phi = 2.0;
	protected int best_n = 10;
	protected final String goodness_of_fit = "fraction";

	protected double[] probs_uniform = new double[]{.5,.5,.5,.5};

	protected NumberFormat formatter = new DecimalFormat("#0.000");
	protected int nF1;

	protected final Map<String, double[]> conjPairRFs = 
			new HashMap<String, double[]>();

	protected final Map<String, Set<FileObject>> fileObj = 
			new HashMap<String, Set<FileObject>>();
	protected final Object lock = new Object();

	protected class FileObject {
		protected final String file;
		protected final String[] markers; 
		protected final int[] start_end_position;

		public FileObject(String file, 
				String[] markers,
				int[] start_end_position) {
			this.file = file;
			this.markers = markers;
			this.start_end_position = start_end_position;
		}
	}

	protected class FileExtraction implements Runnable {
		private final String[] files;

		public FileExtraction(String[] files) {
			this.files = files;
		}
		
		@Override
		public void run() {
			// TODO Auto-generated method stub

			try {
				List<String> scaff_all = new ArrayList<String>();
				String[][] markers = null;
				int[][] start_end_position = null;
				int scaff_n = 0;

				InputStreamObj isObj = new InputStreamObj(this.files[0]);
				isObj.getInputStream("SNP");
				BufferedReader br = Utils.getBufferedReader(isObj.is);
				List<List<String>> markers_all = new ArrayList<List<String>>();

				String marker = br.readLine().split("\\s+")[3];
				String scaff_prev = marker.replaceAll("_[0-9]{1,}$", ""),
						scaff;
				scaff_all.add(scaff_prev);
				String line;
				markers_all.add(new ArrayList<String>());
				int n_=0;
				markers_all.get(n_).add(marker);
				while( (line=br.readLine())!=null ) {
					marker = line.split("\\s+")[3];
					scaff = marker.replaceAll("_[0-9]{1,}$", "");
					if(scaff.equals(scaff_prev))
						markers_all.get(n_).add(marker);
					else {
						markers_all.add(new ArrayList<String>());
						n_++;
						markers_all.get(n_).add(marker);
						scaff_prev = scaff;
						scaff_all.add(scaff_prev);
					}
				}
				br.close();
				isObj.close();
				int cuv = 0;
				scaff_n = scaff_all.size();
				markers = new String[scaff_n][];
				start_end_position = new int[scaff_n][2];
				for(int i=0; i<scaff_n; i++) {
					markers[i] = new String[markers_all.get(i).size()];
					markers_all.get(i).toArray(markers[i]);
					int s = Integer.parseInt(markers[i][0].
							replaceAll(".*[^\\d](\\d+).*", "$1")),
							e = Integer.parseInt(markers[i][1].
									replaceAll(".*[^\\d](\\d+).*", "$1"));
					if(s<=e) {
						start_end_position[i][0] = cuv;
						cuv += markers[i].length;
						start_end_position[i][1] = cuv-1;
					} else {
						start_end_position[i][1] = cuv;
						cuv += markers[i].length;
						start_end_position[i][0] = cuv-1;
					}
				}
				
				// conjugate pairs
				if(scaff_n>1) {
					char[][][] haps = this.readConjHaps(start_end_position);
					int n = this.files.length;
					for(int i=0; i<scaff_n; i++) {
						String scaff_i = scaff_all.get(i);
						for(int j=i+1; j<scaff_n; j++) {
							String scaff_j = scaff_all.get(j);
							double[][] stats = new double[n][4];
							try {
							for(int k=0; k<n; k++) {
								char[][] haps_k = haps[k];
								stats[k][0] = calcRFs(haps_k, i*2,   j*2  );
								stats[k][1] = calcRFs(haps_k, i*2,   j*2+1);
								stats[k][2] = calcRFs(haps_k, i*2+1, j*2  );
								stats[k][3] = calcRFs(haps_k, i*2+1, j*2+1);
							} 
							}catch(java.lang.ArrayIndexOutOfBoundsException e) {
								e.printStackTrace();
							}
							double[] rf = max(stats);
							mininus(1, divide(rf, nF1*Constants._ploidy_H));
							conjPairRFs.put(scaff_i+Constants.scaff_collapsed_str+scaff_j, rf);
							
							myLogger.info("conjugate pair: "+scaff_i+" "+scaff_j);
						}
					}
				}

				for(int i=0; i<this.files.length; i++) {
					for(int j=0; j<scaff_n; j++) {
						scaff = scaff_all.get(j); 
						synchronized(lock) {
							if(!fileObj.containsKey(scaff))
								fileObj.put(scaff, new HashSet<FileObject>());
							fileObj.get(scaff).add(new FileObject(
									this.files[i],
									markers[j],
									start_end_position[j]));
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

		private double[] max(double[][] ds) {
			// TODO Auto-generated method stub
			double tot=Double.NEGATIVE_INFINITY, 
					rf=Double.NEGATIVE_INFINITY;
			double f;
			int r = -1;
			for(int i=0; i<ds.length; i++) 
				if( (f=StatUtils.max(ds[i]))>rf ||
						f==rf && StatUtils.sum(ds[i])>tot ) {
					rf = f;
					tot = StatUtils.sum(ds[i]);
					r = i;
				}
			return ds[r];
		}
		
		private double[] divide(double[] ds, double deno) {
			// TODO Auto-generated method stub
			for(int i=0; i<4; i++) ds[i] /= deno;
			return ds;
		}
		
		private double[] mininus(double d, double[] minus) {
			// TODO Auto-generated method stub
			for(int i=0; i<4; i++) minus[i] = d-minus[i];
			return minus;
		}
	
		private double calcRFs(final char[][] haps, 
				final int i, 
				final int j) {
			// TODO Auto-generated method stub
			int n = nF1*Constants._ploidy_H;
			int z = 0;
			for(int k=0; k<n; k++)
				if(haps[k][i]==haps[k][j])
					z++;
			return z;
		}

		private char[][][] readConjHaps(int[][] start_end_position) {
			// TODO Auto-generated method stub
			int n1 = this.files.length;
			int n2 = nF1*Constants._ploidy_H;
			int n3 = start_end_position.length;
			char[][][] haps = new char[n1][n2][n3*2];
			
			for(int i=0; i<n1; i++) {
				try {
					InputStreamObj isObj = new InputStreamObj(this.files[i]);
					isObj.getInputStream("PHASEDSTATES");
					BufferedReader br = Utils.getBufferedReader(isObj.is);
					
					String line, stateStr;
					String[] s;
					int j = 0;
					while( (line=br.readLine())!=null ) {
						if(!line.startsWith("#")) continue;
						s = line.split("\\s+|:");
						if(Arrays.asList(founder_haps).contains(s[2])) continue;
						stateStr = s[s.length-1];
						for(int k=0; k<n3; k++) {
							haps[i][j][k*2]   = stateStr.charAt(start_end_position[k][0]);
							haps[i][j][k*2+1] = stateStr.charAt(start_end_position[k][1]);
						}
						j++;
					}
					br.close();
					isObj.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			
			
			return haps;
		}
	}

	protected abstract class FileLoader implements Runnable {
		protected final String id;
		protected final FileObject[] files;
		protected final int i;

		public FileLoader(String id, FileObject[] files, int i) {
			this.id = id;
			this.files = files;
			this.i = i;
		}
		
		@Override
		public void run() {
			// TODO Auto-generated method stub

			try {
				double[] ll = new double[this.files.length];
				for(int k=0; k<ll.length; k++) 
					ll[k]=Double.NEGATIVE_INFINITY;
				int marker_n = this.files[0].markers.length;

				if( marker_n<2 ) {
					myLogger.warn("warning: "+
							this.id +" #marker is less than 2.");
					return;
				}

				for(int k=0; k<this.files.length; k++) {
					String file = this.files[k].file;
					try {
						InputStreamObj isObj = new InputStreamObj(file);
						if( !isObj.getInputStream("PHASEDSTATES") ) {
							myLogger.warn("warning: "+file+
									" exsits, but phased states do not.");
							continue;
						} 

						BufferedReader br = Utils.getBufferedReader(isObj.is);
						br.readLine();
						String marker_str = br.readLine();
						br.close();
						isObj.closeIs();

						if( marker_str==null ) {
							myLogger.warn("warning: "+file+
									" exists, but phased states are NULL.");
							continue;
						}
						double frac = marker_n/Double.parseDouble(marker_str);

						String line;
						if( isObj.getInputStream("STDERR") ) {
							BufferedReader br2 = Utils.getBufferedReader(isObj.is);
							String lprob=null;
							while( (line=br2.readLine()) !=null) {
								if(line.startsWith("log"))
									lprob = line;
							}
							br2.close();
							isObj.closeIs();
							if( lprob!=null ) {
								String[] s0 = lprob.split("\\s+");
								ll[k] = Double.parseDouble(s0[3])*frac;
							}
						} else {
							isObj.getInputStream("PHASEDSTATES");
							BufferedReader br2 = Utils.getBufferedReader(isObj.is);
							String lprob = br2.readLine();
							br2.close();
							isObj.closeIs();
							if( lprob!=null ) 
								ll[k] = Double.parseDouble(lprob)*frac;
							myLogger.info(ll[k]);
						}
						isObj.closeIn();
					} catch (IOException e) {
						e.printStackTrace();
						System.exit(1);
					}
				}

				int[] maxN = maxN(ll);
				StringBuilder oos = new StringBuilder();
				boolean[] drop = new boolean[maxN.length];
				int dropped = 0;
				for(int k=0; k<maxN.length; k++) {

					int[] haps_observed = readHaplotypes(maxN[k]);

					long[] observed;
					double p;
					switch(goodness_of_fit) {
					case "fraction":
						double[] phases = new double[Constants._ploidy_H*2];
						for(int z=0; z<phases.length; z++) 
							phases[z] = (double) haps_observed[z];
						double expected = StatUtils.sum(phases)/Constants._ploidy_H/2;
						double maf = StatUtils.max(phases)/expected, 
								mif = StatUtils.min(phases)/expected;
						if( maf>skew_phi || mif<1/skew_phi) {
							myLogger.info(this.files[maxN[k]].file+
									" was dropped due to large haploptype frequency variance. (" +
									cat(phases, ",") +")");
							drop[k] = true;
						}
						if(drop[k]) dropped++;
						oos.append("["+(drop[k]?"drop](maf,":"keep](maf,")+maf+";mif,"+mif+") "+
								cat(haps_observed,",")+"\t"+
								this.files[maxN[k]].file+"\n");
						break;
					case "chisq":
						observed = new long[Constants._haplotype_z];
						for(int z=0; z<observed.length; z++) 
							observed[z] = (long) haps_observed[z];
						p = new ChiSquareTest().chiSquareTest(probs_uniform, observed);
						if(p<skew_phi) drop[k] = true;
						if(drop[k]) dropped++;
						oos.append("["+(drop[k]?"drop](p,":"keep](p,")+formatter.format(p)+") "+
								cat(haps_observed,",")+"\t"+
								this.files[maxN[k]].file+"\n");
						break;
					case "gtest":
						observed = new long[Constants._haplotype_z];
						for(int z=0; z<observed.length; z++) 
							observed[z] = (long) haps_observed[z];
						p = new GTest().gTest(probs_uniform, observed);
						if(p<skew_phi) drop[k] = true;
						if(drop[k]) dropped++;
						oos.append("["+(drop[k]?"drop](p,":"keep](p,")+formatter.format(p)+") "+
								cat(haps_observed,",")+"\t"+
								this.files[maxN[k]].file+"\n");
						break;
					default:
						throw new RuntimeException("Goodness-of-fit test should be fraction, "
								+ "chisq or gTest.");
					}
				}
				myLogger.info(oos.toString());
				myLogger.info(this.id+" - dropped "+dropped);
				if( drop.length-dropped<drop_thres ) {
					myLogger.info("Scaffold "+this.id+" dropped.");
				} else {
					int z=0;
					for(int k=0; k<drop.length; k++) {
						if(!drop[k]) {
							FileObject fobj = this.files[maxN[k]];
							this.collectData(i, z, fobj);
							z++;
						}
						if(z>=best_n) break;
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
		
		abstract protected void collectData(final int i, final int z, final FileObject fobj);
		
		protected int[] readHaplotypes(final int i) {
			// TODO Auto-generated method stub
			try {
				InputStreamObj isObj = new InputStreamObj(this.files[i].file);
				isObj.getInputStream(getPhaseFile());
				BufferedReader br = Utils.getBufferedReader(isObj.is);
				String line, stateStr;
				String[] s;
				int[] haps_observed = new int[Constants._haplotype_z];
				while( (line=br.readLine())!=null ) {
					if(!line.startsWith("#")) continue;
					s = line.split("\\s+|:");
					stateStr = s[s.length-1];
					for(char h : stateStr.toCharArray())
						haps_observed[h>'9'?(h-'a'+9):(h-'1')]++;
				}
				br.close();
				isObj.close();
				return haps_observed;
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				System.exit(1);
			}
			return null;
		}

		protected abstract String getPhaseFile();

		protected int[] maxN(double[] ll) {
			double[] ll0 = Arrays.copyOf(ll, ll.length);
			int n = ll.length;//best_n_phases[0].length;
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

	protected final Map<String, double[][]> mapCalc = 
			new ConcurrentHashMap<String, double[][]>();

	abstract protected class MapCalculator implements Runnable {

		protected final int i;

		public MapCalculator(int i) {
			this.i = i;
		}
		
		protected void calcGDsAll(String phasedStates, int ploidy, 
				String[] parents, int nF1, int[] start_end,  
				double[][] rfAll, int s) {
			// TODO Auto-generated method stub
			char[][] h = readHaplotypes(phasedStates, 
					ploidy, parents, nF1);
			int c = 0;
			if(start_end[0]<=start_end[1]) {
				for(int i=start_end[0]; i<=start_end[1]; i++) { 
					for(int j=i+1; j<=start_end[1]; j++) {
						double r = 0;
						for(int k=0; k<h.length; k++) 
							r += h[k][i]==h[k][j] ? 0 : 1;
						rfAll[c++][s] = r/h.length;
					}
				}
			} else {
				for(int i=start_end[0]; i>=start_end[1]; i--) { 
					for(int j=i-1; j>=start_end[1]; j--) {
						double r = 0;
						for(int k=0; k<h.length; k++) 
							r += h[k][i]==h[k][j] ? 0 : 1;
						rfAll[c++][s] = r/h.length;
					}
				}
			}
		}

		protected double[] calcGDs(String phasedStates, int ploidy, 
				String[] parents, int nF1, int[] start_end) {
			// TODO Auto-generated method stub
			char[][] h = readHaplotypes(phasedStates, ploidy, parents, nF1);
			if(start_end[0]<=start_end[1]) {
				double[] d = new double[start_end[1]-start_end[0]+1];
				for(int i=start_end[0]; i<start_end[1]; i++) {
					double c = 0;
					for(int j=0; j<h.length; j++) 
						c += h[j][i]==h[j][i+1] ? 0 : 1;
					//d[i] = geneticDistance( c/h.length, mapFunc);
					d[i-start_end[0]] = c/h.length;
				}
				return d;
			} else {
				double[] d = new double[start_end[0]-start_end[1]+1];
				for(int i=start_end[0]; i>start_end[1]; i--) {
					double c = 0;
					for(int j=0; j<h.length; j++) 
						c += h[j][i]==h[j][i-1] ? 0 : 1;
					//d[i] = geneticDistance( c/h.length, mapFunc);
					d[start_end[0]-i] = c/h.length;
				}
				return d;
			}
		}

		protected double calcGD(String phasedStates, int ploidy,
				String[] parents, int nF1, int[] start_end) {
			char[][] h = readHaplotypes(phasedStates, ploidy, parents, nF1);
			double c = 0;
			for(int i=0; i<h.length; i++)
				c += h[i][start_end[0]]==h[i][start_end[1]] ? 0 : 1;
			//return geneticDistance( c/h.length, mapFunc);
			return c/h.length;	
		}

	}
	
	protected void writeMap() {
		// TODO Auto-generated method stub
		BufferedWriter bw = Utils.getBufferedWriter(out_prefix+".map");
		try {
			for(String scaff : mapCalc.keySet()) {
				double[][] rfs = mapCalc.get(scaff);
				if(rfs[0]==null) continue;
				bw.write("*");
				bw.write(scaff);
				bw.write("\n");
				for(double[] rf : rfs) {
					if(rf==null) break;
					bw.write(cat(rf,","));
					bw.write("\n");
				}
			}
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
	}

	protected double median(double[] ds) {
		// TODO Auto-generated method stub
		double[] ds0 = removeNEG(ds);
		if(ds0==null) return -1;
		Arrays.sort(ds0);
		int n = ds0.length;
		if (n % 2 == 0)
			return (ds0[n/2] + ds0[n/2-1])/2;
		else
			return ds0[n/2];
	}

	protected double[] removeNEG(double[] ds) {
		List<Double> ds0 = new ArrayList<Double>();
		for(double d : ds)
			if(d>=0)
				ds0.add(d);
		if(ds0.isEmpty()) return null;
		return ArrayUtils.toPrimitive(
				ds0.toArray(new Double[ds0.size()]));

	}

	protected void fill(double[][] dss, 
			double d) {
		// TODO Auto-generated method stub
		for(double[] ds : dss) 
			Arrays.fill(ds, d);
	}
	
	protected char[][] readHaplotypes(String phasedStates, int ploidy,
			String[] parents, int nF1) {
		// TODO Auto-generated method stub
		try {
			InputStreamObj isObj = new InputStreamObj(phasedStates);
			isObj.getInputStream("PHASEDSTATES");
			BufferedReader br = Utils.getBufferedReader(isObj.is);
			br.readLine();
			int m = Integer.parseInt(br.readLine());
			char[][] h = new char[nF1*ploidy][m];
			String line, stateStr;
			String[] s;
			int c = 0;
			while( (line=br.readLine())!=null ) {
				if(!line.startsWith("#")) continue;
				//if(skip++<2) continue;
				s = line.split("\\s+|:");
				if(Arrays.asList(parents).contains(s[2])) continue;
				stateStr = s[s.length-1];
				h[c++] = stateStr.toCharArray();
			}
			br.close();
			isObj.close();
			return h;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.exit(1);
		}

		return null;
	}

	protected static double geneticDistance(double r, String mapFunc) {
		// TODO Auto-generated method stub
		switch(mapFunc.toUpperCase()) {
		case "KOSAMBI":
			return .25*Math.log((1+2*r)/(1-2*r));
		case "HALDANE":
			return -.5*Math.log(1-2*r);	
		default:
			throw new RuntimeException("Undefined genetic mapping function.");
		}
	}

	abstract protected class RfCalculator implements Runnable {
		protected final int i;
		protected final int j;

		public RfCalculator(int i, int j) {
			this.i = i;
			this.j = j;
		}
		
	}
	
	abstract protected class PhasedDataCollection {
		protected final String file;
		protected final String[] markers;
		protected final int[] start_end_position;
		
		protected PhasedDataCollection(String file,
				String[] markers, 
				int[] start_end) {
			// TODO Auto-generated constructor stub
			this.file = file;
			this.markers = markers;
			if(markers.length!=
					Math.abs(start_end[0]-start_end[1])+1)
				throw new RuntimeException("!!!");
			this.start_end_position = start_end;
		}
	}
	
	protected class InputStreamObj { 
		protected ZipFile in = null;
		protected InputStream is = null;
		
		public InputStreamObj(String root) {
			try {
				this.in = new ZipFile(root);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		public void closeIn() {
			// TODO Auto-generated method stub
			try {
				if(this.in != null)
					this.in.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		public void closeIs() {
			// TODO Auto-generated method stub
			try {
				if(this.is != null)
					this.is.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		public void close() {
			// TODO Auto-generated method stub
			this.closeIs();
			this.closeIn();
		}
		
		public boolean getInputStream(String file) {
			// TODO Auto-generated method stub
			try {
				String target = null;
				switch(file.toUpperCase()) {
				case "PHASEDSTATES":
					target = "phasedStates/"+expr_id+".txt";
					break;
				case "STDERR":
					target = "stderr_true";
					break;
				case "EMISS":
					target = "results_hmm/emissionModel.txt";
					break;
				case "TRANS":
					target = "results_hmm/transitionModel.txt";
					break;
				case "SNP":
					target = "snp_"+expr_id+".txt";
					break;
				case "RESAMPLING":
					target = "resampling/"+expr_id+".txt";
				}
				if(in.getEntry(target)==null) {
					in.close();
					return false;
				}
				this.is = in.getInputStream(in.getEntry(target));
				return true;
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			return false;
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
		nF1 = 0;

		for(File file:listFiles) {
			String name = file.getName();
			if( name.startsWith(expr_id) ) {
				if(nF1<1) {
					try {
						InputStreamObj isObj = new InputStreamObj(
								file.getAbsolutePath());
						if( isObj.getInputStream("phasedStates") ) {
							BufferedReader br = Utils.getBufferedReader(isObj.is);
							int n = 0;
							String l;

							while( (l=br.readLine())!=null ) 
								if(l.startsWith("#")) n++;

							nF1 = (n/Constants._ploidy_H)-2;
							myLogger.info(nF1+" F1 samples in the experiment.");
							br.close();
						}
						isObj.close();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			if(nF1>0) break;
		}

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
					map.get(s[1]).add(
							file.getAbsolutePath());
				}
			}
		}
		
		String[] keys = new String[map.keySet().size()];
		map.keySet().toArray(keys);

		this.initial_thread_pool();
		for(int i=0; i<keys.length; i++) {
			List<String> files = map.get(keys[i]);
			executor.submit(new FileExtraction(
					files.toArray(new String[files.size()])));
		}
		this.waitFor();
		
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
}
