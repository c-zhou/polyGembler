package cz1.test;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPOutputStream;

import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.TreeBidiMap;
import org.renjin.eval.Context;
import org.renjin.primitives.io.serialization.RDataWriter;
import org.renjin.primitives.matrix.DoubleMatrixBuilder;
import org.renjin.sexp.DoubleArrayVector;
import org.renjin.sexp.DoubleVector;
import org.renjin.sexp.ListVector;
import org.renjin.sexp.StringArrayVector;
import org.renjin.sexp.StringVector;
import org.renjin.sexp.Vector;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

public class RenjinLOD extends Executor {

	private String[] ds_in;
	private String txt_in;
	private int n_hap;
	private String r_out;

	public static void main(String[] args) {
		RenjinLOD renjinLOD = new RenjinLOD();
		renjinLOD.setParameters(args);
		renjinLOD.run();
	}
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--ds-file                Input distance files.\n"
						+ " -o/--prefix                 Output file.\n"
						+ " -n/--haplotype              Number of haplotypes.\n"
						+ " -g/--txt-in                 Genetic linkage map log file.\n"
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
			myArgsEngine.add( "-i", "--ds-file", true);
			myArgsEngine.add( "-o", "--prefix", true);
			myArgsEngine.add( "-n", "--haplotype", true);
			myArgsEngine.add( "-g", "--txt-in", true);
			myArgsEngine.parse(args);
		}
		
		if(myArgsEngine.getBoolean("-i")) {
			ds_in = myArgsEngine.getString("-i").split(",");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your input zip file.");
		}

		if(myArgsEngine.getBoolean("-o")) {
			r_out = myArgsEngine.getString("-o");
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify your output file prefix.");
		}

		if(myArgsEngine.getBoolean("-n")) {
			n_hap = Integer.parseInt(myArgsEngine.getString("-n"));
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your input zip file.");
		}

		if(myArgsEngine.getBoolean("-g")) {
			txt_in = myArgsEngine.getString("-g");
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify your output file prefix.");
		}
	}

	public void run2() {
		// TODO Auto-generated method stub
		final BidiMap<Integer, String> scaffs = new TreeBidiMap<Integer, String>();
		final List<Double> boundary_lg = new ArrayList<Double>();
		try {
			BufferedReader br = Utils.getBufferedReader(txt_in);
			String line = br.readLine();
			String[] s;
			
			while( line!=null ) {
				if(line.startsWith("$order")) {
					int i = 0;
					while( (line=br.readLine())!=null 
							&& line.length()>0 ) {
						s = line.split("\\s+")[1].split("-");
						for(String scf : s) scaffs.put(i++, 
								scf.replaceAll("^scaffold", ""));
						boundary_lg.add((double) i);
					}
					break;
				}
				line=br.readLine();
			}
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		double[][] recomf = new double[scaffs.size()][scaffs.size()];
		double[][] lodscr = new double[scaffs.size()][scaffs.size()];
		final double constr = Math.log10(.5);
		try {
			String line;
			String[] s;
			Integer scf_i1, scf_i2;
			double f, l, NR, R;
			for(String in : ds_in) {
				BufferedReader br = Utils.getBufferedReader(in);
				while( (line=br.readLine())!=null ) {
					s = line.split("\\s+");
					scf_i1 = scaffs.getKey(s[6]);
					scf_i2 = scaffs.getKey(s[7]);
					if(scf_i1==null || scf_i2==null)
						continue;
					f = Double.parseDouble(s[1]);
					recomf[scf_i1][scf_i2] = f;
					recomf[scf_i2][scf_i1] = f;
					if(f==0) {
						l = 100.0;
					} else {
						NR = (1-f);
						R = f;
						l = Math.min(100.0, 
								n_hap*(Math.log10(NR)*NR+Math.log10(R)*R-constr));
					}
					lodscr[scf_i1][scf_i2] = l;
					lodscr[scf_i2][scf_i1] = l;
				}
				br.close();
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		StringVector scf = new StringArrayVector(scaffs.values());
	
	}
	
	private static class HSVColorMap {
		private final double[] scale;
		private double skew = 0;
		private double range = 0;
		private float saturation = 1.0f;
		private float brightness = 0.8f;
		
		public HSVColorMap(double[] scale) {
			this.scale = scale;
		}
		
		public void setSkew(double skew) {
			this.skew = skew;
			this.range = this.sigmod(0.5);
		}
		
		public double sigmod(double val) {
			return 2/(1+Math.exp(-this.skew*val))-1.0;
		}
		
		public void setSaturation(float saturation) {
			this.saturation = saturation;
		}
		
		public void setBrightness(float brightness) {
			this.brightness = brightness;
		}
		
		public Color get(double val) {
			if(val<this.scale[0] || val>this.scale[1])
				throw new RuntimeException("!!!");
			double scaled_val = val/(this.scale[1]-
					this.scale[0])-.5;
			double hue = 120.0*(1+
					(this.range<1e-12 ? 
					scaled_val+0.5 :
					this.sigmod(scaled_val)/this.range) );
			return Color.getHSBColor(
					(float) Math.min(Math.max(0, hue),240),
					this.saturation, 
					this.brightness);
		}
	}
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		final BidiMap<Integer, String> scaffs = new TreeBidiMap<Integer, String>();
		final List<Double> boundary_lg = new ArrayList<Double>();
		try {
			BufferedReader br = Utils.getBufferedReader(txt_in);
			String line = br.readLine();
			String[] s;
			
			while( line!=null ) {
				if(line.startsWith("$order")) {
					int i = 0;
					while( (line=br.readLine())!=null 
							&& line.length()>0 ) {
						s = line.split("\\s+")[1].split("-");
						for(String scf : s) scaffs.put(i++, 
								scf.replaceAll("^scaffold", ""));
						boundary_lg.add((double) i);
					}
					break;
				}
				line=br.readLine();
			}
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		int n = scaffs.size();
		double[][] recomfJ = new double[n][n];
		double[][] lodscrJ = new double[n][n];
		final double constr = Math.log10(.5);
		try {
			String line;
			String[] s;
			Integer scf_i1, scf_i2;
			double f, l, NR, R;
			for(String in : ds_in) {
				BufferedReader br = Utils.getBufferedReader(in);
				while( (line=br.readLine())!=null ) {
					s = line.split("\\s+");
					scf_i1 = scaffs.getKey(s[6]);
					scf_i2 = scaffs.getKey(s[7]);
					if(scf_i1==null || scf_i2==null)
						continue;
					f = Double.parseDouble(s[1]);
					recomfJ[scf_i1][scf_i2] = f;
					recomfJ[scf_i2][scf_i1] = f;
					if(f==0) {
						l = 100.0;
					} else {
						NR = (1-f);
						R = f;
						l = Math.min(100.0, 
								n_hap*(Math.log10(NR)*NR+Math.log10(R)*R-constr));
					}
					lodscrJ[scf_i1][scf_i2] = l;
					lodscrJ[scf_i2][scf_i1] = l;
				}
				br.close();
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		DoubleMatrixBuilder reccoo = new DoubleMatrixBuilder(n*(n-1), 5);
		int w = 0;
		for(int i=1; i<n; i++) {
			for(int j=0; j<i; j++) {
				reccoo.set(w, 0, j);
				reccoo.set(w, 1, i);
				reccoo.set(w, 2, j+1);
				reccoo.set(w, 3, i+1);
				reccoo.set(w, 4, lodscrJ[i][j]);
				w++;
				reccoo.set(w, 0, i);
				reccoo.set(w, 1, j);
				reccoo.set(w, 2, i+1);
				reccoo.set(w, 3, j+1);
				reccoo.set(w, 4, recomfJ[i][j]);
				w++;
			}
		}
		
		DoubleMatrixBuilder recomf = new DoubleMatrixBuilder(n,n);
		DoubleMatrixBuilder lodscr = new DoubleMatrixBuilder(n,n);
		for(int i=0; i<n; i++) {
			for(int j=0; j<n; j++) {
				recomf.set(i, j, recomfJ[i][j]);
				lodscr.set(i, j, lodscrJ[i][j]);
			}
		}
		
		ScriptEngineManager manager = new ScriptEngineManager();
		ScriptEngine engine = manager.getEngineByName("Renjin"); 
		if(engine == null) { 
			throw new RuntimeException("Renjin not found!!!"); 
		}
		StringVector scf = new StringArrayVector(scaffs.values());
		
		recomf.setRowNames(scf);
		recomf.setColNames(scf);
		lodscr.setRowNames(scf);
		lodscr.setColNames(scf);
		
		try {
			Context context = Context.newTopLevelContext();
			FileOutputStream fos  = new FileOutputStream(r_out);
			GZIPOutputStream zos = new GZIPOutputStream(fos);
			RDataWriter writer = new RDataWriter(context, zos);
			
			ListVector.NamedBuilder dat = new ListVector.NamedBuilder();
			dat.add("scaffs", scf);
			dat.add("boundary_lg", new DoubleArrayVector(boundary_lg));
			dat.add("recomf", recomf.build());
			dat.add("lodscr", lodscr.build());
			dat.add("reccoo", reccoo.build());
			
			writer.save(dat.build());
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private class LowerTriangularDoubleMatrix {
		private final double[] vals; 
		
		public LowerTriangularDoubleMatrix(int d) {
			this.vals = new double[d*(d+1)/2];
		}
		
		public double get(int row, int col) {
			if(row<col) {
				int tmp = row;
				row = col;
				col = tmp;
			}
			
			return this.vals[row*(row+1)/2+col];
		}
		
		public void set(int row, int col, int val) {
			if(row<col) {
				int tmp = row;
				row = col;
				col = tmp;
			}
			this.vals[row*(row+1)/2+col] = val;
		}
	}
}
