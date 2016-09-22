package cz1.model;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.math3.stat.StatUtils;

import cz1.util.IO;

public class GenomeAssemblyComparison {

	public static void main(String[] args) throws IOException, InterruptedException {
		CommandLineParser parser = new PosixParser();
		//gac.colinear("C:\\Users\\chenxi.zhou\\Desktop\\putty\\rdot_samples",
		//		"C:\\Users\\chenxi.zhou\\Desktop\\putty\\chromosome.sizes.txt",
		//		"C:\\Users\\chenxi.zhou\\Desktop\\putty\\itr.bwa2m.b30.nnj.ge3.txt");

		// create the Options
		Options options = new Options();
		options.addOption( "r", "rdot-plot", true, "rdot plot file directory." );
		options.addOption( "c", "chrom-size", true, "chromosomes size file." );
		options.addOption( "g", "genetic-map", true, "genetic map file" );
		options.addOption( "o", "output", true, "output file.");
		options.addOption( "f", "rf",true, "recombination fraction estimation file." );
		options.addOption( "t", "thread", true, "number threads.");
		String rdot=null, chromS=null, gm=null, output=null, 
				rf=null;
		int t=1;
		try {
			// parse the command line arguments
			CommandLine line = parser.parse( options, args );
			if( line.hasOption("r") ) {
				rdot = line.getOptionValue('r');
			}
			if(line.hasOption("c")) {
				chromS = line.getOptionValue("c");
			}
			if(line.hasOption("o")) {
				output = line.getOptionValue("o");
			}
			if(line.hasOption("g")) {
				gm = line.getOptionValue("g");
			}
			if(line.hasOption("f")) {
				rf = line.getOptionValue("f");
			}
			if(line.hasOption("t")) {
				t = Integer.parseInt(line.getOptionValue("t"));
			}
		}
		catch( ParseException exp ) {
			System.out.println( "Unexpected exception:" + exp.getMessage() );
		}

		
		GenomeAssemblyComparison gac = new GenomeAssemblyComparison(t);
		//gac.test();
		//gac.colinear("C:\\Users\\chenxi.zhou\\Desktop\\putty\\rdot_samples",
		//		"C:\\Users\\chenxi.zhou\\Desktop\\putty\\chromosome.sizes.txt",
		//		"C:\\Users\\chenxi.zhou\\Desktop\\putty\\itr.bwa2m.b30.nnj.ge3.txt");
		gac.colinear(rdot, chromS, gm, rf, output);
	}

	private void test() {
		// TODO Auto-generated method stub
		TwoPointSegment tps1 = new TwoPointSegment(
				new double[]{2.0,2.0},
				new double[]{3.0,3.0}),
				tps2 = new TwoPointSegment(
						new double[]{6.0,5.0},
						new double[]{5.0,4.0});
		System.out.println(tps1.distance(tps2));
	}

	private final String refAssembly;
	private final String altAssembly;
	private static BufferedWriter cpWriter;
	private static BufferedWriter dsWriter;
	private static BufferedWriter clWriter;
	private static int THREADS = 1;
	
	public GenomeAssemblyComparison() {
		this.refAssembly = "Trifida";
		this.altAssembly = "Itr";
		THREADS = 1;
	}

	public GenomeAssemblyComparison(int t) {
		this.refAssembly = "Trifida";
		this.altAssembly = "Itr";
		THREADS = t;
	}
	
	private class Colinear implements Runnable {
		private final String ref;
		
		public Colinear(String ref) {
			this.ref = ref;
		}
		
		@Override
		public void run() {
			// TODO Auto-generated method stub
			BufferedReader br;
			String line;
			String[] s;

			Set<String> alts = compare.get(ref);
			List<TwoPointSegment> altRefMap = new ArrayList<TwoPointSegment>();
			for(String alt : alts) {
				try {
					br = IO.getBufferedReader(inputDir+"/"+getFile(ref, alt));
					br.readLine();
					List<TwoPointSegment> segs = new ArrayList<TwoPointSegment>();
					while( (line=br.readLine())!=null ) {
						s = line.split("\\s+");
						double[] p1 = new double[]{Double.parseDouble(s[0]),
								Double.parseDouble(s[1])};
						line = br.readLine();
						s = line.split("\\s+");
						double[] p2 = new double[]{Double.parseDouble(s[0]),
								Double.parseDouble(s[1])};
						br.readLine();
						if(segs.size()==0) {
							segs.add(new TwoPointSegment(p1, p2));
							continue;
						}
						TwoPointSegment tps = new TwoPointSegment(p1, p2),
								tpsL = segs.get(segs.size()-1);
						if( tpsL.direction(tps) &&
								tpsL.distance(tps)<altAssSize.get(alt))
							segs.add(tps);
						else break;
					}
					br.close();

					if(segs.size()==0) continue;
					double l = 0.0, L = altAssSize.get(alt);
					for(int i=0; i<segs.size(); i++) 
						l += segs.get(i).L2();
					double[] refP = segs.get(0).direction1==1 ? 
							new double[] {segs.get(0).p1[0], segs.get(segs.size()-1).p2[0]} : 
								new double[] {segs.get(segs.size()-1).p2[0], segs.get(0).p1[0]};
					double[] altP = segs.get(0).direction1==1 ? 
							new double[] {segs.get(0).p1[1], segs.get(segs.size()-1).p2[1]} : 
								new double[] {segs.get(segs.size()-1).p2[1], segs.get(0).p1[1]};
					if(refP[0]<=1000 || refP[1]>=refAssSize.get(ref)-1000)
							L = Math.max(refP[1]-refP[0], L/2);
					if( l/L>=.7) 
					//System.out.println(ref+"\t"+alt);
						altRefMap.add( new TwoPointSegment(new double[]{refP[0], altP[0]},
							new double[]{refP[1], altP[1]}, alt) );
				} catch (IOException e) {
					e.printStackTrace();
					System.exit(1);
				}
			}
			Collections.sort(altRefMap);
			//System.out.println("########## "+ref);
			StringBuilder cl = new StringBuilder("########## "+ref+"\n");
			for(int i=0; i<altRefMap.size(); i++) {
				cl.append(altRefMap.get(i).os());
				String a = allAltRefMap.get(altRefMap.get(i).alt);
				if(a==null) a="";
				else a=a+"/";
				allAltRefMap.put(altRefMap.get(i).alt, 
						a+ref+"["+(i+1)+"]");
			}
			
			StringBuilder ds = new StringBuilder();
			TwoPointSegment ii, jj;
			for(int i=0; i<altRefMap.size(); i++) {
				ii = altRefMap.get(i);
				for(int j=i+1; j<altRefMap.size(); j++) {
					jj = altRefMap.get(j);
					ds.append(ref+
							"\t"+ii.alt+
							"\t"+jj.alt+
							"\t"+ii.p1[0]+"\t"+ii.p2[0]+
							"\t"+jj.p1[0]+"\t"+jj.p2[0]+
							"\t"+Math.abs(jj.p1[0]-ii.p2[0])+
							"\t"+rfMap.get(ii.alt+"_"+jj.alt)+
							"\n");
				}
			}
			
			try {
				dsWriter.write(ds.toString());
				clWriter.write(cl.toString());
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		}
	}
	
	private static Map<String, Set<String>> compare;
	private static Map<String, Integer> refAssSize;
	private static Map<String, Integer> altAssSize;
	private static String inputDir;
	private static Map<String, String> allAltRefMap;
	private static Map<String, Double> rfMap;
	private static Set<String> contigs;
	private static ExecutorService executor = 
			Executors.newFixedThreadPool(THREADS);
	
	private static void reset() {
		executor = Executors.newFixedThreadPool(THREADS);
	}
	
	public boolean colinear(String input,
			String sizes,
			String linkagemap,
			String rf,
			String output) throws IOException, InterruptedException {
		
		inputDir = input;
		String line;
		String[] s, s0;
		BufferedReader br0 = IO.getBufferedReader(sizes);
		refAssSize = new HashMap<String, Integer>();
		altAssSize = new HashMap<String, Integer>();
		while( (line=br0.readLine())!=null ) {
			s = line.split("\\s+");
			s0 = s[0].split("\\.");
			if(s0[0].equals(refAssembly))
				refAssSize.put(s0[1], Integer.parseInt(s[1]));
			if(s0[0].equals(altAssembly))
				altAssSize.put(s0[1], Integer.parseInt(s[1]));
		}
		br0.close();

		BufferedReader br = IO.getBufferedReader(rf);
		rfMap = new ConcurrentHashMap<String, Double>();
		while( (line=br.readLine())!=null ) {
			s = line.split("\\s+");
			rfMap.put(s[6]+"_"+s[7], Double.parseDouble(s[1]));
			rfMap.put(s[7]+"_"+s[6], Double.parseDouble(s[1]));
		}
		br.close();
		
		contigs = new HashSet<String>();
		br = IO.getBufferedReader(linkagemap);
		while( (line=br.readLine())!=null ) {
			if(line.equals("$contigs")) {
				while( (line=br.readLine())!=null &&
						line.length()!=0)
					contigs.add(line);
				break;
			}
		}
		br.close();
		
		compare = new HashMap<String, Set<String>>();
		File in = new File(inputDir);
		for(File f : in.listFiles()) {
			String n = f.getName();
			if(n.endsWith("chain.rdotplot")) {
				s = n.split("\\.");
				if(!contigs.contains(s[5])) continue;
				if(compare.get(s[1])==null) {
					Set<String> set = new HashSet<String>();
					set.add(s[5]);
					compare.put(s[1], set);
				} else
					compare.get(s[1]).add(s[5]);
			}
		}

		allAltRefMap = new ConcurrentHashMap<String, String>();
		
		dsWriter = IO.getBufferedWriter(output+".ds");
		clWriter = IO.getBufferedWriter(output+".cl");
		reset();
		for(String ref : compare.keySet()) 
			executor.submit(new Colinear(ref));
		executor.shutdown();
		executor.awaitTermination(365, TimeUnit.DAYS);
		dsWriter.close();
		clWriter.close();
		
		for(String alt : allAltRefMap.keySet())
			System.out.println(alt+"\t"+allAltRefMap.get(alt));
		br = IO.getBufferedReader(linkagemap);
		List<String> lg = new ArrayList<String>();
		while( (line=br.readLine())!=null ) {
			if(line.equals("$order")) {
				while( (line=br.readLine())!=null &&
						line.length()!=0) {
					lg.add(line.split("\\s+")[1]);
				}
				break;
			}
		}
		br.close();
		
		cpWriter = IO.getBufferedWriter(output+".lg");
		for(int i=0; i<lg.size(); i++) {
			s = lg.get(i).split("-");
			String a = allAltRefMap.get(s[0]);
			String newLg = (a==null ? s[0] : a);
			for(int j=1; j<s.length; j++) {
				a = allAltRefMap.get(s[j]);
				newLg += "-"+ (a==null ? s[j] : a);
			}
			cpWriter.write("LG"+(i+1)+"\t"+newLg+"\n");
		}
		cpWriter.close();
		
		return true;
	}

	private String getFile(String ref, String alt) {
		// TODO Auto-generated method stub
		return refAssembly+"."+ref+".fa.vs."+altAssembly+"."+alt+".fa.chain.rdotplot";
	}

	private class TwoPointSegment implements Comparable {
		private final double[] p1, p2;
		private final int direction1, direction2;
		private final String ref, alt;
		
		public TwoPointSegment(double[] p1, 
				double[] p2) {
			this.p1 = p1;
			this.p2 = p2; 
			this.direction1 = this.p1[0]<this.p2[0] ? 
					1 : -1;
			this.direction2 = this.p1[1]<this.p2[1] ? 
					1 : -1;
			this.ref = null;
			this.alt = null;
		}
		
		public void print() {
			// TODO Auto-generated method stub
			System.out.println(alt+"\t"+p1[0]+"\t"+p2[0]+"\t"+p1[1]+"\t"+p2[1]);
		}
		
		public String os() {
			// TODO Auto-generated method stub
			return alt+"\t"+p1[0]+"\t"+p2[0]+"\t"+p1[1]+"\t"+p2[1]+"\n";
		}

		public TwoPointSegment(double[] p1, 
				double[] p2, String alt) {
			this.p1 = p1;
			this.p2 = p2; 
			this.direction1 = this.p1[0]<this.p2[0] ? 
					1 : -1;
			this.direction2 = this.p1[1]<this.p2[1] ? 
					1 : -1;
			this.ref = null;
			this.alt = alt;
		}

		public TwoPointSegment(int[] p11, int[] p22) {
			double[] p1 = new double[2], 
					p2 = new double[2];
			for(int i=0; i<2; i++) {
				p1[i] = p11[i];
				p2[i] = p22[i];
			}
			this.p1 = p1;
			this.p2 = p2;
			this.direction1 = this.p1[0]<this.p2[0] ? 
					1 : -1;
			this.direction2 = this.p1[1]<this.p2[1] ? 
					1 : -1;
			this.ref = null;
			this.alt = null;
		}

		public double L1() {
			return Math.abs(
					this.p1[0]-this.p2[0]);
		}

		public double L2() {
			return Math.abs(
					this.p1[1]-this.p2[1]);
		}

		public double slope() {
			return this.L1()/this.L2();
		}

		public boolean direction(TwoPointSegment tpsL) {
			// TODO Auto-generated method stub
			if( this.direction1 != tpsL.direction1 || 
					this.direction2 != tpsL.direction2 ||
					(tpsL.p1[0]-this.p2[0])*this.direction1<0 ||
					(tpsL.p1[1]-this.p2[1])*this.direction2<0 )
				return false;
			return true;
		}

		public double distance(TwoPointSegment tps) {
			if(this.intersect(tps)) return 0;
			double[] allD = new double[4];
			allD[0] = this.distance(tps.p1);
			allD[1] = this.distance(tps.p2);
			allD[2] = tps.distance(this.p1);
			allD[3] = tps.distance(this.p2);
			return StatUtils.min(allD);
		}

		private boolean intersect(TwoPointSegment tps) {
			// TODO Auto-generated method stub
			double dx = this.p2[0]-this.p1[0],
					dy = this.p2[1]-this.p1[1],
					da = tps.p2[0]-tps.p1[0],
					db = tps.p2[1]-tps.p1[1];
			if(da*dy-db*dx==0) return false;
			double s = (dx*(tps.p1[1]-this.p1[1])+dy*(this.p1[0]-tps.p1[0]))/(da*dy-db*dx),
					t = (da *(this.p1[1]-tps.p1[1])+db*(tps.p1[0]-this.p1[0]))/(db*dx-da*dy);	
			return s>=0&&s<=1&&t>=0&&t<=1;
		}

		private double distance(double[] point) {
			if(this.isdot())
				return distance(this.p1, point);
			double px = point[0], py = point[1],
					X1 = this.p1[0], Y1 = this.p1[1],
					X2 = this.p2[0], Y2 = this.p2[1],
					dx = X2-X1, dy = Y2-Y1;
			double t = ((px-X1)*dx+(py-Y1)*dy)/(dx*dx+dy*dy);
			if(t<0) {
				dx = px-X1;
				dy = py-Y1;
			} else if(t>1) {
				dx = px-X2;
				dy = py-Y2;
			} else {
				dx = px-X1-t*dx;
				dy = py-Y1-t*dy;
			}
			return Math.sqrt(dx*dx+dy*dy);
		}

		private double distance(double[] p1, double[] p2) {
			return Math.sqrt(
					Math.pow(p1[0]-p2[0],2)+
					Math.pow(p1[1]-p2[1],2));
		}

		private boolean isdot() {
			return this.p1[0]==this.p2[0] &&
					this.p1[1]==this.p2[1];
		}
		
		@Override
		public int compareTo(Object obj) {
			// TODO Auto-generated method stub
			return (this.p1[0]-((TwoPointSegment)obj).p1[0])>0 ? 1 : -1;
		}
	}
}
