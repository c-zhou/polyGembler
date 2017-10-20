package cz1.ngs.assembly;

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

import org.apache.commons.math3.stat.StatUtils;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

public class GenomeAssemblyTiling extends Executor {

	private String rdot_file = null;
	private String size_file = null;
	private String rf_file = null;
	private double cl_thres = 0.7;
	private String out_prefix = "out";
	private String refAssembly = null;
	private String qryAssembly = null;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -d/--rdot-plot				Rdot plot file directory.\n"
						+ " -s/--size                   Chromosome/contig/scaffold size file.\n"
						+ " -f/--rf                     Recombination frequency file.\n"
						+ " -r/--reference              Reference genome assembly ID.\n"
						+ " -q/--query                  Query genome assembly ID.\n"
						+ " -l/--colinear-thres         Threshold of the colinear part (Default 0.7). \n"
						+ " -t/--threads                Threads (default 1).\n"
						+ " -o/--out-prefix             Output file prefix (default \'out\').\n\n"
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
			myArgsEngine.add( "-d", "--rdot-plot", true );
			myArgsEngine.add( "-s", "--size", true );
			myArgsEngine.add( "-f", "--rf",true );
			myArgsEngine.add( "-r", "--reference",true );
			myArgsEngine.add( "-q", "--query",true );
			myArgsEngine.add( "-l", "--colinear-thres", true );
			myArgsEngine.add( "-t", "--threads", true );
			myArgsEngine.add( "-o", "--out-prefix", true );
			myArgsEngine.parse(args);
		}
		
		if(myArgsEngine.getBoolean("-d")) {
			rdot_file = myArgsEngine.getString("-d");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the location of the Rdor-plot files.");
		}
		
		if(myArgsEngine.getBoolean("-s")) {
			size_file = myArgsEngine.getString("-s");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the chromosome/scaffold/contig size file.");
		}
		
		if(myArgsEngine.getBoolean("-f")) {
			rf_file = myArgsEngine.getString("-f");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the recombination frequency file.");
		}
		
		if(myArgsEngine.getBoolean("-r")) {
			refAssembly = myArgsEngine.getString("-r");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the reference id.");
		}
		
		if(myArgsEngine.getBoolean("-q")) {
			qryAssembly = myArgsEngine.getString("-q");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the query id.");
		}
		
		if(myArgsEngine.getBoolean("-l")) {
			cl_thres = Double.parseDouble(myArgsEngine.getString("-l"));
		}
		
		if(myArgsEngine.getBoolean("-o")) {
			out_prefix = myArgsEngine.getString("-o");
		}
		
		if(myArgsEngine.getBoolean("-t")) {
			THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}
	}

	private static Map<String, Set<String>> compare;
	private static Map<String, Integer> refAssSize;
	private static Map<String, Integer> qryAssSize;
	private static Map<String, String> allQryRefMap;
	private static Map<String, Double> rfMap;
	private static Set<String> contigs;
    private static BufferedWriter dsWriter;
    private static BufferedWriter clWriter;

	@Override
	public void run() {
		// TODO Auto-generated method stub
        String line;
        String[] s, s0;
        BufferedReader br0 = Utils.getBufferedReader(size_file);
        refAssSize = new HashMap<String, Integer>();
        qryAssSize = new HashMap<String, Integer>();
        rfMap = new ConcurrentHashMap<String, Double>();
        contigs = new HashSet<String>();
        compare = new HashMap<String, Set<String>>();
        allQryRefMap = new ConcurrentHashMap<String, String>();
        
        try {
        	while( (line=br0.readLine())!=null ) {
        		s = line.split("\\s+");
        		s0 = s[0].split("\\.");
        		if(s0[0].equals(refAssembly))
        			refAssSize.put(s0[1], Integer.parseInt(s[1]));
        		if(s0[0].equals(qryAssembly))
        			qryAssSize.put(s0[1], Integer.parseInt(s[1]));
        	}
        	br0.close();

        	BufferedReader br = Utils.getBufferedReader(rf_file);
        	while( (line=br.readLine())!=null ) {
        		if(line.startsWith("#")) continue;
        		s = line.split("\\s+");
        		contigs.add(s[5]);
        		contigs.add(s[6]);
        		rfMap.put(s[5]+"_"+s[6], Double.parseDouble(s[0]));
        		rfMap.put(s[6]+"_"+s[5], Double.parseDouble(s[0]));
        	}
        	br.close();
        	
        	File in = new File(rdot_file);
            for(File f : in.listFiles()) {
                String n = f.getName();
                if(n.endsWith("chain.rdotplot") ||
                    n.endsWith("chain.rdotplot.gz")) {
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
            
            dsWriter = Utils.getBufferedWriter(out_prefix+".ds");
            clWriter = Utils.getBufferedWriter(out_prefix+".cl");
            this.initial_thread_pool();
            for(String ref : compare.keySet())
                executor.submit(new Colinear(ref));
            this.waitFor();            
            dsWriter.close();
            clWriter.close();
        } catch (IOException e) {
        	e.printStackTrace();
        	System.exit(1);
        }
	}
	
    private static final Object lock = new Object();
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

            Set<String> qrys = compare.get(ref);
            List<TwoPointSegment> qryRefMap = new ArrayList<TwoPointSegment>();
            for(String qry : qrys) {
                try {
                    br = Utils.getBufferedReader(getFile(ref, qry, rdot_file));
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
                                tpsL.distance(tps)<qryAssSize.get(qry))
                            segs.add(tps);
                        else break;
                    }
                    br.close();

                    if(segs.size()==0) continue;
                    double l = 0.0, L = qryAssSize.get(qry);
                    for(int i=0; i<segs.size(); i++)
                        l += segs.get(i).L2();
                    double[] refP = segs.get(0).direction1==1 ?
                            new double[] {segs.get(0).p1[0], segs.get(segs.size()-1).p2[0]} :
                                new double[] {segs.get(segs.size()-1).p2[0], segs.get(0).p1[0]};
                    double[] qryP = segs.get(0).direction1==1 ?
                            new double[] {segs.get(0).p1[1], segs.get(segs.size()-1).p2[1]} :
                                new double[] {segs.get(segs.size()-1).p2[1], segs.get(0).p1[1]};
                    if(refP[0]<=1000 || refP[1]>=refAssSize.get(ref)-1000)
                            L = Math.max(refP[1]-refP[0], L/2);
                    if( l/L>=cl_thres)
                    //System.out.println(ref+"\t"+qry);
                        qryRefMap.add( new TwoPointSegment(new double[]{refP[0], qryP[0]},
                            new double[]{refP[1], qryP[1]}, ref, qry, l/L) );
                } catch (IOException e) {
                    e.printStackTrace();
                    System.exit(1);
                }
            }
            Collections.sort(qryRefMap);
            //System.out.println("########## "+ref);
            synchronized(lock) {
                StringBuilder cl = new StringBuilder("########## "+ref+"\n");
                for(int i=0; i<qryRefMap.size(); i++) {
                    cl.append(qryRefMap.get(i).os());
                    String a = allQryRefMap.get(qryRefMap.get(i).qry);
                    if(a==null) a="";
                    else a=a+"/";
                    allQryRefMap.put(qryRefMap.get(i).qry,
                            a+ref+"["+(i+1)+"]");
                }

                StringBuilder ds = new StringBuilder();
                TwoPointSegment ii, jj;
                for(int i=0; i<qryRefMap.size(); i++) {
                    ii = qryRefMap.get(i);
                    for(int j=i+1; j<qryRefMap.size(); j++) {
                        jj = qryRefMap.get(j);
                        ds.append(ref+
                            "\t"+ii.qry+
                            "\t"+jj.qry+
                            "\t"+ii.p1[0]+"\t"+ii.p2[0]+
                            "\t"+jj.p1[0]+"\t"+jj.p2[0]+
                            "\t"+Math.abs(jj.p1[0]-ii.p2[0])+
                            "\t"+rfMap.get(ii.qry+"_"+jj.qry)+
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
    }
    
    private String getFile(String ref, String qry, String dir) {
        // TODO Auto-generated method stub
        String f = dir+"/"+refAssembly+"."+ref+".fa.vs."+qryAssembly+"."+qry+".fa.chain.rdotplot";
        return new File(f).exists() ? f :
            dir+"/"+refAssembly+"."+ref+".fa.vs."+qryAssembly+"."+qry+".fa.chain.rdotplot.gz";
    }

    private class TwoPointSegment implements Comparable<TwoPointSegment> {
        private final double[] p1, p2;
        private final int direction1, direction2;
        private final String ref, qry;
        private final double frac;

        public TwoPointSegment(double[] p1,
                double[] p2) {
            this.p1 = p1;
            this.p2 = p2;
            this.direction1 = this.p1[0]<this.p2[0] ?
                    1 : -1;
            this.direction2 = this.p1[1]<this.p2[1] ?
                    1 : -1;
            this.ref = null;
            this.qry = null;

            this.frac = -1;
        }

        public void print() {
            // TODO Auto-generated method stub
            System.out.println(ref+"\t"+qry+"\t"+p1[0]+"\t"+p2[0]+"\t"+p1[1]+"\t"+p2[1]);
        }

        public String os() {
            // TODO Auto-generated method stub
            return ref+"\t"+qry+"\t"+p1[0]+"\t"+p2[0]+"\t"+p1[1]+"\t"+p2[1]+"\t"+frac+"\n";
        }

        public TwoPointSegment(double[] p1,
                double[] p2, String qry) {
            this.p1 = p1;
            this.p2 = p2;
            this.direction1 = this.p1[0]<this.p2[0] ?
                    1 : -1;
            this.direction2 = this.p1[1]<this.p2[1] ?
                    1 : -1;
            this.ref = null;
            this.qry = qry;

            this.frac = -1;
        }

        public TwoPointSegment(double[] p1,
                                double[] p2, String qry, double frac) {
                        this.p1 = p1;
                        this.p2 = p2;
                        this.direction1 = this.p1[0]<this.p2[0] ?
                                        1 : -1;
                        this.direction2 = this.p1[1]<this.p2[1] ?
                                        1 : -1;
                        this.ref = null;
                        this.qry = qry;

                        this.frac = frac;
        }
        
        public TwoPointSegment(double[] p1,
        		double[] p2, String ref, 
        		String qry, double frac) {
        	this.p1 = p1;
        	this.p2 = p2;
        	this.direction1 = this.p1[0]<this.p2[0] ?
        			1 : -1;
        	this.direction2 = this.p1[1]<this.p2[1] ?
        			1 : -1;
        	this.ref = ref;
        	this.qry = qry;

        	this.frac = frac;
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
            this.qry = null;

            this.frac = -1;
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
        public int compareTo(TwoPointSegment obj) {
            // TODO Auto-generated method stub
            return (this.p1[0]-obj.p1[0])>0 ? 1 : -1;
        }
    }
}
