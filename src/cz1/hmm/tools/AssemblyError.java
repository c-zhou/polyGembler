package cz1.hmm.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.apache.commons.math3.stat.inference.GTest;

import cz1.hmm.tools.RFUtils.FileExtraction;
import cz1.hmm.tools.RFUtils.FileLoader;
import cz1.hmm.tools.RFUtils.FileObject;
import cz1.util.Algebra;
import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Utils;

public class AssemblyError extends RFEstimator {
	private final static double breakage_thres = 0.05;
	
	public AssemblyError (String in_haps, 
			String out_prefix,
			String expr_id, 
			int ploidy, 
			String[] founder_haps,
			int threads,
			double skew_phi,
			int drop_thres,
			int best_n) {
		this.in_haps = in_haps;
		this.out_prefix = out_prefix;
		this.expr_id = expr_id;
		Constants.ploidy(ploidy);
		this.founder_haps = founder_haps;
		this.THREADS = threads;
		this.drop_thres = drop_thres;
		this.skew_phi = skew_phi;
		this.probs_uniform = new double[ploidy*2];
		Arrays.fill(this.probs_uniform, .5/ploidy);
		this.best_n = best_n;
	}
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--hap-file				Directory with input haplotype files.\n"
						+ " -o/--prefix					Output file prefix.\n"
						+ " -ex/--experiment-id			Common prefix of haplotype files for this experiment.\n"
						+ " -f/--parent					Parent samples (seperated by a \":\").\n"
						+ " -p/--ploidy					Ploidy of genome (default 2).\n"
						+ " -nb/--best					The most likely nb haplotypes will be used (default 10).\n"
						+ " -phi/--skew-phi				For a haplotype inference, the frequencies of parental \n"
						+ "								haplotypes need to be in the interval [1/phi, phi], \n"
						+ "								otherwise will be discared (default 2).\n"
						+ " -nd/--drop					At least nd haplotype inferences are required for \n"
						+ "								a contig/scaffold to be analysed (default 1).\n"
						+ " -t/--threads				Threads (default 1).\n"	
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
			myArgsEngine.add( "-f", "--parent", true);
			myArgsEngine.add( "-p", "--ploidy", true);
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

		if(myArgsEngine.getBoolean("-f")) {
			founder_haps = myArgsEngine.getString("-f").split(":");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the parent samples (seperated by a \":\").");
		}
		
		if(myArgsEngine.getBoolean("-p")) {
			int ploidy = Integer.parseInt(myArgsEngine.getString("-p"));
			Constants.ploidy(ploidy);
			Constants._haplotype_z = ploidy*2;
			probs_uniform = new double[ploidy*2];
			Arrays.fill(probs_uniform, .5/ploidy);
		}
		
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

	@Override
	public void run() {
		// TODO Auto-generated method stub
		this.initialise();
		
		this.initial_thread_pool();
		
		for(int i=0; i<dc.length; i++) 
			executor.submit(new MyMapCalculator(i));
		this.waitFor();

		writeMap();
		
		render();
		
		myLogger.info("["+Utils.getSystemTime()+"] DONE.");
	}
	
	final private Map<String, int[][]> errs = new HashMap<String, int[][]>();
	
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
	
	public Map<String, int[][]> errs() {
		return this.errs;
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
}
