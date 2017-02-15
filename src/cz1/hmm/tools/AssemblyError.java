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

import cz1.hmm.tools.Copy_11_of_SingleCNVHapRF.PhasedDataCollection;
import cz1.hmm.tools.RFUtils.FileExtraction;
import cz1.hmm.tools.RFUtils.FileLoader;
import cz1.hmm.tools.RFUtils.FileObject;
import cz1.util.Algebra;
import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Utils;

public class AssemblyError extends RFUtils {
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
			myArgsEngine.add( "-ex", "--experiment-id", true);
			myArgsEngine.add( "-i", "--hap-file", true);
			myArgsEngine.add( "-o", "--prefix", true);
			myArgsEngine.add( "-f", "--parent", true);
			myArgsEngine.add( "-p", "--ploidy", true);
			myArgsEngine.add( "-nb", "--best", true);
			myArgsEngine.add( "-t", "--thread", true);
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
			Constants._founder_haps = myArgsEngine.getString("-f");
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

	private String guessExperimentId() {
		// TODO Auto-generated method stub
		File in_dir = new File(in_haps);
		File[] haps = in_dir.listFiles();
		Map<String, Integer> stats = new HashMap<String, Integer>();
		for(File hap : haps) {
			String h = hap.getName().split(".")[0];
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

	@Override
	public void run() {
		// TODO Auto-generated method stub
		File folder = new File(in_haps);
		File[] listFiles = folder.listFiles();
		nF1 = 0;

		for(File file:listFiles) {
			String name = file.getName();
			if( name.startsWith(expr_id) ) {
				if(nF1<1) {
					try {
						InputStream is = this.getInputStream(
								file.getAbsolutePath(),
								"phasedStates");
						if( is!=null ) {
							BufferedReader br = Utils.getBufferedReader(is);
							int n = 0;
							String l;

							while( (l=br.readLine())!=null ) 
								if(l.startsWith("#")) n++;

							nF1 = (n/Constants._ploidy_H)-2;
							myLogger.info(nF1+" F1 samples in the experiment.");
							br.close();
						}
						is.close();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			if(nF1>0) break;
		}

		Map<String, List<File>> map = new HashMap<String, List<File>>();
		List<File> list;
		String[] s;
		for(File file:listFiles) {
			String name = file.getName();
			if( name.startsWith(expr_id) ) {
				name = name.replace(expr_id,"experiment");
				s = name.split("\\.");

				if(map.get(s[1])==null) {
					list = new ArrayList<File>();
					list.add(file);
					map.put(s[1], list);
				} else{
					map.get(s[1]).add(file);
				}
			}
		}

		String[] keys = new String[map.keySet().size()];
		map.keySet().toArray(keys);

		this.initial_thread_pool();
		for(int i=0; i<keys.length; i++) {
			List<File> files = map.get(keys[i]);
			executor.submit(new FileExtraction(
					files.toArray(new File[files.size()])));
		}
		this.waitFor();
		
		this.initial_thread_pool();
		String[] scaff_all = new String[fileObj.keySet().size()];
		fileObj.keySet().toArray(scaff_all);
		dc = new PhasedDataCollection[scaff_all.length][best_n];
		for(int i=0; i<scaff_all.length; i++) {
			Set<FileObject> files = fileObj.get(scaff_all[i]);
			executor.submit(new FileLoader(scaff_all[i],
					files.toArray(new FileObject[files.size()]),
					i));
		}
		this.waitFor();
		
		myLogger.info(map.keySet().size());
		myLogger.info("["+Utils.getSystemTime()+"] LOADING FILES DONE.");
		myLogger.info("["+Utils.getSystemTime()+"] READING LOG LIKELIHOOD DONE.");
		
		this.initial_thread_pool();
		for(int i=0; i<dc.length; i++) 
			executor.submit(new mapCalculator(i));
		this.waitFor();

		write();
		
		render();
		
		myLogger.info("["+Utils.getSystemTime()+"] DONE.");
	}

	private void write() {
		// TODO Auto-generated method stub
		BufferedWriter bw = Utils.getBufferedWriter(out_prefix);
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
	
	final private Map<String, double[]> errs = new HashMap<String, double[]>();
	
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
			final List<Double> breakage_pos = new ArrayList<Double>();
			for(int j=0; j<rf.length; j++) {
 				rf[j] /= bound;
 				if(rf[j]>=breakage_thres) {
 					String[] s = markers[j].split("_");
 					double x = Double.parseDouble(s[s.length-1]);
 					s = markers[j+1].split("_");
 					x += Double.parseDouble(s[s.length-1]);
 					breakage_pos.add(x/2); 
 				}
 			}
			if(breakage_pos.size()==0) continue;
			Double[] ls = new Double[breakage_pos.size()];
			breakage_pos.toArray(ls);
			errs.put(scaff, ArrayUtils.toPrimitive(ls));
		}
		return;
	}
	
	public Map<String, double[]> errs() {
		return this.errs;
	}

	public Set<String> split(String in_vcf) {
		// TODO Auto-generated method stub
		final Set<String> scaff_breakge = new HashSet<String>();
		try {
			String command = "mv "+in_vcf+" "+in_vcf+".bak";
			Process p = this.bash(command);
			p.waitFor();

			BufferedReader br = Utils.getBufferedReader(in_vcf+".bak");
			BufferedWriter bw = Utils.getBufferedWriter(in_vcf);
			 
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
					scaff_breakge.add(scaff);
					double[] tmp = this.errs.get(scaff);
					double[] breakage_pos = new double[tmp.length+1];
					System.arraycopy(tmp, 0, breakage_pos, 0, tmp.length);
					breakage_pos[tmp.length] = Double.POSITIVE_INFINITY;
					int sub = 1;
					while( line!=null ) {
						s = line.split("\\s+");
						if( !scaff.equals(s[0]) ) break;
						if( Double.parseDouble(s[1])>breakage_pos[sub-1] )
							sub++;
						bw.write(line.replaceAll("$"+scaff, scaff+"_"+sub));
						line = br.readLine();
					}
				}
			}
			br.close();
			bw.close();
		} catch (InterruptedException | IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return scaff_breakge;
	}
}
