package cz1.hmm.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
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
import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import org.apache.commons.math3.stat.inference.GTest;

import cz1.util.Algebra;
import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Utils;

public class AssemblyError extends RFUtils {
	
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
		
	}

	
	public AssemblyError (String in_haps, 
			String expr_id, 
			int ploidy, 
			String[] founder_haps,
			int threads,
			double skew_phi,
			int drop_thres,
			int best_n) {
		this.in_haps = in_haps;
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
	
	private void calcMisScaffold(String outputFilePath) throws IOException, InterruptedException {
		
		String[] s;
		BufferedReader br;
		File folder = new File(in_haps);
		//File folder = new File("C:\\Users\\chenxi.zhou\\Desktop\\console out");
		File[] listFiles = folder.listFiles();
		nF1 = 0;

		for(File file:listFiles) {
			String name = file.getName();
			if(file.isDirectory() &&
					name.startsWith(expr_id)) {
				if(nF1<1) {
					String phasedStates = file.getAbsolutePath()+
							"/phasedStates/"+expr_id+".txt";
					if(new File(phasedStates).exists()) {
						br = Utils.getBufferedReader(phasedStates);
						int n = 0;
						String l;
						while( (l=br.readLine())!=null ) 
							if(l.startsWith("#")) n++;
						nF1 = (n/Constants._ploidy_H)-2;
						System.out.println(nF1+" F1 samples in the experiment.");
						br.close();
					}
				}
			}
			if(nF1>0) break;
		}

		Map<String, List<File>> map = new HashMap<String, 
				List<File>>();
		List<File> list;
		for(File file:listFiles) {
			String name = file.getName();
			if(file.isDirectory() &&
					name.startsWith(expr_id)) {
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

		dc = new String[map.keySet().size()][best_n];
		int ii = 0;
		for(String contig : map.keySet()) {
			List<File> files = map.get(contig);
			executor.submit(new FileLoader(contig,
					files.toArray(new File[files.size()]),
					ii++));
		}
		executor.shutdown();
		executor.awaitTermination(365, TimeUnit.DAYS);

		System.out.println(map.keySet().size());
		System.out.println("["+Utils.getSystemTime()+"] LOADING FILES DONE.");
		System.out.println("["+Utils.getSystemTime()+"] READING LOG LIKELIHOOD DONE.");
		
		this.initial_thread_pool();
		mapWriter = Utils.getBufferedWriter(outputFilePath+".map");
		for(int i=0; i<dc.length; i++) 
			executor.submit(new mapCalculator(i));
		executor.shutdown();
		executor.awaitTermination(365, TimeUnit.DAYS);
		mapWriter.close();

		System.out.println("["+Utils.getSystemTime()+"] DONE.");
	}
}
