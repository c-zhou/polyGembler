package cz1.hmm.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import cz1.hmm.model.ModelReader;
import cz1.util.ArgsEngine;
import cz1.util.Utils;

public class AssemblyError extends RFUtils {

	private final static Logger myLogger = Logger.getLogger(AssemblyError.class);
	
	private String out_prefix;
	private double rf_thresh = 0.1;
	private int wbp = 30000;
	private int wnm = 30;
	private String in_vcf = null;
	private String out_vcf;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--hap-file               Directory with input haplotype files.\n"
						+ " -o/--prefix                 Output file prefix.\n"
						+ " -vi/--vcf-in                Input VCF file. If provided, will generate a VCF file contains variants from split scaffolds.\n"
						+ " -vo/--vcf-out               Output VCF file. If not provided, will use output file prefix (-o) option. \n"
						+ " -r/--rf-thresh              Recombination frequency threshold for assembly error detection (default 0.1).\n"
						+ " -wbp/-windows-bp            Window size (#basepairs) for assembly error dectection (default 30000).\n"
						+ " -wnm/-windows-nm            Window size (#markers) for assembly error dectection (default 30).\n"
						+ " -ex/--experiment-id         Common prefix of haplotype files for this experiment.\n"
						+ " -nb/--best                  The most likely nb haplotypes will be used (default 30).\n"
						+ " -phi/--skew-phi             For a haplotype inference, the frequencies of parental \n"
						+ "                             haplotypes need to be in the interval [1/phi, phi], \n"
						+ "                             otherwise will be discared (default 2).\n"
						+ " -nd/--drop                  At least nd haplotype inferences are required for \n"
						+ "                             a contig/scaffold to be analysed (default 1).\n"
						+ " -t/--threads                #threads (default 1).\n"	
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
			myArgsEngine.add( "-i", "--hap-file", true);
			myArgsEngine.add( "-o", "--prefix", true);
			myArgsEngine.add( "-vi", "--vcf-in", true);
			myArgsEngine.add( "-vo", "--vcf-out", true);
			myArgsEngine.add( "-r", "--rf_thresh", true);
			myArgsEngine.add( "-wbp", "--windows-bp", true);
			myArgsEngine.add( "-wnm", "--windows-nm", true);
			myArgsEngine.add( "-ex", "--experiment-id", true);
			myArgsEngine.add( "-nb", "--best", true);
			myArgsEngine.add( "-phi", "--skew-phi", true);
			myArgsEngine.add( "-nd", "--drop", true);
			myArgsEngine.add( "-t", "--threads", true);
			myArgsEngine.parse(args);
		}
		
		if(myArgsEngine.getBoolean("-i")) {
			in_haps = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your location for haplotype files.");
		}

		if(myArgsEngine.getBoolean("-o")) {
			out_prefix = myArgsEngine.getString("-o");
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify your output file prefix.");
		}
		
		if(myArgsEngine.getBoolean("-vi")) {
			in_vcf = myArgsEngine.getString("-vi");
			if(myArgsEngine.getBoolean("-vo")) {
				out_vcf = myArgsEngine.getString("-vo");
			} else {
				out_vcf = out_prefix+".vcf";
			}
		}
		
		if(myArgsEngine.getBoolean("-ex")) {
			expr_id = myArgsEngine.getString("-ex");
		}  else {
			expr_id = guessExperimentId();
			myLogger.warn("No experiment prefix provided, I guess it's "+expr_id+". Please\n"
					+ "specify it with -ex/--experiment-id option if it's incorrect.");
		}
		
		best_n = 30; // use as much as possible up to 30
		if(myArgsEngine.getBoolean("-nb")) {
			best_n = Integer.parseInt(myArgsEngine.getString("-nb"));
		}
		
		if(myArgsEngine.getBoolean("-r")) {
			rf_thresh = Double.parseDouble(myArgsEngine.getString("-r"));
		}
		
		if(myArgsEngine.getBoolean("-wbp")) {
			wbp = Integer.parseInt(myArgsEngine.getString("-wbp"));
		}
		
		if(myArgsEngine.getBoolean("-wnm")) {
			wnm = Integer.parseInt(myArgsEngine.getString("-wnm"));
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
	
	public AssemblyError (String in_haps, 
			String out_prefix,
			String expr_id, 
			int threads,
			double skew_phi,
			int drop_thres,
			int best_n) { 
		this.in_haps = in_haps;
		this.out_prefix = out_prefix;
		this.expr_id = expr_id;
		THREADS = threads;
		this.skew_phi = skew_phi;
		this.drop_thres = drop_thres;
		this.best_n = best_n;
	}
	
	public AssemblyError() {
		// TODO Auto-generated constructor stub
		super();
	}

	BufferedWriter rfWriter;
	Map<String, int[][]> errs = new HashMap<>();
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		super.initialise();
		
		rfWriter = Utils.getBufferedWriter(this.out_prefix+".err");
		this.initial_thread_pool();
		for(String scaff : fileObj.keySet()) 
			executor.submit(new RfCalculator(scaff));
		this.waitFor();
		try {
			rfWriter.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		if(in_vcf!=null) split(in_vcf, out_vcf);
		
		myLogger.info("["+Utils.getSystemTime()+"] DONE.");
	}
	
	private class RfCalculator implements Runnable {
		private final String scaff;
		
		public RfCalculator(String scaff) {
			// TODO Auto-generated constructor stub
			this.scaff = scaff;
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub
			try {
				List<FileObject> objs = fileObj.get(scaff);
				final int n = objs.size();
				ModelReader modelReader = new ModelReader(objs.get(0).file);
				int[] pos = modelReader.getSnpPosition(scaff);
				modelReader.close();
				int m = pos.length;
				int[] dists = new int[m-1];
				for(int i=0; i<dists.length; i++) dists[i] = pos[i+1]-pos[i];
				
				final List<Map<String, char[][]>> haplotypes = new ArrayList<>();
				double[][] jumps = new double[n][];
				for(int i=0; i<n; i++) {
					FileObject obj = objs.get(i);
					modelReader = new ModelReader(obj.file);
					Map<String, char[][]> haps = modelReader.getHaplotypeByPositionRange(obj.position, ploidy);
					modelReader.close();
					for(String f : parents) haps.remove(f);
					haplotypes.add(haps);
					
					double[] jump = new double[m-1];
					char[] h;
					int hapn = 0;
					for(char[][] hap : haps.values()) {
						if(hap[0][0]=='*') continue;
						for(int j=0; j<ploidy; j++) {
							h = hap[j];
							for(int k=0; k<m-1; k++)
								if(h[k]!=h[k+1]) 
									jump[k] += 1;
						}
						++hapn;
					}
					divide(jump, hapn*ploidy);
					jumps[i] = jump;
				}
				
				double[] means = new double[m-1];
				for(int i=0; i<n; i++)
					sum(means, jumps[i], means);
				divide(means, n);
				
				int dist, kb_upstream, kb_downstream;
				Map<String, char[][]> haps;
				boolean iserr;
				double min;
				List<int[]> err = new ArrayList<int[]>();
				for(int j=0; j<m-1; j++) {
					if(means[j]<=rf_thresh) continue;
					
					kb_upstream = Math.max(0, j-wnm);
					dist = 0;
					for(int k=j; k>kb_upstream&&dist<wbp; k--)
						dist += dists[k-1];
					while(dist<wbp&&kb_upstream>0) {
						dist += dists[kb_upstream-1];
						--kb_upstream;
					}
					kb_downstream = Math.min(m-1, j+wnm);
					dist = 0;
					for(int k=j; k<kb_downstream&&dist<wbp; k++)
						dist += dists[k];
					while(dist<wbp&&kb_downstream<m-1) {
						dist += dists[kb_downstream];
						++kb_downstream;
					}
					
					iserr = true;
					min = Double.MAX_VALUE;
					outerloop:
						for(int u=j; u>=kb_upstream; u--) {
							for(int v=j+1; v<=kb_downstream; v++) {
								double rf = 0;
								for(int i=0; i<n; i++) {
									haps = haplotypes.get(i);
									double jump = 0;
									int hapn = 0;
									for(char[][] hap : haps.values()) {
										if(hap[0][0]=='*') continue;
										for(int p=0; p<ploidy; p++) {
											if(hap[p][u]!=hap[p][v]) 
												++jump;
										}
										++hapn;
									}
									rf += jump/hapn/ploidy;
								}
								rf /= n;
								
								if(rf<min) min = rf;
								if(rf<rf_thresh) {
									iserr = false;
									break outerloop;
								}
							}
						}
					if(iserr) {
						// this is an assembly error
						synchronized(lock) {
							rfWriter.write(scaff+" "+pos[j]+" "+pos[j+1]+" "+means[j]+"\n");
						}
						err.add(new int[]{pos[j], pos[j+1]});
					} else {
						myLogger.info(scaff+" "+pos[j]+" "+pos[j+1]+" "+means[j]+" "+min);
					}
				}
				
				if(err.size()>0) {
					int[][] arr = new int[err.size()][];
					for(int i=0; i<arr.length; i++)
						arr[i] = err.get(i);
					errs.put(scaff, arr);
				}
			} catch (Exception e) {
				Thread t = Thread.currentThread();
				t.getUncaughtExceptionHandler().uncaughtException(t, e);
				e.printStackTrace();
				executor.shutdown();
				System.exit(1);
			}
		}
	}
	
	public Map<String, int[][]> errs() {
		// TODO Auto-generated method stub
		return errs;
	}
	
	public Set<String> split(String in_vcf, String out_vcf) {
		// TODO Auto-generated method stub
		Set<String> scaffs = new HashSet<>();
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
					line = br.readLine();
				} else {
					String scaff = s[0];
					int[][] err = this.errs.get(scaff);
					int ec = err.length;
					int sub = 1;
					scaffs.add(scaff+"_"+sub);
					while( line!=null ) {
						s = line.split("\\s+");
						if( !scaff.equals(s[0]) ) break;
						if( sub<=ec&&Double.parseDouble(s[1])>err[sub-1][0] ) {
							++sub;
						}
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
		
		return scaffs;
	}
}
