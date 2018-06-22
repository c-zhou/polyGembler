package cz1.tenx.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;
import org.apache.commons.lang3.ArrayUtils;

import cz1.tenx.model.HiddenMarkovModel;
import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Executor;
import cz1.util.JohnsonTrotter;
import cz1.util.Utils;

public class Phaser extends Executor {

	private int max_iter = 1000;
	private String expr_id = null;
	private String out_prefix = null;
	private int block = 100000;
	private int overlap = 10000;
	private String dat_file;
	private String vcf_file;
	private String file_location;
	private String input_file;
	private String rangeChr;
	private int rangeLowerBound;
	private int rangeUpperBound;
	private int repeat = 10;
	private boolean simulated_annealing = false;
	private boolean clustering = true;
	private String ground_truth = null;
	private int minD = 1;
	
	private static enum Task {makeb, phase, pileup, zzz, eval};
	private Task task = Task.zzz;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		switch(this.task) {
		case zzz:
			myLogger.info(
					"\n\nChoose a task from listed:\n"
							+ " makeb                   Detect breaks and make blocks.\n"
							+ " phase                   Phase a block. \n"
							+ " pileup                  Pileup blocks.\n"
							+ " evaluate                Evaluate the phasing results with ground truth hapltypes.\n\n"
							+ "\n");
			break;
		case makeb:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+" -dat/--data-file             Molecule data file.\n"
							+" -vcf/--vcf-file              VCF file.\n"
							+" -b/--block-size              Block size (default 100000).\n"
							+" -olap/--overlap-size         Overlap size (default 10000).\n"
							+" -r/--range                   Data range (chr/chr:start-end).\n"
							+" -d/--depth                   Minimum depth to make a block (default 1).\n"
							+" -o/--prefix                  Output file location.\n\n"
					);
			break;
		case phase:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+" -dat/--data-file             Molecule data file.\n"
							+" -vcf/--vcf-file              VCF file.\n"
							+" -p/--ploidy                  Ploidy (default 2).\n"
							+" -r/--range                   Data range (chr/chr:start-end).\n"
							+" -a/--repeat                  Repeats (default 10).\n"
							+" -sa/--simulated-annealing    Run simulated annealing.\n"
							+" -no-clus/--no-clustering     Do not run pre-clustering.\n"
							+" -x/--max-iter                Maxmium rounds for EM optimization (default 1000).\n"
							+" -ex/--experiment-id          Common prefix of haplotype files for this experiment.\n"
							+" -@/--num-threads             Number of threads to use (default 1).\n"
							+" -rs/--random-seed            Random seed.\n"
							+" -o/--prefix                  Output file location.\n\n"
					);
			break;
		case pileup:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+" -f/--file-location           Location for the phasing output files.\n"
							+" -ex/--experiment-id          Common prefix of haplotype files for this experiment.\n"
							+" -p/--ploidy                  Ploidy (default 2).\n"
							+" -a/--repeat                  Repeats (default 10).\n"
							+" -@/--num-threads             Number of threads to use (default 1).\n"
							+" -o/--prefix                  Output file location.\n\n"
					);
			break;
		case eval:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+" -i/--input-file              Input file for phasing result.\n"
							+" -p/--ploidy                  Ploidy (default 2).\n"
							+" -hap/--haplotypes            Ground truth hapltypes for evaluation.\n"
					);
			break;
		default:
			throw new RuntimeException("!!!");
		}
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		// create the command line parser

		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}
		
		switch(args[0].toLowerCase()) {
		case "makeb":
			this.task = Task.makeb;
			break;
		case "phase":
			this.task = Task.phase;
			break;
		case "pileup":
			this.task = Task.pileup;
			break;
		case "evaluate":
			this.task = Task.eval;
			break;
		default:
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");	
		}
		
		String[] args2 = new String[args.length-1];
		System.arraycopy(args, 1, args2, 0, args2.length);

		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-dat", "--data-file", true);
			myArgsEngine.add("-vcf", "--vcf-file", true);
			myArgsEngine.add("-b", "--block-size", true);
			myArgsEngine.add("-olap", "--overlap-size", true);
			myArgsEngine.add("-p", "--ploidy", true);
			myArgsEngine.add("-r", "--range", true);
			myArgsEngine.add("-d", "--depth", true);
			myArgsEngine.add("-a", "--repeat", true);
			myArgsEngine.add("-sa", "--simulated-annealing", false);
			myArgsEngine.add("-no-clus", "--no-clustering", false);
			myArgsEngine.add("-x", "--max-iter", true);
			myArgsEngine.add("-ex", "--experiment-id", true);
			myArgsEngine.add("-@", "--num-threads", true);
			myArgsEngine.add("-rs", "--random-seed", true);
			myArgsEngine.add("-hap", "--haplotypes", true);
			myArgsEngine.add("-f", "--file-location", true);
			myArgsEngine.add("-i", "--input-file", true);
			myArgsEngine.add("-o", "--prefix", true);
			myArgsEngine.parse(args2);
		}

		switch(this.task) {
		case makeb:
			if(myArgsEngine.getBoolean("-dat")) {
				this.dat_file = myArgsEngine.getString("-dat");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify your input molecule data file.");
			}

			if(myArgsEngine.getBoolean("-vcf")) {
				this.vcf_file = myArgsEngine.getString("-vcf");
			}  else {
				printUsage();
				throw new IllegalArgumentException("Please specify your input VCF file.");
			}

			if(myArgsEngine.getBoolean("-b")) {
				this.block = Integer.parseInt(myArgsEngine.getString("-b"));
			}
			
			if(myArgsEngine.getBoolean("-d")) {
				this.minD = Integer.parseInt(myArgsEngine.getString("-d"));
			}

			if(myArgsEngine.getBoolean("-olap")) {
				this.overlap = Integer.parseInt(myArgsEngine.getString("-olap"));
			}

			if(myArgsEngine.getBoolean("-r")) {
				this.setDataRange(myArgsEngine.getString("-r"));
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify your data range.");
			}

			if(myArgsEngine.getBoolean("-o")) {
				this.out_prefix = myArgsEngine.getString("-o");
			}  else {
				printUsage();
				throw new IllegalArgumentException("Please specify your output file prefix.");
			}
			break;
		case phase:
			if(myArgsEngine.getBoolean("-dat")) {
				this.dat_file = myArgsEngine.getString("-dat");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify your input molecule data file.");
			}

			if(myArgsEngine.getBoolean("-vcf")) {
				this.vcf_file = myArgsEngine.getString("-vcf");
			}  else {
				printUsage();
				throw new IllegalArgumentException("Please specify your input VCF file.");
			}

			if(myArgsEngine.getBoolean("-p")) {
				Constants.ploidy(Integer.parseInt(myArgsEngine.getString("-p")));
			}

			if(myArgsEngine.getBoolean("-r")) {
				this.setDataRange(myArgsEngine.getString("-r"));
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify your data range.");
			}

			if(myArgsEngine.getBoolean("-a")) {
				this.repeat = Integer.parseInt(myArgsEngine.getString("-a"));
			}

			if(myArgsEngine.getBoolean("-sa")) {
				this.simulated_annealing = true;
			}
			
			if(myArgsEngine.getBoolean("-no-clus")) {
				this.clustering = false;
			}

			if(myArgsEngine.getBoolean("-x")) {
				this.max_iter = Integer.parseInt(myArgsEngine.getString("-x"));
			}

			if(myArgsEngine.getBoolean("-ex")) {
				this.expr_id = myArgsEngine.getString("-ex");
			}  else {
				printUsage();
				throw new IllegalArgumentException("Please specify your experiment identifier.");	
			}

			if(myArgsEngine.getBoolean("-@")) {
				this.THREADS = Integer.parseInt(myArgsEngine.getString("-@"));
			}
			
			if(myArgsEngine.getBoolean("-rs")) {
				Constants.seeding(Long.parseLong(myArgsEngine.getString("-rs")));
			}
			
			if(myArgsEngine.getBoolean("-o")) {
				this.out_prefix = myArgsEngine.getString("-o");
			}  else {
				printUsage();
				throw new IllegalArgumentException("Please specify your output file prefix.");
			}

			break;
		case pileup:
			if(myArgsEngine.getBoolean("-ex")) {
				this.expr_id = myArgsEngine.getString("-ex");
			}  else {
				printUsage();
				throw new IllegalArgumentException("Please specify your experiment identifier.");	
			}
			
			if(myArgsEngine.getBoolean("-f")) {
				this.file_location = myArgsEngine.getString("-f");
			}  else {
				printUsage();
				throw new IllegalArgumentException("Please specify your file location.");	
			}

			if(myArgsEngine.getBoolean("-p")) {
				Constants.ploidy(Integer.parseInt(myArgsEngine.getString("-p")));
			}
			
			if(myArgsEngine.getBoolean("-a")) {
				this.repeat = Integer.parseInt(myArgsEngine.getString("-a"));
			}

			if(myArgsEngine.getBoolean("-@")) {
				this.THREADS = Integer.parseInt(myArgsEngine.getString("-@"));
			}
			
			if(myArgsEngine.getBoolean("-o")) {
				this.out_prefix = myArgsEngine.getString("-o");
			}  else {
				printUsage();
				throw new IllegalArgumentException("Please specify your output file prefix.");
			}

			break;
		case eval:
			if(myArgsEngine.getBoolean("-i")) {
				this.input_file = myArgsEngine.getString("-i");
			}  else {
				printUsage();
				throw new IllegalArgumentException("Please specify your input file.");
			}
			
			if(myArgsEngine.getBoolean("-p")) {
				Constants.ploidy(Integer.parseInt(myArgsEngine.getString("-p")));
			}
			
			if(myArgsEngine.getBoolean("-hap")) {
				this.ground_truth = myArgsEngine.getString("-hap");
			}   else {
				printUsage();
				throw new IllegalArgumentException("Please specify your ground truth haplotype file.");
			}
			break;
		default:
			throw new RuntimeException("!!!");
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub

		switch(this.task) {
		case makeb:
			this.runMakeb();
			break;
		case phase:
			this.runPhase();
			break;
		case pileup:
			this.runPileup();
			break;
		case eval:
			this.runEval();
			break;
		default:
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");	
		}
	}
	
	private void runMakeb() {
		// TODO Auto-generated method stub
		try {
			// read VCF file
			BufferedReader br_vcf = Utils.getBufferedReader(vcf_file);
			String line;
			String[] s;
			int position;
			int var_start = Integer.MAX_VALUE, var_end = Integer.MIN_VALUE;
			int ind = 0;
			final BidiMap<Integer, Integer> idxmap = new DualHashBidiMap<Integer, Integer>();
			
			while( (line=br_vcf.readLine())!=null ){
				if(line.startsWith("#")) continue;
				++ind;
				s = line.split("\\s+");
				if(!s[0].equals(rangeChr)) continue;
				position = Integer.parseInt(s[1]);
				if(position<rangeLowerBound) continue;
				if(position>rangeUpperBound) break;
				if(var_start>ind) var_start = ind;
				if(var_end  <ind)   var_end   = ind;
				idxmap.put(ind, position);
			}
			br_vcf.close();

			Map<Integer, Integer> depth = new HashMap<Integer, Integer>();
			for(int i=var_start; i<var_end; i++) depth.put(i, 0);
			
			// read DAT file
			BufferedReader br_dat = Utils.getBufferedReader(dat_file);
			int n, vstart, vend;
			while( (line=br_dat.readLine())!=null ) {
				s = line.split("\\s+");
				if(!s[1].split(":")[0].equals(rangeChr)) continue;
				n = s.length;
				vstart = Integer.parseInt(s[5]);
				vend   = Integer.parseInt(s[n-3])+s[n-2].length()-1;
				if(vend  <var_start) continue;
				if(vstart>var_end  ) break;
				if(vstart<var_start) vstart = var_start;
				if(vend  >var_end  ) vend   = var_end  ;
				for(int i=vstart; i<vend; i++) depth.put(i, depth.get(i)+1);
			}
			br_dat.close();
			
			/***
			int cursor = rangeLowerBound;
			final Map<Integer, Integer> depth = new HashMap<Integer, Integer>();
			for(final int[] mol : mol_range) {
				vstart = mol[0];
				vend   = mol[1];
				for(int i=vstart; i<=vend; i++) {
					if(vstart<cursor) continue;
					if(depth.containsKey(i))
						depth.put(i, depth.get(i)+1);
					else depth.put(i, 1);
				}
				for(int i=cursor; ;i++) {
					if(depth.containsKey(i)&&depth.get(i)>=minD) {
						depth.remove(i);
						++cursor;
					} else break;
				}
				if(cursor<vstart) cursor = vstart;
			}
			**/
			// get positions to break
			final List<Integer> lowCovR = new ArrayList<Integer>();
			for(final int i : depth.keySet()) 
				if(depth.get(i)<minD) lowCovR.add(i);
			Collections.sort(lowCovR);
			
			n = lowCovR.size()+1;
			int[][] ranges = new int[n][2];
			ranges[0  ][0] = rangeLowerBound;
			ranges[n-1][1] = rangeUpperBound;
			for(int i=1; i<n; i++) {
				ranges[i-1][1] = idxmap.get(lowCovR.get(i-1)  );
				ranges[i  ][0] = idxmap.get(lowCovR.get(i-1)+1);
			}

			BufferedWriter bw = Utils.getBufferedWriter(this.out_prefix+".b");
			int lb, ub;
			for(int i=0; i<n; i++) {
				lb = ranges[i][0];
				ub = ranges[i][1];
				int nB = (int) Math.ceil((double)(ub-lb+1-overlap)/(block-overlap));
				int from, to;
				if(nB<1) nB = 1;
				for(int j=0; j<nB; j++) {
					from = lb+j*(block-overlap);
					to   = Math.min(ub, from+block);
					bw.write(rangeChr+":"+from+"-"+to+"\n");
				}
			}
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return;
		
	}
		
	private void runPhase() {
		// TODO Auto-generated method stub
		myLogger.info("Random seed - "+Constants.seed);
		this.initial_thread_pool();

		for(int j=0; j<this.repeat; j++) {
			this.executor.submit(new Runnable() {
				private int r;
				
				@Override
				public void run() {
					// TODO Auto-generated method stub
					try {
						HiddenMarkovModel hmm = new HiddenMarkovModel(vcf_file, dat_file, rangeChr,
								rangeLowerBound, rangeUpperBound, clustering, simulated_annealing);
						if(hmm.isNullModel()) return;
						hmm.print();
						hmm.train();
						double loglik = hmm.loglik(), loglik1 = 0;
						for(int i=0; i<max_iter; i++) {
							hmm.train();
							loglik1 = hmm.loglik();
							if(loglik1>=loglik && 
									(loglik-loglik1)/loglik<1e-6) 
								break;
							loglik = loglik1;
						}
						hmm.write(out_prefix+"/"+expr_id+"_"+rangeChr+"~"+rangeLowerBound+"-"+rangeUpperBound+"_"+r+".tmp");
						hmm.print();
					} catch (Exception e) {
						myLogger.error("##########ERROR MESSAGE##########");
						myLogger.error(vcf_file);
						myLogger.error(dat_file);
						myLogger.error(rangeChr+":"+rangeLowerBound+"-"+rangeUpperBound);
						myLogger.error("########END ERROR MESSAGE########");
						Thread t = Thread.currentThread();
						t.getUncaughtExceptionHandler().uncaughtException(t, e);
						e.printStackTrace();
						executor.shutdown();
						System.exit(1);
					}
				}

				public Runnable init(final int i) {
					// TODO Auto-generated method stub
					this.r = i;
					return this;
				}
			}.init(j));
		}
		this.waitFor();
	}
	
	private static class Block implements Comparable<Block> {
		private final String chr;
		private final int start;
		private final int end;
		
		public Block(final String chr,
				final int start,
				final int end) {
			this.chr   = chr;
			this.start = start;
			this.end   = end;
		}
		
		@Override    
	    public boolean equals(Object obj) {        
			if (this == obj) return true;        
	        if (obj == null || getClass() != obj.getClass()) 
	        	return false;
	        Block b = (Block) obj;
	        return this.chr.equals(b.chr) &&
	        		this.start == b.start &&
	        		this.end == b.end;
	    }    
	    @Override    
	    public int hashCode() {
	    	return this.toString().hashCode();
	    }    
		
		@Override
		public String toString() {
			return chr+"~"+start+"-"+end;
		}

		@Override
		public int compareTo(Block b) {
			// TODO Auto-generated method stub
			return this.chr.equals(b.chr)?(this.start==b.start?
					Integer.compare(this.end, b.end):
				Integer.compare(this.start, b.start)):
				this.chr.compareToIgnoreCase(b.chr);
		}

		public static boolean overlap(Block b1, Block b2) {
			// TODO Auto-generated method stub
			if(!b1.chr.equals(b2.chr)) return false;
			return b1.end>b2.start;
		}
	}
	
	private void runPileup() {
		// TODO Auto-generated method stub
		if(new File(out_prefix+"/"+expr_id+".phase").exists()) 
			throw new RuntimeException(out_prefix+"/"+expr_id+".phase file exists");
		
		File[] listOfFiles = new File(this.file_location).listFiles();
		final Set<Block> blockSet = new HashSet<Block>();
		String fname;
		String[] s;
		for(File f : listOfFiles) {
			fname = f.getName();
			if(!fname.startsWith(expr_id)) continue;
			fname = fname.replaceAll("^"+expr_id+"_", "");
			fname = fname.split("_")[0];
			s = fname.split("~|-");
			blockSet.add(new Block(s[0], 
					Integer.parseInt(s[1]),
					Integer.parseInt(s[2])));
		}
		final List<Block> blocks = new ArrayList<Block>();
		blocks.addAll(blockSet);
		Collections.sort(blocks);
		
		if(blocks.isEmpty()) throw new RuntimeException("no phasing file found!!!");
		final List<int[][]> phasedData = new ArrayList<int[][]>();
		Iterator<Block> iter1 = blocks.iterator();
		Block tmp_block = iter1.next(), block;
		String block_chr;
		int block_start, block_end;
		
		try {
			BufferedWriter bw = Utils.getBufferedWriter(out_prefix+"/"+expr_id+".phase");

			while(tmp_block!=null) {
				phasedData.clear();
				block = tmp_block;
				phasedData.add(ballot(loadPhasedData(block.toString())));	
				block_chr = block.chr;
				block_start = block.start;
				block_end = block.end;
				while( (tmp_block=(iter1.hasNext()?iter1.next():null))!=null && 
						Block.overlap(block, tmp_block) ) {
					block = tmp_block;
					phasedData.add(ballot(loadPhasedData(block.toString())));
					if(block.end>block_end) block_end = block.end;
				}

				final List<int[][]> pileupData = new ArrayList<int[][]>();
				pileupData.add(phasedData.get(0));
				// TODO parallelise this, say run it in a hierarchical process
				for(int i=1; i<phasedData.size(); i++) {
					if(pileupData.get(i-1)==null||phasedData.get(i)==null) continue;

					final int[][] lastPhase = pileupData.get(i-1);
					final Set<Integer> lastPos = new HashSet<Integer>();
					for(int j : lastPhase[Constants._ploidy_H])
						lastPos.add(j);
					final int[][] thisPhase = phasedData.get(i);
					final Set<Integer> thisPos = new HashSet<Integer>();
					for(int j : thisPhase[Constants._ploidy_H])
						thisPos.add(j);
					thisPos.retainAll(lastPos);
					final int z = thisPos.size();
					if(z<10) {
						//TODO fix this problem
						bw.close();
						throw new RuntimeException("Overlap too small!!!"); 
					}
					final int[][] lastOlapPhase = new int[Constants._ploidy_H+1][z];
					int w, v;
					w = 0;
					for(int j=0; j<lastPhase[0].length; j++) {
						if(thisPos.contains(lastPhase[Constants._ploidy_H][j])) {
							for(int k=0; k<Constants._ploidy_H+1; k++)
								lastOlapPhase[k][w] = lastPhase[k][j];
							++w;
						}
					}
					final int[][] thisOlapPhase = new int[Constants._ploidy_H+1][z];
					final int[][] thisNonOlapPhase = new int[Constants._ploidy_H+1]
							[thisPhase[0].length-z];
					w = 0;
					v = 0;
					for(int j=0; j<thisPhase[0].length; j++) {
						if(thisPos.contains(thisPhase[Constants._ploidy_H][j])) {
							for(int k=0; k<Constants._ploidy_H+1; k++)
								thisOlapPhase[k][w] = thisPhase[k][j];
							++w;
						} else {
							for(int k=0; k<Constants._ploidy_H+1; k++)
								thisNonOlapPhase[k][v] = thisPhase[k][j];
							++v;
						}
					}
					final int[] ballot = ballot(lastOlapPhase, thisOlapPhase);
					final int[][] thisSwapPhase = new int[thisNonOlapPhase.length][];
					for(int j=0; j<Constants._ploidy_H; j++) {
						if(ballot[Constants._ploidy_H+1]==0) {
							thisSwapPhase[j] = thisNonOlapPhase[ballot[j]];
						} else {
							thisSwapPhase[j] = flip(thisNonOlapPhase[ballot[j]]);
						}
					}
					thisSwapPhase[Constants._ploidy_H] = thisNonOlapPhase[Constants._ploidy_H];

					pileupData.add(thisSwapPhase);
				}

				bw.write("BLOCK: "+block_chr+"\t"+block_start+"\t"+block_end+"\n");
				for(int i=0; i<pileupData.size(); i++) {
					int[][] phase = pileupData.get(i);
					for(int j=0; j<phase[0].length; j++) {
						bw.write(block_chr+"\t"+
								String.format("%1$12s", phase[Constants._ploidy_H][j])+"\t"+
								variants.get(block_chr).get(phase[Constants._ploidy_H][j]));
						for(int k=0; k<Constants._ploidy_H; k++) 
							bw.write("\t"+(phase[k][j]==-1?"-":phase[k][j]));
						bw.write("\n");	
					}
				}
				bw.write("********\n");
			}
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return;	
	}

	private void runEval() {
		// TODO Auto-generated method stub
		try {
			
			BufferedReader br = Utils.getBufferedReader(this.ground_truth);
			Map<String, Map<Integer, int[]>> true_haps = new HashMap<String, Map<Integer, int[]>>();
			String[] s;
			int pos;
			String line;
			while( (line=br.readLine())!=null ) {
				if(line.startsWith("#")) continue;
				s = line.trim().split("\\s+");
				final int[] var = new int[Constants._ploidy_H];
				for(int i=0; i<Constants._ploidy_H; i++) 
					var[i] = Integer.parseInt(s[5+i]);
				if(!true_haps.containsKey(s[0]))
					true_haps.put(s[0], new HashMap<Integer, int[]>());
				true_haps.get(s[0]).put(Integer.parseInt(s[1]), var);
			}
			br.close();
			
			long sites = 0, mismatches = 0;
			
			br = Utils.getBufferedReader(this.input_file);
			line = br.readLine();
			while(line!=null) {
				if(!line.startsWith("BLOCK")) throw new RuntimeException("truncated phasing file!!!");

				final List<int[]> true_hap = new ArrayList<int[]>();
				final List<int[]> phas_hap = new ArrayList<int[]>();

				while( (line=br.readLine())!=null && !line.equals("********")) {
					s = line.trim().split("\\s+");
					pos = Integer.parseInt(s[1]);
					if(true_haps.containsKey(s[0])&&
							true_haps.get(s[0]).containsKey(pos)) {
						true_hap.add(true_haps.get(s[0]).get(pos));
						final int[] phas = new int[Constants._ploidy_H];
						for(int i=0; i<Constants._ploidy_H; i++) 
							phas[i] = s[4+i].equals("-")?-1:Integer.parseInt(s[4+i]);
						phas_hap.add(phas);
					}					
				}
				final int n = true_hap.size();
				final int[][] groundTruthSel = new int[Constants._ploidy_H][n];
				final int[][] phasedHaploSel = new int[Constants._ploidy_H][n];

				int[] true_sel, phas_sel;
				for(int i=0; i<n; i++) {
					true_sel = true_hap.get(i);
					phas_sel = phas_hap.get(i);
					for(int j=0; j<Constants._ploidy_H; j++) {
						groundTruthSel[j][i] = true_sel[j];
						phasedHaploSel[j][i] = phas_sel[j];
					}
				}
				final int[] ballot = ballot(groundTruthSel, phasedHaploSel);
				sites += n;
				mismatches += ballot[Constants._ploidy_H];
				
				line = br.readLine();
			}
			br.close();
			myLogger.info("####error rate: "+(double)mismatches/Constants._ploidy_H/sites+
					"("+mismatches+"/"+Constants._ploidy_H*sites+")");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private int[] flip(int[] is) {
		// TODO Auto-generated method stub
		for(int i=0; i<is.length; i++)
			is[i] = 1-is[i];
		return is;
	}

	private int[][] ballot(final int[][][] phasedData) {
		// TODO Auto-generated method stub
		if(phasedData==null) return null;
		int[][] pivot = phasedData[0];
		for(int i=1; i<phasedData.length; i++) {
			int[][] exterior = phasedData[i];
			ballot(pivot, exterior);
		}
		
		final int[][] ballot = new int[pivot.length][pivot[0].length];
		final int[] count = new int[2];
		for(int i=0; i<ballot.length-1; i++) {
			for(int j=0; j<ballot[i].length; j++) {
				Arrays.fill(count, 0);
				for(int k=0; k<phasedData.length; k++) {
					if(phasedData[k][i][j]!=-1)
						++count[phasedData[k][i][j]];
				}
				ballot[i][j]=(count[0]==0&&count[1]==0)?-1:(count[0]>=count[1]?0:1);
			}
		}
		ballot[Constants._ploidy_H] = pivot[Constants._ploidy_H];
		return ballot;
	}

	private int[] ballot(final int[][] pivot, final int[][] exterior) {
		// TODO Auto-generated method stub
		// return [1...P,#diff,flip]
		int[][] jtPerm = JohnsonTrotter.perm(Constants._ploidy_H);
		int[][] dist0 = new int[jtPerm.length][Constants._ploidy_H]; //original
		int[][] dist1 = new int[jtPerm.length][Constants._ploidy_H]; //switch
		for(int i=0; i<Constants._ploidy_H; i++) {
			dist0[0][i] = hammingDist(pivot[i], exterior[i], false);
			dist1[0][i] = hammingDist(pivot[i], exterior[i], true);
		}
		int[] tmp;
		int tmp_int;
		for(int i=1; i<jtPerm.length; i++) {
			tmp = exterior[jtPerm[i][0]];
			exterior[jtPerm[i][0]] = exterior[jtPerm[i][1]];
			exterior[jtPerm[i][1]] = tmp;
			System.arraycopy(dist0[i-1], 0, dist0[i], 0, Constants._ploidy_H);
			dist0[i][jtPerm[i][0]] = hammingDist(pivot[jtPerm[i][0]],exterior[jtPerm[i][0]],false);
			dist0[i][jtPerm[i][1]] = hammingDist(pivot[jtPerm[i][1]],exterior[jtPerm[i][1]],false);
			System.arraycopy(dist1[i-1], 0, dist1[i], 0, Constants._ploidy_H);
			dist1[i][jtPerm[i][0]] = hammingDist(pivot[jtPerm[i][0]],exterior[jtPerm[i][0]],true);
			dist1[i][jtPerm[i][1]] = hammingDist(pivot[jtPerm[i][1]],exterior[jtPerm[i][1]],true);
		}

		// swap (0, 1) to reset exterior
		tmp = exterior[0];
		exterior[0] = exterior[1];
		exterior[1] = tmp;

		// find minimum hamming distance
		int minD = Integer.MAX_VALUE, k = -1, tmpD;
		for(int i=0; i<jtPerm.length; i++) {
			if( (tmpD=sumInts(dist0[i]))<minD ) {
				k = i;
				minD = tmpD;
			}
		}
		int flip = 0;
		for(int i=0; i<jtPerm.length; i++) {
			if( (tmpD=sumInts(dist1[i]))<minD ) {
				k = i;
				minD = tmpD;
				flip = 1;
			}
		}

		int[] ballot = new int[Constants._ploidy_H+2];
		for(int i=0; i<Constants._ploidy_H; i++)
			ballot[i] = i;
		
		// swapping exterior phased data
		for(int i=1; i<=k; i++) {
			tmp = exterior[jtPerm[i][0]];
			exterior[jtPerm[i][0]] = exterior[jtPerm[i][1]];
			exterior[jtPerm[i][1]] = tmp;
			tmp_int = ballot[jtPerm[i][0]];
			ballot[jtPerm[i][0]] = ballot[jtPerm[i][1]];
			ballot[jtPerm[i][1]] = tmp_int;
		}

		// re-checking
		int d = 0;
		for(int i=0; i<Constants._ploidy_H; i++) 
			d += hammingDist(pivot[i], exterior[i], flip==0?false:true);

		if(minD!=d) throw new RuntimeException("!!!");
		ballot[Constants._ploidy_H  ] = minD;
		ballot[Constants._ploidy_H+1] = flip;
		
		myLogger.info("ballot hamming dist, #"+minD);
		
		return ballot;
	}
	
	private int sumInts(int[] is) {
		// TODO Auto-generated method stub
		int sum = 0;
		for(final int i : is) sum += i;
		return sum;
	}

	private int hammingDist(int[] is, int[] is2, boolean flip) {
		// TODO Auto-generated method stub
		int d = 0;
		if(flip) {
			for(int i=0; i<is.length; i++) {
				if(is[i]==-1||is2[i]==-1)
					continue;
				if(is[i]==is2[i]) ++d;
			}
		} else {
			for(int i=0; i<is.length; i++) {
				if(is[i]==-1||is2[i]==-1)
					continue;
				if(is[i]!=is2[i]) ++d;
			}
		}
		return d;
	}

	private final Map<String, Map<Integer, String>> variants = new HashMap<String, Map<Integer, String>>();
	
	private int[][][] loadPhasedData(final String i) {
		// TODO Auto-generated method stub
		// load variants data
		try {
			BufferedReader br = Utils.getBufferedReader(out_prefix+"/"+expr_id+"_"+i+"_0.tmp");
			String line;
			String[] s;
			while( (line=br.readLine())!=null ) {
				if(line.startsWith("#")) continue;
				s = line.trim().split("\\s+");
				if(!variants.containsKey(s[0])) variants.put(s[0], new HashMap<Integer, String>());
				variants.get(s[0]).put(Integer.parseInt(s[1]), s[2]+"\t"+s[3]);
			}
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		int[][][] data = null;
		for(int j=0; j<repeat; j++) {
			try {
				if(new File(out_prefix+"/"+expr_id+"_"+i+"_"+j+".tmp").exists()) {
					BufferedReader br = Utils.getBufferedReader(out_prefix+"/"+expr_id+"_"+i+"_"+j+".tmp");
					String line;
					String[] s;
					int n = 0;
					
					if(data==null) {
						line = br.readLine();
						line = line.replaceAll("^##", "");
						s = line.split("\\s+");
						data = new int[repeat][Constants._ploidy_H+1][Integer.parseInt(s[0])];
					}
					
					while( (line=br.readLine())!=null ) {
						if(line.startsWith("#")) continue;
						s = line.trim().split("\\s+");
						data[j][Constants._ploidy_H][n] = Integer.parseInt(s[1]);
						for(int k=0; k<Constants._ploidy_H; k++) {
							if(s[4+k].equals("-"))
								data[j][k][n] = -1;
							else
								data[j][k][n] = Integer.parseInt(s[4+k]);
						}
						++n;
					}
					br.close();
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return data;
	}

	private void setDataRange(String dat_range) {
		// TODO Auto-generated method stub
		String[] s = dat_range.split(":");
		rangeChr = s[0];
		if(s.length==2) {
			s = s[1].split("-");
			rangeLowerBound = Integer.parseInt(s[0]);
			rangeUpperBound = Integer.parseInt(s[1]);
			if(rangeLowerBound>rangeUpperBound) 
				throw new RuntimeException("invalid data range!!!");
		} else {
			rangeLowerBound = Integer.MAX_VALUE;
			rangeUpperBound = Integer.MIN_VALUE;
			try {
				BufferedReader br = Utils.getBufferedReader(this.vcf_file);
				String line;
				int i;
				while( (line=br.readLine())!=null ){
					if(line.startsWith("#")) continue;
					s = line.split("\\s+");
					if(!s[0].equals(rangeChr)) continue;
					if( (i=Integer.parseInt(s[1]))<rangeLowerBound )
						rangeLowerBound = i;
					if( (i=Integer.parseInt(s[1]))>rangeUpperBound )
						rangeUpperBound = i;
				}
				br.close();		
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
}
