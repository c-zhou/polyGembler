package cz1.tenx.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

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
	private int block = 1000000;
	private int overlap = 100000;
	private String dat_file;
	private String vcf_file;
	private String rangeChr;
	private int rangeLowerBound;
	private int rangeUpperBound;
	private int repeat = 30;
	private boolean simulated_annealing = false;
	private boolean clustering = true;
	private String ground_truth = null;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+" -dat/--data-file             Molecule data file.\n"
						+" -vcf/--vcf-file              VCF file.\n"
						+" -b/--block-size              Block size (default 1000000).\n"
						+" -olap/--overlap-size         Overlap size (default 100000).\n"
						+" -p/--ploidy                  Ploidy (default 2).\n"
						+" -r/--range                   Data range (chr/chr:start-end).\n"
						+" -f/--repeat                  Repeats (default 30).\n"
						+" -sa/--simulated-annealing    Run simulated annealing.\n"
						+" -no-clus/--no-clustering     Do not run pre-clustering.\n"
						+" -x/--max-iter                Maxmium rounds for EM optimization (default 1000).\n"
						+" -ex/--experiment-id          Common prefix of haplotype files for this experiment.\n"
						+" -@/--num-threads             Number of threads to use (default 1).\n"
						+" -rs/--random-seed            Random seed.\n"
						+" -hap/--haplotypes            Ground truth hapltypes for evaluation.\n"
						+" -o/--prefix                  Output file location.\n\n"
				);
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		// create the command line parser

		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}

		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-dat", "--data-file", true);
			myArgsEngine.add("-vcf", "--vcf-file", true);
			myArgsEngine.add("-b", "--block-size", true);
			myArgsEngine.add("-olap", "--overlap-size", true);
			myArgsEngine.add("-p", "--ploidy", true);
			myArgsEngine.add("-r", "--range", true);
			myArgsEngine.add("-f", "--repeat", true);
			myArgsEngine.add("-sa", "--simulated-annealing", false);
			myArgsEngine.add("-no-clus", "--no-clustering", false);
			myArgsEngine.add("-x", "--max-iter", true);
			myArgsEngine.add("-ex", "--experiment-id", true);
			myArgsEngine.add("-@", "--num-threads", true);
			myArgsEngine.add("-rs", "--random-seed", true);
			myArgsEngine.add("-hap", "--haplotypes", true);
			myArgsEngine.add("-o", "--prefix", true);
			myArgsEngine.parse(args);
		}

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

		if(myArgsEngine.getBoolean("-olap")) {
			this.overlap = Integer.parseInt(myArgsEngine.getString("-olap"));
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

		if(myArgsEngine.getBoolean("-f")) {
			this.repeat = Integer.parseInt(myArgsEngine.getString("-f"));
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
		
		if(myArgsEngine.getBoolean("-hap")) {
			this.ground_truth = myArgsEngine.getString("-hap");
		}

		if(myArgsEngine.getBoolean("-o")) {
			this.out_prefix = myArgsEngine.getString("-o");
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify your output file prefix.");
		}
	}

	private int[][] dataB;
	private double[][] logLik;
	private int[] noLoci;
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		myLogger.info("Random seed - "+Constants.seed);
		this.makeDataBlock();
		this.initial_thread_pool();
		
		for(int i=0; i<this.dataB.length; i++) {
			for(int j=0; j<this.repeat; j++) {
				this.executor.submit(new Runnable() {
					private int b;
					private int r;

					@Override
					public void run() {
						// TODO Auto-generated method stub
						try {
							HiddenMarkovModel hmm = new HiddenMarkovModel(vcf_file, dat_file, rangeChr, dataB[b][0], dataB[b][1], clustering, simulated_annealing);
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
							logLik[b][r] = loglik1;
							hmm.write(out_prefix+"/"+expr_id+"_"+b+"_"+r+".tmp");
							hmm.print();
							noLoci[b] = hmm.noLoci();
						} catch (Exception e) {
							myLogger.error("##########ERROR MESSAGE##########");
							myLogger.error(vcf_file);
							myLogger.error(dat_file);
							myLogger.error(rangeChr+":"+dataB[b][0]+"-"+dataB[b][1]);
							myLogger.error("########END ERROR MESSAGE########");
							Thread t = Thread.currentThread();
							t.getUncaughtExceptionHandler().uncaughtException(t, e);
							e.printStackTrace();
							executor.shutdown();
							System.exit(1);
						}
					}

					public Runnable init(final int b, final int r) {
						// TODO Auto-generated method stub
						this.b = b;
						this.r = r;
						return this;
					}

				}.init(i, j));
			}
		}
		this.waitFor();
		this.pileup();
		if(this.ground_truth!=null) this.evaluate();
	}
	
	private void evaluate() {
		// TODO Auto-generated method stub
		try {
			BufferedReader br = Utils.getBufferedReader(out_prefix+"/"+expr_id+".phase");
			final Set<Integer> positions = new HashSet<Integer>();
			String line;
			while( (line=br.readLine())!=null ) {
				if(line.startsWith("#")) continue;
				positions.add(Integer.parseInt(line.trim().split("\\s+")[0]));
			}
			br.close();
			
			final Set<Integer> posSel = new HashSet<Integer>();
			br = Utils.getBufferedReader(this.ground_truth);
			final List<int[]> groundTruthSelList = new ArrayList<int[]>();
			String[] s;
			int pos;
			while( (line=br.readLine())!=null ) {
				if(line.startsWith("#")) continue;
				s = line.trim().split("\\s+");
				if(positions.contains(pos=Integer.parseInt(s[1]))) {
					posSel.add(pos);
					final int[] var = new int[Constants._ploidy_H];
					for(int i=0; i<Constants._ploidy_H; i++) 
						var[i] = Integer.parseInt(s[5+i]);
					groundTruthSelList.add(var);	
				}
			}
			br.close();
			final int[][] groundTruthSel = new int[Constants._ploidy_H][groundTruthSelList.size()];
			int[] phase;
			for(int i=0; i<groundTruthSelList.size(); i++) {
				phase = groundTruthSelList.get(i);
				for(int j=0; j<Constants._ploidy_H; j++)
					groundTruthSel[j][i] = phase[j]; 
			}
			final int[][] phasedHapSel = new int[Constants._ploidy_H][groundTruthSelList.size()];
			br = Utils.getBufferedReader(out_prefix+"/"+expr_id+".phase");
			int n = 0;
			while( (line=br.readLine())!=null ) {
				if(line.startsWith("#")) continue;
				s = line.trim().split("\\s+");
				if(posSel.contains(pos=Integer.parseInt(s[0]))) {
					for(int i=0; i<Constants._ploidy_H; i++) 
						phasedHapSel[i][n] = Integer.parseInt(s[1+i]);
					++n;
				}
			}
			br.close();
			
			final int[] ballot = ballot(groundTruthSel, phasedHapSel);
			myLogger.info("####error rate: "+(double)ballot[Constants._ploidy_H]/Constants._ploidy_H/n+
					"("+ballot[Constants._ploidy_H]+"/"+Constants._ploidy_H*n+")");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void pileup() {
		// TODO Auto-generated method stub
		final List<int[][]> phasedData = new ArrayList<int[][]>();
		for(int i=0; i<dataB.length; i++) 
			phasedData.add(ballot(loadPhasedData(i)));
		final List<int[][]> pileupData = new ArrayList<int[][]>();
		pileupData.add(phasedData.get(0));
		for(int i=1; i<dataB.length; i++) {
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
			if(z<10) throw new RuntimeException("Overlap too small!!!"); 
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
		
		try {
			BufferedWriter bw = Utils.getBufferedWriter(out_prefix+"/"+expr_id+".phase");
			for(int i=0; i<pileupData.size(); i++) {
				int[][] phase = pileupData.get(i);
				for(int j=0; j<phase[0].length; j++) {
					bw.write(String.format("%1$12s", phase[Constants._ploidy_H][j]));
					for(int k=0; k<Constants._ploidy_H; k++) 
						bw.write("\t"+phase[k][j]);
					bw.write("\n");	
				}
			}
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return;
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
					++count[phasedData[k][i][j]];
				}
				ballot[i][j]=(count[0]>=count[1]?0:1);
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
			for(int i=0; i<is.length; i++)
				if(is[i]==is2[i]) ++d;
		} else {
			for(int i=0; i<is.length; i++)
				if(is[i]!=is2[i]) ++d;
		}
		return d;
	}

	private int[][][] loadPhasedData(final int i) {
		// TODO Auto-generated method stub
		int[][][] data = null;
		for(int j=0; j<repeat; j++) {
			try {
				if(new File(out_prefix+"/"+expr_id+"_"+i+"_"+j+".tmp").exists()) {
					if(data==null) data = new int[repeat][Constants._ploidy_H+1][noLoci[i]];
					BufferedReader br = Utils.getBufferedReader(out_prefix+"/"+expr_id+"_"+i+"_"+j+".tmp");
					String line;
					String[] s;
					int n = 0;
					while( (line=br.readLine())!=null ) {
						if(line.startsWith("#")) continue;
						s = line.trim().split("\\s+");
						data[j][Constants._ploidy_H][n] = Integer.parseInt(s[0]);
						for(int k=0; k<Constants._ploidy_H; k++)
							data[j][k][n] = Integer.parseInt(s[3+k]);
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

	private void makeDataBlock() {
		// TODO Auto-generated method stub
		int nB = (int) Math.ceil((double)(rangeUpperBound-rangeLowerBound+1-overlap)/(block-overlap));
		if(nB<1) nB = 1;
		dataB = new int[nB][2];
		for(int i=0; i<nB; i++) {
			int from = rangeLowerBound+i*(block-overlap);
			dataB[i] = new int[]{from, Math.min(rangeUpperBound, from+block)};
		}
		logLik = new double[nB][repeat];
		noLoci = new int[nB];
		for(int i=0; i<nB; i++) 
			myLogger.error("#"+String.format("%04d", i)+" | "+rangeChr+":"+dataB[i][0]+"-"+dataB[i][1]);
		return;
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
