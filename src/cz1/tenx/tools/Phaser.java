package cz1.tenx.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import cz1.tenx.model.HiddenMarkovModel;
import cz1.util.Algebra;
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
	private boolean simulated_annealing = true;

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
						+" -nosa/--no-sa                Do not run simulated annealing.\n"
						+" -x/--max-iter                Maxmium rounds for EM optimization (default 1000).\n"
						+" -ex/--experiment-id          Common prefix of haplotype files for this experiment.\n"
						+" -@/--num-threads             Number of threads to use (default 1).\n"
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
			myArgsEngine.add("-nosa", "--no-sa", false);
			myArgsEngine.add("-x", "--max-iter", true);
			myArgsEngine.add("-ex", "--experiment-id", true);
			myArgsEngine.add("-@", "--num-threads", true);
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

		if(myArgsEngine.getBoolean("-nosa")) {
			this.simulated_annealing = false;
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
							HiddenMarkovModel hmm = new HiddenMarkovModel(vcf_file, dat_file, rangeChr, dataB[b][0], dataB[b][1], simulated_annealing);
							if(hmm.isNullModel()) return;
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
							noLoci[b] = hmm.noLoci();
						} catch (Exception e) {
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
	}

	private void pileup() {
		// TODO Auto-generated method stub
		final List<int[][]> phasedData = new ArrayList<int[][]>();
		for(int i=0; i<dataB.length; i++) 
			phasedData.add(ballot(loadPhasedData(i)));
	}

	private int[][] ballot(final int[][][] phasedData) {
		// TODO Auto-generated method stub
		if(phasedData==null) return null;
		int[][] jtPerm = JohnsonTrotter.perm(Constants._ploidy_H);
		int[][] pivot = phasedData[0];
		for(int i=1; i<phasedData.length; i++) {
			int[][] exterior = phasedData[i];
			int[][] dist = new int[jtPerm.length][Constants._ploidy_H];
			for(int j=0; j<Constants._ploidy_H; j++) 
				dist[0][j] = hammingDist(pivot[j], exterior[j]);
			int[] tmp;
			for(int j=1; j<jtPerm.length; j++) {
				tmp = exterior[jtPerm[j][0]];
				exterior[jtPerm[j][0]] = exterior[jtPerm[j][1]];
				exterior[jtPerm[j][1]] = tmp;
				System.arraycopy(dist[j-1], 0, dist[j], 0, Constants._ploidy_H);
				dist[j][jtPerm[j][0]] = hammingDist(pivot[jtPerm[j][0]],exterior[jtPerm[j][0]]);
				dist[j][jtPerm[j][1]] = hammingDist(pivot[jtPerm[j][1]],exterior[jtPerm[j][1]]);
			}
			
			// swap (0, 1) to reset exterior
			tmp = exterior[0];
			exterior[0] = exterior[1];
			exterior[1] = tmp;
			
			// find minimum hamming distance
			int minD = Integer.MAX_VALUE, k = -1, tmpD;
			for(int j=0; j<jtPerm.length; j++) {
				if( (tmpD=sumInts(dist[j]))<minD ) {
					k = j;
					minD = tmpD;
				}
			}
			
			// swapping exterior phased data
			for(int j=1; j<=k; j++) {
				tmp = exterior[jtPerm[j][0]];
				exterior[jtPerm[j][0]] = exterior[jtPerm[j][1]];
				exterior[jtPerm[j][1]] = tmp;
			}
			
			// re-checking
			int d = 0;
			for(int j=0; j<Constants._ploidy_H; j++) 
				d += hammingDist(pivot[j], exterior[j]);
			
			if(minD!=d) throw new RuntimeException("!!!");
		
			myLogger.info("#### "+minD);
		}
		return null;
		
	}

	private int sumInts(int[] is) {
		// TODO Auto-generated method stub
		int sum = 0;
		for(final int i : is) sum += i;
		return sum;
	}

	private int hammingDist(int[] is, int[] is2) {
		// TODO Auto-generated method stub
		int d = 0;
		for(int i=0; i<is.length; i++)
			if(is[i]!=is2[i]) ++d;
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
		dataB = new int[nB][2];
		for(int i=0; i<nB; i++) {
			int from = rangeLowerBound+i*(block-overlap);
			dataB[i] = new int[]{from, from+block};
		}
		logLik = new double[nB][repeat];
		noLoci = new int[nB];
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
