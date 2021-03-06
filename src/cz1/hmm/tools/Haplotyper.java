package cz1.hmm.tools;

import java.io.File;
import java.util.Arrays;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import cz1.hmm.data.DataCollection;
import cz1.hmm.data.DataEntry;
import cz1.hmm.model.BaumWelchTrainer;
import cz1.hmm.model.ModelTrainer;
import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Constants.Field;
import cz1.util.Executor;

public class Haplotyper extends Executor {
	
	private final static Logger myLogger = LogManager.getLogger(Haplotyper.class);
	
	private final static double minImprov = 1e-4;
	private final static double max_init_seperation = 1e7; 
	
	private String in_zip = null;
	private String out_prefix = null;
	private String[] scaff = null;
	private double[] seperation = null;
	private boolean[] reverse = new boolean[]{false};
	private int max_iter = 1000;
	private Field field = Field.AD;
	private String expr_id = null;
	private int[] start_pos = null;
	private int[] end_pos = null;
	private int ploidy = 2;
	private String[] parents;

	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
							+" -i/--input                   Input zipped file.\n"
							+" -o/--prefix                  Output file location.\n"
							+" -ex/--experiment-id          Common prefix of haplotype files for this experiment.\n"
							+" -c/--scaffold                The scaffold/contig/chromosome id will run.\n"
							+" -cs/--start-position         The start position of the scaffold/contig/chromosome.\n"
							+" -ce/--end-position           The end position of the scaffold/contig/chromosome.\n"
							+" -x/--max-iter                Maxmium rounds for EM optimization (default 1000).\n"
							+" -p/--ploidy                  Ploidy of genome (default 2).\n"
							+" -f/--parent                  Parent samples (separated by a \":\").\n"
							+" -s/--initial-seperation      Initialisations of distances between the adjacent scaffolds \n"
							+"                              if multiple scaffolds will be jointly inferred. The separation \n"
							+"                              could be either physical distances or recombination frequencies, \n"
							+"                              i.e., if all values provided is below 0.5, the \n"
							+"                              program will take them as recombination frequencies. \n"
							+"                              Distances should be separated by \":\".\n"
							+" -r/--reverse                 Take either 'true' or 'false', indicating whetherr the \n"
							+"                              scaffold is reversed before inferring haplotypes. Multiple \n"
							+"                              scaffolds are separated by \":\".\n"
							+" -G/--genotype                Use genotypes to infer haplotypes. Mutually exclusive with \n"
							+"                              option -D/--allele-depth.\n"
							+" -D/--allele-depth            Use allele depth to infer haplotypes. Mutually exclusive \n"
							+"                              with option -G/--genotype.(default)\n"
							+" -S/--random-seed             Random seed for this run.\n"
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
			myArgsEngine.add("-i", "--input", true);
			myArgsEngine.add("-o", "--prefix", true);
			myArgsEngine.add("-ex", "--experiment-id", true);
			myArgsEngine.add("-c", "--scaffold", true);
			myArgsEngine.add("-cs", "--start-position", true);
			myArgsEngine.add("-ce", "--end-position", true);
			myArgsEngine.add("-x", "--max-iter", true);
			myArgsEngine.add("-p", "--ploidy", true);
			myArgsEngine.add("-f", "--parent", true);
			myArgsEngine.add("-s", "--initial-seperation", true);
			myArgsEngine.add("-r", "--reverse", true);
			myArgsEngine.add("-G", "--genotype", false);
			myArgsEngine.add("-D", "--allele-depth", false);
			myArgsEngine.add("-S", "--random-seed", true);
		}
		myArgsEngine.parse(args);
		
		if(myArgsEngine.getBoolean("-i")) {
			in_zip = myArgsEngine.getString("-i");
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
			expr_id = new File(in_zip).getName().
					replaceAll(".zip$", "").
					replace(".", "").
					replace("_", "");
		}
		
		if(myArgsEngine.getBoolean("-c")) {
			scaff = myArgsEngine.getString("-c").split(":");
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify the scaffold(s).");
		}
		
		if( myArgsEngine.getBoolean("-cs") ^ myArgsEngine.getBoolean("-ce") ) {
			printUsage();
			throw new IllegalArgumentException("Need to specify both start position "
					+ "and end position or none of them.");
		}
		
		if(myArgsEngine.getBoolean("-cs")) {
			String[] s = myArgsEngine.getString("-cs").split(":");
			start_pos = new int[scaff.length];
			int n = Math.min(scaff.length, s.length);
			for(int i=0; i!=n; i++) 
				start_pos[i] = Integer.parseInt(s[i]);
			for(int i=n; i<start_pos.length; i++)
				start_pos[i] = Integer.MIN_VALUE;
		}
		
		if(myArgsEngine.getBoolean("-ce")) {
			String[] s = myArgsEngine.getString("-ce").split(":");
			end_pos = new int[scaff.length];
			int n = Math.min(scaff.length, s.length);
			for(int i=0; i!=n; i++) 
				end_pos[i] = Integer.parseInt(s[i]);
			for(int i=n; i<end_pos.length; i++)
				end_pos[i] = Integer.MAX_VALUE;
		}
		
		if(myArgsEngine.getBoolean("-x")) {
			max_iter = Integer.parseInt(myArgsEngine.getString("-x"));
		}
		
		if(myArgsEngine.getBoolean("-p")) {
			this.ploidy = Integer.parseInt(myArgsEngine.getString("-p"));
		}
		
		this.parents = new String[2];
		if(myArgsEngine.getBoolean("-f")) {
			String[] s = myArgsEngine.getString("-f").split(":");
			if(s.length>2) 
				throw new IllegalArgumentException("Please specify at most TWO parent samples. Seperated by a \":\" if necessary.");
			for(int i=0; i<s.length; i++) parents[i] = s[i];
		}
		
		if(scaff.length>1) {
			if(myArgsEngine.getBoolean("-s")) {
				String[] ss = myArgsEngine.getString("-s").split(":");
				if(ss.length<scaff.length-1)
					throw new RuntimeException("Number of scaffolds does not match number of initial seperations!!!");
				seperation = new double[scaff.length-1];
				for(int i=0; i<seperation.length; i++)
					seperation[i] = Double.parseDouble(ss[i]);
			} else {
				seperation = new double[scaff.length-1];
				for(int i=0; i<seperation.length; i++)
					seperation[i] = Math.max(Math.round(
							Constants.rand.nextDouble()*
							max_init_seperation),1);
			}
			boolean isRF = true;
			for(int i=0; i<seperation.length; i++)
				if(seperation[i]>1.0)
					isRF = false;
			if(isRF) {
				for(int i=0; i<seperation.length; i++) {
					if(seperation[i]>=0.5) seperation[i] = 0.4999999;
					if(seperation[i]<=0)   seperation[i] = 0.0000001;
					seperation[i] = -.5*Math.log(1-2*seperation[i])*100000000;
				}
			}
			
			myLogger.info("##using initial seperation: ");
			for(int i=0; i<seperation.length; i++)
				myLogger.info(seperation[i]);
			myLogger.info("####");
		}
		
		if(myArgsEngine.getBoolean("-r")) {
			String[] dd = myArgsEngine.getString("-r").split(":");
			if(dd.length<scaff.length)
				throw new RuntimeException("Number of scaffolds does not match number of reverses!!!");
			reverse = new boolean[scaff.length];
			for(int i=0; i<reverse.length; i++)
				reverse[i] = Boolean.parseBoolean(dd[i]);
		} else {
			reverse = new boolean[scaff.length];
			Arrays.fill(reverse, false);
		}
		
		int i = 0;
		if(myArgsEngine.getBoolean("-G")) {
			field = Field.GT;
			i++;
		}
		
		if(myArgsEngine.getBoolean("-D")) {
			field = Field.AD;
			i++;
		}
		
		if(i>1) throw new RuntimeException("Options -G/--genotype and "
				+ "-D/--allele-depth are mutually exclusive.");
		
		if(myArgsEngine.getBoolean("-S")) {
			Constants.seed = Long.parseLong(myArgsEngine.getString("-S"));
			Constants.setRandomGenerator();
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
        long currentNanoTime1 = System.nanoTime();
        Runtime jr = Runtime.getRuntime();

		myLogger.info("Random seed - "+Constants.seed);

		DataEntry[] de = start_pos==null ? DataCollection.readDataEntry(in_zip, scaff, ploidy) :
			DataCollection.readDataEntry(in_zip, scaff, start_pos, end_pos, ploidy);

		double ll, ll0;
		
		myLogger.info("=> STAGE I. training emission model with no transitions allowed.");
		final ModelTrainer model = new ModelTrainer(de, seperation, reverse, field, ploidy, parents);
		
		if(!model.runnable()) return;
		
		ll0 = Double.NEGATIVE_INFINITY;
		for(int i=0; i<max_iter; i++) {
			model.train();
			ll = model.loglik();
			if( ll==0 || ll0!=Double.NEGATIVE_INFINITY && (ll0-ll)/ll0 < minImprov)
				break;
			ll0 = ll;
			myLogger.info("#iteration "+model.iteration()+": loglik "+ll);
		}

		myLogger.info("=> STAGE II. training emission model with transitions allowed.");
		final BaumWelchTrainer model1 = BaumWelchTrainer.copyOf(model);
		
		if(!model1.runnable()) return;
		
		ll0 = Double.NEGATIVE_INFINITY;
		for(int i=0; i<max_iter; i++) {
			model1.train();
			ll = model1.loglik();
			if( ll==0 || ll0!=Double.NEGATIVE_INFINITY && (ll0-ll)/ll0 < minImprov)
				break;
			ll0 = ll;
			myLogger.info("#iteration "+model1.iteration()+": loglik "+ll);
		}
		
		String scaff_str = scaff[0]+
				(start_pos==null||start_pos[0]==Integer.MIN_VALUE?"":"_"+start_pos[0])+
				(end_pos==null||end_pos[0]==Integer.MAX_VALUE?"":"_"+end_pos[0]);
		for(int i=1; i<scaff.length; i++) {
			if(scaff_str.length()+scaff[i].length()+32<=Constants.MAX_FILE_ID_LENGTH)
				scaff_str += Constants.collapsed_str+scaff[i]+
				(start_pos==null||start_pos[i]==Integer.MIN_VALUE?"":"_"+start_pos[i])+
				(end_pos==null||end_pos[i]==Integer.MAX_VALUE?"":"_"+end_pos[i]);
			else {
				scaff_str += Constants.collapsed_str+"etc"+scaff.length;
				break;
			}
		}
		
		model1.write(out_prefix, expr_id, scaff_str);
		
        long currentNanoTime2 = System.nanoTime();

        jr.gc();
        long totalMemory = jr.totalMemory();
        long freeMemory = jr.freeMemory();

        myLogger.error("Memory Usage: "+(totalMemory-freeMemory)/1024d/1024d+" megabytes.");
        myLogger.error("Time Usage: "+(currentNanoTime2-currentNanoTime1)/1e9+" seconds.");
	}
}
