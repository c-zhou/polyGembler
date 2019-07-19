package cz1.hmm.tools;

import java.io.File;
import java.util.Arrays;

import cz1.hmm.data.DataCollection;
import cz1.hmm.data.DataEntry;
import cz1.hmm.model.BaumWelchTrainer;
import cz1.hmm.model.ModelTrainer;
import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Constants.Field;
import cz1.util.Executor;

public class Haplotyper extends Executor {
	
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
	private int ploidy;
	private String[] parents;
	
	public Haplotyper() {}
	
	public Haplotyper(String in_zip,
			String out_prefix,
			String[] scaff,
			double[] seperation,
			boolean[] reverse,
			int max_iter,
			int ploidy,
			Field field) {
		this.in_zip = in_zip;
		this.out_prefix = out_prefix;
		this.scaff = scaff;
		this.seperation = seperation;
		this.reverse = reverse;
		this.max_iter = max_iter;
		this.ploidy = ploidy;
		this.field = field;
		this.expr_id = new File(in_zip).getName().
				replaceAll(".zip$", "").
				replace(".", "").
				replace("_", "");
	}
	
	public Haplotyper(String in_zip,
			String out_prefix,
			String[] scaff,
			double[] seperation,
			boolean[] reverse,
			int max_iter,
			int ploidy,
			Field field,
			String expr_id) {
		this.in_zip = in_zip;
		this.out_prefix = out_prefix;
		this.scaff = scaff;
		this.seperation = seperation;
		this.reverse = reverse;
		this.max_iter = max_iter;
		this.ploidy = ploidy;
		this.field = field;
		this.expr_id = expr_id;
	}
	
	public Haplotyper(String in_zip,
			String out_prefix, 
			String[] scaff,
			int ploidy, 
			Field field) {
		// TODO Auto-generated constructor stub
		this.in_zip = in_zip;
		this.out_prefix = out_prefix;
		this.scaff = scaff;
		this.ploidy = ploidy;
		this.field = field;
		this.expr_id = new File(in_zip).getName().
				replaceAll(".zip$", "").
				replace(".", "").
				replace("_", "");
	}
	
	public Haplotyper(String in_zip,
			String out_prefix, 
			String[] scaff,
			int ploidy, 
			Field field,
			String expr_id,
			int max_iter) {
		// TODO Auto-generated constructor stub
		this.in_zip = in_zip;
		this.out_prefix = out_prefix;
		this.scaff = scaff;
		this.ploidy = ploidy;
		this.field = field;
		this.expr_id = expr_id;
		this.max_iter = max_iter;
	}
	
	public Haplotyper(String in_zip,
			String out_prefix, 
			String scaff,
			String seperation,
			String reverse,
			int ploidy, 
			Field field,
			String expr_id,
			int max_iter) {
		// TODO Auto-generated constructor stub
		this.in_zip = in_zip;
		this.out_prefix = out_prefix;
		this.scaff = scaff.split(":");
		String[] s = seperation.split(":");
		this.seperation = new double[s.length];
		for(int i=0; i<s.length; i++) 
			this.seperation[i] = Double.parseDouble(s[i]);
		s = reverse.split(":");
		this.reverse = new boolean[s.length];
		for(int i=0; i<s.length; i++) 
			this.reverse[i] = Boolean.parseBoolean(s[i]);
		this.ploidy = ploidy;
		this.field = field;
		this.expr_id = expr_id;
		this.max_iter = max_iter;
	}

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
							+" -x/--max-iter                Maxmium rounds for EM optimization (default 100).\n"
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
							+"                              with option -G/--genotype (default).\n"
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
			myArgsEngine.add("-hf", "--hmm-file", true);
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
			myArgsEngine.add("-e", "--train-exp", false);
			myArgsEngine.add("-S", "--random-seed", true);
			myArgsEngine.add("-pp", "--print-plot", false);
			myArgsEngine.add("-sp", "--save-plot", true);
			myArgsEngine.parse(args);
		}
		
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
		
		if(myArgsEngine.getBoolean("-f")) {
			this.parents = myArgsEngine.getString("-f").split(":");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the parent samples (seperated by a \":\").");
		}
		
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
			for(int i=0; i<seperation.length; i++) 
				seperation[i] = -.5*Math.log(1-2*seperation[i])*10000000;
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
		
		if(i>1) throw new RuntimeException("Options -G/--genotype, "
				+ "-D/--allele-depth are mutually exclusive!!!");
		
		if(myArgsEngine.getBoolean("-S")) {
			Constants.seed = Long.parseLong(myArgsEngine.getString("-S"));
			Constants.setRandomGenerator();
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		myLogger.info("Random seed - "+Constants.seed);

		DataEntry[] de = start_pos==null ? DataCollection.readDataEntry(in_zip, scaff, ploidy) :
			DataCollection.readDataEntry(in_zip, scaff, start_pos, end_pos, ploidy);

		double ll, ll0;
		
		myLogger.info("=> STAGE I. training emission model with no transitions allowed.");
		final ModelTrainer model = new ModelTrainer(de, seperation, reverse, field, ploidy, parents);
		ll0 = Double.NEGATIVE_INFINITY;
		for(int i=0; i<max_iter; i++) {
			model.train();
			ll = model.loglik();
			myLogger.info("#iteration "+model.iteration()+": loglik "+ll);
			if(ll<ll0)
				throw new RuntimeException("Fatal error, likelihood decreased!");
			if( ll0!=Double.NEGATIVE_INFINITY && 
					Math.abs((ll-ll0)/ll0) < minImprov)
				break;
			ll0 = ll;
		}

		myLogger.info("=> STAGE II. training emission model with transitions allowed.");
		final BaumWelchTrainer model1 = BaumWelchTrainer.copyOf(model);
		ll0 = Double.NEGATIVE_INFINITY;
		for(int i=0; i<max_iter; i++) {
			model1.train();
			ll = model1.loglik();
			myLogger.info("#iteration "+model1.iteration()+": loglik "+ll);
			if(ll<ll0)
				throw new RuntimeException("Fatal error, likelihood decreased!");
			if( ll0!=Double.NEGATIVE_INFINITY && 
					Math.abs((ll-ll0)/ll0) < minImprov)
				break;
			ll0 = ll;
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
	}
}
