package cz1.hmm.tools;

import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.awt.print.PageFormat;
import java.awt.print.PrinterException;
import java.awt.print.PrinterJob;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;

import javax.imageio.ImageIO;
import javax.swing.JFrame;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.renjin.invoke.codegen.scalars.IntegerType;

import com.itextpdf.awt.PdfGraphics2D;
import com.itextpdf.text.Document;
import com.itextpdf.text.DocumentException;
import com.itextpdf.text.Rectangle;
import com.itextpdf.text.pdf.PdfContentByte;
import com.itextpdf.text.pdf.PdfTemplate;
import com.itextpdf.text.pdf.PdfWriter;

import cz1.hmm.data.DataCollection;
import cz1.hmm.data.DataEntry;
import cz1.hmm.model.HiddenMarkovModel;
import cz1.hmm.model.HiddenMarkovModelBWT;
import cz1.hmm.model.HiddenMarkovModelVBT;
import cz1.hmm.swing.HMMFrame;
import cz1.hmm.swing.HMMPanel;
import cz1.hmm.swing.Printer;
import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Constants.Field;
import cz1.util.Executor;
import cz1.util.Utils;

public class Resampler extends Executor {

	private HMMFrame hmmf = null;
	private HMMPanel hmmp = null;
	private String in_zip = null;
	private String out_prefix = null;
	private String[] scaff = null;
	private double[] seperation = null;
	private boolean trainExp = false;
	private boolean vbt = false;
	private boolean[] reverse = new boolean[]{false};
	private int max_iter = 100;
	private Field field = Field.PL;
	private String plot_pdf = null;
	private String expr_id = null;
	private int[] start_pos = null;
	private int[] end_pos = null;
	private String hmm_file = null;
	private int resampling = 100;
	private double loglik_diff = 0;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
							+" -i/--input                   Input zipped file.\n"
							+" -o/--prefix                  Output file location.\n"
							+" -ex/--experiment-id          Common prefix of haplotype files for this experiment.\n"
							+" -hf/--hmm-file               A zipped HMM file. If provided the initial transition and \n"
							+"                              and emission probabilities will be read from the file instead \n"
							+"                              of randomly selected. \n"
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
							+"                              option -D/--allele-depth and -L/--genetype likelihood.\n"
							+" -D/--allele-depth            Use allele depth to infer haplotypes. Mutually exclusive \n"
							+"                              with option -G/--genotype and -L/--genetype likelihood.\n"
							+" -L/--genotype-likelihood     Use genotype likelihoods to infer haplotypes. Mutually \n"
							+"                              exclusive with option -G/--genotype and -L/--allele-depth \n"
							+"                              (default).\n"
							+" -b/--segmental-kmeans        Use Viterbi training instead of Baum-Welch algorithm.\n"
							+" -ld/--loglike-diff           Log likelihood difference for Veterbi path (default 0)."
							+" -e/--train-exp               Re-estimate transition probabilities between founder/parental \n"
							+"                              haplotypes at each step.\n"
							+" -rs/--resampling             Resampling from the HMM (default 100).\n"
							+" -S/--random-seed             Random seed for this run.\n"
							+" -pp/--print-plot             Plot the hidden Markov model.\n"
							+" -sp/--save-plot              Save the plot as a pdf file. The file name should be provided here.\n"
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
			myArgsEngine.add("-L", "--genotype-likelihood", false);
			myArgsEngine.add("-b", "--segmental-kmeans", false);
			myArgsEngine.add("-ld", "--loglike-diff", true);
			myArgsEngine.add("-e", "--train-exp", false);
			myArgsEngine.add("-S", "--random-seed", true);
			myArgsEngine.add("-rs", "--resampling", true);
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
		
		if(myArgsEngine.getBoolean("-hf")) {
			hmm_file = myArgsEngine.getString("-hf");
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
			Constants.ploidy(Integer.parseInt(myArgsEngine.getString("-p")));
		}
		
		if(myArgsEngine.getBoolean("-f")) {
			Constants._founder_haps = myArgsEngine.getString("-f");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the parent samples (seperated by a \":\").");
		}
		
		if(myArgsEngine.getBoolean("-ld")) {
			loglik_diff = Double.parseDouble(myArgsEngine.getString("-ld"));
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
						Constants._max_initial_seperation),1);
		}
		boolean isRF = Constants.isRF(seperation);
		if(isRF) {
			for(int i=0; i<seperation.length; i++) 
				seperation[i] = Constants.haldane(seperation[i]);
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
		
		if(myArgsEngine.getBoolean("-L")) {
			field = Field.PL;
			i++;
		}
		if(i>1) throw new RuntimeException("Options -G/--genotype, "
				+ "-D/--allele-depth, and -L/--genotype-likelihood "
				+ "are exclusive!!!");
		
		if(myArgsEngine.getBoolean("-b")) {
			vbt = true;
			throw new RuntimeException("Viterbi training not supported yet!!!");
		}
		
		if(myArgsEngine.getBoolean("-S")) {
			Constants.seed = Long.parseLong(myArgsEngine.getString("-S"));
			Constants.setRandomGenerator();
		}
		
		if(myArgsEngine.getBoolean("-rs")) {
			resampling = Integer.parseInt(myArgsEngine.getString("-rs"));
		}
		
		if(myArgsEngine.getBoolean("-e")) {
			trainExp = true;
		}
		
		if(myArgsEngine.getBoolean("-pp")) {
			Constants.plot(true);
		}
		
		if(myArgsEngine.getBoolean("-sp")) {
			plot_pdf = myArgsEngine.getString("-sp");
			Constants.plot(true);
			Constants.printPlots(true);
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		myLogger.info("Random seed - "+Constants.seed);

		
		DataEntry[] de = start_pos==null ?
			DataCollection.readDataEntry(in_zip, scaff) :
			DataCollection.readDataEntry(in_zip, scaff, start_pos, end_pos);

				// DataEntry[] de = DataCollection.readDataEntry(in_zip, Constants._ploidy_H);
		final HiddenMarkovModel hmm = vbt ? 
			new HiddenMarkovModelVBT(de, seperation, reverse, trainExp, field):
			(hmm_file==null ?
			new HiddenMarkovModelBWT(de, seperation, reverse, trainExp, field):
			new HiddenMarkovModelBWT(de, seperation, reverse, trainExp, field, hmm_file));
		
		double ll, ll0 = hmm.loglik();
		
		String scaff_str = scaff[0]+
				(start_pos==null||start_pos[0]==Integer.MIN_VALUE?"":"_"+start_pos[0])+
				(end_pos==null||end_pos[0]==Integer.MAX_VALUE?"":"_"+end_pos[0]);
		for(int i=1; i<scaff.length; i++) {
			if(scaff_str.length()+scaff[i].length()+32<=Constants.MAX_FILE_ID_LENGTH)
				scaff_str += Constants.scaff_collapsed_str+scaff[i]+
				(start_pos==null||start_pos[i]==Integer.MIN_VALUE?"":"_"+start_pos[i])+
				(end_pos==null||end_pos[i]==Integer.MAX_VALUE?"":"_"+end_pos[i]);
			else {
				scaff_str += Constants.scaff_collapsed_str+"etc"+scaff.length;
				break;
			}
		}
		
		hmm.write(out_prefix, expr_id, scaff_str, resampling, loglik_diff);
	}
}
