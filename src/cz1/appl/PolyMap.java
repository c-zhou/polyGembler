package cz1.appl;

import java.util.Arrays;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

import cz1.data.DataCollection;
import cz1.data.DataEntry;
import cz1.model.HiddenMarkovModelDupSHEXA;
import cz1.model.HiddenMarkovModelDupSLHEX;
import cz1.model.HiddenMarkovModel;
import cz1.util.Constants;

public class PolyMap {

	public static void main(String[] args) {

		// create the command line parser
		CommandLineParser parser = new PosixParser();
		
		// create the Options
		Options options = new Options();
		options.addOption( "e", "experiment", true, "experiment name." );
		options.addOption( "w", "workspace", true, "directory contains input files." );
		options.addOption( "o", "output", true, "output directory name." );
		options.addOption( "c", "contig", true, "contig id.");
		options.addOption( "i", "max-iteration", true, "maximum iterations.");
		options.addOption( "d", "duplicate-hs", false, "duplicate hidden states.");
		options.addOption( "p", "ploidy", true, "ploidy of genome.");
		options.addOption( "P", "parents", true, "parental haplotypes.");
		options.addOption( "s", "seed", true, "random seed.");
		options.addOption( "S", "initial-seperation", true, "initial seperation.");
		options.addOption( "R", "direction", true, "direction.");
		options.addOption( "t", "train-exp", false, "train the exp parameter for RF.");
		options.addOption( "u", "unif-rand", false, "sample initial RF from a uniform distribution.");
		
		String experiment = null, workspace = null, output = null;
		String[] contig = null;
		double[] seperation = null;
		boolean dupS = false, trainExp = false, unifR = false;
		boolean[] reverse = null;
		int max_iter = 30;
		int ploidy = 2;
		
		try {
			// parse the command line arguments
			CommandLine line = parser.parse( options, args );
			if( line.hasOption("e") ) {
				experiment = line.getOptionValue('e');
			} else
				throw new RuntimeException("!!!");
			if(line.hasOption("w")) {
				workspace = line.getOptionValue("w");
			} else
				throw new RuntimeException("!!!");
			if(line.hasOption("o")) {
				output = line.getOptionValue("o");
			} else {
				output = workspace;
			}
			if(line.hasOption("c")) {
				contig = line.getOptionValue("c").split(":");
			} else
				throw new RuntimeException("!!!"); 
			if(line.hasOption("S")) {
				String[] ss = line.getOptionValue("S").split(":");
				if(ss.length<contig.length-1)
					throw new RuntimeException("!!!");
				seperation = new double[contig.length-1];
				for(int i=0; i<seperation.length; i++)
					seperation[i] = Double.parseDouble(ss[i]);
			} else {
				seperation = new double[contig.length-1];
				for(int i=0; i<seperation.length; i++)
					seperation[i] = Math.max(Math.round(
							Constants.rand.nextDouble()*
							Constants._max_initial_seperation),1);
			}
			if(line.hasOption("R")) {
				String[] dd = line.getOptionValue("R").split(":");
				if(dd.length<contig.length)
					throw new RuntimeException("!!!");
				reverse = new boolean[contig.length];
				for(int i=0; i<reverse.length; i++)
					reverse[i] = Boolean.parseBoolean(dd[i]);
			} else {
				reverse = new boolean[contig.length];
				Arrays.fill(reverse, false);
			}
			if(line.hasOption("i")) {
				max_iter = Integer.parseInt(line.getOptionValue("i"));
			}
			if(line.hasOption("d")) {
				dupS = true;
			}
			if(line.hasOption("t")) {
				trainExp = true;
			}
			if(line.hasOption("u")) {
				unifR = true;
			}
			if(line.hasOption("p")) {
				ploidy = Integer.parseInt(line.getOptionValue("p"));
				Constants._ploidy_H = ploidy;
				Constants._haplotype_z = ploidy*2;
			}
			if(line.hasOption("P")) {
				Constants._founder_haps = line.getOptionValue("P");
			} else
				throw new RuntimeException("!!!");
			if(line.hasOption("s")) {
				Constants.seed = Long.parseLong(line.getOptionValue("s"));
				Constants.setRandomGenerator();
			}
			boolean isRF = isRF(seperation);
			if(isRF && unifR)
				for(int i=0; i<seperation.length; i++) 
					seperation[i] = haldane(seperation[i]*
							Constants.rand.nextDouble());
			else if(isRF)
				for(int i=0; i<seperation.length; i++) 
					seperation[i] = haldane(seperation[i]);
		} catch( ParseException exp ) {
			System.out.println( "Unexpected exception:" + exp.getMessage() );
		}

		System.out.println(Constants.seed);
		
		DataEntry[] de = DataCollection.readDataEntry(
				workspace+
				Constants.file_sep+
				experiment+".zip", 
				contig);
		
		//de.print();
		//System.exit(1);
		
		/*
		de.get(5);
		*/
		//for(int j=0; j<de.length; j++)
		//	for(int i=24; i<de[j].getSample().length; i++)
		//		de[j].remove(i);
		/**/
		//de.remove(new int[]{6,7,8,9,10,12,13,14});
		
		HiddenMarkovModel hmm;
		if(Constants._ploidy_H<6)
			hmm = new HiddenMarkovModelDupSLHEX(de, seperation, reverse, trainExp);
		else
			hmm = new HiddenMarkovModelDupSHEXA(de, seperation, reverse, trainExp);
		hmm.print(true);
		
		double ll, ll0 = hmm.loglik();
		
		for(int i=0; i<max_iter; i++) {
			//System.out.println(i);
			hmm.train();
			ll = hmm.loglik();
			
			//hmm.print(true);
			
			//if(ll<ll0) {
			//	System.err.println("error!");
				//System.exit(1);
			//	break;
			//}
			
			if( ll0!=Double.NEGATIVE_INFINITY && 
					Math.abs((ll-ll0)/ll0)< 1e-4)
				break;
			ll0 = ll;
		}
		hmm.print(true);
		
		String contig_str = contig[0];
		for(int i=1; i<contig.length; i++) {
			if(contig_str.length()+contig[i].length()+8
					<=Constants.MAX_FILE_ID_LENGTH)
				contig_str += Constants.scaff_collapsed_str+contig[i];
			else {
				contig_str += Constants.scaff_collapsed_str+"etc"+contig.length;
				break;
			}
		}
		hmm.write(output, experiment.
				replace(".", "").
				replace("_", ""), contig_str);
	}

	private static double haldane(double r) {
		// TODO Auto-generated method stub
		return -.5*Math.log(1-2*r)/Constants._con_base_r;
	}

	private static boolean isRF(double[] rf) {
		// TODO Auto-generated method stub
		for(int i=0; i<rf.length; i++)
			if(rf[i]>1.0) return false;
		return true;
	}
}
