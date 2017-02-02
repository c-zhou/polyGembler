package cz1.appl;

import java.io.File;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;

import cz1.hmm.data.DataCollection;
import cz1.hmm.tools.PolyGembler;
import cz1.simulation.tools.GBSSimulator;
import cz1.simulation.tools.PopulationSimulator;

public class Main {
	
	protected final static Logger myLogger = 
			Logger.getLogger(Main.class);
	static {
		BasicConfigurator.configure();
	}
	
	public static void main(String[] args) {
	
		if(args.length<1) {
			printUsage();
			throw new RuntimeException("Undefined tool!!!");
		}
		String[] args2 = new String[args.length-1];
		System.arraycopy(args, 1, args2, 0, args2.length);
		switch(args[0].toLowerCase()) {
		case "popsimulation":
			PopulationSimulator popsimulator = new PopulationSimulator();
			popsimulator.setParameters(args2);
			popsimulator.run();
			break;
		case "gbssimulation":
			GBSSimulator gbssimulator = new GBSSimulator();
			gbssimulator.setParameters(args2);
			gbssimulator.run();
			break;
		case "datapreparation":
			String vcf = args[1];
			String zip = new File(vcf).getName().
					replaceAll(".vcf.gz$", "").
					replaceAll(".vcf$", "");
			String out_dir = args[2];
			DataCollection.zip(out_dir, zip, vcf);
			break;
		case "haplotypephasing":
			PolyGembler.main(args2);
			break;
		default:
			printUsage();
			throw new RuntimeException("Undefined tool!!!");
		}
	}
	
	private static void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " popsimulaion 		Simulate a full-sib mapping population. \n"
						+ " gbssimulation 		Simulate GBS data. \n"
						+ " datapreparation  	Prepare data for haplotype phasing. \n"
						+ " haplotypephasing  	Contig/scaffold haplotype construction from a mapping population.\n");
	}
}
