package cz1.appl;

import java.io.File;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;

import cz1.gbs.tools.GBSpileup;
import cz1.hmm.data.DataCollection;
import cz1.hmm.tools.DataPreparation;
import cz1.hmm.tools.Gembler;
import cz1.simulation.tools.GBSSimulator;
import cz1.simulation.tools.PopulationSimulator;

public class PolyGembler {
	
	protected final static Logger myLogger = 
			Logger.getLogger(PolyGembler.class);
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
		case "gbspileup":
			GBSpileup gbspileup = new GBSpileup();
			gbspileup.setParameters(args2);
			gbspileup.run();
			break;
		case "datapreparation":
			DataPreparation datapreparation = new DataPreparation();
			datapreparation.run();
			break;
		case "haplotypephasing":
			Gembler.main(args2);
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
						+ " gbspileup 			Variant calling from GBS data. \n"
						+ " datapreparation  	Prepare data for haplotype phasing. \n"
						+ " haplotypephasing  	Contig/scaffold haplotype construction from a mapping population.\n");
	}
}
