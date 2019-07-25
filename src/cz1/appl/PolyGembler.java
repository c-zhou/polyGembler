package cz1.appl;

import org.apache.log4j.Logger;

import cz1.hmm.tools.AssemblyError;
import cz1.hmm.tools.DataPreparation;
import cz1.hmm.tools.Gembler;
import cz1.hmm.tools.Haplotyper;
import cz1.hmm.tools.MappingAnalysis;
import cz1.hmm.tools.NNsuperscaffold;
import cz1.hmm.tools.Pseudomolecule;
import cz1.hmm.tools.SinglePointAnalysis;
import cz1.hmm.tools.TwoPointAnalysis;
import cz1.simulation.tools.GBSSimulator;
import cz1.simulation.tools.PopulationSimulator;

public class PolyGembler {
	
	protected final static Logger myLogger = Logger.getLogger(PolyGembler.class);
	
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
			DataPreparation datapreparation = new DataPreparation();
			datapreparation.setParameters(args2);
			datapreparation.run();
			break;
		case "haplotyper":
			Haplotyper haplotyper = new Haplotyper();
			haplotyper.setParameters(args2);
			haplotyper.run();
			break;
        case "asmerr":
        	AssemblyError asmErr = new AssemblyError();
        	asmErr.setParameters(args2);
        	asmErr.run();
        	break;
        case "singlepoint":
            SinglePointAnalysis singlePoint = new SinglePointAnalysis();
            singlePoint.setParameters(args2);
            singlePoint.run();
            break;
        case "twopoint":
            TwoPointAnalysis twoPoint = new TwoPointAnalysis();
            twoPoint.setParameters(args2);
            twoPoint.run();
            break;
        case "map":
            MappingAnalysis mapping = new MappingAnalysis();
            mapping.setParameters(args2);
            mapping.run();
            break;
        case "superscaffold":
        	NNsuperscaffold superscaffold = new NNsuperscaffold();
        	superscaffold.setParameters(args2);
        	superscaffold.run();
        	break;
        case "chromosomer":
        	Pseudomolecule pseudoM = new Pseudomolecule();
        	pseudoM.setParameters(args2);
        	pseudoM.run();
        	break;
        case "gembler":
			Gembler gembler = new Gembler();
			gembler.setParameters(args2);
			gembler.run();
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
						+ " popsimulation       Simulate a full-sib mapping population.\n"
						+ " gbssimulation       Simulate GBS data.\n"
						+ " datapreparation     Prepare data for haplotype phasing.\n"
						+ " haplotyper          Haplotype phasing from a mapping population.\n"
						+ " asmerr              Correct assembly errors.\n"
                        + " singlepoint         Signle-point analysis: estimate recombination fraction between markers within contigs/scaffolds.\n"
                        + " twopoint            Two-point analysis: estimate pairwise recombination fraction between contigs/scaffolds.\n"
                        + " map                 Construct linkage maps.\n"
						+ " superscaffold       Construct superscaffold using nearest neighbour joining.\n"
						+ " chromosomer         Construct pseudo chromosomes. \n"
						+ " gembler             Run PolyGembler pipeline to construct genetic linkage maps/pseudomolecules.\n\n");
	}
}
