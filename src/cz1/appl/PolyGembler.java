package cz1.appl;

import org.apache.log4j.Logger;

import cz1.gbs.tools.GBSpileup;
import cz1.hmm.tools.DataPreparation;
import cz1.hmm.tools.Gembler;
import cz1.hmm.tools.Haplotyper;
import cz1.hmm.tools.LinkageMapConstructor;
import cz1.hmm.tools.RFEstimatorML;
import cz1.hmm.tools.RFEstimatorRS2;
import cz1.hmm.tools.RFEstimatorRS3;
import cz1.hmm.tools.Resampler;
import cz1.hmm.tools.SuperScaffoldConstructor;
import cz1.hmm.tools.VCFResampling;
import cz1.ngs.assembly.GenomeAssemblyTiling;
import cz1.ngs.tools.Consensus;
import cz1.ngs.tools.HetCorr;
import cz1.simulation.tools.GBSSimulator;
import cz1.simulation.tools.PopulationSimulator;

public class PolyGembler {
	
	protected final static Logger myLogger = 
			Logger.getLogger(PolyGembler.class);
	
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
			datapreparation.setParameters(args2);
			datapreparation.run();
			break;
		case "haplotyper":
			Haplotyper haplotyper = new Haplotyper();
			haplotyper.setParameters(args2);
			haplotyper.run();
			break;
		case "rfestimator":
			RFEstimatorML rfEstimator = new RFEstimatorML();
			rfEstimator.setParameters(args2);
			rfEstimator.run();
			break;
		case "rfestimator2":
			RFEstimatorRS3 rfEstimator2 = new RFEstimatorRS3();
			rfEstimator2.setParameters(args2);
			rfEstimator2.run();
			break;
		case "resample":
			Resampler resampler = new Resampler();
			resampler.setParameters(args2);
			resampler.run();
			break;
		case "gembler":
			Gembler gembler = new Gembler();
			gembler.setParameters(args2);
			gembler.run();
			break;
		case "tiling":
			GenomeAssemblyTiling tiling = new GenomeAssemblyTiling();
			tiling.setParameters(args2);
			tiling.run();
			break;
		case "nnss":
			SuperScaffoldConstructor nnss = new SuperScaffoldConstructor();
			nnss.setParameters(args2);
			nnss.run();
			break;
		case "map":
			LinkageMapConstructor lgc = new LinkageMapConstructor();
			lgc.setParameters(args2);
			lgc.run();
			break;
		case "dataresampling":
			VCFResampling dataResampler = new VCFResampling();
			dataResampler.setParameters(args2);
			dataResampler.run();
			break;
		case "datacorrection":
			HetCorr hetCorr = new HetCorr();
			hetCorr.setParameters(args2);
			hetCorr.run();
			break;
		case "consensus":
			Consensus consensus = new Consensus();
			consensus.setParameters(args2);
			consensus.run();
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
						+ " popsimulation       Simulate a full-sib mapping population. \n"
						+ " gbssimulation       Simulate GBS data. \n"
						+ " gbspileup           Variant calling from GBS data. \n"
						+ " datapreparation     Prepare data for haplotype phasing. \n"
						+ " haplotyper          Contig/scaffold haplotype construction from a mapping population.\n"
						+ " rfestimator         Recombination fraction estimation.\n"
						+ " resample            Haplotype resampling.\n"
						+ " tiling              Attempts to construct a tiling path out of the query assembly as \n"
						+ "                     mapped to the reference assembly using the LAST alignment results.\n"
						+ " nnss                Make 1-nearest neighbour superscaffolds.\n"
						+ " map                 Contruct linkage maps.\n"
						+ " dataresampling      Resample from ZIP data.\n"
						+ " datacorrection      Heterozygosity and sequencing errors correction.\n"
						+ " consensus           Consensus with paired reads link analysis.\n"
						+ " gembler             Run PolyGembler pipeline to construct genetic linkage maps/pseudomolecules.\n\n");
	}
}
