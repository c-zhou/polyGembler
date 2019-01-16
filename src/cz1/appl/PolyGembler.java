package cz1.appl;

import org.apache.log4j.Logger;

import cz1.gbs.tools.GBSpileup;
import cz1.hmm.tools.DataPreparation;
import cz1.hmm.tools.Gembler;
import cz1.hmm.tools.Haplotyper;
import cz1.hmm.tools.LinkageMapConstructor;
import cz1.hmm.tools.RFEstimatorML;
import cz1.hmm.tools.Resampler;
import cz1.hmm.tools.SuperScaffoldConstructor;
import cz1.hmm.tools.TwoPointConstructor;
import cz1.hmm.tools.VCFResampling;
import cz1.ngs.assembly.GenomeAssemblyTiling;
import cz1.ngs.tools.AStats;
import cz1.ngs.tools.Anchor;
import cz1.ngs.tools.Consensus;
import cz1.ngs.tools.Contigger;
import cz1.ngs.tools.Graphmap;
import cz1.ngs.tools.HetCorr;
import cz1.ngs.tools.Redundas;
import cz1.ngs.tools.Scaffolder;
import cz1.ngs.tools.Stitcher;
import cz1.simulation.tools.GBSSimulator;
import cz1.simulation.tools.PopulationSimulator;
import cz1.tenx.tools.TenXArcs;

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
		case "twopoint":
			TwoPointConstructor twopoint = new TwoPointConstructor();
			twopoint.setParameters(args2);
			twopoint.run();
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
		case "hetcorr":
			HetCorr hetCorr = new HetCorr();
			hetCorr.setParameters(args2);
			hetCorr.run();
			break;
		case "anchor":
			Anchor anchor = new Anchor();
			anchor.setParameters(args2);
			anchor.run();
			break;
		case "graphmap":
			Graphmap graphmap = new Graphmap();
			graphmap.setParameters(args2);
			graphmap.run();
			break;	
		case "contigger":
			Contigger contigger = new Contigger();
			contigger.setParameters(args2);
			contigger.run();
			break;	
		case "consensus":
			Consensus consensus = new Consensus();
			consensus.setParameters(args2);
			consensus.run();
			break;
		case "redundas":
			Redundas redundas = new Redundas();
			redundas.setParameters(args2);
			redundas.run();
			break;	
		case "a-stats":
			AStats astats = new AStats();
			astats.setParameters(args2);
			astats.run();
			break;
		case "scaffolder":
			Scaffolder scaffolder = new Scaffolder();
			scaffolder.setParameters(args2);
			scaffolder.run();
			break;
		case "tenxarcs":
			TenXArcs tenxarcs = new TenXArcs();
			tenxarcs.setParameters(args2);
			tenxarcs.run();
			break;
		case "stitcher":
			Stitcher stitcher = new Stitcher();
			stitcher.setParameters(args2);
			stitcher.run();
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
						+ " twopoint            Make two-point superscaffolds.\n"
						+ " map                 Contruct linkage maps.\n"
						+ " dataresampling      Resample from ZIP data.\n"
						+ " hetcorr             Heterozygosity and sequencing errors correction.\n"
						+ " anchor              Anchor contigs/scaffolds against a reference.\n"
						+ " graphmap            Map noisy long reads (Nanopore and Pacbio) to an assembly graph.\n"
						+ " contigger           Contigging using paired-end and/or long reads.\n"
						+ " scaffolder          Scaffolding using paired-end and/or long reads.\n"
						+ " tenxarcs            Scaffolding using 10x genomics' linked reads.\n"
						+ " stitcher            Scaffolding using a reference.\n"
                        + " consensus           Consensus with paired reads link analysis.\n"
						+ " redundas            Remove redundancies in the genome assembly.\n"
						+ " a-stats             Compute Myers' a-statistic for a set of contigs using read alignments.\n"
						+ " gembler             Run PolyGembler pipeline to construct genetic linkage maps/pseudomolecules.\n\n");
	}
}
