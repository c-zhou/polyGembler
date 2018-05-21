package cz1.appl;

import org.apache.log4j.Logger;

import cz1.tenx.tools.BAMstats;
import cz1.tenx.tools.GenomeComparison;
import cz1.tenx.tools.TenXMoleculeStats;
import cz1.tenx.tools.TenXMoleculeStatsZ;
import cz1.tenx.tools.GenomeComparisonZZ;
import cz1.tenx.tools.GenomeComparisonZ;
import cz1.tenx.tools.TenXSamtools;
import cz1.tenx.tools.TenXVcftools;
import cz1.tenx.tools.TenxMoleculeTools;


public class TenX {

	protected final static Logger myLogger = 
			Logger.getLogger(TenX.class);

	public static void main(String[] args) {

		if(args.length<1) {
			printUsage();
			throw new RuntimeException("Undefined tool!!!");
		}
		String[] args2 = new String[args.length-1];
		System.arraycopy(args, 1, args2, 0, args2.length);
		switch(args[0].toLowerCase()) {
		case "samtools":
			TenXSamtools tenXSamtools = new TenXSamtools();
			tenXSamtools.setParameters(args2);
			tenXSamtools.run();
			break;
		case "stats":
			TenXMoleculeStats tenXMStats = new TenXMoleculeStats();
			tenXMStats.setParameters(args2);
			tenXMStats.run();
			break;
		case "statsz":
			TenXMoleculeStatsZ tenXMStatsZ = new TenXMoleculeStatsZ();
			tenXMStatsZ.setParameters(args2);
			tenXMStatsZ.run();
			break;
		case "bamstats":
			BAMstats bamStats = new BAMstats();
			bamStats.setParameters(args2);
			bamStats.run();
			break;
		case "compare":
			GenomeComparison genomeComparison = new GenomeComparison();
			genomeComparison.setParameters(args2);
			genomeComparison.run();
			break;
		case "comparez":
			GenomeComparisonZ genomeComparisonZ = new GenomeComparisonZ();
			genomeComparisonZ.setParameters(args2);
			genomeComparisonZ.run();
			break;
		case "comparezz":
			GenomeComparisonZZ genomeComparisonZZ = new GenomeComparisonZZ();
			genomeComparisonZZ.setParameters(args2);
			genomeComparisonZZ.run();
			break;
		case "merge":
			TenxMoleculeTools tools = new TenxMoleculeTools();
			tools.setParameters(args2);
			tools.run();
			break;
		case "vcftools":
			TenXVcftools vcftools = new TenXVcftools();
			vcftools.setParameters(args2);
			vcftools.run();
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
						+ " samtools              TenX Samtools.\n"
						+ " stats                 TenX BAM file molecular statistics.\n"
						+ " compare               Compare two TenX data.\n"
						+ " statsz                TenX BAM file molecular statistics.\n"
						+ " comparez              Compare two TenX data.\n"
						+ " comparezz             Compare two TenX data.\n"
						+ " merge                 Merge molecule file.\n"
						+ " vcftools              VCF tools.\n"
						+ " bamstats              TenX BAM file statistics.\n\n");
	}
}
