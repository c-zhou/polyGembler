package cz1.appl;

import org.apache.log4j.Logger;

import cz1.tenx.tools.TenXMoleculeStats;
import cz1.tenx.tools.TenXSamtools;


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
						+ " stats                 TenX BAM file molecular statistics.\n\n");
	}
}
