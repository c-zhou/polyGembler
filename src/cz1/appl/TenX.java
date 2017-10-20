package cz1.appl;

import org.apache.log4j.Logger;

import cz1.tenx.tools.GenomeTools;
import cz1.tenx.tools.TenXBamFilter;
import cz1.tenx.tools.TenXBamTools;
import cz1.tenx.tools.TenXFastqTools;


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
		case "fastqrename":
			TenXFastqTools tenXFastqTools = new TenXFastqTools();
			tenXFastqTools.setParameters(args2);
			tenXFastqTools.run();
			break;
		case "addreadbc":
			TenXBamTools tenXBamTools = new TenXBamTools();
			tenXBamTools.setParameters(args2);
			tenXBamTools.run();
			break;
		case "bamfilter":
			TenXBamFilter tenXBamFilter = new TenXBamFilter();
			tenXBamFilter.setParameters(args2);
			tenXBamFilter.run();
			break;
		case "uniqueregion":
			GenomeTools genomeTools = new GenomeTools();
			genomeTools.setParameters(args2);
			genomeTools.run();
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
						+ " fastqrename           Rename barcoded fastq files.\n"
						+ " addreadbc             Add barcode to aligned reads.\n"
						+ " bamfilter             Bam file filter.\n"
						+ " uniqueregion          Unique region in the reference genome.\n\n");
	}
}
