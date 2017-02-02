package cz1.appl;

import java.io.File;

import cz1.data.DataCollection;
import cz1.simulation.tools.GBSSimulator;
import cz1.simulation.tools.PopulationSimulator;
import cz1.tools.PolyGembler;

public class Main {

	public static void main(String[] args) {
		String[] args2 = new String[args.length-1];
		System.arraycopy(args, 1, args2, 0, args2.length);
		switch(args[0].toLowerCase()) {
		case "popsimulaion":
			PopulationSimulator popsimulator = new PopulationSimulator();
			popsimulator.setParameters(args2);
			popsimulator.run();
			break;
		case "gbssimulaion":
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
			throw new RuntimeException("Undefined tool!!!");
		}
	}
}
