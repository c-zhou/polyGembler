package cz1.test;

import cz1.data.DataCollection;

public class DataCollectionTest {
	
	public static void main(String[] args) {
		//String wd = ConstantsTest.wd;
		
		/**
		DataCollection.reformatSOAPOutputFile(
				"C:/Users/chenxi.zhou/Desktop/3_tetra_trifida.gbs.tassel.m50.vcf", 
				"C:/Users/chenxi.zhou/Desktop/3_tetra_trifida.gbs.tassel.m50.SOAPpopIndel",
				"C:/Users/chenxi.zhou/Desktop/3_tetra_trifida.gbs.tassel.m50.SOAPpopIndel.vcf", 
				"C:/Users/chenxi.zhou/Desktop/3_tetra_trifida.gbs.tassel.m50.SOAPpopIndel.reformat.vcf");
		*/
		
		//String wd = "C:\\Users\\chenxi.zhou\\Documents\\Sweetpotato\\polyhp\\data\\";
		//DataCollection.zip(wd, "simulation.gatk", 
		//		wd+"simulation.gatk.vcf");
		
		String wd = "C:\\Users\\chenxi.zhou\\Documents\\Sweetpotato\\polyhp\\data\\";
				DataCollection.zip(wd, "itr.fb.m50", 
						wd+"itr.fb.m50.vcf");
		
		//DataCollection.writeSOAPInputFileAll(wd+"3_tetra_trifida.gbs.tassel.m50.zip", 
		//		wd+"3_tetra_trifida.gbs.tassel.m50.SOAPpopIndel");
		/**
		String wd = args[0];
		
		String vcfFilePath = "potato.10m.vcf";
		String experimentId = vcfFilePath.replace(".vcf", "");
		DataCollection.zip(wd, experimentId, wd+vcfFilePath);

		vcfFilePath = "potato.5m.vcf";
		experimentId = vcfFilePath.replace(".vcf", "");
		DataCollection.zip(wd, experimentId, wd+vcfFilePath);
		

		vcfFilePath = "simulation.gbs.tassel.m50.vcf.gz";
		experimentId = vcfFilePath.replace(".vcf.gz", "");
		DataCollection.zip(wd, experimentId, wd+vcfFilePath);
		

		vcfFilePath = "trifida.gbs.tassel.itr10.m50.vcf.gz";
		experimentId = vcfFilePath.replace(".vcf.gz", "");
		DataCollection.zip(wd, experimentId, wd+vcfFilePath);
		*/
	}
	
}
