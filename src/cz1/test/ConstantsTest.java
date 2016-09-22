package cz1.test;

import cz1.data.DataCollection;

public class ConstantsTest {
	
	public final static String wd = "C:/Users/chenxi.zhou/Documents/Sweetpotato/polyhp/data/";
	public final static String zip = "hexa.simulation.pgsim.f574.m20.zip";
	//public final static String zip = "simulation.gbs.tassel.m50.zip";
	//public final static String zip = "tetra.simulation.zip";
	//public final static String contig = "chr2_24149038_29137011";
	public final static String contig = DataCollection.getContigIDFromIndex(wd+zip, 1);
	//public final static String contig = "5364";
}
