package cz1.util;

import java.util.Random;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

public class Constants {

	public static long seed = System.nanoTime();
	public static Random rand = new Random(seed);
	public static RandomGenerator rg = new Well19937c(seed);
	
	public final static double MIN_EXP_DOUBLE = Math.log(Double.MIN_VALUE);
	public final static double MAX_EXP_DOUBLE = Math.log(Double.MAX_VALUE);
	public final static double threshMax = 1e-100d;
	public final static double threshMin = 1e-200d;
	public final static double logThreshMax = Math.log(threshMax);
	public final static double logThreshMin = Math.log(threshMin);
	public static enum Field { PL, AD, GT, GL }
	public final static String collapsed_str = "____";
	public final static int MAX_FILE_ID_LENGTH = 128;

	public static void seeding(long s) {
		seed = s;
		setRandomGenerator();
	}
	
	public static void seeding() {
		seeding(System.nanoTime());
	}
		
	public static void setRandomGenerator() {
		// TODO Auto-generated method stub
		rand = new Random(seed);
		rg = new Well19937c(seed);
	}

	public static void throwRuntimeException(String message) {
		// TODO Auto-generated method stub
		throw new RuntimeException(message);
	}
}



