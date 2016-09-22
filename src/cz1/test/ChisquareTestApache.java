package cz1.test;

import java.util.Arrays;

import org.apache.commons.math3.stat.inference.ChiSquareTest;

public class ChisquareTestApache {

	public static void main(String[] args) {
		final double[] probs_uniform = new double[8];
		Arrays.fill(probs_uniform, .125);
		long[] observed = new long[] {87,70,68,90,77,60,68,110};
		double p = new ChiSquareTest().chiSquareTest(probs_uniform, observed);
		System.err.println(p);
	}
}
