package cz1.test;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.stat.StatUtils;

import cz1.util.Constants;
import cz1.util.IO;

public class BetaDistributionTest {
	public static void main(String[] args) {
		
		for(int i=0; i<100; i++) {
			double Jm = new BetaDistribution(Constants.rg, 1, 2).sample();
			System.err.println(Jm);
		}
		
		IO.print(StatUtils.mode(new double[]{0.1,0.2,0.3}));
	}

}
