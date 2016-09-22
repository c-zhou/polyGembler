package cz1.test;

import java.util.Random;

import org.apache.commons.math3.stat.StatUtils;

import cz1.util.Algebra;
import cz1.util.IO;

public class SumExpsTest {
	
	public static void main(String[] args) {
		Random random = new Random(12345678);
		double[] e = new double[10000000];
		for(int i=0; i<e.length; i++)
			e[i] = Double.NEGATIVE_INFINITY;
		
		//double[] e = new double[]{Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY};
		//IO.print(e);
		
		long s1 = System.nanoTime();
		//double a1 = StatUtils.sum(e);
		double a1=0, b=0;
		for(int i=0; i<e.length; i++)
			//if(b!=0)
				//a1 = Math.pow(Math.exp(e[i]/2),2);
				a1 = Math.log(b);
		System.out.println(a1);
		long e1 = System.nanoTime();
		
		
		
		System.out.println(Double.MAX_EXPONENT);
		System.out.println(Double.MIN_EXPONENT);
		System.out.println(Double.MIN_VALUE);
		double m = Math.log(Double.MIN_VALUE);
		System.out.println(m);
		System.out.println(Math.exp(m-.1));

		
		

		
		


		long s3 = System.nanoTime();
		double a3 = Algebra.sumExps(e);
		long e3 = System.nanoTime();
		
		long s4 = System.nanoTime();
		double a4 = Algebra.sumExps(e);
		long e4 = System.nanoTime();
		
		
		for(int i=0; i<e.length; i++) e[i]*=-1;
		long s2 = System.nanoTime();
		double a2 = Algebra.sumExps(e);
		long e2 = System.nanoTime();
		
	

		
		System.out.println((e1-s1)+",2="+(e2-s2)+",3="+(e3-s3)+",4="+(e4-s4)+","+a2+","+a3);
	}

}
