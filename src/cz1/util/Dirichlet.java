package cz1.util;

import java.util.Arrays;

import org.apache.commons.math3.stat.StatUtils;

import cern.jet.random.Gamma;
import cern.jet.random.engine.RandomEngine;

public class Dirichlet extends Sampler {
	/**
	 * 
	 */
	private static final long serialVersionUID = -297762133965462712L;
	Gamma[] g;
	
	public Dirichlet(double[] dist, double u){
		super(dist, u);
		this.g = new Gamma[dist.length];
		set(u);
	}

	public Dirichlet(Double[] dist, double u){
		super(dist, u);
		this.g = new Gamma[dist.length];
		set(u);
	}
	
	public Dirichlet(int length, double u) {
		this(getUniform(length),u);
	}

	private static double[] getUniform(int length) {
		double[] res = new double[length];
		Arrays.fill(res,1.0/(double)length);
		return res;
	}

	public void set(double u){
		if(u==Double.POSITIVE_INFINITY){
			// throw new RuntimeException("should not set u to pos inf");
			// return;
		}
		for(int i=0; i<g.length; i++){
			if(dist[i]>0){
				if(g[i]!=null) g[i].setState(dist[i]*u, 1);
				g[i] = new Gamma(dist[i]*u, 1, re) ;
			}
		}
	}
	
	RandomEngine re = new RandomEngine(){
		private static final long serialVersionUID = 
				-4676280456745648481L;

		public int nextInt() {
			return Constants.rand.nextInt();
		}
	};
	
	public Double[] sample(double s){
		if(u==Double.POSITIVE_INFINITY) return dist;
		Double[] res = new Double[dist.length];
		double sum=0;
		for(int i=0; i<res.length; i++){
			res[i] = g[i]==null ? 0 : 
				Math.max(Constants.eps,g[i].nextDouble());
			sum+=res[i];
		}
		for(int i=0; i<res.length; i++){
			res[i] = res[i]/sum*s;
		}
		return res;
	}

	public Double[] sample(){
		return sample(1.0);
	}
	
	public double u() {
		return u;
	}
	
	public static double logDensity(double[] alpha, double[] x) {
		if(StatUtils.sum(x)!=1) return 0;
		double A = org.apache.commons.math3.special.
				Gamma.logGamma(StatUtils.sum(alpha));
		double B = 0;
		for(int i=0; i<alpha.length; i++) 
			B += org.apache.commons.math3.special.
				Gamma.logGamma(alpha[i]);
		double C = 0;
		for(int i=0; i<x.length; i++)
			C += (alpha[i]-1)*Math.log(x[i]);
		
		System.out.println(A);
		System.out.println(B);
		System.out.println(C);
		return A-B+C;
	}
}

