package cz1.util;

import java.io.Serializable;

import org.apache.commons.lang3.ArrayUtils;

public abstract class Sampler implements Serializable {

	private static final long serialVersionUID = -2023477953560276351L;	
	protected Double[] dist;
	protected final double u;
	
	public Sampler(double[] dist, double u){
		this.dist = ArrayUtils.toObject(dist);
		this.u = u;
	}
	
	public Sampler(float[] dist, double u){
		this.dist = new Double[dist.length];
		for(int i=0; i<dist.length; i++){
			this.dist[i] = (double) dist[i];
		}
		this.u = u;
	}

	public Sampler(Double[] dist, double u){
		this.dist = dist;
		this.u = u;
	}

	public abstract Double[] sample();
}
