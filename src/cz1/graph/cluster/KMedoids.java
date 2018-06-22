package cz1.graph.cluster;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class KMedoids {
	private final double[][] mat;
	private final int k;
	private final int n;
	private final int maxIter = 10000;
	private final Set<Integer> medoidS;
	private final int[] medoids;
	private final List<Set<Integer>> cluster;
	private double cost = 0;
	private int optscale = 1;
	
	public KMedoids(final double[][] mat,
			final int k) {
		this(mat, k, false);
	}
	
	public KMedoids(final double[][] mat,
			final int k,
			final boolean maximise) {
		if(true) throw new RuntimeException("Yet to implemented!!!");
		this.mat = mat;
		this.k   = k;
		this.n   = mat.length;
		this.medoidS = new HashSet<Integer>();
		this.medoids = new int[k];
		this.cluster = new ArrayList<Set<Integer>>();
		for(int i=0; i<k; i++) cluster.add(new HashSet<Integer>());
		if(maximise) this.optscale = -1;
		this.cluster();
	}

	private void cluster() {
		// TODO Auto-generated method stub
		final List<Integer> datP = new ArrayList<Integer>();
		for(int i=0; i<n; i++) datP.add(i);
		Collections.shuffle(datP);
		for(int i=0; i<k ;i++) {
			medoidS.add(datP.get(i));
			medoids[i] = datP.get(i);
		}
		assign(this);
	
		for(int iter=0; iter<maxIter; iter++) {
			KMedoids best = this;
			for(int i=0; i<k; i++) {
				int medoid = medoids[i];
				for(int j=0; j<n; j++) {
					if(medoidS.contains(medoid)) continue;
					KMedoids copyOfThis = this.copy();
					copyOfThis.medoidS.remove(medoid);
					copyOfThis.medoidS.add(j);
					copyOfThis.medoids[i] = j;
					assign(copyOfThis);
					if(copyOfThis.cost<best.cost) best = copyOfThis;
				}
			}
			
			if(best.cost>=this.cost) break;
			this.copy(best);
		}
	}
	
	private void copy(KMedoids instance) {
		// TODO Auto-generated method stub
		this.medoidS.clear();
		this.medoidS.addAll(instance.medoidS);
		System.arraycopy(instance.medoids, 0, 
				this.medoids, 0, k);
		for(int i=0; i<k; i++) {
			Set<Integer> clus = this.cluster.get(i);
			clus.clear();
			clus.addAll(instance.cluster.get(i));
		}
		this.cost = instance.cost;
	}

	private KMedoids copy() {
		// TODO Auto-generated method stub
		KMedoids copyOfThis = new KMedoids(this.mat, 
				this.k, this.optscale==-1);
		copyOfThis.medoidS.addAll(this.medoidS);
		System.arraycopy(this.medoids, 0,
				copyOfThis.medoids, 0, k);
		return copyOfThis;
	}

	private static void assign(KMedoids instance) {
		// TODO Auto-generated method stub
		for(int i=0; i<instance.n; i++) {
			if(instance.medoidS.contains(i)) continue;
			double c = Double.MAX_VALUE;
			double d;
			int assign = -1;
			for(int j=0; j<instance.k; j++) {
				d = instance.optscale*instance.mat[i][instance.medoids[j]];
				if(d<c) {
					c = d;
					assign = j;
				}
			}
			instance.cost += c;
			instance.cluster.get(assign).add(i);
		}
	}
	
	public int[] getMedoids() {
		return this.medoids;
	}
	
	public List<Set<Integer>> getClusters() {
		List<Set<Integer>> clusters = new ArrayList<Set<Integer>>();
		for(int i=0; i<k; i++) {
			Set<Integer> clus = new HashSet<Integer>();
			clus.add(this.medoids[i]);
			clus.addAll(this.cluster.get(i));
		}
		return clusters;
	}
	
	public double getCost() {
		return optscale*this.cost;
	}
}
