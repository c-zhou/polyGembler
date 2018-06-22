package cz1.graph.cluster;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;

import cz1.util.Utils;

public class SpectralClustering {
	private final double[][] mat;
	private final int k;
	private final int n;
	private final List<Set<Integer>> cluster;
	
	public SpectralClustering(final double[][] mat,
			final int k) {
		if(true) throw new RuntimeException("Yet to implemented!!!");
		this.mat = mat;
		this.k   = k;
		this.n   = mat.length;
		this.cluster = new ArrayList<Set<Integer>>();
		this.cluster();
	}

	private void cluster() {
		// TODO Auto-generated method stub
		double[][] diag = new double[n][n];
		for(int i=0; i<n; i++) {
			double s = 0;
			for(int j=0; j<n; j++) s+=mat[i][j];
			diag[i][i] = s;
		}
		double[][] laplacian = new double[n][n];
		for(int i=0; i<n; i++) 
			for(int j=0; j<n; j++)
				laplacian[i][j] = diag[i][j]-mat[i][j];
		EigenDecomposition eigenD = new EigenDecomposition(new Array2DRowRealMatrix(laplacian));
		double[][] eigen = new double[n][k];
		for(int i=n-k; i<n; i++) {
			double[] e = eigenD.getEigenvector(i).toArray();
			for(int j=0; j<n; j++)
				eigen[j][i-n+k] = e[j];
		}
		KMedoids kMedoids = new KMedoids(eigen, k);
		final List<Set<Integer>> clusters = kMedoids.getClusters();
		for(int i=0; i<k; i++)
			this.cluster.get(i).addAll(clusters.get(i));
		return;
	}
	
	public List<Set<Integer>> getClusters() {
		return this.cluster;
	}
}
