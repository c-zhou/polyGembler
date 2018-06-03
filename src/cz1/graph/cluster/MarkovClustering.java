package cz1.graph.cluster;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.ArrayUtils;

import cz1.algebra.matrix.SparseMatrix;

public class MarkovClustering {
	private final SparseMatrix mat;
	private final double maxResidual;
	private final double pGamma;
	private final double loopGain;
	private final double maxZero;
	
	public MarkovClustering(SparseMatrix mat, 
			double maxResidual, 
			double pGamma, 
			double loopGain, 
			double maxZero) {
		this.mat = mat;
		this.maxResidual = maxResidual;
		this.pGamma = pGamma;
		this.loopGain = loopGain;
		this.maxZero = maxZero;
	}
	
	private double residual = 1.;
	private int step = 0;
	
	public void run() { 
        // add cycles 
		mat.setDiag(loopGain); 
        mat.normalise(); 
        
        while (residual > maxResidual) { 
        	++step;
        	mat.expand(); 
            residual = mat.inflate(pGamma, maxZero);
        }
    }
	
	public String progress() {
		return "step "+step+", residual "+residual;
	}

	final List<Set<Integer>> clusters = new ArrayList<Set<Integer>>();
	
	public void parse() {
		// TODO Auto-generated method stub	
		int n = mat.getCols();
		for(int i=0; i<n; i++) {
			if(mat.get(i, i)>0) {
				Set<Integer> clus = mat.getColumn(i);
				boolean add = true;
				for(int j=0; j<clusters.size(); j++) {
					if( clus.containsAll(clusters.get(j))&&
							clusters.get(j).containsAll(clus) ) {
						add = false;
						break;
					}
				}
				if(add) clusters.add(clus);
			}
		}
	}

	public List<Set<Integer>> getCluster() {
		// TODO Auto-generated method stub
		return this.clusters;
	}
}
