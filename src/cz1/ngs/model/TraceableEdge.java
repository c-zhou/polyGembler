package cz1.ngs.model;

import org.jgrapht.graph.DefaultWeightedEdge;

public class TraceableEdge extends DefaultWeightedEdge {

	/**
	 * 
	 */
	private static final long serialVersionUID = -9185940962140705345L;

	private double penalty = Double.POSITIVE_INFINITY;
	private double score   = Double.NEGATIVE_INFINITY;
	
	public void setPenalty(double penalty) {
		this.penalty = penalty;
	}
	
	public void setScore(double score) {
		this.score = score;
	}
	
	public double getPenalty() {
		return this.penalty;
	}
	
	public double getScore() {
		return this.score;
	}
}
