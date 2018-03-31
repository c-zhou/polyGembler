package cz1.ngs.model;

import org.jgrapht.graph.DirectedWeightedPseudograph;

public class TraceableDirectedWeightedPseudograph<V extends Comparable<V>> extends DirectedWeightedPseudograph<TraceableVertex<V>, TraceableEdge> {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1101335273597375221L;

	public TraceableDirectedWeightedPseudograph(Class<? extends TraceableEdge> edgeClass) {
		super(edgeClass);
		// TODO Auto-generated constructor stub
	}
	
	public void changeVertexStatus(TraceableVertex<V> vertex) {
		vertex.changeStatus();
	}
	
	public void setVertexStatus(TraceableVertex<V> vertex, boolean status) {
		vertex.setStatus(status);
	}
	
	public boolean getVertexStatus(TraceableVertex<V> vertex) {
		return vertex.getStatus();
	}
	
	public void setVertexBackTrace(TraceableVertex<V> vertex, TraceableVertex<V> backtrace) {
		vertex.setBackTrace(backtrace);
	}
	
	public void setVertexScore(TraceableVertex<V> vertex, double score) {
		vertex.setScore(score);
	}
	
	public void setVertexPenalty(TraceableVertex<V> vertex, double penalty) {
		vertex.setPenalty(penalty);
	}

	public void setVertexSAMSegment(TraceableVertex<V> vertex, SAMSegment segment) {
		vertex.setSAMSegment(segment);
	}
	
	public TraceableVertex<V> getVertexBackTrace(TraceableVertex<V> vertex) {
		return vertex.getBackTrace();
	}
	
	public double getVertexPenalty(TraceableVertex<V> vertex) {
		return vertex.getPenalty();
	}
	
	public double getVertexScore(TraceableVertex<V> vertex) {
		return vertex.getScore();
	}
	
	public V getVertexId(TraceableVertex<V> vertex) {
		return vertex.getId();
	}
	
	public SAMSegment getVertexSAMSegment(TraceableVertex<V> vertex) {
		return vertex.getSAMSegment();
	}

}
