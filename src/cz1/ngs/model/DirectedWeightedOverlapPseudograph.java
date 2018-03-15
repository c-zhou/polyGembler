package cz1.ngs.model;

import org.jgrapht.graph.DirectedWeightedPseudograph;

public class DirectedWeightedOverlapPseudograph<V> extends DirectedWeightedPseudograph<V, OverlapEdge> {
	/**
	 * 
	 */
	private static final long serialVersionUID = -4494520686889061388L;
	
	public DirectedWeightedOverlapPseudograph(Class<? extends OverlapEdge> edgeClass) {
		// TODO Auto-generated constructor stub
		super(edgeClass);
	}
	
	public void setEdgeOverlap(OverlapEdge edge, int olap) {
		edge.olap = olap;
	}
	
	public void setEdgeCigar(OverlapEdge edge, String cigar) {
		edge.cigar = cigar;
	}
	
	public int getEdgeOverlap(OverlapEdge edge) {
		return edge.olap;
	}
	
	public String getEdgeCigar(OverlapEdge edge) {
		return edge.cigar;
	}
}