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
	
	public void setEdgeOverlapF(OverlapEdge edge, double olap) {
		edge.olapF = olap;
	}
	
	public void setEdgeOverlapR(OverlapEdge edge, double olap) {
		edge.olapR = olap;
	}
	
	public void setEdgeOverlap(OverlapEdge edge, double olap) {
		edge.olap = olap;
	}
	
	public void setEdgeCigar(OverlapEdge edge, String cigar) {
		edge.cigar = cigar;
	}
	
	public void setEdgeOverlapInfo(OverlapEdge edge, OverlapResult olapInfo) {
		edge.olapInfo = olapInfo;
	}
	
	public void setEdgeRealigned(OverlapEdge e) {
		// TODO Auto-generated method stub
		e.setRealigned();
	}
	
	public double getEdgeOverlapF(OverlapEdge edge) {
		return edge.olapF;
	}
	
	public double getEdgeOverlapR(OverlapEdge edge) {
		return edge.olapR;
	}
	
	public double getEdgeOverlap(OverlapEdge edge) {
		return edge.olap;
	}
	
	public String getEdgeCigar(OverlapEdge edge) {
		return edge.cigar;
	}
	
	public OverlapResult getEdgeOverlapInfo(OverlapEdge edge) {
		return edge.olapInfo;
	}

	public double getEdgeOverlap(V source, V target) {
		// TODO Auto-generated method stub
		return getEdge(source, target).olap;
	}

	public void setEdgeOverlap(V source, V target, double olap) {
		// TODO Auto-generated method stub
		getEdge(source, target).setOlap(olap);
	}
	
	/***
	public V getVertexId(SequenceVertex<V> vertex) {
		return vertex.id;
	}
	
	public String getVertexSeq(SequenceVertex<V> vertex) {
		return vertex.seq;
	}
	
	public int getVertexLen(SequenceVertex<V> vertex) {
		return vertex.len;
	}
	
	public boolean getVertexRev(SequenceVertex<V> vertex) {
		return vertex.rev;
	}
	
	public void setVertexSeq(SequenceVertex<V> vertex, String seq) {
		vertex.seq = seq;
		vertex.len = seq.length();
	}
	
	public void setVertexRev(SequenceVertex<V> vertex, boolean rev) {
		vertex.rev = rev;
	}
	**/
}