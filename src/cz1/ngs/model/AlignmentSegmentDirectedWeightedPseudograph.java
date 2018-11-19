package cz1.ngs.model;

import java.util.HashMap;
import java.util.Map;

import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.DirectedWeightedPseudograph;

public class AlignmentSegmentDirectedWeightedPseudograph<V extends Comparable<V>> extends DirectedWeightedPseudograph<V, DefaultWeightedEdge> {

	private final Map<V, CompoundAlignmentSegment> innerMap = new HashMap<V, CompoundAlignmentSegment>();
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1101335273597375221L;
	
	public AlignmentSegmentDirectedWeightedPseudograph(Class<? extends DefaultWeightedEdge> edgeClass) {
		// TODO Auto-generated constructor stub
		super(edgeClass);
	}
	
	public void addVertex(V vertex, CompoundAlignmentSegment seg) {
		if(!this.containsVertex(vertex)&&!this.addVertex(vertex)) 
			throw new RuntimeException("!!!");
		this.innerMap.put(vertex, seg);
	}
	
	public CompoundAlignmentSegment getAlignmentSegmentByVertex(V vertex) {
		return this.innerMap.get(vertex);
	}
}
