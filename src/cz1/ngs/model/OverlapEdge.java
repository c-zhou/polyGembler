package cz1.ngs.model;

import org.jgrapht.graph.DefaultWeightedEdge;

public class OverlapEdge extends DefaultWeightedEdge {
	/**
	 * 
	 */
	private static final long serialVersionUID = -5051004486663077338L;
	
	protected int olap = 0;
	protected String cigar = null;
	
	public int olap() {
		return this.olap;
	}
	
	public String cigar() {
		return this.cigar;
	}
	
	public void setOlap(int olap) {
		this.olap = olap;
	}
	
	public void setCigar(String cigar) {
		this.cigar = cigar;
	}
}
