package cz1.ngs.model;

import org.jgrapht.graph.DefaultWeightedEdge;

public class OverlapEdge extends DefaultWeightedEdge {
	/**
	 * 
	 */
	private static final long serialVersionUID = -5051004486663077338L;
	
	protected double olap = Double.NaN;
	protected String cigar = null;
	protected OverlapResult olapInfo = null;
	
	public double olap() {
		return this.olap;
	}
	
	public String cigar() {
		return this.cigar;
	}
	
	public OverlapResult olapInfo() {
		return this.olapInfo;
	}
	
	public void setOlap(double olap) {
		this.olap = olap;
	}
	
	public void setCigar(String cigar) {
		this.cigar = cigar;
	}
	
	public void setOlapInfo(OverlapResult olapInfo) {
		this.olapInfo = olapInfo;
	}
}
