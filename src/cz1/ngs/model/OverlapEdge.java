package cz1.ngs.model;

import org.jgrapht.graph.DefaultWeightedEdge;

public class OverlapEdge extends DefaultWeightedEdge {
	/**
	 * 
	 */
	private static final long serialVersionUID = -5051004486663077338L;
	
	protected double olap  = Double.NaN;
	protected double olapF = Double.NaN;
	protected double olapR = Double.NaN;
	protected String cigar = null;
	protected OverlapResult olapInfo = null;
	
	protected boolean realigned = false;
	
	public double olap() {
		return this.olap;
	}
	
	public double olapF() {
		return this.olapF;
	}
	
	public double olapR() {
		return this.olapR;
	}
	
	public String cigar() {
		return this.cigar;
	}
	
	public OverlapResult olapInfo() {
		return this.olapInfo;
	}
	
	public void setOlap(double olap) {
		// if(olap<0) throw new RuntimeException("!!!");
		this.olap = olap;
	}
	
	public void setOlapF(double olap) {
		this.olapF = olap;
	}
	
	public void setOlapR(double olap) {
		this.olapR = olap;
	}
	
	public void setCigar(String cigar) {
		this.cigar = cigar;
	}
	
	public void setOlapInfo(OverlapResult olapInfo) {
		this.olapInfo = olapInfo;
	}
	
	public void setRealigned() {
		this.realigned = true;
	}
	
	public void setRealigned(boolean realigned) {
		this.realigned = realigned;
	}
	
	public boolean isRealigned() {
		return this.realigned;
	}
}
