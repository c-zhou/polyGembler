package cz1.ngs.model;

import com.google.common.collect.DiscreteDomain;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

public class TraceableVertex<V extends Comparable<V>> implements Comparable<TraceableVertex<V>> {
	private final V id;
	private TraceableVertex<V> backtrace = null;
	private double penalty = Double.POSITIVE_INFINITY;
	private double score   = Double.NEGATIVE_INFINITY;
	private RangeSet<Integer> sinterval = null;
	private SAMSegment segment;
	// specify the status of the vertex
	// e.g., is visited or not during traversal
	private boolean status = false;
	
	public TraceableVertex(V id) {
		this.id = id;
	}
	
	public void changeStatus() {
		this.status  = !this.status;
	}
	
	public void setStatus(boolean status) {
		this.status  = status;
	}
	
	public void setSInterval() {
		this.sinterval = TreeRangeSet.create();
	}
	
	public void setSInterval(RangeSet<Integer> sinterval) {
		this.sinterval = sinterval;
	}
	
	public RangeSet<Integer> getSInterval() {
		return this.sinterval;
	}
	
	public void addSInterval(Range<Integer> range) {
		this.sinterval.add(range);
	}
	

	public void setSInterval(int lowerEndPoint, int upperEndPoint) {
		// TODO Auto-generated method stub
		this.setSInterval();
		this.addSInterval(lowerEndPoint, upperEndPoint);
	}
	
	public void addSInterval(int lowerEndPoint, int upperEndPoint) {
		this.sinterval.add(Range.closed(lowerEndPoint, upperEndPoint).
				canonical(DiscreteDomain.integers()));
	}
	
	public RangeSet<Integer> getSIntervalCopy() {
		RangeSet<Integer> rangeS = TreeRangeSet.create();
		rangeS.addAll(this.sinterval);
		return rangeS;
	}
	
	public boolean getStatus() {
		return this.status;
	}
	
	public void setBackTrace(TraceableVertex<V> backtrace) {
		this.backtrace = backtrace;
	}
	
	public void setScore(double score) {
		this.score = score;
	}
	
	public void setPenalty(double penalty) {
		this.penalty = penalty;
	}

	public void setSAMSegment(SAMSegment segment) {
		this.segment = segment;
	}
	
	public TraceableVertex<V> getBackTrace() {
		return this.backtrace;
	}
	
	public double getPenalty() {
		return this.penalty;
	}
	
	public double getScore() {
		return this.score;
	}
	
	public V getId() {
		return this.id;
	}
	
	public SAMSegment getSAMSegment() {
		return this.segment;
	}
	
	@Override
	public int hashCode() {
		return this.id.hashCode();
	}
	
	@Override
	public boolean equals(Object o) {
		if (o == this) return true;
		if (!(o instanceof TraceableVertex)) {
            return false;
        }
		@SuppressWarnings("unchecked")
		TraceableVertex<V> vertex = (TraceableVertex<V>) o;
		return this.id.equals(vertex.id);
	}
	
	@Override
	public String toString() {
		return this.id.toString();
	}

	@Override
	public int compareTo(TraceableVertex<V> v) {
		// TODO Auto-generated method stub
		return this.id.compareTo(v.id);
	}
}
