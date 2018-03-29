package cz1.ngs.model;

public class TraceableVertex<V> {
	private final V id;
	private TraceableVertex<V> backtrace = null;
	private double score = Double.NEGATIVE_INFINITY;
	private SAMSegment segment;
	
	public TraceableVertex(V id) {
		this.id = id;
	}
	
	public TraceableVertex(V id, TraceableVertex<V> backtrace, double score, SAMSegment segment) {
		this.id = id;
		this.backtrace = backtrace;
		this.score = score;
		this.segment = segment;
	}
	
	public void setBackTrace(TraceableVertex<V> backtrace) {
		this.backtrace = backtrace;
	}
	
	public void setScore(double score) {
		this.score = score;
	}

	public void setSAMSegment(SAMSegment segment) {
		this.segment = segment;
	}
	
	public TraceableVertex<V> getBackTrace() {
		return this.backtrace;
	}
	
	public double getScore() {
		return this.score;
	}
	
	public V getVertexId() {
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
}
