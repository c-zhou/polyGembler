package cz1.ngs.model;

import java.io.Serializable;

public class SequenceVertex<V> implements Serializable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 4754142281268346991L;
	
	protected final V id;               // vertex id
	protected String seq = null;        // sequence
	protected int len = -1;             // sequence length
	protected boolean rev = false;      // is reverse?
	
	public SequenceVertex (V id) {
		this.id = id;
	}
	
	public SequenceVertex (V id, String seq) {
		this.id = id;
		this.seq = seq;
		this.len = seq.length();
	}
	
	public SequenceVertex (V id, String seq, boolean rev) {
		this.id = id;
		this.seq = seq;
		this.len = seq.length();
		this.rev = rev;
	}
	
	public SequenceVertex (V id, int len) {
		this.id = id;
		this.len = len;
	}
	
	public SequenceVertex (V id, int len, boolean rev) {
		this.id = id;
		this.len = len;
		this.rev = rev;
	}
	
	public SequenceVertex (V id, boolean rev) {
		this.id = id;
		this.rev = rev;
	}
	
	public V id() {
		return this.id;
	}
	
	public String seq() {
		return this.seq;
	}
	
	public int len() {
		return this.len;
	}
	
	public boolean rev() {
		return this.rev;
	}
	
	public void setSeq(String seq) {
		this.seq = seq;
		this.len = seq.length();
	}
	
	public void setRev(boolean rev) {
		this.rev = rev;
	}
	
	@Override
	public int hashCode() {
		return this.id.hashCode();
	}
}
