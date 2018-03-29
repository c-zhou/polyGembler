package cz1.ngs.model;

public class CompoundBlast6Segment extends Blast6Segment {
	
	private int score;
	private int nSeg;
	private int[][] segs;
	
	public CompoundBlast6Segment(
			final String qseqid,   // query (e.g., gene) sequence id
			final String sseqid,   // subject (e.g., reference genome) sequence id
			final double pident,   // percentage of identical matches
			final int length,      // alignment length
			final int mismatch,    // number of mismatches
			final int gapopen,     // number of gap openings
			final int qstart,      // start of alignment in query
			final int qend,        // end of alignment in query
			final int sstart,      // start of alignment in subject
			final int send,        // end of alignment in subject
			final double evalue,   // expect value
			final double bitscore  // bit score
			) {
		// TODO Auto-generated method stub
		super(qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore);
	}
	
	public CompoundBlast6Segment(
			final String qseqid,   // query (e.g., gene) sequence id
			final String sseqid,   // subject (e.g., reference genome) sequence id
			final double pident,   // percentage of identical matches
			final int length,      // alignment length
			final int mismatch,    // number of mismatches
			final int gapopen,     // number of gap openings
			final int qstart,      // start of alignment in query
			final int qend,        // end of alignment in query
			final int sstart,      // start of alignment in subject
			final int send,        // end of alignment in subject
			final double evalue,   // expect value
			final double bitscore, // bit score
			final int[][] segs,    // segments [qstart, qend, sstart, send]
			final int score        // segments weighted score  
			) {
		// TODO Auto-generated method stub
		super(qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore);
		this.score = score;
		this.nSeg = segs.length;
		this.segs = segs;
	}
	
	public static CompoundBlast6Segment compoundBlast6Segment(String b6Record) {
		// TODO Auto-generated method stub
		if(b6Record==null) return null;
		String[] s = b6Record.split("\\s+");
		return new CompoundBlast6Segment(
				s[0],
				s[1],
				Double.parseDouble(s[2]),
				Integer.parseInt(s[3]),
				Integer.parseInt(s[4]),
				Integer.parseInt(s[5]),
				Integer.parseInt(s[6]),
				Integer.parseInt(s[7]),
				Integer.parseInt(s[8]),
				Integer.parseInt(s[9]),
				Double.parseDouble(s[10]),
				Double.parseDouble(s[11]) );
	}
	
	public void setSegs(final int[][] segs) {
		this.segs = segs;
		this.nSeg = segs.length;
	}
	
	public void setScore(final int score) {
		this.score = score;
	}
	
	public int getNSeg() {
		return this.nSeg;
	}
	
	public int[][] getSegs() {
		return this.segs;
	}
	
	public int getScore() {
		return this.score;
	}
}
