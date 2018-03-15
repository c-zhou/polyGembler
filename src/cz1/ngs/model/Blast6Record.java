package cz1.ngs.model;

import java.util.Comparator;

public class Blast6Record extends AlignmentRecord {

	// blastn outfmt 6
	private final double pident;   // percentage of identical matches
	private final int length;      // alignment length
	private final int mismatch;    // number of mismatches
	private final int gapopen;     // number of gap openings
	private final double evalue;   // expect value
	private final double bitscore;    // bit score
	
	public Blast6Record(
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
		super(qseqid, sseqid, qstart, qend, sstart, send);
		this.pident = pident;
		this.length = length;
		this.mismatch = mismatch;
		this.gapopen = gapopen;
		this.evalue = evalue;
		this.bitscore = bitscore;
	}

	public static Blast6Record blast6Record(String b6Record) {
		// TODO Auto-generated method stub
		if(b6Record==null) return null;
		String[] s = b6Record.split("\\s+");
		return new Blast6Record(
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

	public double pident() {
		return this.pident;
	}
	
	public int length() {
		return this.length;
	}
	
	public int mismatch() {
		return this.mismatch;
	}
	
	public int gapopen() {
		return this.gapopen;
	}
	
	public double evalue() {
		return this.evalue;
	}
	
	public double bitscore() {
		return this.bitscore;
	}
	
	public static class SubjectCoordinationComparator 
		implements Comparator<Blast6Record> {
		
		@Override
		public int compare(Blast6Record b1, Blast6Record b2) {
			// TODO Auto-generated method stub
			// check position on the subject sequence
			int b1_sstart = (b1.sstart<b1.send?b1.sstart:b1.send);
			int b2_sstart = (b2.sstart<b2.send?b2.sstart:b2.send);
			if(b1_sstart!=b2_sstart) return b1_sstart-b2_sstart;
			// larger match length ranks higher
			if(b1.length!=b2.length) return b2.length-b1.length;
			// larger identity ranks higher
			return Double.compare(b2.pident, b1.pident);
		}
	}
	
	public static class MatchIndentityComparator 
		implements Comparator<Blast6Record> {

		@Override
		public int compare(Blast6Record b1, Blast6Record b2) {
			// TODO Auto-generated method stub
			// larger match length ranks higher
			if(b1.length!=b2.length) return b2.length-b1.length;
			// larger identity ranks higher
			return Double.compare(b2.pident, b1.pident);
		}

	}
	
	public static class SegmentSizeComparator 
		implements Comparator<Blast6Record> {
		@Override
		public int compare(Blast6Record record1, Blast6Record record2) {
			// TODO Auto-generated method stub
			// large segment size ranks higher
			return record2.length-record1.length;
		}
		
	}
	
	@Override
	public String toString() {
		return this.qseqid+"\t"+
				this.sseqid+"\t"+
				this.pident+"\t"+
				this.length+"\t"+
				this.mismatch+"\t"+
				this.gapopen+"\t"+
				this.qstart+"\t"+
				this.qend+"\t"+
				this.sstart+"\t"+
				this.send+"\t"+
				this.evalue+"\t"+
				this.bitscore+"\t";
	}
	
	public static Blast6Record collinear(final Blast6Record record1, final Blast6Record record2, final double max_shift) {
		// TODO Auto-generated method stub
		
		if(AlignmentRecord.reverse(record1, record2) ||
				AlignmentRecord.sdistance(record1, record2)>max_shift ||
				AlignmentRecord.qdistance(record1, record2)>max_shift ||
				AlignmentRecord.pdistance(record1, record2)>max_shift) {
			return null;
		}

		// merge collinear alignment segments
		int qstart = Math.min(record1.true_qstart(), 
				record2.true_qstart());
		int qend = Math.max(record1.true_qend(),
				record2.true_qend());
		int sstart = Math.min(record1.true_sstart(), 
				record2.true_sstart());
		int send = Math.max(record1.true_send(),
				record2.true_send());
		double pident = Math.max(record1.pident(), record2.pident());
		int length = qend-qstart;
		
		return record1.reverse() ? 
				new Blast6Record(record1.qseqid(),record1.sseqid(),pident,length,-1,-1,qstart,qend,send,sstart,-1,-1):
				new Blast6Record(record1.qseqid(),record1.sseqid(),pident,length,-1,-1,qstart,qend,sstart,send,-1,-1);
	}
}
