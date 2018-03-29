package cz1.ngs.model;

import cz1.util.Constants;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class SAMSegment extends AlignmentSegment {
	
	private final String cigar;
	private final int qual;
	private final int nm; // edit distance to the subject sequence
	private final int as; // alignment score
	
	public SAMSegment(String qseqid,
			String sseqid, 
			int qstart, 
			int qend, 
			int sstart, 
			int send) {
		// TODO Auto-generated constructor stub
		super(qseqid, sseqid, qstart, qend, sstart, send);
		cigar = null;
		qual = -1;
		nm = -1;
		as = -1;
	}
	
	public SAMSegment(String qseqid,
			String sseqid, 
			int qstart, 
			int qend, 
			int sstart, 
			int send,
			String cigar,
			int qual,
			int nm,
			int as) {
		// TODO Auto-generated constructor stub
		super(qseqid, sseqid, qstart, qend, sstart, send);
		this.cigar = qstart<=qend ? cigar : Constants.cgRev(cigar);
		this.qual = qual;
		this.nm = nm;
		this.as = as;
	}
	
	public String cigar() {
		return this.cigar;
	}
	
	public int qual() {
		return this.qual;
	}
	
	public int nm() {
		return this.nm;
	}
	
	public int as() {
		return this.as;
	}

	public static SAMSegment samRecord(SAMRecord sam_rc) {
		// TODO Auto-generated method stub
		if(sam_rc==null) return null;
		return new SAMSegment(sam_rc.getReadName(), 
				sam_rc.getReferenceName(),
				sam_rc.getReadPositionAtReferencePosition(sam_rc.getAlignmentStart()),
				sam_rc.getReadPositionAtReferencePosition(sam_rc.getAlignmentEnd()),
				sam_rc.getAlignmentStart(),
				sam_rc.getAlignmentEnd(),
				sam_rc.getCigarString(),
				sam_rc.getMappingQuality(),
				sam_rc.getIntegerAttribute("NM"),
				sam_rc.getIntegerAttribute("AS"));
	}
	
	public static SAMSegment samRecord(SAMRecord sam_rc, boolean rev, int qry_ln) {
		// TODO Auto-generated method stub
		
		if(sam_rc==null) return null;
		
		CigarElement c = sam_rc.getCigar().getFirstCigarElement();
		int hc = c.getOperator()==CigarOperator.HARD_CLIP ? c.getLength() : 0;
		
		String qseqid = sam_rc.getReadName();
		int qstart = sam_rc.getReadPositionAtReferencePosition(sam_rc.getAlignmentStart())+hc;
		int qend   = sam_rc.getReadPositionAtReferencePosition(sam_rc.getAlignmentEnd())  +hc;
		int sstart = sam_rc.getAlignmentStart();
		int send   = sam_rc.getAlignmentEnd();
		String cigar = sam_rc.getCigarString();
				
		if( sam_rc.getReadNegativeStrandFlag() ) {
			qseqid  = qseqid+"'";
			int tmp = qstart;
			qstart  = qry_ln-qend+1;
			qend    = qry_ln-tmp +1;
			cigar   = Constants.cgRev(cigar);
		}
		
		return new SAMSegment(qseqid, 
				sam_rc.getReferenceName(),
				qstart,
				qend,
				sstart,
				send,
				cigar,
				sam_rc.getMappingQuality(),
				sam_rc.getIntegerAttribute("NM"),
				sam_rc.getIntegerAttribute("AS"));
	}
}
