package cz1.ngs.model;

import java.util.Comparator;
import java.util.List;

import cz1.util.Constants;
import htsjdk.samtools.Cigar;
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
		if(!sam_rc.getReadNegativeStrandFlag())
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
		
		List<CigarElement> cigars = sam_rc.getCigar().getCigarElements();
		int nL = 0;
		for(CigarElement cigar : cigars) {
			switch(cigar.getOperator()) {
			case H:
			case S:
			case M:
			case I:
				nL += cigar.getLength();
				break;
			default:
				break;
			}
		}
		
		return new SAMSegment(sam_rc.getReadName(), 
				sam_rc.getReferenceName(),
				nL-sam_rc.getReadPositionAtReferencePosition(sam_rc.getAlignmentStart())+1,
				nL-sam_rc.getReadPositionAtReferencePosition(sam_rc.getAlignmentEnd())+1,
				sam_rc.getAlignmentStart(),
				sam_rc.getAlignmentEnd(),
				sam_rc.getCigarString(),
				sam_rc.getMappingQuality(),
				sam_rc.getIntegerAttribute("NM"),
				sam_rc.getIntegerAttribute("AS"));
	}
	
	public static SAMSegment samRecord(SAMRecord sam_rc, boolean rev, int seq_ln) {
		// TODO Auto-generated method stub
		
		if(sam_rc==null) return null;

		if(rev) {
			// reverse query sequence for reverse complementary alignment
			
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
				qstart  = seq_ln-qend+1;
				qend    = seq_ln-tmp +1;
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
		} else {
			// reverse reference sequence for reverse complementary alignment
			
			CigarElement c = sam_rc.getCigar().getFirstCigarElement();
			int hc = c.getOperator()==CigarOperator.HARD_CLIP ? c.getLength() : 0;

			String sseqid = sam_rc.getReferenceName();
			int sstart = sam_rc.getAlignmentStart();
			int send   = sam_rc.getAlignmentEnd();
			int qstart = sam_rc.getReadPositionAtReferencePosition(sstart)+hc;
			int qend   = sam_rc.getReadPositionAtReferencePosition(send)  +hc;
			String cigar = sam_rc.getCigarString();

			if( sam_rc.getReadNegativeStrandFlag() ) {
				sseqid  = sseqid+"'";
				int tmp = sstart;
				sstart  = seq_ln-send+1;
				send    = seq_ln-tmp +1;
				cigar   = Constants.cgRev(cigar);
			}

			return new SAMSegment(sam_rc.getReadName(),
					sseqid,
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
	
	public static class SubjectCoordinationComparator 
	implements Comparator<SAMSegment> {

		@Override
		public int compare(SAMSegment b1, SAMSegment b2) {
			// TODO Auto-generated method stub
			// check position on the subject sequence
			int b1_sstart = (b1.sstart<b1.send?b1.sstart:b1.send);
			int b2_sstart = (b2.sstart<b2.send?b2.sstart:b2.send);
			return b1_sstart-b2_sstart;
		}
	}
	
	public static class SegmentSizeComparator 
	implements Comparator<SAMSegment> {
		@Override
		public int compare(SAMSegment record1, SAMSegment record2) {
			// TODO Auto-generated method stub
			// large segment size ranks higher
			return record2.qlength()-record1.qlength();
		}	
	}
	
	public static SAMSegment collinear(final SAMSegment record1, final SAMSegment record2, final double max_shift) {
		// TODO Auto-generated method stub
		
		if(AlignmentSegment.reverse(record1, record2) ||
				AlignmentSegment.sdistance(record1, record2)>max_shift ||
				AlignmentSegment.qdistance(record1, record2)>max_shift ||
				AlignmentSegment.pdistance(record1, record2)>max_shift) {
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
		
		return record1.reverse() ? 
				new SAMSegment(record1.qseqid(),record1.sseqid(),qstart,qend,send,sstart):
				new SAMSegment(record1.qseqid(),record1.sseqid(),qstart,qend,sstart,send);
	}
}
