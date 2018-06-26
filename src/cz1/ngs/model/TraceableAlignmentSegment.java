package cz1.ngs.model;

public class TraceableAlignmentSegment extends AlignmentSegment {

	public final static int match_score = 2;
	public final static int mismatch_penalty = -4;
	public final static int gap_open = -20;
	public final static int gap_extension = -1;
	
	private TraceableAlignmentSegment next = null;
	private TraceableAlignmentSegment previous = null;
	private double objective = Double.NEGATIVE_INFINITY;
	private int mer_count = 0;
	private boolean to_query_start = false;
	private boolean to_query_end   = false;
	private boolean to_subject_start = false;
	private boolean to_subject_end   = false;
	private boolean to_start = false;
	private boolean to_end   = false;
	private boolean end_to_end = false; 
	private int sstart_clip = 0;
	private int send_clip = 0;
	private int qstart_clip = 0;
	private int qend_clip = 0;
	private int score = 0;

	public TraceableAlignmentSegment(String qseqid, String sseqid, 
			int qstart, int qend, int sstart, int send) {
		// TODO Auto-generated constructor stub
		super(qseqid, sseqid, qstart, qend, sstart, send, true);
	}

	public void setClip(int qstart_clip, int qend_clip, int sstart_clip, int send_clip) {
		// TODO Auto-generated method stub
		this.qstart_clip = qstart_clip;
		this.qend_clip = qend_clip;
		this.sstart_clip = sstart_clip;
		this.send_clip = send_clip;
	}

	public void setQueryStartClip(int qstart_clip) {
		// TODO Auto-generated method stub
		this.qstart_clip = qstart_clip;
	}

	public void setQueryEndClip(int qend_clip) {
		// TODO Auto-generated method stub
		this.qend_clip = qend_clip;
	}

	public int getQueryStartClip() {
		// TODO Auto-generated method stub
		return this.qstart_clip;
	}

	public int getQueryEndClip() {
		// TODO Auto-generated method stub
		return this.qend_clip;
	}

	public void setSubjectStartClip(int sstart_clip) {
		// TODO Auto-generated method stub
		this.sstart_clip = sstart_clip;
	}

	public void setSubjectEndClip(int send_clip) {
		// TODO Auto-generated method stub
		this.send_clip = send_clip;
	}

	public int getSubjectStartClip() {
		// TODO Auto-generated method stub
		return this.sstart_clip;
	}

	public int getSubjectEndClip() {
		// TODO Auto-generated method stub
		return this.send_clip;
	}

	public int getClip() {
		// TODO Auto-generated method stub
		return this.sstart_clip+this.send_clip;
	}

	public void setEndToEnd(boolean to_query_start, boolean to_query_end,
			boolean to_subject_start, boolean to_subject_end) {
		// TODO Auto-generated method stub
		this.to_query_start = to_query_start;
		this.to_query_end = to_query_end;
		this.to_subject_start = to_subject_start;
		this.to_subject_end = to_subject_end;
		this.to_start = to_query_start||to_subject_start;
		this.to_end   = to_query_end||to_subject_end;
		this.end_to_end = to_start&&to_end;
	}

	public void setToQueryStart(boolean to_query_start) {
		// TODO Auto-generated method stub
		this.to_query_start = to_query_start;
	}

	public void setToQueryEnd(boolean to_query_end) {
		// TODO Auto-generated method stub
		this.to_query_end = to_query_end;
	}

	public void setToSubjectStart(boolean to_subject_start) {
		// TODO Auto-generated method stub
		this.to_subject_start = to_subject_start;
	}

	public void setToSubjectEnd(boolean to_subject_end) {
		// TODO Auto-generated method stub
		this.to_subject_end = to_subject_end;
	}

	public boolean getToQueryStart() {
		// TODO Auto-generated method stub
		return this.to_query_start;
	}

	public boolean getToQueryEnd() {
		// TODO Auto-generated method stub
		return this.to_query_end;
	}

	public boolean getToSubjectStart() {
		// TODO Auto-generated method stub
		return this.to_subject_start;
	}

	public boolean getToSubjectEnd() {
		// TODO Auto-generated method stub
		return this.to_subject_end;
	}

	public void setEndToEnd(boolean end_to_end) {
		// TODO Auto-generated method stub
		this.end_to_end = end_to_end;
	}

	public void setToStart(boolean to_start) {
		// TODO Auto-generated method stub
		this.to_start = to_start;
	}

	public void setToEnd(boolean to_end) {
		// TODO Auto-generated method stub
		this.to_end = to_end;
	}

	public boolean getToStart() {
		// TODO Auto-generated method stub
		return this.to_start;
	}


	public boolean getToEnd() {
		// TODO Auto-generated method stub
		return this.to_end;
	}

	public boolean getEndToEnd() {
		// TODO Auto-generated method stub
		return this.end_to_end;
	}

	public TraceableAlignmentSegment(String qseqid, String sseqid, 
			int qstart, int qend, int sstart, int send, int mer_count) {
		// TODO Auto-generated constructor stub
		super(qseqid, sseqid, qstart, qend, sstart, send, true);
		this.mer_count = mer_count;
	}

	public void setTraceForward(TraceableAlignmentSegment next) {
		this.next = next;
	}

	public void setTraceBackward(TraceableAlignmentSegment previous) {
		this.previous = previous;
	}

	public void calcScore() {
		this.score = (this.qend-this.qstart+1)*match_score+
				Math.min(qstart_clip, sstart_clip)*gap_extension+
				Math.min(qend_clip, send_clip)*gap_extension;
	}

	public void setScore(int score) {
		this.score = score;
	}

	public int getScore() {
		return this.score;
	}

	public void setObjective(double objective) {
		this.objective = objective;
	}

	public void addObjective(double dobj) {
		this.objective += dobj;
	}

	public void setMerCount(int mer_count) {
		this.mer_count = mer_count;
	}

	public TraceableAlignmentSegment getTraceForward() {
		return this.next;
	}

	public TraceableAlignmentSegment getTraceBackward() {
		return this.previous;
	}

	public double getObjective() {
		return this.objective;
	}

	public int getMerCount() {
		return this.mer_count;
	}

	public static TraceableAlignmentSegment collinear(final TraceableAlignmentSegment record1, 
			final TraceableAlignmentSegment record2, final double max_shift) {
		// TODO Auto-generated method stub

		// if(TraceableAlignmentSegment.sdistance(record1, record2)>max_shift ||
		// 		TraceableAlignmentSegment.qdistance(record1, record2)>max_shift ||
		//		TraceableAlignmentSegment.pdistance(record1, record2)>max_shift) {
		//	return null;
		// }

		// how about we only check the pdistance?
		if(TraceableAlignmentSegment.pdistance(record1, record2)>max_shift)
			return null;

		// merge collinear alignment segments
		int qstart = Math.min(record1.qstart(), record2.qstart());
		int qend = Math.max(record1.qend(), record2.qend());
		int sstart = Math.min(record1.sstart(), record2.sstart());
		int send = Math.max(record1.send(), record2.send());

		return new TraceableAlignmentSegment(record1.qseqid, record1.sseqid,qstart,qend,sstart,send);
	}
}
