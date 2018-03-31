package cz1.ngs.model;

import java.util.Comparator;

import org.apache.commons.math3.stat.StatUtils;

public class AlignmentSegment {
	
	protected final String qseqid;   // query (e.g., gene) sequence id
	protected final String sseqid;   // subject (e.g., reference genome) sequence id
	protected final int qstart;      // start of alignment in query
	protected final int qend;        // end of alignment in query
	protected final int sstart;      // start of alignment in subject
	protected final int send;        // end of alignment in subject
	
	protected final double sintercept;
	protected final double qintercept;
	
	public AlignmentSegment(
			final String qseqid,   // query (e.g., gene) sequence id
			final String sseqid,   // subject (e.g., reference genome) sequence id
			final int qstart,      // start of alignment in query
			final int qend,        // end of alignment in query
			final int sstart,      // start of alignment in subject
			final int send         // end of alignment in subject
			) {
		this.qseqid = qseqid;
		this.sseqid = sseqid;
		if(qstart<=qend) {
			this.qstart = qstart;
			this.qend = qend;
			this.sstart = sstart;
			this.send = send;
		} else {
			this.qstart = qend;
			this.qend = qstart;
			this.sstart = send;
			this.send = sstart;
		}
		this.sintercept = this.calc_sintercept();
		this.qintercept = this.calc_qintercept();
	}

	@Override
    public boolean equals(Object obj) {

        if (obj == this) return true;
        if (!(obj instanceof AlignmentSegment)) {
            return false;
        }

        AlignmentSegment b = (AlignmentSegment) obj;

        return b.qseqid.equals(qseqid) &&
        		b.sseqid.equals(sseqid) &&
        		b.qstart==qstart &&
        		b.qend==qend &&
        		b.sstart==sstart;
    }

    @Override
    public int hashCode() {
        int hc = 17;
        hc = 31 * hc + qseqid.hashCode();
        hc = 31 * hc + sseqid.hashCode();
        hc = 31 * hc + qstart;
        hc = 31 * hc + qend;
        hc = 31 * hc + sstart;
        hc = 31 * hc + send;
        return hc;
    }
    
    public String qseqid() {
    	return this.qseqid;
    }
    
    public String sseqid() {
    	return this.sseqid;
    }
    
    public int qstart() {
    	return this.qstart;
    }
    
    public int qend() {
    	return this.qend;
    }
    
    public int sstart() {
    	return this.sstart;
    }
    
    public int send() {
    	return this.send;
    }
    
    public static class SubjectCoordinationComparator 
    	implements Comparator<AlignmentSegment> {
		@Override
		public int compare(AlignmentSegment b1, AlignmentSegment b2) {
			// TODO Auto-generated method stub
			return b1.true_sstart()-b2.true_sstart();
		}
	}
    
    public static class SInterceptComparator 
    	implements Comparator<AlignmentSegment> {
		@Override
		public int compare(AlignmentSegment record1, AlignmentSegment record2) {
			// TODO Auto-generated method stub
			return Double.compare(record1.sintercept, record2.sintercept);
		}
	}
    
    public static class QInterceptComparator 
    	implements Comparator<AlignmentSegment> {
    	@Override
    	public int compare(AlignmentSegment record1, AlignmentSegment record2) {
    		// TODO Auto-generated method stub
    		return Double.compare(record1.qintercept, record2.qintercept);
    	}
    }
    
    private double calc_qintercept() {
		// TODO Auto-generated method stub
		return ((double)sstart*qend-(double)send*qstart)/(sstart-send);
	}

	public int true_sstart() {
		// TODO Auto-generated method stub
		return Math.min(this.sstart, this.send);
	}

	public int true_qstart() {
		// TODO Auto-generated method stub
		return Math.min(this.qstart, this.qend);
	}

	public int true_send() {
		// TODO Auto-generated method stub
		return Math.max(this.sstart, this.send);
	}

	public int true_qend() {
		// TODO Auto-generated method stub
		return Math.max(this.qstart, this.qend);
	}
	
	private double calc_sintercept() {
		// TODO Auto-generated method stub
		return ((double)send*qstart-(double)sstart*qend)/(qstart-qend);
	}

	public double sL() {
		return Math.abs(
				this.sstart-this.send);
	}

	public double qL() {
		return Math.abs(
				this.qstart-this.qend);
	}

	public double slope() {
		return (double)(this.qstart-this.qend)/
				(this.sstart-this.send);
	}

	public double abs_slop() {
		return Math.abs(this.slope());
	}

	public double pdistance(AlignmentSegment record) {
		// TODO Auto-generated method stub
		// perpendicular distance
		if(this.intersect(record)) return 0;
		double[] allD = new double[4];
		allD[0] = this.distance(new double[]{record.sstart, record.qstart});
		allD[1] = this.distance(new double[]{record.send, record.qend});
		allD[2] = record.distance(new double[]{this.sstart, this.qstart});
		allD[3] = record.distance(new double[]{this.send, this.qend});
		return StatUtils.min(allD);
	}
	
	public static double pdistance(AlignmentSegment record1, 
			AlignmentSegment record2) {
		return record1.pdistance(record2);
	}

	private boolean intersect(AlignmentSegment record) {
		// TODO Auto-generated method stub
		double s1_x = this.send-this.sstart,
				s1_y = this.qend-this.qstart,
				s2_x = record.send-record.sstart,
				s2_y = record.qend-record.qstart;
		
		if(s1_x*s2_y==s2_x*s1_y) return false;
		double s = (-s1_y * (this.sstart - record.sstart) + s1_x * (this.qstart - record.qstart)) / (-s2_x * s1_y + s1_x * s2_y);
	    double t = ( s2_x * (this.qstart - record.qstart) - s2_y * (this.sstart - record.sstart)) / (-s2_x * s1_y + s1_x * s2_y);
		
	    return s>=0&&s<=1&&t>=0&&t<=1;
	}
	
	public static boolean intersect(AlignmentSegment record1, 
			AlignmentSegment record2) {
		// TODO Auto-generated method stub
		return record1.intersect(record2);
	}

	private double distance(double[] point) {
		if(this.isdot())
			return distance(new double[]{this.sstart, this.qstart}, 
					point);

		double px = point[0], py = point[1],
				X1 = this.sstart, Y1 = this.qstart,
				X2 = this.send, Y2 = this.qend,
				dx = X2-X1, dy = Y2-Y1;
		double t = ((px-X1)*dx+(py-Y1)*dy)/(dx*dx+dy*dy);
		if(t<0) {
			dx = px-X1;
			dy = py-Y1;
		} else if(t>1) {
			dx = px-X2;
			dy = py-Y2;
		} else {
			dx = px-X1-t*dx;
			dy = py-Y1-t*dy;
		}
		return Math.sqrt(dx*dx+dy*dy);
	}

	public static double distance(AlignmentSegment record, 
			double[] point) {
		return record.distance(point);
	}
		
	public static double distance(double[] p1, double[] p2) {
		return Math.sqrt(
				Math.pow(p1[0]-p2[0],2)+
				Math.pow(p1[1]-p2[1],2));
	}

	private boolean isdot() {
		return this.sstart==this.send &&
				this.qstart==this.qend;
	}
	
	public int sdirect() {
		return this.sstart<this.send ? 1 : -1;
	}
	
	public int qdirect() {
		return this.qstart<this.qend ? 1 : -1;
	}
		
	public double sintercept() {
		return this.sintercept;
	}

	public double qintercept() {
		return this.qintercept;
	}
	
	public static boolean forward(AlignmentSegment record1, 
			AlignmentSegment record2) {
		return record1.forward()&&record2.forward() ||
				record1.reverse()&&record2.reverse();
	}
	
	public static boolean reverse(AlignmentSegment record1, 
			AlignmentSegment record2) {
		return !forward(record1, record2);
	}
	
	public boolean forward(AlignmentSegment record) {
		return this.forward()&&record.forward() ||
				this.reverse()&&record.reverse();
	}
	
	public boolean reverse(AlignmentSegment record) {
		return !this.forward(record);
	}
	
	public boolean forward() {
		return this.sstart<=this.send&&this.qstart<=this.qend ||
				this.sstart>=this.send&&this.qstart>=this.qend;
	}
	
	public boolean reverse() {
		return !this.forward();
	}
	
	@Override
	public String toString() {
		return this.qseqid+"\t"+
				this.sseqid+"\t"+
				this.qstart+"\t"+
				this.qend+"\t"+
				this.sstart+"\t"+
				this.send+"\t";
	}
	
	public void print() {
		// TODO Auto-generated method stub
		System.out.println(this.toString());
	}

	public int sdistance(AlignmentSegment record) {
		// TODO Auto-generated method stub
		if(this.soverlap(record)) return 0;
		if(!this.sseqid.equals(record.sseqid)) return Integer.MAX_VALUE;
		return Math.max(this.true_sstart(), record.true_sstart())-
				Math.min(this.true_send(), record.true_send());
	}
	
	public static int sdistance(AlignmentSegment record1,
			AlignmentSegment record2) {
		// TODO Auto-generated method stub
		if(soverlap(record1, record2)) return 0;
		if(!record1.sseqid.equals(record2.sseqid)) return Integer.MAX_VALUE;
		return Math.max(record1.true_sstart(), record2.true_sstart())-
				Math.min(record1.true_send(), record2.true_send());
	}
	
	public int qdistance(AlignmentSegment record) {
		// TODO Auto-generated method stub
		if(this.qoverlap(record)) return 0;
		if(!this.qseqid.equals(record.qseqid)) return Integer.MAX_VALUE;
		return Math.max(this.true_qstart(), record.true_qstart())-
				Math.min(this.true_qend(), record.true_qend());
	}
	
	public static int qdistance(AlignmentSegment record1,
			AlignmentSegment record2) {
		// TODO Auto-generated method stub
		if(qoverlap(record1, record2)) return 0;
		if(!record1.qseqid.equals(record2.qseqid)) return Integer.MAX_VALUE;
		return Math.max(record1.true_qstart(), record2.true_qstart())-
				Math.min(record1.true_qend(), record2.true_qend());
	}
	
	public boolean soverlap(AlignmentSegment record) {
		// TODO Auto-generated method stub
		if(!this.sseqid.equals(record.sseqid)) return false;
		return ((double)(this.true_sstart()-record.true_send()))*
				(record.true_sstart()-this.true_send())>=0;
	}

	public boolean qoverlap(AlignmentSegment record) {
		// TODO Auto-generated method stub
		if(!this.qseqid.equals(record.qseqid)) return false;
		return ((double)(this.true_qstart()-record.true_qend()))*
				(record.true_qstart()-this.true_qend())>=0;
	}

	public static boolean soverlap(AlignmentSegment record1,
			AlignmentSegment record2) {
		// TODO Auto-generated method stub
		if(!record1.sseqid.equals(record2.sseqid)) return false;
		return ((double)(record1.true_sstart()-record2.true_send()))*
				(record2.true_sstart()-record1.true_send())>=0;
	}

	public static boolean qoverlap(AlignmentSegment record1,
			AlignmentSegment record2) {
		// TODO Auto-generated method stub
		if(!record1.qseqid.equals(record2.qseqid)) return false;
		return ((double)(record1.true_qstart()-record2.true_qend()))*
				(record2.true_qstart()-record1.true_qend())>=0;
	}
	
	public static AlignmentSegment collinear(final AlignmentSegment record1, final AlignmentSegment record2, final double max_shift) {
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
				new AlignmentSegment(record1.qseqid(),record1.sseqid(),qstart,qend,send,sstart):
				new AlignmentSegment(record1.qseqid(),record1.sseqid(),qstart,qend,sstart,send);
	}

	public int qlength() {
		// TODO Auto-generated method stub
		return Math.abs(this.qend-this.qstart+1);
	}
	
	public int slength() {
		// TODO Auto-generated method stub
		return Math.abs(this.send-this.sstart+1);
	}
}
