package cz1.ngs.model;

import java.util.Comparator;

import org.apache.commons.math3.stat.StatUtils;

public class BlastRecord {
	
	protected final String qseqid;   // query (e.g., gene) sequence id
	protected final String sseqid;   // subject (e.g., reference genome) sequence id
	protected final int qstart;      // start of alignment in query
	protected final int qend;        // end of alignment in query
	protected final int sstart;      // start of alignment in subject
	protected final int send;        // end of alignment in subject
	
	protected final double sintercept;
	protected final double qintercept;
	
	public BlastRecord(
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
        if (!(obj instanceof BlastRecord)) {
            return false;
        }

        BlastRecord b = (BlastRecord) obj;

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
    	implements Comparator<BlastRecord> {
		@Override
		public int compare(BlastRecord b1, BlastRecord b2) {
			// TODO Auto-generated method stub
			return b1.true_sstart()-b2.true_sstart();
		}
	}
    
    public static class SInterceptComparator 
    	implements Comparator<BlastRecord> {
		@Override
		public int compare(BlastRecord record1, BlastRecord record2) {
			// TODO Auto-generated method stub
			return Double.compare(record1.sintercept, record2.sintercept);
		}
	}
    
    public static class QInterceptComparator 
    	implements Comparator<BlastRecord> {
    	@Override
    	public int compare(BlastRecord record1, BlastRecord record2) {
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

	public double distance(BlastRecord record) {
		if(this.intersect(record)) return 0;
		double[] allD = new double[4];
		allD[0] = this.distance(new double[]{record.sstart, record.qstart});
		allD[1] = this.distance(new double[]{record.send, record.qend});
		allD[2] = record.distance(new double[]{this.sstart, this.qstart});
		allD[3] = record.distance(new double[]{this.send, this.qend});
		return StatUtils.min(allD);
	}
	
	public static double distance(BlastRecord record1, 
			BlastRecord record2) {
		return record1.distance(record2);
	}

	private boolean intersect(BlastRecord record) {
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
	
	public static boolean intersect(BlastRecord record1, 
			BlastRecord record2) {
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

	public static double distance(BlastRecord record, 
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
	
	public static boolean forward(BlastRecord record1, 
			BlastRecord record2) {
		return record1.forward()&&record2.forward() ||
				record1.reverse()&&record2.reverse();
	}
	
	public static boolean reverse(BlastRecord record1, 
			BlastRecord record2) {
		return !forward(record1, record2);
	}
	
	public boolean forward(BlastRecord record) {
		return this.forward()&&record.forward() ||
				this.reverse()&&record.reverse();
	}
	
	public boolean reverse(BlastRecord record) {
		return !this.forward(record);
	}
	
	public boolean forward() {
		return this.sstart<=this.send&&this.qstart<=this.qend ||
				this.sstart>=this.send&&this.qstart>=this.qend;
	}
	
	public boolean reverse() {
		return !this.forward();
	}
}
