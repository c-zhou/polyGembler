package cz1.ngs.model;

public class FASTQRead {
	private final String id;
	private final String str;
	private final String qual;
	
	private final static String qstr = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
	
	public FASTQRead(final String id, final String str, final String qual) {
		this.id   = id;
		this.str  = str;
		this.qual = qual;
	}
	
	public static double getBaseError(char qual) {
		//TODO: distinguish between different coding system
		// use default Phred+33
		int i = qstr.indexOf(qual);
		if(i<0) 
			throw new RuntimeException("invalid base quality score!!!");
		// here we make a maximum of 0.5 since e=1.0 is controversial
		return Math.min(0.5, Math.pow(10, -i/10.0));
	}
}
