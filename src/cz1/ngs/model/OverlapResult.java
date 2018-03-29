package cz1.ngs.model;

public final class OverlapResult {
	private final String fromId;
	private final String toId;
	private final int aLen;
	private final int bLen;
	
	private int a1;
	private int a2;
	private int b1;
	private int b2;
	private double score;
	private double rawScore;
	
	private boolean realigned = false;
	
	protected OverlapResult(String fromId, String toId, double score, double rawScore, int a1, int a2, int aLen, int b1, int b2, int bLen, boolean realigned) {
		this.fromId = fromId;
		this.toId = toId;
		this.score = score;
		this.rawScore = rawScore;
		this.a1 = a1;
		this.a2 = a2;
		this.aLen = aLen;
		this.b1 = b1;
		this.b2 = b2;
		this.bLen = bLen;
		this.realigned = realigned;
	}
	
	protected OverlapResult(String fromId, String toId, double score, double rawScore, int a1, int a2, int aLen, int b1, int b2, int bLen) {
		this.fromId = fromId;
		this.toId = toId;
		this.score = score;
		this.rawScore = rawScore;
		this.a1 = a1;
		this.a2 = a2;
		this.aLen = aLen;
		this.b1 = b1;
		this.b2 = b2;
		this.bLen = bLen;
	}

	/**
	 * @return the fromId
	 */
	public String getFromId()
	{
		return this.fromId;
	}

	/**
	 * @return the toId
	 */
	public String getToId()
	{
		return this.toId;
	}
	
	/**
	 * @return from start
	 */
	public int getFromStart()
	{
		return this.a1;
	}
	
	/**
	 * @return from end
	 */
	public int getFromEnd()
	{
		return this.a2;
	}

	/**
	 * @return from length
	 */
	public int getFromLen()
	{
		return this.aLen;
	}
	
	/**
	 * @return to start
	 */
	public int getToStart()
	{
		return this.b1;
	}
	
	/**
	 * @return to end
	 */
	public int getToEnd()
	{
		return this.b2;
	}

	/**
	 * @return to length
	 */
	public int getToLen()
	{
		return this.bLen;
	}

	/**
	 * @return the score
	 */
	public double getScore()
	{
		return this.score;
	}
	
	/**
	 * @return the rawScore
	 */
	public double getRawScore()
	{
		return this.rawScore;
	}
	
	public void setFromStart(int a1) {
		this.a1 = a1;
	}
	
	public void setFromEnd(int a2) {
		this.a2 = a2;
	}
	
	public void setToStart(int b1) {
		this.b1 = b1;
	}
	
	public void setToEnd(int b2) {
		this.b2 = b2;
	}
	
	public void setScore(double score) {
		this.score = score;
	}

	public void setRawScore(double rawScore) {
		this.rawScore = rawScore;
	}
	
	public void setRealigned() {
		this.realigned = true;
	}
	
	public void setRealigned(boolean realigned) {
		this.realigned = realigned;
	}
	
	public boolean isRealigned() {
		return this.realigned;
	}
	
	public String toString()
	{
		return String.format("%s %s %.6f %.6f %d %d %d %d %d %d",
				this.fromId, 
				this.toId,
				this.score,
				this.rawScore,
				this.a1,
				this.a2,
				this.aLen,
				this.b1,
				this.b2,
				this.bLen);
	}


}
