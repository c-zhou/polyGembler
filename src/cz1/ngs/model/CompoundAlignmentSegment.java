package cz1.ngs.model;

import java.util.List;

import com.google.common.collect.RangeSet;

import htsjdk.samtools.SAMRecord;

public class CompoundAlignmentSegment extends AlignmentSegment {

	private String key = null;
	private final RangeSet<Integer> subCov;
	private final RangeSet<Integer> qryCov;
	private final int slen;
	private final int qlen;
	private double score = 0;
	private List<SAMRecord> samRecords = null;
	
	public CompoundAlignmentSegment(
			final String key,      // key
			final String qseqid,   // query (e.g., gene) sequence id
			final String sseqid,   // subject (e.g., reference genome) sequence id
			final int qstart,      // start of alignment in query
			final int qend,        // end of alignment in query
			final int sstart,      // start of alignment in subject
			final int send,        // end of alignment in subject
			final RangeSet<Integer> subCov,
			final RangeSet<Integer> qryCov,
			final int slen,
			final int qlen,
			final List<SAMRecord> samRecords) {
	
		super(qseqid, sseqid, qstart,qend,sstart,send, true);
		this.key      = key;
		this.subCov   = subCov;
		this.qryCov   = qryCov;
		this.slen     = slen;
		this.qlen     = qlen;
		this.samRecords = samRecords;
	}
	
	public CompoundAlignmentSegment(
			final String key,      // key
			final String qseqid,   // query (e.g., gene) sequence id
			final String sseqid,   // subject (e.g., reference genome) sequence id
			final int qstart,      // start of alignment in query
			final int qend,        // end of alignment in query
			final int sstart,      // start of alignment in subject
			final int send,        // end of alignment in subject
			final RangeSet<Integer> subCov,
			final RangeSet<Integer> qryCov,
			final int slen,
			final int qlen) {
	
		super(qseqid, sseqid, qstart,qend,sstart,send, true);
		this.key      = key;
		this.subCov   = subCov;
		this.qryCov   = qryCov;
		this.slen     = slen;
		this.qlen     = qlen;
	}
	
	public CompoundAlignmentSegment(
			final String qseqid,   // query (e.g., gene) sequence id
			final String sseqid,   // subject (e.g., reference genome) sequence id
			final int qstart,      // start of alignment in query
			final int qend,        // end of alignment in query
			final int sstart,      // start of alignment in subject
			final int send,         // end of alignment in subject
			final RangeSet<Integer> subCov,
			final RangeSet<Integer> qryCov,
			final int slen,
			final int qlen) {
	
		super(qseqid, sseqid, qstart,qend,sstart,send, true);
		this.subCov   = subCov;
		this.qryCov   = qryCov;
		this.slen     = slen;
		this.qlen     = qlen;
	}
	
	public void setSamRecords(final List<SAMRecord> samRecords) {
		this.samRecords = samRecords;
	}
	
	public List<SAMRecord> getSamRecords() {
		return this.samRecords;
	}
	
	public void setScore(double score) {
		this.score = score;
	}
	
	public double getScore() {
		return this.score;
	}
	
	public String getKey() {
		return this.key;
	}
	
	public int getSubLen() {
		return this.slen;
	}
	
	public int getQryLen() {
		return this.qlen;
	}
	
	public RangeSet<Integer> getSubCov() {
		return this.subCov;
	}
	
	public RangeSet<Integer> getQryCov() {
		return this.qryCov;
	}
}
