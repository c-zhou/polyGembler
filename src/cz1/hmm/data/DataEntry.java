package cz1.hmm.data;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;

public class DataEntry {
	private String id = null;	// contig name
	private String[] marker_id = null;
	private double[] position = null; // positions along contig
	private List<String[]> allele = null; // alleles
	private List<List<int[]>> ad = null; // allele depth
	private List<List<double[]>> gl = null; // genotype likelihood
	private List<List<String[]>> gt = null; //genotypes
	private String[] sample = null; //samples
	
	public DataEntry(String id, 
			double[] position,
			List<String[]> allele,
			List<List<int[]>> ad,
			List<List<double[]>> gl,
			List<List<String[]>> gt,
			List<String> sample) {
		this.id = id;
		this.position = position;
		this.markers();
		this.allele = allele;
		this.ad = ad;
		this.gl = gl;
		this.gt = gt;
		this.sample = sample.toArray(new String[sample.size()]);
	}
	
	private void markers() {
		// TODO Auto-generated method stub
		this.marker_id = new String[this.position.length];
		for(int i=0; i<marker_id.length; i++)
			marker_id[i] = this.id+"_"+(int) this.position[i];
	}
	
	public double[] getPosition() {
		return this.position;
	}
	
	public int[] getIntegerPosition() {
		int[] d = new int[this.position.length];
		for(int i=0; i<d.length; i++)
			d[i] = (int) this.position[i];
		return d;
	}
	
	public List<String[]> getAllele() {
		return this.allele;
	}
	
	public List<List<int[]>> getAlleleDepth() {
		return this.ad;
	}
	
	public List<List<double[]>> getGenotypeLikelihood() {
		return this.gl;
	}
	
	public List<List<String[]>> getGenotype() {
		return this.gt;
	}
	
	public String[] getMarker() {
		return this.marker_id;
	}
	
	public String[] getSample() {
		// TODO Auto-generated method stub
		return this.sample;
	}

	public void get(int from) {
		// TODO Auto-generated method stub
		int n = this.sample.length;
		for(int i=0; i<this.allele.size(); i++) {
			this.ad.get(i).subList(from, n).clear();
			this.gl.get(i).subList(from, n).clear();
			this.gt.get(i).subList(from, n).clear();
		}
		this.sample = Arrays.copyOfRange(this.sample, 0, from);
	}
	
	public void remove(int index) {
		// TODO Auto-generated method stub
		for(int i=0; i<this.allele.size(); i++) {
			this.ad.get(i).remove(index);
			this.gl.get(i).remove(index);
			this.gt.get(i).remove(index);
		}
		this.sample = ArrayUtils.remove(this.sample, index);
	}

	public void remove(int[] is) {
		// TODO Auto-generated method stub
		for(int i=0; i<is.length; i++) {
			this.position = ArrayUtils.
					remove(this.position,is[i]-i);
			this.marker_id = ArrayUtils.
					remove(this.marker_id,is[i]-i);
			this.allele.remove(is[i]-i);
			this.ad.remove(is[i]-i);
			this.gl.remove(is[i]-i);
			this.gt.remove(is[i]-i);
		}
	}

	public void reverse() {
		// TODO Auto-generated method stub
		this.id = this.id+"#R";
		reverse(this.position);
		reverse(this.marker_id);
		// reverse(this.baf);
		Collections.reverse(this.allele);
		Collections.reverse(this.ad);
		Collections.reverse(this.gl);
		Collections.reverse(this.gt);
		rescale();
	}

	private void rescale() {
		// TODO Auto-generated method stub
		double p = this.position[0], c;
		this.position[0] = 0;
		for(int i=1; i<this.position.length; i++) {
			c = this.position[i];
			this.position[i] = this.position[i-1]+
					Math.abs(c-p);
			p = c;
		}
	}

	private void reverse(double[][] mat) {
		// TODO Auto-generated method stub
		double temp;
		for (int start=0,end=mat.length-1;
				start<=end;
				start++,end--) {
			for(int i=0; 
					i<mat[start].length; 
					i++) {
				temp = mat[start][i];
				mat[start][i] = mat[end][i];
				mat[end][i] = temp;
			}
		}
	}

	private void reverse(double[] arr) {
		// TODO Auto-generated method stub
		double temp;
		for (int start=0,end=arr.length-1;
				start<=end;
				start++,end--) {
			temp = arr[start];
			arr[start] = arr[end];
			arr[end] = temp;
		}
	}
	
	private void reverse(String[] arr) {
		// TODO Auto-generated method stub
		String temp;
		for (int start=0,end=arr.length-1;
				start<=end;
				start++,end--) {
			temp = arr[start];
			arr[start] = arr[end];
			arr[end] = temp;
		}
	}

	public int modelLength() {
		// TODO Auto-generated method stub
		return this.allele.size();
	}

	public void addAll(DataEntry dataEntry, 
			double distance) {
		// TODO Auto-generated method stub
		this.id = this.id+":"+dataEntry.id;
		add(dataEntry.position, distance);
		add(dataEntry.marker_id);
		// add(dataEntry.baf);
		this.allele.addAll(dataEntry.allele);
		if(this.ad!=null) this.ad.addAll(dataEntry.ad);
		if(this.gl!=null) this.gl.addAll(dataEntry.gl);
		if(this.gt!=null) this.gt.addAll(dataEntry.gt);
	}

	private void add(double[] position2, double distance) {
		// TODO Auto-generated method stub
		int n = this.modelLength();
		double[] position = new double[n+position2.length];
		System.arraycopy(this.position, 0, position, 0, n);
		double start = this.position[n-1]+distance-position2[0];
		for(int i=this.position.length; i<position.length; i++)
			position[i] = position2[i-this.position.length]+start;
		this.position = position;
	}
	
	private void add(String[] marker_id) {
		int n = this.modelLength(), n2 = marker_id.length;
		String[] markers = new String[n+n2];
		System.arraycopy(this.marker_id, 0, markers, 0, n);
		System.arraycopy(marker_id, 0, markers, n, n2);
		this.marker_id = markers;
	}

	public List<Character> getAlleleA() {
		// TODO Auto-generated method stub
		return null;
	}

	public List<Character> getAlleleB() {
		// TODO Auto-generated method stub
		return null;
	}

	public String getId() {
		// TODO Auto-generated method stub
		return this.id;
	}
}
