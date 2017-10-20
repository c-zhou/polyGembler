package cz1.hmm.data;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.inference.ChiSquareTest;

import cz1.util.Algebra;
import cz1.util.Combination;
import cz1.util.Constants;
import cz1.util.Utils;

public class DataEntry {
	private String id = null;	// contig name
	private String[] marker_id = null;
	private double[] position = null; // positions along contig
	private List<String[]> allele = null; // alleles
	private List<List<int[]>> ad = null; // allele depth
	private List<List<double[]>> gl = null; // genotype likelihood
	private List<List<String[]>> gt = null; //genotypes
	private String[] sample = null; //samples
	private double[][] configuration = null;
	// private double[][] baf;
	
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
		this.configuration = this.configuration();
		// this.baf();
		// this.filter();
	}
	
	private void markers() {
		// TODO Auto-generated method stub
		this.marker_id = new String[this.position.length];
		for(int i=0; i<marker_id.length; i++)
			marker_id[i] = this.id+"_"+(int) this.position[i];
	}

	/**
	private void baf() {
		// TODO Auto-generated method stub
		//if(this.pl.get(0).get(0).length>3)
		//	throw new RuntimeException("need to change for polyploidy!!!");
		int cp = this.pl.get(0).get(0).length;
		this.baf = new double[this.allele.size()][2];
		for(int i=0; i<this.baf.length; i++) {
			double f0 = 0.0, f1 = 0.0;
			int n = 0;
			for(int j=0; j<this.pl.get(i).size(); j++) {
				double[] ll = this.pl.get(i).get(j);
				if(ll[0]==-1) 
					for(int k=0; k<ll.length; k++)
						f0 += 1.0/cp*k;
				else {
					for(int k=1; k<ll.length; k++) {
						f0 += ll[k]*k;
						f1 += ll[k]*k;
					}
					n++;
				}
			}
			this.baf[i][0] = f0/this.pl.get(i).size()
					/(cp-1);
			this.baf[i][1] = f1/n/(cp-1);
		}
	}
	**/

	private void filter() {
		// TODO Auto-generated method stub
		List<Integer> rl = new ArrayList<Integer>();
		int h = Constants._ploidy_H;
		double thresh = 1.0/h/4;
		for(int i=0; i<this.allele.size(); i++) {
			if(this.allele.get(i).length != 2) {
				rl.add(i);
				continue;
			}
			double[] ds = new double[2];
			for(int j=0; j<this.gl.get(i).size(); j++) {
				double[] d = this.gl.get(i).get(j);
				for(int k=0; k<d.length; k++) 
					if(d[k]>0) {
						ds[0] += (h-k)*d[k];
						ds[1] += k*d[k];
					}
			}
			if(StatUtils.min(ds)/StatUtils.sum(ds)<thresh) 
				rl.add(i);
		}
		this.remove(ArrayUtils.toPrimitive(
				rl.toArray(new Integer[rl.size()])));
		return;
	}
	
	private void filter1() {
		// TODO Auto-generated method stub
		List<Integer> rl = new ArrayList<Integer>();
		for(int i=0; i<this.allele.size(); i++) {
			if(this.allele.get(i).length != 2) {
				rl.add(i);
				continue;
			}
			double[] ds = new double[Constants._ploidy_H+1];
			for(int j=0; j<this.gl.get(i).size(); j++) {
				double[] d = this.gl.get(i).get(j);
				for(int k=0; k<d.length; k++)
					ds[k] += d[k];
			}
			if((ds[0]+ds[Constants._ploidy_H])/
					StatUtils.sum(ds)>0.9) 
				rl.add(i);
		}
		this.remove(ArrayUtils.toPrimitive(rl.toArray(new Integer[rl.size()])));
		return;
	}
	
	private void filter2() {
		// TODO Auto-generated method stub
		List<Integer> rl = new ArrayList<Integer>();
		for(int i=0; i<this.allele.size(); i++) {
			if(this.allele.get(i).length != 2) {
				rl.add(i);
				continue;
			}
			double[] ds = new double[this.configuration[0].length];
			for(int j=0; j<this.gl.get(i).size(); j++) {
				double[] d = this.gl.get(i).get(j);
				for(int k=0; k<d.length; k++)
					ds[k] += d[k]>0 ? d[k] : 0;
			}
			System.out.println(distortion(ds));
			Utils.print(ds);
			if(distortion(ds)<0.01) rl.add(i);
		}
		this.remove(ArrayUtils.toPrimitive(rl.toArray(new Integer[rl.size()])));
		return;
	}
	
	private double distortion(double[] ds) {
		// TODO Auto-generated method stub
		double[] dt = new double[this.configuration.length];
		long[] dsL = new long[ds.length];
		double s = StatUtils.sum(ds);
		for(int i=0; i<dsL.length; i++)
			dsL[i] = Math.round(ds[i]);
		for(int i=0; i<dt.length; i++) {
			//for(int j=0; j<ds.length; j++)
			//	dt[i] += Math.pow(
			//			ds[j]-this.configuration[i][j], 
			//			2);
			double[] conf = this.configuration[i];
			List<Long> ds_i = new ArrayList<Long>();
			List<Double> conf_i = new ArrayList<Double>();
			boolean b = false;
			for(int j=0; j<conf.length; j++) 
				if(conf[j]!=0) {
					ds_i.add(dsL[j]);
					conf_i.add(conf[j]);		
				} else if (ds[j]/s>0.3) {
					//System.out.println(ds[j]/s);
					b = true;
					break;
				}
			if(b) continue;
			dt[i] = new ChiSquareTest().chiSquareTest(
					ArrayUtils.toPrimitive(conf_i.
							toArray(new Double[conf_i.size()])), 
					ArrayUtils.toPrimitive(ds_i.
							toArray(new Long[ds_i.size()])));
		}
		return Math.sqrt(StatUtils.max(dt));
	}

	private double[][] configuration() {
		int h = Constants._ploidy_H, g = h+1;
		char[][] genotypes = new char[g][h];
		for(int i=0; i<g; i++) {
			Arrays.fill(genotypes[i], 0, h-i, 'A');
			Arrays.fill(genotypes[i], h-i, h, 'B');
		}
		List<double[]> configs = new ArrayList<double[]>();
		for(int i=0; i<g; i++) {
			int[] gameteA = this.gamete(genotypes[i]);
			for(int j=i; j<g; j++)
				if( i==0&&j==0 || 
					i==(g-1)&&j==(g-1) ||
					i==0&&j==(g-1) )
					//two homozygotes
					continue;
				else {
					int[] gameteB = this.gamete(genotypes[j]);
					double[] os = new double[g];
					for(int k=0; k<gameteA.length; k++) 
						for(int l=0; l<gameteB.length; l++) 
							os[gameteA[k]+gameteB[l]] += 1.0;
					addConf(configs, Algebra.normalize(os));
				}
		}
		double[][] conf = new double[configs.size()][g];
		for(int i=0; i<configs.size(); i++) conf[i] = configs.get(i);
		return conf;
	}

	private void addConf(List<double[]> configs, 
			double[] norm) {
		// TODO Auto-generated method stub
		boolean exist = false;
		for(int i=0; i<configs.size(); i++) {
			double[] a = configs.get(i);
			boolean y = true;
			for(int j=0; j<a.length; j++) {
				y = y && a[j]==norm[j];
			}
			exist = y;
		}
		if(!exist) configs.add(norm);
	}

	private int[] gamete(char[] genotype) {
		// TODO Auto-generated method stub
		List<List<Character>> gametes = 
				Combination.combination(
						ArrayUtils.toObject(genotype), 
						genotype.length/2);
		int[] a = new int[gametes.size()];
		for(int j=0; j<a.length; j++)
			for(int k=0; k<gametes.get(j).size(); k++)
				if(gametes.get(j).get(k)=='A')
					a[j] += 1.0;
		return a;
	}

	public String getId() {
		return this.id;
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
	
	public void print() {
		this.printId();
		this.printPosition();
		this.printAllele();
		this.printAlleleDepth();
		this.printGenotypeLikelihood();
		this.printGenotype();
	}

	private void printGenotype() {
		// TODO Auto-generated method stub
		Utils.println("---Genotype---");
		int m = 1;
		Utils.print("Sample "+"\t\t|\t");
		for(int i=1; i<=this.sample.length; i++)
			Utils.print(String.format("%04d", i)+"\t");
		Utils.print("\n");
		for(List<String[]> list : this.gt) {
			Utils.print("Marker "+m+++"\t|\t");
			for(String[] chars : list) 
				Utils.print(StringUtils.join(chars, ',')+";\t");
			Utils.println();
		}
	}

	private void printGenotypeLikelihood() {
		// TODO Auto-generated method stub
		Utils.println("---Likelihood---");
		int m = 1;
		for(List<double[]> list : this.gl) {
			Utils.print("Marker "+m+++"\t|\t");
			for(double[] doubles : list) 
				Utils.print(StringUtils.join(doubles, ',')+"; ");
			Utils.println();
		}
	}

	private void printAlleleDepth() {
		// TODO Auto-generated method stub
		Utils.println("---Allele Depth---");
		int m = 1;
		for(List<int[]> list : this.ad) {
			Utils.print("Marker "+m+++"\t|\t");
			for(int[] ints : list) 
				Utils.print(StringUtils.join(ints, ',')+"; ");
			Utils.println();
		}
	}

	private void printAllele() {
		// TODO Auto-generated method stub
		Utils.println("---Allele---");
		for(String[] allele : this.allele)
			Utils.print(StringUtils.join(allele, ',')+"\t");
		Utils.println();
	}

	private void printPosition() {
		// TODO Auto-generated method stub
		Utils.println("---Position---");
		Utils.print(this.position);
	}

	private void printId() {
		// TODO Auto-generated method stub
		Utils.println("---Block Id---");
		Utils.println(this.id);
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

	/**
	private void add(double[][] baf2) {
		// TODO Auto-generated method stub
		int n = this.baf.length, m = baf2.length;
		double[][] baf = new double[n+m][];
		for(int i=0; i<n; i++)
			baf[i] = this.baf[i];
		for(int i=0; i<m; i++)
			baf[i+n] = baf2[i].clone();
		this.baf = baf;
	}
	**/

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
}
