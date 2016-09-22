package cz1.util;

import java.util.Arrays;
import java.util.List;
import java.util.Collections;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.math.NumberUtils;
import org.apache.commons.math3.stat.StatUtils;

public class Algebra {
	public static int largestElementIndex(double[] array) {
		// TODO Auto-generated method stub
		double largest = Double.NEGATIVE_INFINITY;
		int index = -1;
		for(int i=0; i<array.length; i++)
			if(array[i] > largest) {
				index = i;
				largest = array[i];
			}
		return index;
	}
	
	/*** sum of logs
	 * @param array
	 * @return
	 * calculate products of probabilities considering under-flow problem
	 */
	public static double sumLogs(double[] array) {
		if(allNegativeInfinity(array)) return Double.NEGATIVE_INFINITY;
		if(array.length>1) 
			return array[0]+Math.log(1+Math.exp(
					sumLogs(Arrays.copyOfRange(array, 1, array.length))-array[0]));
		else
			return array[0];
	}

	public static double sumLogs(Double[] array) {
		return sumLogs(ArrayUtils.toPrimitive(array));
	}
	
	/*** sum of exponentials
	 * @param array
	 * @return
	 * similar to sum of logs
	 */
	public static double sumExps(double[] array) {
		//if(allNegativeInfinity(array)) return Double.NEGATIVE_INFINITY;
		double a = NumberUtils.max(array), s = 0;
		if(a==Double.NEGATIVE_INFINITY) 
			return Double.NEGATIVE_INFINITY;
		double r;
		for(int i=0; i<array.length; i++) {
			r = array[i]-a;
			if(r>Constants.MIN_EXP_DOUBLE)
				s += Math.exp(r);
		}
		return a+Math.log(s);
	}

	public static double sumExps(Double[] array) {
		return sumExps(ArrayUtils.toPrimitive(array));
	}
	
	private static boolean allNegativeInfinity(double[] array) {
		// TODO Auto-generated method stub
		for(int i=0; i<array.length; i++)
			if(array[i] != Double.NEGATIVE_INFINITY) return false;
		return true;
	}

	public static int sum12n(int i) {
		// TODO Auto-generated method stub
		return i*(i+1)/2;
	}

	public static int distance(List<Character> list0, List<Character> list1) {
		// TODO Auto-generated method stub
		if(list0.size() != list1.size()) {
			System.err.println("two lists have different size. program halt.");
			System.exit(1);
		}
		int d = 0;
		for(int i=0; i<list0.size(); i++)
			d += list0.get(i)==list1.get(i) ? 0 : 1;
		return d;
	}

	public static double[] normalize(double[] array) {
		// TODO Auto-generated method stub
		double s = StatUtils.sum(array);
		if(s==0) return array;
		for(int i=0; i<array.length; i++) array[i]/=s;
		return array;
	}
	
	public static double[] normalizedLog(double[] array) {
		// TODO Auto-generated method stub
		double s = StatUtils.sum(array);
		for(int i=0; i<array.length; i++) 
			array[i] = Math.log(array[i]/s);
		return array;
	}

	public static int maxIndex(double[] array) {
		// TODO Auto-generated method stub
		if(array.length<1) return -1;
		double d = array[0];
		int index = 0;
		for(int i=1; i<array.length; i++) {
			if(array[i]>d) {
				d = array[i];
				index = i;
			}
		}
		return index;
	}
	
	public static int minIndex(double[] array) {
		// TODO Auto-generated method stub
		if(array.length<1) return -1;
		double d = array[0];
		int index = 0;
		for(int i=1; i<array.length; i++) {
			if(array[i]<d) {
				d = array[i];
				index = i;
			}
		}
		return index;
	}
	
	public static double min(double[][] array) {
		double min = Double.POSITIVE_INFINITY;
		for(int i=0; i<array.length; i++) {
			double m = StatUtils.min(array[i]);
			if(m<min) min = m;
		}
		return min;
	}

	public static double max(double[][] array) {
		// TODO Auto-generated method stub
		double max = Double.NEGATIVE_INFINITY;
		for(int i=0; i<array.length; i++) {
			double m = StatUtils.max(array[i]);
			if(m>max) max = m;
		}
		return max;
	}
	
	public static double sum(Double[] array) {
		// TODO Auto-generated method stub
		double s = 0;
		for(Double a : array) s += a;
		return s;
	}
	
	public static double sum(double[] array) {
		// TODO Auto-generated method stub
		double s = 0;
		for(double a : array) s += a;
		return s;
	}

	public static double[][] transpose(double[][] mat) {
		int a = mat.length, b = mat[0].length;
		double[][] t = new double[b][a];
		for(int i=0; i<a; i++)
			for(int j=0; j<b; j++)
				t[j][i] = mat[i][j];
		return t;
	}
	
	public static double[] sum(double[][] mat, 
			boolean byrow, 
			boolean log) {
		if(!byrow) mat = transpose(mat);
		double[] s = new double[mat.length];
		for(int i=0; i<s.length; i++)
			s[i] = log ? sumExps(mat[i]) : sum(mat[i]);
		if(!byrow) mat = transpose(mat);
		return s;
	}
	
	public static double[] sum(double[][] mat, 
			boolean byrow) {
		return sum(mat, byrow, false);
	}
	
	public static double sum(double[][] array) {
		// TODO Auto-generated method stub
		double s = 0;
		for(int i=0; i<array.length; i++)
			s += sum(array[i]);
		return s;
	}
	
    public static int[] sort(int[] array, boolean descending) {
       List<Integer> list = Arrays.asList(
               ArrayUtils.toObject(array));
       Collections.sort(list);
       if(descending) Collections.reverse(list);
       return ArrayUtils.toPrimitive(list.toArray(
                   new Integer[list.size()]));
    }

	public static int min(int[] array) {
		// TODO Auto-generated method stub
		int minimum = Integer.MAX_VALUE;
		for(int i=0; i<array.length; i++)
			if(array[i]<minimum)
				minimum = array[i];
		return minimum;
	}

	public static int positiveMinimum(int[] array) {
		// TODO Auto-generated method stub
		int minimum = Integer.MAX_VALUE;
		boolean b = true;
		for(int i=0; i<array.length; i++)
			if(array[i]>0) {
				b = false;
				break;
			}
		if(b) return -1;
		for(int i=0; i<array.length; i++)
			if(array[i]>0 && array[i]<minimum)
				minimum = array[i];
		return minimum;
	}

	public static int positiveMaximum(int[] array) {
		// TODO Auto-generated method stub
		int maximum = Integer.MIN_VALUE;
		boolean b = true;
		for(int i=0; i<array.length; i++)
			if(array[i]>0) {
				b = false;
				break;
			}
		if(b) return -1;
		for(int i=0; i<array.length; i++)
			if(array[i]>0 && array[i]>maximum)
				maximum = array[i];
		return maximum;
	}
	
	public static double positiveMinimum(double[][] array) {
		// TODO Auto-generated method stub
		double minimum = Double.POSITIVE_INFINITY;
		for(int i=0; i<array.length; i++) 
			for(int j=0; j<array[i].length; j++)
				minimum = array[i][j]>0 ? (array[i][j]<minimum ? 
						array[i][j] : minimum) : minimum;
		return minimum;
	}

	public static double positiveMaximum(double[][] array) {
		// TODO Auto-generated method stub
		double maximum = Double.NEGATIVE_INFINITY;
		for(int i=0; i<array.length; i++) 
			for(int j=0; j<array[i].length; j++)
				maximum = array[i][j]>0 ? (array[i][j]>maximum ? 
						array[i][j] : maximum) : maximum;
		return maximum;
	}

	public static double positiveMinimum(double[] ds) {
		// TODO Auto-generated method stub
		double minimum = Double.POSITIVE_INFINITY;
		for(int i=0; i<ds.length; i++) 
			minimum = ds[i]>0 ? (ds[i]<minimum ? 
						ds[i] : minimum) : minimum;
		return minimum;
	}
	
	public static double positiveMaximum(double[] ds) {
		// TODO Auto-generated method stub
		double maximum = Double.NEGATIVE_INFINITY;
		for(int i=0; i<ds.length; i++)
			maximum = ds[i]>0 ? (ds[i]>maximum ? 
				ds[i] : maximum) : maximum;
		return maximum;
	}
}
