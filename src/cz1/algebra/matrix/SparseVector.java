package cz1.algebra.matrix;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

public class SparseVector extends HashMap<Integer, Double> {

	/**
	 * 
	 */
	private static final long serialVersionUID = -1081854300031972017L;
	private final int length;
	
	public SparseVector(int length) { 
        super();
        this.length = length;
    } 
 	
	@Override 
	public Double get(Object key) { 
		Double b = super.get(key); 
		if (b == null) { 
			return 0.; 
		} 
		return b; 
	}
 
    @Override 
    public Double put(Integer key, Double value) {
        if (value == 0) return remove(key);
        return super.put(key, value); 
    }
    
    public void normalise() { 
    	if(super.isEmpty()) return;
    	double invsum = 1. / sum(); 
        for (int i : keySet()) { 
        	put(i, invsum*get(i)); 
        }
    }
    
    private double sum() { 
        double sum = 0; 
        for (double a : values()) { 
        	sum += a; 
        } 
        return sum; 
    }
    
    public final int getLength() { 
        return length; 
    }
    
    public void prune(double threshold) { 
        for (Iterator<Integer> it = keySet().iterator(); it.hasNext();) { 
            int key = it.next(); 
            if (Math.abs(get(key)) < threshold) { 
                it.remove(); 
            } 
        } 
    }

	public double[] getDense() {
		// TODO Auto-generated method stub
		double[] vector = new double[length];
		for(Map.Entry<Integer, Double> entry : entrySet())
			vector[entry.getKey()] = entry.getValue();
		return vector;
	}

	public void hadamardPower(double s) {
		// TODO Auto-generated method stub
		final Set<Integer> keys = new HashSet<Integer>();
		keys.addAll(this.keySet());
		for(int i : keys) { 
			put(i, Math.pow(get(i), s)); 
		}
	}

	public double max() {
		// TODO Auto-generated method stub
		if(this.isEmpty()) return 0.;
		double max = Double.NEGATIVE_INFINITY;
		for(double d : this.values())
			if(d>max) max = d;
		return max;
	}

	public double powSum(double p) {
		// TODO Auto-generated method stub
		double sum = 0.;
		for(double d : this.values())
			sum += Math.pow(d, p);
		return sum;
	}
}
