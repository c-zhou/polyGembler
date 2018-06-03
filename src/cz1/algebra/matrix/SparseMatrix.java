package cz1.algebra.matrix;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import cz1.util.Utils;

public class SparseMatrix extends ArrayList<SparseVector> {

	/**
	 * 
	 */
	private static final long serialVersionUID = 8363795995066777942L;
	private final int rows;
	private final int cols;
	
	public SparseMatrix(int rows, int cols) { 
        this.rows = rows;
        this.cols = cols;
        for(int i=0; i<rows; i++)
        	super.add(new SparseVector(cols));
    }
	
	public double[][] getDense() { 
        double[][] mat = new double[rows][]; 
        for (int i = 0; i < size(); i++) { 
            mat[i] = get(i).getDense(); 
        } 
        return mat; 
    } 
	
	public Set<Integer> getColumn(int j) {
		if(j>=cols) throw new ArrayIndexOutOfBoundsException();
		final Set<Integer> colNonZero = new HashSet<Integer>();
		for(int i=0; i<rows; i++) {
			if(get(i, j)>0) colNonZero.add(i);
		}
		return colNonZero;
	}
	
	public double get(int i, int j) { 
        checkBound(i, j);
		return get(i).get(j);
    }
	
	private void checkBound(int i, int j) {
		// TODO Auto-generated method stub
		if (i>=rows||j>=cols) 
        	throw new ArrayIndexOutOfBoundsException();
	}

	public void set(int i, int j, double a) { 
        checkBound(i, j);
		get(i).put(j, a); 
    }
	
	public void normalise() { 
		for (SparseVector vec : this)
			vec.normalise(); 
	}

	public void setDiag(double diag) {
		// TODO Auto-generated method stub
		int n = Math.min(rows, cols);
		for(int i=0; i<n; i++) set(i, i, diag);
	}

	public void expand() {
		// TODO Auto-generated method stub
		this.square();
	}

	private void square() {
		// TODO Auto-generated method stub
		SparseMatrix mat = this.deepCopy();
		for(int i=0; i<rows; i++) { 
            SparseVector v = mat.get(i);
			for(int j=0; j<cols; j++) { 
				double a = 0;
				for(int k : v.keySet()) {
					a += v.get(k)*mat.get(k, j);
				}
				this.set(i, j, a);
			}
        }		
	}

	public SparseMatrix deepCopy() {
		// TODO Auto-generated method stub
		SparseMatrix mat = new SparseMatrix(rows, cols);
		for(int i=0; i<rows; i++) {
			SparseVector a = this.get(i);
			SparseVector b = mat.get(i) ;
			for(Map.Entry<Integer, Double> entry : a.entrySet()) 
				b.put(entry.getKey(), entry.getValue());
		}
		return mat;
	}

	public double inflate(double pGamma, double maxZero) {
		// TODO Auto-generated method stub
		hadamardPower(pGamma);
		prune(maxZero);
		normalise();
		
		double res = 0.;
		for (int i=0; i<rows; i++) { 
            SparseVector row = get(i); 
            if(row.isEmpty()) continue;
            double max = row.max(); 
            double sumsq = row.powSum(2.); 
            res = Math.max(res, max - sumsq); 
        }
		
        return res;
	}

	public void hadamardPower(double s) { 
        for (int i=0; i<rows; i++) { 
            get(i).hadamardPower(s); 
        } 
    }
	
	public void prune(double threshold) { 
        for (int i=0; i<rows; i++) 
            get(i).prune(threshold); 
    }

	public int getRows() {
		// TODO Auto-generated method stub
		return rows;
	}
	
	public int getCols() {
		// TODO Auto-generated method stub
		return cols;
	}
}
