package cz1.util;

import java.io.PrintStream;
import java.io.OutputStream;

/******************************************************************************
 *  Compilation:  javac JohnsonTrotter.java
 *  Execution:    java JohnsonTrotter n
 *  
 *  Generate permutations by transposing adjacent elements using the
 *  Johnson-Trotter algorithm.
 *
 *  This program is a Java version based on the program SJT.c
 *  writen by Frank Ruskey,
 *  modified by Chenxi Zhou.
 *  
 *     http://theory.cs.uvic.ca/inf/perm/PermInfo.html
 * 
 *  % java JohnsonTrotter 3
 *  012   (2 1)
 *  021   (1 0)
 *  201   (2 1)
 *  210   (0 1)
 *  120   (1 2)
 *  102   (0 1)
 *
 ******************************************************************************/


public class JohnsonTrotter {

	private final static PrintStream StdOut = //System.out;
			new PrintStream(new OutputStream(){
				public void write(int b) {
				}
			});
	private static int cursor = 1;
	
	public static char[][] perm(char[] cs) {
		cursor = 1;
		int[][] p = perm(cs.length);
		char[][] ps = new char[p.length][2];
		for(int i=0; i<p.length; i++)
			for(int j=0; j<2; j++)
				ps[i][j] = cs[p[i][j]];
		return ps;
	}
	
    public static int[][] perm(int n) {
    	cursor = 1;
    	int k = perm_count(n);
    	int[][] swap = new int[k][2];
        int[] p   = new int[n];     // permutation
        int[] pi  = new int[n];     // inverse permutation
        int[] dir = new int[n];     // direction = +1 or -1
        for (int i = 0; i < n; i++) {
            dir[i] = -1;
            p[i]  = i;
            pi[i] = i;
        }
        perm(0, p, pi, dir, swap);
        StdOut.printf("   (0 1)\n");
        return swap;
    }

    private static int perm_count(int n) {
		// TODO Auto-generated method stub
		if(n==0) return 1;
    	return n*perm_count(n-1);
	}

	public static void perm(int n, int[] p, int[] pi, 
			int[] dir, int[][] swap) { 

        // base case - print out permutation
        if (n >= p.length) {
            for (int i = 0; i < p.length; i++)
                StdOut.print(p[i]);
            return;
        }

        perm(n+1, p, pi, dir, swap);
        for (int i = 0; i <= n-1; i++) {

            // swap 
            StdOut.printf("   (%d %d)\n", pi[n], pi[n] + dir[n]);
            swap[cursor++] = new int[]{p[pi[n]], p[pi[n] + dir[n]]};
            int z = p[pi[n] + dir[n]];
            p[pi[n]] = z;
            p[pi[n] + dir[n]] = n;
            pi[z] = pi[n];
            pi[n] = pi[n] + dir[n];  

            perm(n+1, p, pi, dir, swap); 
        }
        dir[n] = -dir[n];
    }

    public static void main(String[] args) {
        int[][] swap = perm(6);
        for(int i=0; i<swap.length; i++)
        	StdOut.printf("(%d %d)\n", swap[i][0], swap[i][1]);
    }
}