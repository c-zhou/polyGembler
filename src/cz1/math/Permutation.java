package cz1.math;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.apache.log4j.Logger;
import org.paukov.combinatorics.Factory;
import org.paukov.combinatorics.Generator;
import org.paukov.combinatorics.ICombinatoricsVector;

public class Permutation {
	private final static Logger myLogger = Logger.getLogger(Permutation.class);
	
	public static <T> List<List<T>> permutation(T[] elements){
		// Create the initial vector
		ICombinatoricsVector<T> initialVector = Factory.createVector(elements);
		// Create a simple combination generator to generate 3-combinations of the initial vector
		Generator<T> gen = Factory.createPermutationGenerator(initialVector);
		// Locate all possible combinations
		List<List<T>> permutations = new ArrayList<List<T>>();
		for (ICombinatoricsVector<T> combination : gen) {
			permutations.add(combination.getVector());
		}
		return permutations;
	}
	
	@SuppressWarnings("unchecked")
	public static <T> List<List<T>> permutation(List<T> elements){
		if(elements.isEmpty()) {
			List<List<T>> perms = new ArrayList<List<T>>();
			perms.add(new ArrayList<T>());
			return perms;
		}
		int N = elements.size();
		Object[] elements_array = new Object[N];
		elements.toArray(elements_array);
		return permutation((T[])elements_array);
	}

	@SuppressWarnings("unchecked")
	public static <T> List<List<T>> permutation(Set<T> elements){
		if(elements.isEmpty()) {
			List<List<T>> perms = new ArrayList<List<T>>();
			perms.add(new ArrayList<T>());
			return perms;
		}
		int N = elements.size();
		Object[] elements_array = new Object[N];
		elements.toArray(elements_array);
		return permutation((T[])elements_array);
	}
	
	public static List<List<Integer>> permutation(int ele) {
		// TODO Auto-generated method stub
		Integer[] elements = new Integer[ele];
		for(int i=0; i<ele; i++)
			elements[i] = i;
		return permutation(elements);
	}
	
	public static <T> List<List<T>> multiPermutation(T[] elements, int K) {
		// TODO Auto-generated method stub
		// Create the initial vector
		ICombinatoricsVector<T> initialVector = Factory.createVector(elements);
		// Create a simple combination generator to generate 3-combinations of the initial vector
		Generator<T> gen = Factory.createSimpleCombinationGenerator(initialVector, K);
		// Locate all possible combinations
		List<List<T>> combinations = new ArrayList<List<T>>();
		for (ICombinatoricsVector<T> combination : gen) {
			combinations.add(combination.getVector());
		}
		List<List<T>> permutations = new ArrayList<List<T>>();
		for(List<T> combination : combinations) {
			List<List<T>> perm = permutation(combination);
			permutations.addAll(perm);
		}
		return permutations;
	}
	
	public static long factorial(int n){
		if(n<0) {
			myLogger.error("Can NOT factorial a negative number. Program halted.");
			System.exit(1);
		}
		long f=1;
		for(int k=1; k<=n; k++) {
			if(Long.MAX_VALUE/f<n) { 
				myLogger.error("Factorial out of range (larger than "
						+ Long.MAX_VALUE+"). Program halted.");
				System.exit(1);
			}
			f*=k;
		}
		return f;
	}
}
