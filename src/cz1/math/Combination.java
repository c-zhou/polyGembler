package cz1.math;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.paukov.combinatorics.Factory;
import org.paukov.combinatorics.Generator;
import org.paukov.combinatorics.ICombinatoricsVector;

public class Combination {

	private final static Logger myLogger = LogManager.getLogger(Combination.class);
	
	public static <T> List<List<T>> combination(T[] elements, int K){
		if(K > elements.length){
			myLogger.error("Invalid input, K > N");
			System.exit(1);
		}
		// Create the initial vector
		ICombinatoricsVector<T> initialVector = Factory.createVector(elements);
		// Create a simple combination generator to generate 3-combinations of the initial vector
		Generator<T> gen = Factory.createSimpleCombinationGenerator(initialVector, K);
		// Locate all possible combinations
		List<List<T>> combinations = new ArrayList<List<T>>();
		for (ICombinatoricsVector<T> combination : gen) {
			combinations.add(combination.getVector());
		}
		return combinations;
	}
	
	public static <T> List<List<T>> multiCombination(T[] elements, int K){
		// Create the initial vector
		ICombinatoricsVector<T> initialVector = Factory.createVector(elements);
		// Create a simple combination generator to generate 3-combinations of the initial vector
		Generator<T> gen = Factory.createMultiCombinationGenerator(initialVector, K);
		// Locate all possible combinations
		List<List<T>> combinations = new ArrayList<List<T>>();
		for (ICombinatoricsVector<T> combination : gen) {
			combinations.add(combination.getVector());
		}
		return combinations;
	}
	
	@SuppressWarnings("unchecked")
	public static <T> List<List<T>> combination(List<T> elements, int K){
		if(elements.isEmpty()) {
			List<List<T>> combs = new ArrayList<List<T>>();
			combs.add(new ArrayList<T>());
			return combs;
		}
		int N = elements.size();
		Object[] elements_array = new Object[N];
		elements.toArray(elements_array);
		return combination((T[])elements_array, K);
	}
	
	@SuppressWarnings("unchecked")
	public static <T> List<List<T>> combination(Set<T> elements, int K){
		if(elements.isEmpty()) {
			List<List<T>> combs = new ArrayList<List<T>>();
			combs.add(new ArrayList<T>());
			return combs;
		}
		int N = elements.size();
		Object[] elements_array = new Object[N];
		elements.toArray(elements_array);
		return combination((T[])elements_array, K);
	}
	
	@SuppressWarnings("unchecked")
	public static <T> List<List<T>> multiCombination(List<T> elements, int K){
		if(elements.isEmpty()) {
			List<List<T>> combs = new ArrayList<List<T>>();
			combs.add(new ArrayList<T>());
			return combs;
		}
		int N = elements.size();
		Object[] elements_array = new Object[N];
		elements.toArray(elements_array);
		return multiCombination((T[])elements_array, K);
	}
	
	@SuppressWarnings("unchecked")
	public static <T> List<List<T>> multiCombination(Set<T> elements, int K){
		if(elements.isEmpty()) {
			List<List<T>> combs = new ArrayList<List<T>>();
			combs.add(new ArrayList<T>());
			return combs;
		}
		int N = elements.size();
		Object[] elements_array = new Object[N];
		elements.toArray(elements_array);
		return multiCombination((T[])elements_array, K);
	}
	
	public static int nchoosek(int n, int k) {
		int nCk = 1;
		for (int i = 0; i < k; i++) {
            nCk = nCk * (n-i) / (i+1);
        }
		return nCk;
	}
	
	public static int nmultichoosek(int n, int k) {
		return nchoosek(n+k-1, k);
	}
	
	// all Stirling partitions of the second kind
	public static <T> List<List<List<T>>> stirlingS2Partitions(List<T> elements) {
		
		List<List<List<T>>> partitions = new ArrayList<>();
		if(elements.isEmpty()) {
			List<List<T>> empty = new ArrayList<>();
			partitions.add(empty);
			return partitions;
		}

		int limit = 1 << (elements.size() - 1);

        for (int j = 0; j < limit; ++j) {
            List<List<T>> parts = new ArrayList<>();
            List<T> part1 = new ArrayList<>();
            List<T> part2 = new ArrayList<>();
            parts.add(part1);
            parts.add(part2);
            int i = j;
            for (T item : elements) {
                parts.get(i&1).add(item);
                i >>= 1;
            }
            for (List<List<T>> b : stirlingS2Partitions(part2)) {
                List<List<T>> holder = new ArrayList<>();
                holder.add(part1);
                holder.addAll(b);
                partitions.add(holder);
            }
        }
        return partitions;
    }
}
