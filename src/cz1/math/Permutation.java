package cz1.math;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.log4j.Logger;
import org.paukov.combinatorics.Factory;
import org.paukov.combinatorics.Generator;
import org.paukov.combinatorics.ICombinatoricsVector;

public class Permutation {
	private final static Logger myLogger = Logger.getLogger(Permutation.class);
	
	public static ArrayList<List<Character>> permutation(Character[]  elements){
		// Create the initial vector
		ICombinatoricsVector<Character> initialVector = Factory.createVector(elements);
		// Create a simple combination generator to generate 3-combinations of the initial vector
		Generator<Character> gen = Factory.createPermutationGenerator(initialVector);
		// Locate all possible combinations
		ArrayList<List<Character>> permutations = new ArrayList<List<Character>>();
		for (ICombinatoricsVector<Character> combination : gen) {
			permutations.add(combination.getVector());
		}
		return permutations;
	}
	
	public static ArrayList<List<Integer>> permutation(Integer[]  elements){
		// Create the initial vector
		ICombinatoricsVector<Integer> initialVector = Factory.createVector(elements);
		// Create a simple combination generator to generate 3-combinations of the initial vector
		Generator<Integer> gen = Factory.createPermutationGenerator(initialVector);
		// Locate all possible combinations
		ArrayList<List<Integer>> permutations = new ArrayList<List<Integer>>();
		for (ICombinatoricsVector<Integer> combination : gen) {
			permutations.add(combination.getVector());
		}
		return permutations;
	}
	
	public static ArrayList<List<Character>> permutation(List<Character>  elements){
		int N = elements.size();
		Character[] elements_array = new Character[N];
		elements.toArray(elements_array);
		return permutation(elements_array);
	}

	public static ArrayList<List<Character>> permutation(char[] elements){
		return permutation(ArrayUtils.toObject(elements));
	}
	
	public static ArrayList<List<Character>> permutation(Set<Character>  elements){
		int N = elements.size();
		Character[] elements_array = new Character[N];
		elements.toArray(elements_array);
		return permutation(elements_array);
	}
	
	public static List<List<String>> permutation(String[] elements) {
		// TODO Auto-generated method stub
		// Create the initial vector
		ICombinatoricsVector<String> initialVector = Factory.createVector(elements);
		// Create a simple combination generator to generate 3-combinations of the initial vector
		Generator<String> gen = Factory.createPermutationGenerator(initialVector);
		// Locate all possible combinations
		ArrayList<List<String>> permutations = new ArrayList<List<String>>();
		for (ICombinatoricsVector<String> combination : gen) {
			permutations.add(combination.getVector());
		}
		return permutations;
	}
	
	public static List<List<String>> multiPermutation(String[] elements, int K) {
		// TODO Auto-generated method stub
		// Create the initial vector
		ICombinatoricsVector<String> initialVector = Factory.createVector(elements);
		// Create a simple combination generator to generate 3-combinations of the initial vector
		Generator<String> gen = Factory.createSimpleCombinationGenerator(initialVector, K);
		// Locate all possible combinations
		ArrayList<List<String>> combinations = new ArrayList<List<String>>();
		for (ICombinatoricsVector<String> combination : gen) {
			combinations.add(combination.getVector());
		}
		ArrayList<List<String>> permutations = new ArrayList<List<String>>();
		for(List<String> combination : combinations) {
			List<List<String>> perm = permutation(
					combination.toArray(new String[combination.size()]));
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

	public static List<List<Integer>> permutation(int[] ele) {
		// TODO Auto-generated method stub
		return permutation(ArrayUtils.toObject(ele));
	}

	public static List<List<Integer>> permutation(int ele) {
		// TODO Auto-generated method stub
		int[] eles = new int[ele];
		for(int i=0; i<ele; i++)
			eles[i] = i;
		return permutation(eles);
	}
}
