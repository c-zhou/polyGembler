package cz1.util;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.log4j.Logger;
import org.paukov.combinatorics.Factory;
import org.paukov.combinatorics.Generator;
import org.paukov.combinatorics.ICombinatoricsVector;

import cz1.tenx.model.HiddenMarkovModel;

public class Combination {

	private final static Logger myLogger = Logger.getLogger(Combination.class);
	
	public static ArrayList<List<Character>> combination(Character[]  elements, int K){
		if(K > elements.length){
			myLogger.error("Invalid input, K > N");
			System.exit(1);
		}
		// Create the initial vector
		ICombinatoricsVector<Character> initialVector = Factory.createVector(elements);
		// Create a simple combination generator to generate 3-combinations of the initial vector
		Generator<Character> gen = Factory.createSimpleCombinationGenerator(initialVector, K);
		// Locate all possible combinations
		ArrayList<List<Character>> combinations = new ArrayList<List<Character>>();
		for (ICombinatoricsVector<Character> combination : gen) {
			combinations.add(combination.getVector());
		}
		return combinations;
	}
	
	public static ArrayList<List<Integer>> combination(Integer[]  elements, int K){
		if(K > elements.length){
			myLogger.error("Invalid input, K > N");
			System.exit(1);
		}
		// Create the initial vector
		ICombinatoricsVector<Integer> initialVector = Factory.createVector(elements);
		// Create a simple combination generator to generate 3-combinations of the initial vector
		Generator<Integer> gen = Factory.createSimpleCombinationGenerator(initialVector, K);
		// Locate all possible combinations
		ArrayList<List<Integer>> combinations = new ArrayList<List<Integer>>();
		for (ICombinatoricsVector<Integer> combination : gen) {
			combinations.add(combination.getVector());
		}
		return combinations;
	}
	
	public static ArrayList<List<String>> combination(String[]  elements, int K){
		if(K > elements.length){
			myLogger.error("Invalid input, K > N");
			System.exit(1);
		}
		// Create the initial vector
		ICombinatoricsVector<String> initialVector = Factory.createVector(elements);
		// Create a simple combination generator to generate 3-combinations of the initial vector
		Generator<String> gen = Factory.createSimpleCombinationGenerator(initialVector, K);
		// Locate all possible combinations
		ArrayList<List<String>> combinations = new ArrayList<List<String>>();
		for (ICombinatoricsVector<String> combination : gen) {
			combinations.add(combination.getVector());
		}
		return combinations;
	}
	
	public static ArrayList<List<Character>> multiCombination(Character[]  elements, int K){
		// Create the initial vector
		ICombinatoricsVector<Character> initialVector = Factory.createVector(elements);
		// Create a simple combination generator to generate 3-combinations of the initial vector
		Generator<Character> gen = Factory.createMultiCombinationGenerator(initialVector, K);
		// Locate all possible combinations
		ArrayList<List<Character>> combinations = new ArrayList<List<Character>>();
		for (ICombinatoricsVector<Character> combination : gen) {
			combinations.add(combination.getVector());
		}
		return combinations;
	}
	
	public static ArrayList<List<String>> multiCombination(String[]  elements, int K){
		// Create the initial vector
		ICombinatoricsVector<String> initialVector = Factory.createVector(elements);
		// Create a simple combination generator to generate 3-combinations of the initial vector
		Generator<String> gen = Factory.createMultiCombinationGenerator(initialVector, K);
		// Locate all possible combinations
		ArrayList<List<String>> combinations = new ArrayList<List<String>>();
		for (ICombinatoricsVector<String> combination : gen) {
			combinations.add(combination.getVector());
		}
		return combinations;
	}
	
	public static ArrayList<List<Character>> combination(List<Character>  elements, int K){
		int N = elements.size();
		if(K > N){
			myLogger.error("Invalid input, K > N");
			System.exit(1);
		}
		Character[] elements_array = new Character[N];
		elements.toArray(elements_array);
		return combination(elements_array, K);
	}

	public static ArrayList<List<Character>> combination(Set<Character>  elements, int K){
		int N = elements.size();
		if(K > N){
			myLogger.error("Invalid input, K > N");
			System.exit(1);
		}
		Character[] elements_array = new Character[N];
		elements.toArray(elements_array);
		return combination(elements_array, K);
	}
	
	public static ArrayList<List<Character>> multiCombination(List<Character>  elements, int K){
		return multiCombination(elements.toArray(new Character[elements.size()]), K);
	}
	
	public static ArrayList<List<Character>> multiCombination(char[] elements, int K){
		return multiCombination(ArrayUtils.toObject(elements), K);
	}
	
	public static ArrayList<List<Character>> multiCombination(Set<Character>  elements, int K){
		return multiCombination(elements.toArray(new Character[elements.size()]), K);
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

	public static List<List<Integer>> combination(int[] ele, int k) {
		// TODO Auto-generated method stub
		return combination(ArrayUtils.toObject(ele), k);
	}

	public static List<List<Integer>> combination(int ele, int k) {
		// TODO Auto-generated method stub
		int[] eles = new int[ele];
		for(int i=0; i<ele; i++)
			eles[i] = i;
		return combination(eles, k);
	}
}
