package cz1.model;

import java.util.Map;

public class GeneticMap {

	
	private final Scaffold root = 
			new Scaffold("root",null);  
	
	
	
	private class Scaffold {
		private Map<Scaffold, Double> distance;
		private Map<Scaffold, double[]> distance4;
		private String id_str;
		private String newick_str;
		private Scaffold root;
		
		public Scaffold(String id_str, 
				Scaffold root) {
			this.id_str = id_str;
			this.root = root;
		}
		
	}
}
