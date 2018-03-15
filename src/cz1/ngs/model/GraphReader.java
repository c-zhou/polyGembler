package cz1.ngs.model;

public class GraphReader {

	private enum Graph_type {de_bruijn, overlap};

	public DirectedWeightedOverlapPseudograph<String> read(String graph_file) {
		Graph_type graph_type = this.graphType(graph_file);
		
		switch (graph_type) {
		case de_bruijn:
			throw new RuntimeException("need to provide k-mer size for reading de bruijn assembly graph!!!");
		case overlap:
			return this.readOverlapGraph(graph_file); 
		default:
			throw new RuntimeException("errors occured reading the assembly graph file!!!");
		}
	}

	public DirectedWeightedOverlapPseudograph<String> read(String graph_file, int k) {
		Graph_type graph_type = this.graphType(graph_file);
		
		switch (graph_type) {
		case de_bruijn:
			return this.readDeBruijnGraph(graph_file);
		case overlap:
			return this.readOverlapGraph(graph_file);
		default:
			throw new RuntimeException("errors occured reading the assembly graph file!!!");
		}	
	}
	
	private DirectedWeightedOverlapPseudograph<String> readDeBruijnGraph(String graph_file) {
		// TODO Auto-generated method stub
		throw new RuntimeException("de bruijn assembly graph reader yet to implement!!!");
		// return null;
	}

	private DirectedWeightedOverlapPseudograph<String> readOverlapGraph(String graph_file) {
		// TODO Auto-generated method stub
		DirectedWeightedOverlapPseudograph<String> graph = 
				new DirectedWeightedOverlapPseudograph<String>(OverlapEdge.class);
		
		return graph;
	}
	
	private Graph_type graphType(String asm_graph) {
		// TODO Auto-generated method stub

		// TODO distinguish between the overlap and de-bruijn assembly graph
		return Graph_type.overlap;
	}
	
	public static void main(String[] args) {
		
	}
}
