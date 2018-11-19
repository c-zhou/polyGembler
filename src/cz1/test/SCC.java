package cz1.test;

import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.jgrapht.Graph;
import org.jgrapht.alg.CycleDetector;
import org.jgrapht.alg.KosarajuStrongConnectivityInspector;
import org.jgrapht.alg.interfaces.StrongConnectivityAlgorithm;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.DirectedWeightedPseudograph;

public class SCC {

	public static void main(String[] args) {
		final DirectedWeightedPseudograph<String, DefaultWeightedEdge> graph = new DirectedWeightedPseudograph<>(DefaultWeightedEdge.class);
		
		graph.addVertex("A");
		graph.addVertex("B");
		graph.addVertex("C");
		graph.addVertex("D");
		graph.addVertex("E");
		graph.addVertex("F");
		graph.addVertex("G");
		graph.addVertex("H");
		graph.addVertex("I");
		
		graph.addEdge("A", "H");
		graph.addEdge("C", "I");
		graph.addEdge("A", "B");
		graph.addEdge("B", "C");
		graph.addEdge("C", "E");
		graph.addEdge("E", "F");
		graph.addEdge("F", "G");
		graph.addEdge("G", "A");
		graph.addEdge("A", "D");
		graph.addEdge("D", "E");

		final StrongConnectivityAlgorithm<String, DefaultWeightedEdge> sscDetector = new KosarajuStrongConnectivityInspector<>(graph);
		List<Graph<String, DefaultWeightedEdge>> sccs = sscDetector.getStronglyConnectedComponents();
		for(final Graph<String, DefaultWeightedEdge> scc : sccs) {
			System.out.println(scc.vertexSet().size());
		}
	}

}
