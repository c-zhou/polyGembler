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
import org.jgrapht.traverse.TopologicalOrderIterator;

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
		
		graph.addEdge("A", "B");
		graph.addEdge("A", "C");
		graph.addEdge("B", "C");
		graph.addEdge("B", "D");
		graph.addEdge("C", "D");
		graph.addEdge("C", "G");
		graph.addEdge("D", "E");
		graph.addEdge("F", "G");
		graph.addEdge("G", "D");
		graph.addEdge("G", "H");
		graph.addEdge("I", "E");
		graph.addEdge("I", "H");
		
		
		
		final StrongConnectivityAlgorithm<String, DefaultWeightedEdge> sscDetector = new KosarajuStrongConnectivityInspector<>(graph);
		List<Graph<String, DefaultWeightedEdge>> sccs = sscDetector.getStronglyConnectedComponents();
		for(final Graph<String, DefaultWeightedEdge> scc : sccs) {
			System.out.println(scc.vertexSet().size());
		}
		
		final TopologicalOrderIterator<String, DefaultWeightedEdge> topoIter = new TopologicalOrderIterator<>(graph);
		try {
			while(topoIter.hasNext()) {
				System.out.println(topoIter.next());
			}
		} catch (IllegalArgumentException e) {
			System.out.println("Graph is not a DAG");
		}

		System.out.println("DONE.");
	}

}
