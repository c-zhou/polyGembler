package cz1.test;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
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
		
		DirectedWeightedPseudograph<Integer, DefaultWeightedEdge> subgraph = new DirectedWeightedPseudograph<>(DefaultWeightedEdge.class);
		subgraph.addVertex(1);
		subgraph.addVertex(2);
		subgraph.addVertex(3);
		subgraph.addVertex(4);
		subgraph.addVertex(5);
		subgraph.addVertex(6);
		subgraph.addVertex(7);
		subgraph.addVertex(8);
		subgraph.addVertex(9);
		subgraph.addVertex(10);
		
		subgraph.addEdge(1, 2);
		subgraph.addEdge(1, 3);
		subgraph.addEdge(2, 3);
		subgraph.addEdge(2, 4);
		subgraph.addEdge(3, 7);
		subgraph.addEdge(4, 5);
		subgraph.addEdge(4, 8);
		subgraph.addEdge(5, 2);
		subgraph.addEdge(5, 3);
		subgraph.addEdge(5, 6);
		subgraph.addEdge(5, 8);
		subgraph.addEdge(6, 3);
		subgraph.addEdge(6, 7);
		subgraph.addEdge(6, 8);
		subgraph.addEdge(6, 9);
		subgraph.addEdge(6, 10);
		subgraph.addEdge(7, 10);
		
		int source_segix, target_segix;
		int out = Integer.MAX_VALUE;
		for(int v : subgraph.vertexSet()) out = Math.min(v, out);
		
		Map<Integer, Integer> hierarch = new HashMap<>();
		
		int level = 0;
		hierarch.put(out, level);
		Set<Integer> upper_layer = new HashSet<>(), lower_layer = new HashSet<>();
		upper_layer.add(out);
		
		while(true) {
			lower_layer.clear();
			++level;
			for(int v : upper_layer) {
				for(DefaultWeightedEdge edge : subgraph.outgoingEdgesOf(v)) {
					target_segix = subgraph.getEdgeTarget(edge);
					if(!hierarch.containsKey(target_segix)) 
						lower_layer.add(target_segix);
				}
			}
			for(int z : lower_layer) 
				hierarch.put(z, level);
			upper_layer.clear();
			upper_layer.addAll(lower_layer);
			if(upper_layer.isEmpty()) break;
		}
		
		DirectedWeightedPseudograph<Integer, DefaultWeightedEdge> subgraph2 = 
				new DirectedWeightedPseudograph<>(DefaultWeightedEdge.class);
		
		for(int v : hierarch.keySet()) subgraph2.addVertex(v);
		for(DefaultWeightedEdge edge : subgraph.edgeSet()) {
			source_segix = subgraph.getEdgeSource(edge);
			target_segix = subgraph.getEdgeTarget(edge);
			if(hierarch.get(source_segix)<=hierarch.get(target_segix))
				subgraph2.addEdge(source_segix, target_segix);
		}
		
		subgraph = subgraph2;
		
		System.out.println("Graph reconstructed "+"#V "+subgraph.vertexSet().size()+", #E "+subgraph.edgeSet().size());
		
		
		
		
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
