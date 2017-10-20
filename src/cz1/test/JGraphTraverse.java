package cz1.test;

import java.util.HashSet;
import java.util.Set;

import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.jgrapht.traverse.ClosestFirstIterator;

public class JGraphTraverse {
	final static DirectedWeightedPseudograph<Character, DefaultWeightedEdge> assembly_graph = 
			new DirectedWeightedPseudograph<Character, DefaultWeightedEdge>(DefaultWeightedEdge.class);
	
	public static void main(String[] args) {
		for(int i=0; i!=5; i++) assembly_graph.addVertex((char)('A'+i));
		DefaultWeightedEdge e;
		e = assembly_graph.addEdge('A', 'B');
		assembly_graph.setEdgeWeight(e, 1);
		e = assembly_graph.addEdge('A', 'C');
		assembly_graph.setEdgeWeight(e, 2);
		e = assembly_graph.addEdge('C', 'D');
		assembly_graph.setEdgeWeight(e, 7);
		e = assembly_graph.addEdge('C', 'E');
		assembly_graph.setEdgeWeight(e, 3);
		e = assembly_graph.addEdge('D', 'B');
		assembly_graph.setEdgeWeight(e, 3);
		
		final Set<Character> visited = new HashSet<Character>();
		char v = 'A';
		int r = 8;
		ClosestFirstIterator<Character, DefaultWeightedEdge> radius_search = 
				new ClosestFirstIterator<Character, DefaultWeightedEdge>(
						assembly_graph, v, r);
		visited.clear();
		while(radius_search.hasNext()) 
			visited.add(radius_search.next());
		System.out.print("radius search from vtex "+v+", "+visited.size()+
				" vtex within "+r+" radius: ");
		for(Character vv : visited) System.out.print(vv+";");
		System.out.println();
	}
	
}
