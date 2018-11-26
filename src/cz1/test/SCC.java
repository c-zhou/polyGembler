package cz1.test;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Deque;
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

import com.google.common.collect.DiscreteDomain;
import com.google.common.collect.ImmutableRangeSet;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

import cz1.util.Utils;

public class SCC {

	private static int extension(RangeSet<Integer> cov1, RangeSet<Integer> cov2) {
		// TODO Auto-generated method stub
		RangeSet<Integer> ins = ImmutableRangeSet.copyOf(cov1).intersection(cov2);
		RangeSet<Integer> ext = ImmutableRangeSet.copyOf(ins).complement().intersection(cov2);
		int ln = 0;
		for(Range<Integer> r : ext.asRanges())
			ln += r.upperEndpoint()-r.lowerEndpoint();
		return ln;
	}
	
	public static void main(String[] args) throws NumberFormatException, IOException {

		RangeSet<Integer> cov1 = TreeRangeSet.create();
		RangeSet<Integer> cov2 = TreeRangeSet.create();
		cov1.add(Range.closed(1, 3).canonical(DiscreteDomain.integers()));
		cov1.add(Range.closed(6, 8).canonical(DiscreteDomain.integers()));
		cov2.add(Range.closed(6, 16).canonical(DiscreteDomain.integers()));
		cov1.add(Range.closed(18, 28).canonical(DiscreteDomain.integers()));
		System.out.println(extension(cov1, cov2));
		
		DirectedWeightedPseudograph<Integer, DefaultWeightedEdge> subgraph = new DirectedWeightedPseudograph<>(DefaultWeightedEdge.class);
		subgraph.addVertex(1);
		subgraph.addVertex(2);
		subgraph.addVertex(3);
		subgraph.addEdge(1,2);
		subgraph.addEdge(1,3);
		subgraph.addEdge(2,3);
		subgraph.removeAllVertices(new HashSet<Integer>(subgraph.vertexSet()));
		System.out.println(subgraph.vertexSet().size());
		System.out.println(subgraph.edgeSet().size());
		
		
		BufferedReader br1 = Utils.getBufferedReader("/shares/coin/c.zhou/8x8/data/CANU_out_NASPOT_10_O_mOVL050_v1.7_010/blast_out/kkk");
		String dagLine = br1.readLine();
		String[] s;
		int source, target;

		while( (dagLine=br1.readLine())!=null&&!dagLine.startsWith("#") ) {
			s = dagLine.split("\\s+");
			source = Integer.parseInt(s[0]);
			target = Integer.parseInt(s[1]);
			if(!subgraph.containsVertex(source))
				subgraph.addVertex(source);
			if(!subgraph.containsVertex(target))
				subgraph.addVertex(target);
			subgraph.addEdge(source, target);
		}

		List<Integer> outsV = new ArrayList<>();
		for(int v : subgraph.vertexSet()) 
			if(subgraph.incomingEdgesOf(v).isEmpty())
				outsV.add(v);
		if(outsV.isEmpty()) {
			for(int v : subgraph.vertexSet()) outsV.add(v);
		}
		Collections.sort(outsV);

		for(int i=0; i<100; i++) {
			DirectedWeightedPseudograph<Integer, DefaultWeightedEdge> subgraph2 = 
					new DirectedWeightedPseudograph<>(DefaultWeightedEdge.class);

			int out = outsV.get(i);

			Deque<DefaultWeightedEdge> stack = new ArrayDeque<DefaultWeightedEdge>();

			List<DefaultWeightedEdge> outs = new ArrayList<>(subgraph.outgoingEdgesOf(out));
			Collections.sort(outs, new Comparator<DefaultWeightedEdge>() {

				@Override
				public int compare(DefaultWeightedEdge arg0, DefaultWeightedEdge arg1) {
					// TODO Auto-generated method stub
					return subgraph.getEdgeTarget(arg1)-subgraph.getEdgeTarget(arg0);
				}

			});

			for(DefaultWeightedEdge e : outs)
				stack.push(e);

			Map<Integer, Boolean> visited = new HashMap<>();
			for(int v : subgraph.vertexSet()) visited.put(v, false);
			DefaultWeightedEdge edge;

			while(!stack.isEmpty()) {
				edge = stack.pop();
				source = subgraph.getEdgeSource(edge);
				target = subgraph.getEdgeTarget(edge);

				if(visited.get(target)) continue;
				if(!subgraph2.containsVertex(source)) subgraph2.addVertex(source);
				if(!subgraph2.containsVertex(target)) subgraph2.addVertex(target);
				DefaultWeightedEdge ea = subgraph2.addEdge(source, target);
				// System.out.println(ea.toString());

				List<DefaultWeightedEdge> edges = new ArrayList<>(subgraph.outgoingEdgesOf(target));
				Collections.sort(edges, new Comparator<DefaultWeightedEdge>() {

					@Override
					public int compare(DefaultWeightedEdge arg0, DefaultWeightedEdge arg1) {
						// TODO Auto-generated method stub
						return subgraph.getEdgeTarget(arg1)-subgraph.getEdgeTarget(arg0);
					}

				});

				for(DefaultWeightedEdge e : edges) stack.push(e);

				visited.put(source, true);
			}
			
			//CycleDetector<Integer, DefaultWeightedEdge> cycleDetector = new CycleDetector<>(subgraph2);
			//System.out.println(subgraph2.vertexSet().size()+" "+subgraph2.edgeSet().size()+" "+cycleDetector.detectCycles());
			System.out.println(subgraph2.vertexSet().size()+" "+subgraph2.edgeSet().size());
			
		}
		/****




		int source_segix, target_segix;


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
		 ***/
	}

}
