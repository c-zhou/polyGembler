package cz1.test;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.HashSet;
import java.util.Set;

import javax.swing.JFrame;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.TreeBidiMap;
import org.jgraph.JGraph;
import org.jgrapht.ext.JGraphModelAdapter;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.jgrapht.graph.ListenableDirectedWeightedGraph;

import cz1.ngs.model.TraceableVertex;

public class Traceable {
	public static void main(String[] args) {

		final BidiMap<String, Integer> bidi = new TreeBidiMap<String, Integer>();
		bidi.put("a", 1);
		bidi.put("a", 2);
		bidi.put("b", 2);
		System.out.println(bidi.size());
		
		
		
		
		final DirectedWeightedPseudograph<TraceableVertex<String>, DefaultWeightedEdge> razor = 
				new DirectedWeightedPseudograph<TraceableVertex<String>, DefaultWeightedEdge>(DefaultWeightedEdge.class);

		//JGraphModelAdapter<TraceableVertex<String>, DefaultWeightedEdge> jgAdapter = 
		//		new JGraphModelAdapter<TraceableVertex<String>, DefaultWeightedEdge>(razor);
		//JGraph jgraph = new JGraph(jgAdapter);
		
		TraceableVertex<String> a = new TraceableVertex<String>("a");
		TraceableVertex<String> b = new TraceableVertex<String>("b");
		TraceableVertex<String> c = new TraceableVertex<String>("c");
		TraceableVertex<String> d = new TraceableVertex<String>("d");
		
		razor.addVertex(a);
		razor.addVertex(b);
		razor.addVertex(c);
		razor.addVertex(d);

		razor.setEdgeWeight(razor.addEdge(a,a),1);
		razor.setEdgeWeight(razor.addEdge(a,b),1);
		razor.setEdgeWeight(razor.addEdge(a,c),1);
		razor.setEdgeWeight(razor.addEdge(a,d),1);
		razor.setEdgeWeight(razor.addEdge(b,a),1);
		razor.setEdgeWeight(razor.addEdge(b,b),1);
		razor.setEdgeWeight(razor.addEdge(b,c),1);
		razor.setEdgeWeight(razor.addEdge(b,d),1);
		razor.setEdgeWeight(razor.addEdge(c,a),1);
		razor.setEdgeWeight(razor.addEdge(c,b),1);
		razor.setEdgeWeight(razor.addEdge(c,c),1);
		razor.setEdgeWeight(razor.addEdge(c,d),1);
		razor.setEdgeWeight(razor.addEdge(d,a),1);
		razor.setEdgeWeight(razor.addEdge(d,b),1);
		razor.setEdgeWeight(razor.addEdge(d,c),1);
		razor.setEdgeWeight(razor.addEdge(d,d),1);


		//JFrame frame = new JFrame();
		//frame.getContentPane().add(jgraph);
		//frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		//frame.pack();
		//frame.setVisible(true);
	
		Set<TraceableVertex<String>> visited = new HashSet<TraceableVertex<String>>();
        Deque<TraceableVertex<String>> queue = new ArrayDeque<TraceableVertex<String>>();

        queue.push(a);
        visited.add(a);
        a.setScore(10);
        
        Set<DefaultWeightedEdge> out_edges;
		double max_score = Double.NEGATIVE_INFINITY, source_score, target_score, score;
		int source_ln;
		TraceableVertex<String> max_vertex = null;
		boolean isLeaf;
		
		TraceableVertex<String> source_vertex, target_vertex;
		double edge_weight;
		
        while(!queue.isEmpty()) {
        	source_vertex = queue.pop();
        	source_ln = 10;
        	source_score = source_vertex.getScore()-source_ln;
        	isLeaf = true;
        	out_edges = razor.outgoingEdgesOf(source_vertex);
        	for(DefaultWeightedEdge out : out_edges) {
        		target_vertex = razor.getEdgeTarget(out);
        		target_score = target_vertex.getScore();
        		edge_weight = razor.getEdgeWeight(out);
        		score = source_score+edge_weight;
        		
        		if( visited.contains(target_vertex) && 
        				(score<=target_score ||
        				isLoopback(razor, source_vertex, target_vertex)) ) 
        			continue;
        		
        		isLeaf = false;
        		target_vertex.setBackTrace(source_vertex);
        		target_vertex.setScore(score);
        		queue.push(target_vertex);
        		visited.add(target_vertex);
        	}
        	
        	if(isLeaf && source_vertex.getScore()>max_score) {
        		max_score = source_vertex.getScore();
        		max_vertex = source_vertex;
        	}
        }
	
	}
	
	
	private static boolean isLoopback(DirectedWeightedPseudograph<TraceableVertex<String>, DefaultWeightedEdge> graph,
			TraceableVertex<String> source,
			TraceableVertex<String> target) {
		// TODO Auto-generated method stub
		if(source.equals(target)) return true;
		while( (source = source.getBackTrace())!=null )
        	if(source.equals(target)) return true;
		return false;
	}
}
