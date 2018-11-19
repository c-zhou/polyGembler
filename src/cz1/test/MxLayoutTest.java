package cz1.test;

import java.awt.Dimension;
import java.util.Random;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.jgrapht.ListenableGraph;
import org.jgrapht.ext.JGraphXAdapter;
import org.jgrapht.graph.DefaultListenableGraph;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.renjin.repackaged.guava.collect.Range;

import com.mxgraph.layout.mxCircleLayout;
import com.mxgraph.layout.mxCompactTreeLayout;
import com.mxgraph.layout.mxEdgeLabelLayout;
import com.mxgraph.layout.mxFastOrganicLayout;
import com.mxgraph.layout.mxIGraphLayout;
import com.mxgraph.swing.mxGraphComponent;

public class MxLayoutTest {
	
	private static final Dimension DEFAULT_SIZE = new Dimension(800, 800);
	
    private static void createAndShowGui() {
        JFrame frame = new JFrame("DemoGraph");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        
        ListenableGraph<String, DefaultWeightedEdge> g = buildGraph();
        JGraphXAdapter<String, DefaultWeightedEdge> graphAdapter = 
                new JGraphXAdapter<String, DefaultWeightedEdge>(g);
        
        mxIGraphLayout layout = new mxCircleLayout(graphAdapter);
        //mxIGraphLayout layout = new mxFastOrganicLayout(graphAdapter);
        
        layout.execute(graphAdapter.getDefaultParent());
        
        frame.add(new mxGraphComponent(graphAdapter));
        frame.setPreferredSize(DEFAULT_SIZE);
        frame.pack();
        frame.setLocationByPlatform(true);
        frame.setVisible(true);
    }

    public static void main(String[] args) {
    	Range<Integer> a = Range.closed(1, 2);
    	Range<Integer> b = Range.closed(4, 4);
    	System.out.println(a.isConnected(b));
    	Range<Integer> c = a.intersection(b);
    	System.out.println(c.lowerEndpoint());
    	System.out.println(c.upperEndpoint());
    	
        SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                createAndShowGui();
            }
        });
    }

    public static ListenableGraph<String, DefaultWeightedEdge> buildGraph() {
        ListenableGraph<String, DefaultWeightedEdge> g = 
            new DefaultListenableGraph<>(new DirectedWeightedPseudograph<>(DefaultWeightedEdge.class));

        for(int i=0; i<10; i++)
        	g.addVertex("x"+i);
        
        Random r = new Random();
        for(int i=0; i<10; i++)
        	for(int j=0; j<10; j++)
        		if(r.nextDouble()<0.1)
        			g.addEdge("x"+i, "x"+j);
        
        return g;
    }
}