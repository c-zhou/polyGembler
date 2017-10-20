package cz1.tenx.model;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import org.jgrapht.GraphPath;
import org.jgrapht.alg.shortestpath.DijkstraShortestPath;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import org.jgrapht.graph.WeightedPseudograph;
import org.jgrapht.traverse.ClosestFirstIterator;

import cz1.util.Utils;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader.ValidationStringency;

public class AssemblyGraph {

	final DirectedWeightedPseudograph<Integer, DefaultWeightedEdge> assembly_graph = 
			new DirectedWeightedPseudograph<Integer, DefaultWeightedEdge>(DefaultWeightedEdge.class);
	final Map<Integer, Integer> revc_map = new HashMap<Integer, Integer>();
	
	public void make(String fastg_file) {
		try {
			BufferedReader fastg_reader = Utils.getBufferedReader(fastg_file);
			String edge, edge2, c, c2;
			boolean d, d2;
			String[] s, s2, s3;
			int a, a2, l, l2;
			int line_count = 0;
			while( (edge=fastg_reader.readLine())!=null ) {
				edge2 = fastg_reader.readLine();
				
				c = edge.split(":")[0];
				c2 = edge2.split(":")[0];
				d = c.endsWith(",1");
				d2 = c2.endsWith(",1");
				if(d&&d2) {
					s = c.split(",");
					a = Integer.parseInt(s[0]);
					revc_map.put(a, a);
					assembly_graph.addVertex(a);
					s = c2.split(",");
					a = Integer.parseInt(s[0]);
					revc_map.put(a, a);
					assembly_graph.addVertex(a);
					continue;
				}
				s = c.split(",");
				a = Integer.parseInt(s[0]);
				if(d) {
					revc_map.put(a, a);
					revc_map.put(a+1, a);
					assembly_graph.addVertex(a);
				} else {
					revc_map.put(a, a+1);
					revc_map.put(a+1, a+1);
					assembly_graph.addVertex(a+1);
				}
				line_count += 2;
				if(line_count%1000000==0) System.out.println(line_count+" records loaded.");
			}
			fastg_reader.close();
		
			fastg_reader = Utils.getBufferedReader(fastg_file);
			line_count = 0;
			DefaultWeightedEdge e;
			while( (edge=fastg_reader.readLine())!=null ) {
				if(edge.endsWith(":")) continue;
				s = edge.split(":");
				s2 = s[0].split(",");
				a = revc_map.get(Integer.parseInt(s2[0]));
				l = Integer.parseInt(s2[1]);
				s2 = s[1].split(";");
				if(s[0].endsWith(",1")) { //outgoing edges
					for(int i=0; i!=s2.length; i++) {
						s3 = s2[i].split(",");
						a2 = revc_map.get(Integer.parseInt(s3[0]));
						l2 = Integer.parseInt(s3[1]);
						e = assembly_graph.addEdge(a, a2);
						assembly_graph.setEdgeWeight(e, (l+l2)/2.0);
					}
				} else { //incoming edges
					for(int i=0; i!=s2.length; i++) {
						s3 = s2[i].split(",");
						a2 = revc_map.get(Integer.parseInt(s3[0]));
						l2 = Integer.parseInt(s3[1]);
						e = assembly_graph.addEdge(a2, a);
						assembly_graph.setEdgeWeight(e, (l+l2)/2.0);
					}
				}
				line_count ++;
				if(line_count%1000000==0) System.out.println(line_count+" records loaded.");
			}
			fastg_reader.close();
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	public void map(String bam_file, String tenx_bc) {
		final SAMFileReader inputSam = new SAMFileReader(new File(bam_file));
		inputSam.setValidationStringency(ValidationStringency.SILENT);
		SAMRecordIterator iter=inputSam.iterator();
		SAMRecord record = null;
		Set<Integer> edges = new HashSet<Integer>();
		while( iter.hasNext() ) {
			record = iter.next();
			if( record.getReadName().startsWith(tenx_bc) )
				break;
		}
		while( record.getReadName().startsWith(tenx_bc) ) {
			if(!record.getReadUnmappedFlag())
				edges.add(Integer.parseInt(record.getReferenceName()));
			if(!iter.hasNext()) break;
			record = iter.next();
		}
		iter.close();
		inputSam.close();
		
		/*** dijkstra shortest path
		for(Integer e : edges) {
			for(Integer e2 : edges) {
				if(e!=e2) {
					GraphPath<Integer,DefaultWeightedEdge> shortest_path = 
							DijkstraShortestPath.findPathBetween(assembly_graph, e, e2);
                    StringBuilder os = new StringBuilder();
                    os.append(e);
                    os.append("->");
                    os.append(e2);
                    if(shortest_path==null) {
                            os.append(" path not existed.\n");
                            System.out.println(os.toString());
                            continue;
                    }
                    double l = shortest_path.getWeight();
                    List<Integer> v = shortest_path.getVertexList();
                    os.append("(");
                    os.append(l);
                    os.append("): ");
                    os.append(v.get(0));
                    for(Integer i : v) { os.append("-"); os.append(i); }
                    os.append("\n");
                    System.out.println(os.toString());
				}
			}
		}
		***/
		
		System.out.println(tenx_bc+" edge size, "+edges.size());
		Integer[] edge_arr = edges.toArray(new Integer[edges.size()]);
		Set<Integer> edge_copy = new HashSet<Integer>(edges);
		final Set<Integer> visited = new HashSet<Integer>();
		Integer vtex;
		for(int i=0; i!=edge_arr.length; i++) {
			if(!edge_copy.contains(edge_arr[i])) continue;
			edge_copy.remove(edge_arr[i]);
			ClosestFirstIterator<Integer, DefaultWeightedEdge> radius_search = 
					new ClosestFirstIterator<Integer, DefaultWeightedEdge>(assembly_graph,
							edge_arr[i], 50000);
			visited.clear();
			while(radius_search.hasNext()) {
				vtex = radius_search.next();
				if(edges.contains(vtex))  {
					visited.add(vtex);
					edge_copy.remove(vtex);
				}
			}
			System.out.print("radius search from vtex "+edge_arr[i]+", "+visited.size()+
					" vtex within 50Kb radius: ");
			for(Integer v : visited) System.out.print(v+";");
			System.out.println();
		}
	}
	
	public static void main(String[] args) {
		AssemblyGraph graph = new AssemblyGraph();
		//graph.make("C:\\Users\\chenxi.zhou\\Desktop\\10x_igv\\batatas_81mer.fastg.gz");
		graph.make(args[0]);
		DirectedWeightedPseudograph<Integer, DefaultWeightedEdge> assembly_graph = 
				graph.assembly_graph;
		System.out.println("edge set size, "+assembly_graph.edgeSet().size());
		System.out.println("vtex set size, "+assembly_graph.vertexSet().size());
		//graph.map("C:\\Users\\chenxi.zhou\\Desktop\\10x_igv\\barcoded.fastq2.tail.k81.sorted.bam",
		//		"BX:Z:AAACACCAGCTGCGAA-1");
		//graph.map(args[1], args[2]);
		
		Set<Integer> vtex_set = assembly_graph.vertexSet();
		int n = vtex_set.size();
		Integer[] vtex_arr = vtex_set.toArray(new Integer[n]);
		Random rnd = new Random();
		
		Integer v, v2;
		for(int i=0; i!=1000; i++) {
			v = vtex_arr[rnd.nextInt(n)];
			v2 = vtex_arr[rnd.nextInt(n)];
			if(v==v2) {
				i--;
				continue;
			}
			
			GraphPath<Integer,DefaultWeightedEdge> shortest_path = 
					DijkstraShortestPath.findPathBetween(assembly_graph, v, v2);
            StringBuilder os = new StringBuilder();
            os.append(v);
            os.append("->");
            os.append(v2);
            if(shortest_path==null) {
                    os.append(" path not existed.\n");
                    System.out.println(os.toString());
                    continue;
            }
            double l = shortest_path.getWeight();
            List<Integer> pv = shortest_path.getVertexList();
            os.append("(");
            os.append(l);
            os.append("): ");
            os.append(pv.get(0));
            for(Integer vv : pv) { os.append("-"); os.append(vv); }
            os.append("\n");
            System.out.println(os.toString());
		}
		
		final Set<Integer> visited = new HashSet<Integer>();
		for(int i=0; i!=1000; i++) {
			v = vtex_arr[rnd.nextInt(n)];
			ClosestFirstIterator<Integer, DefaultWeightedEdge> radius_search = 
					new ClosestFirstIterator<Integer, DefaultWeightedEdge>(
							assembly_graph, v, 25000);
			visited.clear();
			while(radius_search.hasNext()) 
				visited.add(radius_search.next());
			System.out.println("radius search from vtex "+v+", "+visited.size()+
					" vtex within 25Kb radius.");
		}
	}
}