package cz1.test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;
import org.apache.commons.math3.stat.inference.TestUtils;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.jgrapht.GraphPath;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.DirectedWeightedPseudograph;

import com.google.common.collect.DiscreteDomain;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

import cz1.ngs.model.BFSShortestPath;
import cz1.ngs.model.DirectedWeightedOverlapPseudograph;
import cz1.ngs.model.OverlapEdge;
import cz1.ngs.model.Sequence;
import cz1.util.Constants;
import cz1.util.Dirichlet;
import cz1.util.JohnsonTrotter;
import cz1.util.Utils;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class AssemblyGraph {

	private final static int kmer_size = 81;
	private final static DirectedWeightedPseudograph<Integer, DefaultWeightedEdge> assembly_graph = 
			new DirectedWeightedPseudograph<Integer, DefaultWeightedEdge>(DefaultWeightedEdge.class);
	private final static BidiMap<Integer, String> contig_coordinate = new DualHashBidiMap<Integer, String>();
	private static Map<String, Sequence> qry_seqs;
	
	private static void makeAssemblyGraph(String query_file) {
		// TODO Auto-generated method stub
		List<Sequence> qry_seqs = Sequence.parseFastaFileWithRevCmpAsList(query_file);
		Map<String, Set<Integer>> pref_mer = new HashMap<String, Set<Integer>>();
		Map<String, Set<Integer>> suff_mer = new HashMap<String, Set<Integer>>();
		String seq_str, kmer;
		int seq_ln;
		Sequence seq;
		for(int i=0; i<qry_seqs.size(); i++) {
			seq = qry_seqs.get(i);
			contig_coordinate.put(i, seq.seq_sn());
			assembly_graph.addVertex(i);
			
			seq_str = seq.seq_str();
			seq_ln = seq.seq_ln();
			kmer = seq_str.substring(0, kmer_size);
			if(!pref_mer.containsKey(kmer))
				pref_mer.put(kmer, new HashSet<Integer>());
			pref_mer.get(kmer).add(i);
			kmer = seq_str.substring(seq_ln-kmer_size, seq_ln);
			if(!suff_mer.containsKey(kmer))
				suff_mer.put(kmer, new HashSet<Integer>());
			suff_mer.get(kmer).add(i);
		}
		Set<String> pref_set = pref_mer.keySet();
		Set<String> suff_set = suff_mer.keySet();
		for(String mer : suff_set) 
			// from a suff_mer to a pref_mer
			if(pref_set.contains(mer)) 
				for(Integer i : suff_mer.get(mer)) 
					for(Integer j : pref_mer.get(mer)) 
						assembly_graph.setEdgeWeight(assembly_graph.addEdge(i, j), 1.0);
		
		System.out.println("Assembly graph "+ assembly_graph.vertexSet().size()+" vetices and "+
				assembly_graph.edgeSet().size()+" edges.");

		return;
	}

	private static void validateAssemblyGraph(String query_file, String bam_in) {
		final SamReaderFactory factory =
		          SamReaderFactory.makeDefault()
		              .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
		            		  SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
		              .validationStringency(ValidationStringency.SILENT);
		
		qry_seqs = Sequence.parseFastaFileWithRevCmpAsMap(query_file);
		
		final SamReader in1 = factory.open(new File(bam_in));
		SAMRecordIterator iter1 = in1.iterator();
		final List<SAMRecord> buffered_records = new ArrayList<SAMRecord>();
		SAMRecord tmp_record1 = iter1.hasNext() ? iter1.next() : null;
		String record_id = tmp_record1.getReadName();
		SAMRecord r1, r2;
		int i1, i2;
		StringBuilder oos = new StringBuilder();
		String source, sink;
		
		final BFSShortestPath<Integer, DefaultWeightedEdge> bfs = 
				new BFSShortestPath<Integer, DefaultWeightedEdge>(assembly_graph);
		
		while(tmp_record1!=null) {
			buffered_records.clear();
			buffered_records.add(tmp_record1);
			while( (tmp_record1=iter1.hasNext()?iter1.next():null)!=null && 
					tmp_record1.getReadName().equals(record_id)) {
				buffered_records.add(tmp_record1);
			}
			
			if(buffered_records.size()==1) {
				if(tmp_record1!=null) record_id = tmp_record1.getReadName();
				continue;
			}
			
			// sort records by head clipping size 
			buffered_records.sort(new Comparator<SAMRecord>() {

				@Override
				public int compare(SAMRecord r0, SAMRecord r1) {
					// TODO Auto-generated method stub
					return Integer.compare(head_clip(r0), head_clip(r1));
				}
				
				private int head_clip(SAMRecord r) {	
					CigarElement cigar = r.getReadNegativeStrandFlag() ? 
							r.getCigar().getLastCigarElement() : r.getCigar().getFirstCigarElement();
					return cigar.getOperator().isClipping() ? cigar.getLength() : 0;
				}
			});
			
			// calculate distance for adjacent alignment record pairs
			r1 = buffered_records.get(0);
			// for writing
			oos.setLength(0);
			oos.append(record_id);
			oos.append("[");
			oos.append(buffered_records.size());
			oos.append(", ");
			oos.append(readLength(r1));
			// find reference index
			source = r1.getReferenceName()+(r1.getReadNegativeStrandFlag()?"'":"");
			i1 = contig_coordinate.getKey(source);
			oos.append(", ");
			oos.append(source);
			oos.append("] ");
			
			for(int i=1; i<buffered_records.size(); i++) {
				r2 = buffered_records.get(i);
				sink = r2.getReferenceName()+(r2.getReadNegativeStrandFlag()?"'":"");
				i2 = contig_coordinate.getKey(sink);
				oos.append(", (");
				GraphPath<Integer, DefaultWeightedEdge> e = bfs.getPath( i1, i2 );
				oos.append(sink);
				oos.append(",");
				oos.append(e==null ? "Inf" : e.getLength());
				oos.append(",");
				oos.append(pathLength(e));
				oos.append(")");
				
				i1 = i2;
			}
			oos.append("\n");
			
			System.out.println(oos.toString());
			
			if(tmp_record1!=null) record_id = tmp_record1.getReadName();
		}
		
		try {
			iter1.close();
			in1.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
			
	private static int pathLength(GraphPath<Integer, DefaultWeightedEdge> path) {
		// TODO Auto-generated method stub
		if(path==null) return -1;
		int nL = 0;
		List<Integer> vertices = path.getVertexList();
		for(int k=1; k<vertices.size()-1; k++)
			nL += qry_seqs.get(contig_coordinate.get(vertices.get(k))).seq_ln();
		nL -= (vertices.size()-1)*kmer_size;
		return nL>0 ? nL : 0;
	}

	private static int readLength(SAMRecord record) {
		// TODO Auto-generated method stub
		int nL = record.getReadLength();
		if(!record.getReadUnmappedFlag()) {	
			Cigar cigar = record.getCigar();
			CigarElement cs = cigar.getFirstCigarElement();
			CigarElement ce = cigar.getLastCigarElement();

			nL += (cs.getOperator()==CigarOperator.H?cs.getLength():0)+
					(ce.getOperator()==CigarOperator.H?ce.getLength():0);
		}
		return nL;
	}
	
	private final static double d_max2 = 50d;
	private final static double k_dst2 = 10d;
	private final static double olap_min2 = 30;
	private final static int default_flank = 50;
	private final static int merK = 12;
	private final static double m_clip2 = 0.1;
	
	public static void main(String[] args) {
		long tic, toc;
		tic = System.nanoTime();
		//makeAssemblyGraph(args[0]);
		toc = System.nanoTime();
		System.out.println("Assembly graph construction elapsed time: "+(toc-tic)+"ns.");
	
		tic = System.nanoTime();
		//validateAssemblyGraph(args[0], args[1]);
		toc = System.nanoTime();
		System.out.println("Processing BAM file elapsed time: "+(toc-tic)+"ns.");
		
		RangeSet<Integer> rangeS = TreeRangeSet.create();
		rangeS.add(Range.closed(1, 2).canonical(DiscreteDomain.integers()));
		RangeSet<Integer> rangeS2 = TreeRangeSet.create();
		rangeS2.addAll(rangeS);
		rangeS.add(Range.closed(2, 3).canonical(DiscreteDomain.integers()));
		rangeS.add(Range.closed(0, 1).canonical(DiscreteDomain.integers()));
		rangeS2.add(Range.closed(11, 22).canonical(DiscreteDomain.integers()));
		System.out.println();
		
		System.out.println(String.format("%12s","aa"));
		
		System.out.println(Long.parseLong("11111111111111111111111111111111",2));
		
		
		long a = 1;
		a <<= 32;
		a += 2;
		int i1 = (int) a;
		int j1 = (int) (a>>32);
		System.out.println(i1+" "+j1);
		
		String source = "CCATTAGAATTCCAAATACAAAAAAAGAAAAGAAAAGTCTTTCCCTCAAATAAAACAAAGCAATTATGCACTGCAAAGAATTAGGTTTTAATTGAATTAC"+
"CCAAATACAAACAAATCTTACCAAAAATGTCCTTAAAAGTTTTTGTAAAAAAAAAAAAAAAAAAATTGAATGCAAAAAATATATTTTTAAAATTTGTTTA"+
"AAAACGGAAAGATGATGAAATTAGGCAACTAAAGCCTTAAGTCCAATTACTTACCAATGAGGCAATGATCTACAGTTAATTAAACCGATAAGCAATGAAC"+
"ATGAATAGTCAAGTTTAGGCTTCCCTAAGTTTGTCTTCTGAACTAGTTCAACGGCGGTTGCAGTGAGAGGTTCACCCGATTAATAACGAACACAATAAAA"+
"AAACAGTGATGAACACACAGATCTTGATAACTCAGTTTAGAGATTAAACTCCTACGTCTGGGGGCATCACATGAAATGACTGCCAGATTCTCCATCGATC"+
"TCAAAGGAAACGATATATTGAATCACAATGATGAACACGAAGGACTTCAAGAATGTATAGCTACAACTTAGCTAATCAAGACGAGCTCCTATTTTGAAGA"+
"ATCGATGTGAGATCTCCTGATCACCCAAATTGAGGCTCCCCTCAATTGCAACAGACTAGCGTCAGCTCTATATGATAATCAGCCCGAAAACAACAACGGA"+
"TGAGATGACTGACTTGACTTCAAATCT";
		String target = "ACTTCAAATCTTTATTCAATGGAGAATAATATAAAGTCACTCGAATTAGCTCTCTAAATATACTTAAGGTGTATCAAATCTTCTCTGTAAAAATCTTCATCATCAGGCCCCGAGATAATCCTTTATATATATGTCCCACATAGATCTAATCCTCTATATATATATATATATATATATATATATATATATATATATCCACATAAATCTTAACTGAATAACTGTCTTGAAAATCTGCTCGAACAATCTTGTCAGCTTCTTCACTTCAACAAGTCTAACTGCTCATCAACTTTGAACCATTATATTCTCAATTCCAGGCTCTATTTACCTTAGTACATAAATTGAACTATATGGCTCGAACCCACTCCCTTCCATGTGAGAGTGTAAATTGGGTGCCACTAAATCACAAGGTCGTTGATACATTAACATATTTCTAAACATTCATACTTTTTAAATCATATTAATATTATGCATGAGTGTTTATTGAATTATTGTATTGCATATAATAAGAAAAATGTCACTCTACACAATTGAAGCATCAAATACAAAAATTTGTACAATTGAACCATTTAACATTAAAAAGAAATGTGTGTAATTGACATTTTTATAGAACATATACTACCTTGTAAAGAATCATGAGTATTTAAAAAAAAAATATAGTAGTGAACAAACGGTTGGTCAATTACATAACCGACTAATACATTGAATAACAAACAAAAATTTAACATTTACAACCGACTGAAAATCAATCGCCAATTTTTATCAATTATGCCGGTAACAAAGAAAACAAGAGCATACTATTCAAAATTAATTAGCACTTGCCTAGCTCAGGTGGTATTCTTCGGTTTCTTCCCC";
		
		// so first we find a rough overlap with kmer hits 
		final Set<String> mer_bank = new HashSet<String>();
		final int source_ln = source.length();
		final int target_ln = target.length();
		final int mln = Math.min(source_ln, target_ln);
		for(int i=source_ln-mln; i<=source_ln-merK; i++)
			mer_bank.add(source.substring(i, i+merK));
		int hits = 0, merC = 0;
		for(int i=0; i<=mln-merK; i++) {
			if(mer_bank.contains(target.substring(i, i+merK))) {
				hits = i;
				++merC;
			}
			if(i-hits>d_max2) break;
		}
		
		
		int a1, a2, b1, b2, a1_flank, a2_flank, b1_flank, b2_flank, aLen;
		
		a1 = source_ln-hits;
		a2 = source_ln;
		b1 = 0;
		b2 = hits;
		
		a1_flank = Math.max(0, a1-default_flank);
		a2_flank = source_ln;
		b1_flank = 0;
		b2_flank = Math.min(target_ln, b2+default_flank);
		
		SequencePair<DNASequence, NucleotideCompound> seqPair = Constants.localPairMatcher(source.substring(a1_flank,a2_flank), target.substring(b1_flank,b2_flank));
		// what do we do if we don't find overlap?
		// ignore the link? maybe
		aLen = seqPair.getLength();
		//if(aLen<olap_min2) return;
		
		a1 = seqPair.getIndexInQueryAt(1);
		a2 = seqPair.getIndexInQueryAt(aLen);
		b1 = seqPair.getIndexInTargetAt(1);
		b2 = seqPair.getIndexInTargetAt(aLen);
		
		a1 += a1_flank;
		a2 += a1_flank;
		
		double clip = m_clip2*aLen;
		int a_clip = source_ln-a2,  b_clip = b1-1;
		//if(a_clip>clip||b_clip>clip) return;
		
		int a_match = a2-a1+1, b_match = b2-b1+1;
		int match = Math.min(a_match, b_match);
		int delete = b_clip+Math.max(b_match-match, 0);
		int insert = a_clip+Math.max(a_match-match, 0);
		
		System.out.println(seqPair.toString());
		System.out.println(aLen);
		System.out.println(a1);
		System.out.println(a2);
		System.out.println(b1);
		System.out.println(b2);
		System.out.println(delete+"D"+match+"M"+insert+"I");
		
		System.out.println(Constants.getOlapFromCigar("-100M"));
		System.out.println(Constants.getOlapFromCigar(Constants.cgRevCmp("0M")));
		
		SequencePair<DNASequence, NucleotideCompound> seqPair2 = Constants.localPairMatcher(source, target);
		System.out.println(seqPair2.toString());

		System.out.println(seqPair2.getNumIdenticals());
		System.out.println(seqPair2.getNumSimilars());
		
		
		System.out.println(TestUtils.chiSquareTest(new double[]{2d/6,4d/6}, new long[]{200,390}));
	
		
		int rangeUpperBound = 103;
		int rangeLowerBound = 1;
		int overlap = 1;
		int block = 10;
		int nB = (int) Math.ceil((double)(rangeUpperBound-rangeLowerBound+1-overlap)/(block-overlap));
		int[][] dataB = new int[nB][2];
		for(int i=0; i<nB; i++) {
			int from = rangeLowerBound+i*(block-overlap);
			dataB[i] = new int[]{from, from+block};
		}
		Utils.print(dataB);
		
		final List<String> aaa = new ArrayList<String>();
		aaa.add("aaa");
		aaa.add(null);
		aaa.add("bbb");
		aaa.add("ccc");
		System.out.println(aaa.size());
		
		System.out.println();
		int[][] jtPerm = JohnsonTrotter.perm(6);
		Utils.print(jtPerm);
		
		int[] p = new int[]{1,2,3,4,5,6};
		int tmp;
		for(int i=0; i<jtPerm.length; i++) {
			tmp = p[jtPerm[i][0]];
			p[jtPerm[i][0]] = p[jtPerm[i][1]];
			p[jtPerm[i][1]] = tmp;
			Utils.print(p);
		}
		
		double baf = 0.999;
		Dirichlet diri = new Dirichlet(new double[]{1-baf, baf}, 
				Constants._mu_theta_e);
		System.out.println(diri.sample()[0]);
	
		

		final DirectedWeightedOverlapPseudograph<String> gfa = 
				new DirectedWeightedOverlapPseudograph<String>(OverlapEdge.class);
		gfa.addVertex("a");
		gfa.addVertex("b");
		gfa.addVertex("c");
		gfa.addVertex("d");
		gfa.addVertex("e");
		gfa.addVertex("f");
		
		OverlapEdge edge;
		edge = gfa.addEdge("a","b");
		edge.setOlapF(1);
		edge = gfa.addEdge("a","c");
		edge.setOlapF(2);
		edge = gfa.addEdge("b","d");
		edge.setOlapF(3);
		edge = gfa.addEdge("c","e");
		edge.setOlapF(2);
		edge = gfa.addEdge("c","f");
		edge.setOlapF(3);
		edge = gfa.addEdge("a","f");
		edge.setOlapF(1);
		
		Map<String, Integer> len = new HashMap<String, Integer>();
		len.put("a", 10);
		len.put("b", 99);
		len.put("c", 8);
		len.put("d", 7);
		len.put("e", 10);
		len.put("f", 8);
		
		final Map<String, Map<String, Double>> distMat = new HashMap<String, Map<String, Double>>();
		Map.Entry<Double, Set<String>> nearest;
		Set<String> neighbors;
		final Set<String> visited = new HashSet<String>();
		double distance, d;
		String reached;
		final TreeMap<Double, Set<String>> visitor = new TreeMap<Double, Set<String>>();
		
		final double radius = 3.5;
		
		for(final String sourceV : gfa.vertexSet()) {
			visited.clear();
			
			final Set<String> root = new HashSet<String>();
			root.add(sourceV);
			// offer root node
			visitor.put(.0, root);
			final Map<String, Double> dist = new HashMap<String, Double>();
			
			while(!visitor.isEmpty()) {
				// poll this nearst node
				nearest = visitor.pollFirstEntry();
				distance = nearest.getKey();
				neighbors = nearest.getValue();
				
				for(final String neighbor : neighbors) {
					visited.add(neighbor);
					
					for(final OverlapEdge outEdge : gfa.outgoingEdgesOf(neighbor)) {
						reached = gfa.getEdgeTarget(outEdge);
						if(!visited.contains(reached)) { // not visited yet 
							d  = distance-outEdge.olapF();
							if(d<=radius) dist.put(reached, Math.max(.0, d));
							else continue;
							d += len.get(reached);
							if(visitor.containsKey(d)) {
								visitor.get(d).add(reached);
							} else {
								final Set<String> node = new HashSet<String>();
								node.add(reached);
								visitor.put(d, node);
							}
						}
					}
				}
			}
			distMat.put(sourceV, dist);
		}
		
		System.out.println();
	}
}




