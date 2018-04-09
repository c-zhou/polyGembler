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

import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.DualHashBidiMap;
import org.jgrapht.GraphPath;
import org.jgrapht.graph.DefaultWeightedEdge;
import org.jgrapht.graph.DirectedWeightedPseudograph;

import com.google.common.collect.DiscreteDomain;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

import cz1.ngs.model.BFSShortestPath;
import cz1.ngs.model.Sequence;
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
	}
}




