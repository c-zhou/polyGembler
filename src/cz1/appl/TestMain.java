package cz1.appl;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.DiscreteDomain;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

import cz1.backup.Anchor1;
import cz1.ngs.model.Blast6Record;
import cz1.ngs.model.Sequence;
import cz1.ngs.tools.Redundas;
import cz1.util.Utils;

public class TestMain {

	public static void main(String[] args) {
		
		/**
		List<Sequence> seqs = Sequence.parseFastaFileAsList("C://Users//chenxi.zhou//Desktop//polyassembler//program_development//AP017304.1.fa");
		
		String seq_str = seqs.get(0).seq_str();
		
		String seq_str2 = seq_str.substring(20000);
		String seq_str3 = seq_str2.substring(0,30000)+Sequence.polyN(10000)+seq_str2.substring(31000);
		seq_str3 += Sequence.revCompSeq(seq_str.substring(0, 20000));
		try {
			BufferedWriter bw = Utils.getBufferedWriter("C://Users//chenxi.zhou//Desktop//polyassembler//program_development//AP017304.2.fa");
			bw.write(Sequence.formatOutput("AP017304.2", seq_str3));
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		**/
		/**
		RangeSet<Integer> ranges = TreeRangeSet.create();
		ranges.add(Range.closed(1, 4));
		ranges.add(Range.closed(3, 6));
		System.out.println(ranges.encloses(Range.closed(2, 5)));
		**/
		/**
		for(int i=1; i<=17; i++) {
			Anchor anchor = new Anchor();
			anchor.setParameters(new String[] {
					"-s",
					"C://Users//chenxi.zhou//Desktop//polyassembler//program_development//AP017304.2.fa",
					"-q",
					"C://Users//chenxi.zhou//Desktop//polyassembler//program_development//cp_genome//before_rr//before_rr_S"+i+".fasta",
					"-b",
					"C://Users//chenxi.zhou//Desktop//polyassembler//program_development//cp_genome//blat2//"+i+".txt",
					"-o",
					"C://Users//chenxi.zhou//Desktop//polyassembler//program_development//cp_genome//out2//"+i+"anchored",
			});
			anchor.run();
		}
		**/
		/**
		final List<Integer> list1 = new ArrayList<Integer>();
		list1.add(1);
		list1.add(2);
		final List<Integer> list2 = new ArrayList<Integer>();
		list2.add(3);
		list2.add(4);
		list1.addAll(list2);
		System.out.println(list1.size());
		System.out.println(list2.size());
		list2.clear();
		System.out.println(list1.size());
		System.out.println(list2.size());
		list2.add(5);
		System.out.println(list1.size());
		System.out.println(list2.size());
		**/
	}
}





