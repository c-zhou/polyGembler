package cz1.test;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import cz1.ngs.model.Sequence;
import cz1.util.Utils;

public class ContigUtils {

	public static void main(String[] args) throws IOException {
		final List<Sequence> seqs = Sequence.parseFastaFileAsList(args[0]);
		Collections.sort(seqs, new Comparator<Sequence>() {

			@Override
			public int compare(Sequence arg0, Sequence arg1) {
				// TODO Auto-generated method stub
				return arg1.seq_ln()-arg0.seq_ln();
			}
			
		});
		
		BufferedWriter bw = Utils.getBufferedWriter(args[1]);
		for(final Sequence seq : seqs) bw.write(seq.formatOutput());
		bw.close();
	}
}
