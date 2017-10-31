package cz1.ngs.model;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import cz1.util.Utils;

public class Sequence implements Comparable<Sequence> {
	private final int seq_no; // sequence index 0-
	private final int seq_ln; // sequence length
	private final String seq_sn; // sequence name
	private final String seq_str; //sequence string
	
	public final static char[] nucleotide = new char[]{'A','C','G','T', 'N', 'a','c','g','t', 'n'};
	public final static Map<Character, Character> revcmp = new HashMap<Character, Character>();
	static {
		revcmp.put('A', 'T');
		revcmp.put('C', 'G');
		revcmp.put('G', 'C');
		revcmp.put('T', 'A');
		revcmp.put('N', 'N');
		revcmp.put('a', 't');
		revcmp.put('c', 'g');
		revcmp.put('g', 'c');
		revcmp.put('t', 'a');
		revcmp.put('n', 'n');
	}
	
	public Sequence (final int seq_no,
			final int seq_ln) {
		this.seq_no = seq_no;
		this.seq_ln = seq_ln;
		this.seq_sn = null;
		this.seq_str = null;
	}
	
	public Sequence (final String seq_sn,
			final String seq_str) {
		this.seq_no = -1;
		this.seq_ln = seq_str.length();
		this.seq_sn = seq_sn;
		this.seq_str = seq_str;
	}
	
	public Sequence (final String seq_str) {
		this.seq_no = -1;
		this.seq_ln = seq_str.length();
		this.seq_sn = null;
		this.seq_str = seq_str;
	}

	@Override
	public int compareTo(Sequence contig) {
		// TODO Auto-generated method stub
		return this.seq_ln-contig.seq_ln;
	}
	
	public int seq_no() {
		// sequence index 0-
		return this.seq_no;
	}
	
	public int seq_ln() {
		// sequence length
		return this.seq_ln;
	}
	
	public String seq_sn() {
		// sequence name
		return this.seq_sn;
	}
	
	public String seq_str() {
		//sequence string
		return this.seq_str;
	}

	public static List<Sequence> parseFastaFileAsList(String seq_fa) {
		// TODO Auto-generated method stub
		final List<Sequence> sequences = new ArrayList<Sequence>();
		try {
			BufferedReader br_fa = Utils.getBufferedReader(seq_fa);
			StringBuilder str_buf = new StringBuilder();
			String line = br_fa.readLine();
			String seq_sn = null;
			while(line!=null) {
				if(line.startsWith(">")) 
					seq_sn = line.replace(">","").split("\\s+")[0];

				str_buf.setLength(0);
				while( (line=br_fa.readLine())!=null && !line.startsWith(">") ) 
					str_buf.append(line);

				sequences.add(new Sequence(seq_sn, str_buf.toString()));
			}
			br_fa.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return sequences;
	}
	
	public static Map<String, Sequence> parseFastaFileAsMap(String seq_fa) {
		// TODO Auto-generated method stub
		final List<Sequence> sequences = parseFastaFileAsList(seq_fa);
		final Map<String, Sequence> seq_map = new HashMap<String, Sequence>();
		for(Sequence sequence : sequences) seq_map.put(sequence.seq_sn(), sequence);
		return seq_map;
	}
	
	public String revCompSeq() {
		// TODO Auto-generated method stub
		return revCompSeq(this.seq_str);
	}
	
	public static String revCompSeq(String seq) {
		StringBuilder seq_buf = new StringBuilder(seq);
		int n = seq.length();
		for(int i=0; i!=n; i++) 
			seq_buf.setCharAt(i, revcmp.get(seq_buf.charAt(i)));
		return seq_buf.reverse().toString();
	}

	public String formatOutput() {
		// TODO Auto-generated method stub
		return this.formatOutput(100);
	}
	
	public String formatOutput(int line_width) {
		// TODO Auto-generated method stub
		StringBuilder os = new StringBuilder();
		os.append(">");
		os.append(seq_sn);
		os.append("\n");
		Pattern p = Pattern.compile("(.{" + line_width + "})", Pattern.DOTALL);
	    Matcher m = p.matcher(this.seq_str);
	    os.append(m.replaceAll("$1" + "\n"));
	    if(os.charAt(os.length()-1)!='\n') os.append("\n");
	    return os.toString();
	}
	
	public static String formatOutput(String seq_sn, 
			String seq_str) {
		// TODO Auto-generated method stub
	    return formatOutput(seq_sn, seq_str, 100);
	}
	
	public static String formatOutput(String seq_sn, 
			String seq_str, 
			int line_width) {
		// TODO Auto-generated method stub
		StringBuilder os = new StringBuilder();
		os.append(">");
		os.append(seq_sn);
		os.append("\n");
		Pattern p = Pattern.compile("(.{" + line_width + "})", Pattern.DOTALL);
	    Matcher m = p.matcher(seq_str);
	    os.append(m.replaceAll("$1" + "\n"));
	    if(os.charAt(os.length()-1)!='\n') os.append("\n");
	    return os.toString();
	}

	public static Sequence polyN(int n) {
		// TODO Auto-generated method stub
		StringBuilder s = new StringBuilder();
		for(int i=0; i<n; i++)
			s.append("N");
		return new Sequence("GAP", s.toString());
	}
}
