package cz1.appl;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang3.StringUtils;

import cz1.util.IO;

public class GenomeAssemblyStatistics {

	private final static Set<Sequence> sequences = new HashSet<Sequence>();
	
	public static void main(String[] args) {
		read("C:\\Users\\chenxi.zhou\\Desktop\\putty\\Mx23Hm_81mer_bM.contig");
		System.out.println(count('A'));
		System.out.println(count('T'));
		System.out.println(count('G'));
		System.out.println(count('C'));
		System.out.println(count('N'));
		System.out.println(nucleotide());
		System.out.println(GCcontent());
		System.out.println(sequences.size());
		System.out.println(length());
		System.out.println(((double) length())/
				sequences.size());
		List<Sequence> list_seqs = new ArrayList<Sequence>(sequences);
		Collections.sort(list_seqs);
		System.out.println(list_seqs.get(list_seqs.size()-1).length);
		System.out.println(list_seqs.get(0).length);
		System.out.println(calcWeightedMedianStatistic(list_seqs, .5));
		
		double[][] intervals = new double[][] {
				new double[]{100, 200},
				new double[]{200, 300},
				new double[]{300, 400},
				new double[]{400, 500},
				new double[]{500, 1000},
				new double[]{1000, 2000},
				new double[]{2000, 3000},
				new double[]{3000, 4000},
				new double[]{4000, 5000},
				new double[]{5000, 10000},
				new double[]{10000, 20000},
				new double[]{20000, 30000},
				new double[]{30000, 40000},
				new double[]{40000, 50000},
				new double[]{50000, Double.POSITIVE_INFINITY}
		};
		
		for(int i=0; i<intervals.length; i++)
			interval_stats(intervals[i][0],intervals[i][1]);
	}
	
	private static void read(String fastaFile) {
		try{
			BufferedReader br = IO.getBufferedReader(fastaFile);
			String line = br.readLine();
			while( line!=null ) {
				if(line.startsWith(">")) {
					System.out.println(line);
					StringBuilder seq = new StringBuilder();
					while( (line=br.readLine())!=null 
							&& !line.startsWith(">"))
						seq.append(line);
					sequences.add(new Sequence(seq.toString()));
				}
			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private static void interval_stats(double lower_bound,
			double upper_bound) {
		//semi closed interval, lower_bound included, upper_bound excluded
		int n = 0, bases = 0;
		for(Sequence seq : sequences) 
			if(seq.length>=lower_bound
					&& seq.length<upper_bound) {
				n++;
				bases += seq.length;
			}
		System.out.println(n);
		System.out.println(bases);
		System.out.println(((double)bases)/n);
	}
	
	private static double GCcontent() {
		return ((double) (count('G')+count('C')))/
				nucleotide();
	}

	private static int nucleotide() {
		return count('A')+count('C')+
				count('G')+count('T');
	}
	
	private static int length() {
		int n = 0;
		for(Sequence seq : sequences) 
			n += seq.length;
		return n;
	}
	
	private static int count(char b) {
		int n = 0;
		for(Sequence seq : sequences) 
			n += seq.count(b);
		return n;
	}
	
    private static int calcWeightedMedianStatistic (
    		List<Sequence> contigSortedAscending, 
    		double p) {
        int L = 0;
        for(int i=0; i<contigSortedAscending.size(); i++) {
            L += contigSortedAscending.get(i).length;
        }

        int N = (int) Math.ceil(p*L);
        int C = 0, B = 0, R = 0;
        while(C<N) {
            R = N-C;
            C+=contigSortedAscending.get(B++).length;
        }
        B--;

        if(L%2==0 && R==0) {
            return (contigSortedAscending.get(B).length+
            		contigSortedAscending.get(B-1).length)/2;
        } else if (L%2!=0 && R==0) {
            return contigSortedAscending.get(B-1).length;
        } else {
            return contigSortedAscending.get(B).length;
        }
    }
	
	private static class Sequence implements Comparable<Sequence> {
		private final String sequence;
		private final int length;
		private final int A;
		private final int C;
		private final int G;
		private final int T;
		private final int N;
		
		public Sequence(String sequence) {
			this.sequence = sequence;
			this.length = sequence.length();
			this.A = count('A');
			this.C = count('C');
			this.G = count('G');
			this.T = count('T');
			this.N = count('N');
		}

		private int count(char b) {
			// TODO Auto-generated method stub
			return StringUtils.countMatches(this.sequence, b);
		}

		@Override
		public int compareTo(Sequence seq) {
			// TODO Auto-generated method stub
			return this.length-seq.length;
		}
	}
}
