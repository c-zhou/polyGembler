package cz1.hmm.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.TreeBidiMap;
import org.apache.commons.lang3.ArrayUtils;

import cz1.external.Queue;
import cz1.hmm.model.LinkageMap;
import cz1.util.IO;

public class NearestNeighbourJoining {
	
	public static void main2(String[] args) {
		// create the command line parser
		CommandLineParser parser = new PosixParser();

		// create the Options
		Options options = new Options();
		options.addOption( "f", "rf-file", true, "recombination frequency file." );
		options.addOption( "c", "contig-file", true, "contig file." );
		options.addOption( "n", "minimum-snp", true, "minimum number of SNPs to analyse." );
		options.addOption( "o", "output", true, "output file." );
		options.addOption( "t", "rf-thresh", true, "recombination frequency threshold." );
		
		String rfFile=null, contigFile=null, outputFile=null;
		int minimumSNP=0;
		double thresh=1;
		try {
			// parse the command line arguments
			CommandLine line = parser.parse( options, args );
			if( line.hasOption("f") ) {
				rfFile = line.getOptionValue('f');
			} else
				throw new RuntimeException("!!!");
			if( line.hasOption("c") ) {
				contigFile = line.getOptionValue('c');
			} else
				throw new RuntimeException("!!!");
			if( line.hasOption("n") ) 
				minimumSNP = Integer.parseInt(line.getOptionValue('n'));
			if( line.hasOption("o") ) {
				outputFile = line.getOptionValue('o');
			} else {
				outputFile = System.nanoTime()+".gmo";
			}
			if( line.hasOption("t") )
				thresh = Double.parseDouble(line.getOptionValue('t'));
		} catch( ParseException exp ) {
			throw new RuntimeException( "Unexpected exception:" + exp.getMessage() );
		}

		NearestNeighbourJoining nnj = new NearestNeighbourJoining();
		nnj.initialise(rfFile, contigFile, minimumSNP);
		nnj.construct(outputFile, thresh);
	}
	
	
	private void construct(String outputFile, double thresh) {
		// TODO Auto-generated method stub
		boolean[][] adj_matrix = new boolean[MaximumContigIndex+1][MaximumContigIndex+1];
		for(int i=0; i<=MaximumContigIndex; i++) {
			int k = -1;
			double d = Double.POSITIVE_INFINITY;
			
			outerloop:
				for(int j=0; j<=MaximumContigIndex; j++) {
					if(j==i) continue;
					if(DistanceMAT[i][j]==null) break outerloop;
					if(DistanceMAT[i][j].distance<d) {
						d = DistanceMAT[i][j].distance;
						k = j;
					}
				}
			
			if(d<=thresh) {
				adj_matrix[i][k] = true;
				adj_matrix[k][i] = true;
			}
		}
		
		int[] clans = traverse(adj_matrix);
		
		/**
		BufferedWriter bw;
		try {
			bw = IO.getBufferedWriter(outputFile);
			for(int i=0; i<clans.length; i++) {
				bw.write(clans[i]);
				bw.write("\n");
			}
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		*/
		
		int maxC = 0;
		for(int i=0; i<=MaximumContigIndex; i++) 
			if(clans[i]<maxC) maxC=clans[i];
		for(int i=0; i<=maxC; i++) {
			Set<Integer> tmpS = new HashSet<Integer>();
			for(int j=0; j<=MaximumContigIndex; j++)
				if(clans[j]==i) tmpS.add(j);
			Integer[] tmpI = new Integer[tmpS.size()];
			tmpS.toArray(tmpI);
			int[] clani = ArrayUtils.toPrimitive(tmpI);
			LinkageMap lg;
			int n = clani.length;
			int j1, j2, k1, k2;
			double[] distance4;
			double max_d = Double.NEGATIVE_INFINITY;
			
			if(n==1) lg = new LinkageMap(clani);
			else {
				double[][] distanceMat = new double[n*2][];
				for(int j=0; j<n; j++) {
					j1 = j*2;
					j2 = j*2+1;
					for(int k=j; k<n; k++) {
						k1 = k*2;
						k2 = k*2+1;
						if(j==k) {
							distanceMat[j1][k1] = distanceMat[j2][k2] = 0;
							distanceMat[j1][k2] = distanceMat[j2][k1] = Double.NEGATIVE_INFINITY;
						} else {
							distance4 = this.DistanceMAT[i][j].distance4;
							distanceMat[j1][k1] = distanceMat[j2][k2] = 0;
						}
					}
				}
				
			}
		}
	}

	private int[] traverse(boolean[][] adj_matrix) {
		// TODO Auto-generated method stub
		boolean[] visited = new boolean[MaximumContigIndex+1];
		int[] clans = new int[MaximumContigIndex+1];
		int cn = -1;
		for(int i=0; i<=MaximumContigIndex; i++) {
			if(visited[i]) continue;
			cn++;
			Queue<Integer> q = new Queue<Integer>();
			q.enqueue(i);
			while(!q.isEmpty()) {
				int n = q.dequeue();
				if(visited[n]) continue;
				clans[n] = cn;
				visited[n] = true;
				for(int j=0; j<=MaximumContigIndex; j++)
					if(adj_matrix[i][j]) q.enqueue(j);
			}
		}
		return clans;
	}

	private class LinkageGroupPair {
		private final double distance;
		private final double[] distance4;

		public LinkageGroupPair(double distance,
				double[] distance4) {
			this.distance = distance;
			this.distance4 = distance4;
		}
	}

	private final BidiMap<Short, String> PrimitiveContigIndexMap = 
			new TreeBidiMap<Short, String>();
	private final Map<Short, Double> PrimitiveContigSizeMap = new HashMap<Short, Double>();
	private static short MaximumContigIndex = -1;
	private LinkageGroupPair[][] DistanceMAT = null;
	
	private void initialise(String rfFile, String contigFile, int minimumSNP) {
		// TODO Auto-generated method stub
		try {
			String line;
			String[] s;
			BufferedReader br_contig = IO.getBufferedReader(contigFile);
			while( (line=br_contig.readLine())!=null ) {
				s = line.split("\\s+");
				if(Integer.parseInt(s[1])>=minimumSNP) {
					MaximumContigIndex++;
					PrimitiveContigIndexMap.put(MaximumContigIndex, s[0]);
					try {
						PrimitiveContigSizeMap.put(MaximumContigIndex, Double.parseDouble(s[2]));
					} catch (NumberFormatException e) {
						PrimitiveContigSizeMap.put(MaximumContigIndex, 0.0);
					}
				}
			}
			br_contig.close();

			DistanceMAT = new LinkageGroupPair[MaximumContigIndex+1][MaximumContigIndex+1];
			
			BufferedReader br_rf = IO.getBufferedReader(rfFile);
			double distance;
			short c1, c2;
			while( (line=br_rf.readLine())!=null ) {
				s = line.split("\\s+");
				if(!PrimitiveContigIndexMap.containsValue(s[6]) || 
						!PrimitiveContigIndexMap.containsValue(s[7]))
					continue;
				distance = Double.parseDouble(s[1]);
				c1 = PrimitiveContigIndexMap.getKey(s[6]);
				c2 = PrimitiveContigIndexMap.getKey(s[7]);

				DistanceMAT[c1][c2] = new LinkageGroupPair(
						distance,
						new double[]{
								Double.parseDouble(s[2]),
								Double.parseDouble(s[3]),
								Double.parseDouble(s[4]),
								Double.parseDouble(s[5])
						});
				DistanceMAT[c2][c1] = new LinkageGroupPair(
						distance,
						new double[]{
								Double.parseDouble(s[2]),
								Double.parseDouble(s[4]),
								Double.parseDouble(s[3]),
								Double.parseDouble(s[5])
						});
			}
			br_rf.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
