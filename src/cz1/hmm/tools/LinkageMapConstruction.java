package cz1.hmm.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.collections4.BidiMap;
import org.apache.commons.collections4.bidimap.TreeBidiMap;
import org.apache.commons.math3.stat.StatUtils;

import cz1.external.RedBlackBST;
import cz1.util.IO;

public class LinkageMapConstruction {

	public static void main(String[] args) {
		// create the command line parser
		CommandLineParser parser = new PosixParser();

		// create the Options
		Options options = new Options();
		options.addOption( "f", "rf-file", true, "recombination frequency file." );
		options.addOption( "c", "contig-file", true, "contig file." );
		options.addOption( "n", "minimum-snp", true, "minimum number of SNPs to analyse." );
		options.addOption( "o", "output", true, "output file." );
		options.addOption( "F", "fold", false, "construted from folded scaffolds");
		options.addOption( "t", "rf-thresh", true, "recombination frequency threshold." );
		
		String rfFile=null, contigFile=null, outputFile=null;
		boolean fold = false;
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
			if( line.hasOption("F") )
				fold = true;
		} catch( ParseException exp ) {
			throw new RuntimeException( "Unexpected exception:" + exp.getMessage() );
		}

		LinkageMapConstruction lmc = new LinkageMapConstruction();
		lmc.initialise(rfFile, contigFile, minimumSNP, fold);
		lmc.construct(outputFile, thresh);
	}

	private class LinkageGroupPair {
		private final short[] id;
		private final double[] distance4;

		public LinkageGroupPair(short[] id, 
				double[] distance4) {
			this.id = id;
			this.distance4 = distance4;
		}
	}

	private class NodeBST {
		private final Set<LinkageGroupPair> queue = 
				new HashSet<LinkageGroupPair>();

		public NodeBST(LinkageGroupPair lp1) {
			this.add(lp1);
		}

		public void add(LinkageGroupPair item) {
			this.queue.add(item);
		}

		public void remove(LinkageGroupPair item) {
			this.queue.remove(item);
		}

		public LinkageGroupPair retrieve() {
			LinkageGroupPair item = 
					this.queue.iterator().next();
			this.queue.remove(item);
			return item;
		}

		public void remove(short contig) {
			Set<LinkageGroupPair> rm = new HashSet<LinkageGroupPair>();
			for(LinkageGroupPair lp1 : this.queue )
				if( lp1.id[0]==contig ||
				lp1.id[1]==contig )
					rm.add(lp1);
			this.queue.removeAll(rm);
		}

		public int size() {
			return this.queue.size();
		}
	}

	private final RedBlackBST<Double, NodeBST> DistanceMAT = 
			new RedBlackBST<Double, NodeBST>();
	private final BidiMap<Short, String> PrimitiveContigIndexMap = 
			new TreeBidiMap<Short, String>();
	private final Map<Short, Double> PrimitiveContigSizeMap = new HashMap<Short, Double>();

	private static short MaximumContigIndex = -1;
	private final BidiMap<Short, String> ContigIndexMap = new TreeBidiMap<Short, String>();
	private final Map<Integer, LinkageGroupPair> ContigPairMap = 
			new HashMap<Integer, LinkageGroupPair>();
	private final Map<Short, Set<Double>> ContigBSTMap = new HashMap<Short, Set<Double>>();

	private void initialise(String rfFile, String contigFile, int minimumSNP, boolean fold) {
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

			BufferedReader br_rf = IO.getBufferedReader(rfFile);
			double key;
			short c1, c2;
			int int_key;
			while( (line=br_rf.readLine())!=null ) {
				s = line.split("\\s+");
				if(!PrimitiveContigIndexMap.containsValue(s[6]))
					if(fold && !primitiveContig(s[6])) {
						MaximumContigIndex++;
						PrimitiveContigIndexMap.put(MaximumContigIndex, s[6]);
					}
					else continue;
					
				if(!PrimitiveContigIndexMap.containsValue(s[7]))
					if(fold && !primitiveContig(s[7])) {
						MaximumContigIndex++;
						PrimitiveContigIndexMap.put(MaximumContigIndex, s[7]);
					}
					else continue;
				key = Double.parseDouble(s[1]);
				c1 = PrimitiveContigIndexMap.getKey(s[6]);
				c2 = PrimitiveContigIndexMap.getKey(s[7]);

				ContigIndexMap.put(c1, ""+c1);
				ContigIndexMap.put(c2, ""+c2);
				LinkageGroupPair pair = new LinkageGroupPair(
						new short[]{c1, c2},
						new double[]{
								Double.parseDouble(s[2]),
								Double.parseDouble(s[3]),
								Double.parseDouble(s[4]),
								Double.parseDouble(s[5])
						});
				int_key = (((0<<16)+c1)<<16)+c2;

				ContigPairMap.put(int_key, pair);
				NodeBST node = DistanceMAT.get(key);
				if(node==null) {
					DistanceMAT.put(key, new NodeBST(pair));
				} else {
					node.add(pair);
				}
				for(int i=6; i<=7; i++) {
					Short key_i = PrimitiveContigIndexMap.getKey(s[i]);
					if( !ContigBSTMap.containsKey(key_i) )
						ContigBSTMap.put(key_i, new HashSet<Double>());
					ContigBSTMap.get(key_i).add(key);
				}
			}
			br_rf.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private boolean primitiveContig(String scaffold) {
		// TODO Auto-generated method stub
		return scaffold.indexOf(':')>0;
	}

	private class Cluster {
		private final Set<Integer> group;
		private final double distance;

		public Cluster(Set group, double distance) {
			this.group = group;
			this.distance = distance;
		}
	}
	private final List<Cluster> MetaCluster = new ArrayList<Cluster>();
	private BufferedWriter oos = null;

	private void construct(String outputFile, double thresh) {
		try {
			oos = IO.getBufferedWriter(outputFile);
			int step = 0;
			oos(step++, -1);
			while(this.DistanceMAT.min()<=thresh && 
					this.DistanceMAT.size()>1) {
				double key = DistanceMAT.min();
				NodeBST val = DistanceMAT.get(key);
				LinkageGroupPair pair = val.retrieve();
				if(val.size()==0) DistanceMAT.delete(key);
				Short c1 = pair.id[0],
						c2 = pair.id[1];
				double[] distance4 = pair.distance4;
				int concat = 0;
				while( distance4[concat]!=key ) concat++;

				String str1 = ContigIndexMap.get(c1),
						str2 = ContigIndexMap.get(c2);
				if(str1==null || str2==null) 
					System.out.println();

				Short new_c = ++MaximumContigIndex;
				ContigBSTMap.put(new_c, new HashSet<Double>());
				int[] rJ = null;
				String new_str = null;
				switch(concat) {
				case 0:
					rJ = new int[]{2,3,2,3};
					new_str = strReverse(str1)+"-"+str2;
					break;
				case 1:
					rJ = new int[]{2,3,0,1};
					new_str = strReverse(str1)+"-"+strReverse(str2);
					break;
				case 2:
					rJ = new int[]{0,1,2,3};
					new_str = str1+"-"+str2;
					break;
				case 3:
					rJ = new int[]{0,1,0,1};
					new_str = str1+"-"+strReverse(str2);
					break;
				default:
					throw new RuntimeException("!!!");
				}

				//System.err.println("step: "+new_str);
				ContigIndexMap.remove(c1);
				ContigIndexMap.remove(c2);
				trimmingBST(c1);
				trimmingBST(c2);

				Set<Short> contigs = ContigIndexMap.keySet();
				int key_cc;
				double[] distance41, distance42;
				LinkageGroupPair lp1;
				for(Short c : contigs) {
					distance41 = distance4(c1, c);
					distance42 = distance4(c2, c);
					double[] new_distance4 = new double[4];
					new_distance4[0] = distance41[rJ[0]];
					new_distance4[1] = distance41[rJ[1]];
					new_distance4[2] = distance42[rJ[2]];
					new_distance4[3] = distance42[rJ[3]];
					short[] new_id = new short[]{new_c, c};
					lp1 = new LinkageGroupPair(new_id, new_distance4);
					key_cc = (((0<<16)+new_c)<<16)+c;
					ContigPairMap.put(key_cc, lp1);
					double min_d = StatUtils.min(new_distance4);
					if(this.DistanceMAT.contains(min_d)) 
						this.DistanceMAT.get(min_d).add(lp1);
					else
						this.DistanceMAT.put(min_d, new NodeBST(lp1));
					ContigBSTMap.get(new_c).add(min_d);
					ContigBSTMap.get(c).add(min_d);
				}
				ContigIndexMap.put(new_c, new_str);
				oos(step++, key);
			}
			oos.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void oos(int step, double key) {
		// TODO Auto-generated method stub
		try {
			oos.write("step "+step+++"@["+key+"]:");
			for(Short key_out : ContigIndexMap.keySet()) {
				oos.write(" ");
				String str_out = ContigIndexMap.get(key_out);
				String[] s = str_out.split("-");
				oos.write(PrimitiveContigIndexMap.get(Short.parseShort(s[0])));
				for(int i=1; i<s.length; i++) {
					oos.write("-");
					oos.write(PrimitiveContigIndexMap.get(Short.parseShort(s[i])));
				}
			}
			oos.write("\n");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}


	private void check(short c1) {
		// TODO Auto-generated method stub
		for(Double key : this.DistanceMAT.keys()) {
			Set<LinkageGroupPair> lp1s = this.DistanceMAT.get(key).queue;
			for(LinkageGroupPair lp1 : lp1s) {
				short[] cs = lp1.id;
				if(c1==cs[0] || c1==cs[1]) {
					throw new RuntimeException(key+" "+c1+" !!!");
				}
			}
		}
	}

	private void trimmingBST(Short c1) {
		// TODO Auto-generated method stub
		Set<Double> ds = ContigBSTMap.get(c1);
		for(Double d : ds) {
			NodeBST bst = this.DistanceMAT.get(d);
			if(bst==null) continue;
			bst.remove(c1);
			if(bst.size()==0) this.DistanceMAT.delete(d);
		}
		ContigBSTMap.remove(c1);
	}

	private double[] distance4(Short c2, Short c) {
		// TODO Auto-generated method stub		
		LinkageGroupPair lp1 = 
				ContigPairMap.get((((0<<16)+c2)<<16)+c);
		double[] distance4;
		if(lp1==null) { 
			lp1 = ContigPairMap.get((((0<<16)+c)<<16)+c2);
			distance4 = new double[4];
			System.arraycopy(lp1.distance4, 0, distance4, 0, 4);
			double tmp = distance4[1];
			distance4[1] = distance4[2];
			distance4[2] = tmp;
		} else {
			distance4 = lp1.distance4;
		}
		return distance4;
	}

	private String strReverse(String str) {
		// TODO Auto-generated method stub
		return strReverse(str, "-");
	}

	private String strReverse(String str, String delimeter) {
		// TODO Auto-generated method stub
		String[] s = str.split(delimeter);
		StringBuilder os = new StringBuilder();
		for(int i=s.length-1; i>0; i--) {
			os.append(s[i]);
			os.append(delimeter);
		}
		os.append(s[0]);
		return os.toString();
	}
}












