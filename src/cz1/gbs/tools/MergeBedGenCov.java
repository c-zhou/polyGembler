package cz1.gbs.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.lang3.builder.EqualsBuilder;
import org.apache.commons.lang3.builder.HashCodeBuilder;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

public class MergeBedGenCov extends Executor {

	private String input_covs = null;
	private String output_cov = "./merged.cov.gz";
	private boolean rev_zero = false;
	private long maxCov = Long.MAX_VALUE;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--input-dir      Input sam/bam file.\n"
						+ " -0/--reverse-zero   Chr00 is the maximum.\n"
						+ " -x/--max-coverage   A position is omitted if the coverage is "
						+ "                     greater than this number (default no limit).\n"
						+ " -o/--prefix         Output file. \n\n");
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}

		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-i", "--input-dir", true);
			myArgsEngine.add("-0", "--reverse-zero", true);
			myArgsEngine.add("-x", "--max-coverage", true);
			myArgsEngine.add("-o", "--prefix", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-i")) {
			input_covs = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the folder of coverage files.");
		}

		if (myArgsEngine.getBoolean("-o")) {
			output_cov = myArgsEngine.getString("-o");
		}
		
		if (myArgsEngine.getBoolean("-0")) {
			rev_zero = true;
		}
		
		if (myArgsEngine.getBoolean("-x")) {
			maxCov = Long.parseLong(myArgsEngine.getString("-x"));
		}
	}

	private class RefPosObj implements Comparator<RefPosObj>, Comparable<RefPosObj> {
		private final String chr_id;
		private final long position;
		
		public RefPosObj(String chr_id, long position) {
			this.chr_id = chr_id;
			this.position = position;
		}
		
		@Override
	    public int hashCode() {
	        return new HashCodeBuilder(17, 31). // two randomly chosen prime numbers
	            // if deriving: appendSuper(super.hashCode()).
	            append(chr_id).
	            append(position).
	            toHashCode();
	    }

	    @Override
	    public boolean equals(Object obj) {
	       if (!(obj instanceof RefPosObj))
	            return false;
	        if (obj == this)
	            return true;

	        RefPosObj rhs = (RefPosObj) obj;
	        return new EqualsBuilder().
	            // if deriving: appendSuper(super.equals(obj)).
	            append(chr_id, rhs.chr_id).
	            append(position, rhs.position).
	            isEquals();
	    }

		@Override
		public int compare(RefPosObj obj, RefPosObj obj1) {
			// TODO Auto-generated method stub
			if(!obj.chr_id.equals("Chr00") &&
					!obj1.chr_id.equals("Chr00")) {
				return obj.chr_id.equals(obj1.chr_id) ? 
						(int)(obj.position-obj1.position) : 
							obj.chr_id.compareTo(obj1.chr_id);
			}
			if(obj.chr_id.equals("Chr00") &&
					!obj1.chr_id.equals("Chr00"))
				return rev_zero ? 1 : -1;
			if(!obj.chr_id.equals("Chr00") &&
					obj1.chr_id.equals("Chr00"))
				return rev_zero ? -1 : 1;
			return (int)(obj.position-obj1.position);
		}

		@Override
		public int compareTo(RefPosObj obj) {
			// TODO Auto-generated method stub
			return this.compare(this, obj);
		}
	}
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		File folder = new File(input_covs);
		File[] listOfFiles = folder.listFiles();
		int cov_n = listOfFiles.length;
		
		try {
			BufferedReader[] cov_br = new BufferedReader[cov_n];
			BufferedWriter out_cov = output_cov.endsWith(".gz") ? 
					Utils.getGZIPBufferedWriter(output_cov) :
						Utils.getBufferedWriter(output_cov) ;
			boolean[] hasMoreLine = new boolean[cov_n];
			List<TreeMap<RefPosObj, Integer>> cov_pool = new ArrayList<TreeMap<RefPosObj, Integer>>();
			
			String[] s;
			long a, b;
			int c;
			RefPosObj min_RefPosObj = null, tmp_RefPosObj;
			for (int i = 0; i < cov_n; i++) {
				cov_br[i] = new BufferedReader(new FileReader(listOfFiles[i]));
				hasMoreLine[i] = true;
				cov_pool.add(new TreeMap<RefPosObj, Integer>());
				s = cov_br[i].readLine().split("\\s+");
				a = Long.parseLong(s[1]);
				b = Long.parseLong(s[2]);
				c = Integer.parseInt(s[3]);
				tmp_RefPosObj = new RefPosObj(s[0], a);
				
				if(min_RefPosObj==null || 
						tmp_RefPosObj.compareTo(min_RefPosObj)<0)
					min_RefPosObj = tmp_RefPosObj;
				for(long j=a; j<b; j++)
					cov_pool.get(i).put(new RefPosObj(s[0], j), c);
			}
			
			int cov;
			TreeMap<RefPosObj, Integer> tmp_cov;
			String line;
			while(true) {
				cov = 0;
				for(int i=0; i<cov_n; i++) {
					tmp_cov = cov_pool.get(i);
					if(tmp_cov.isEmpty()) continue;
					if(min_RefPosObj.compareTo(tmp_cov.firstKey())==0) {
						cov += tmp_cov.get(min_RefPosObj);
						tmp_cov.remove(min_RefPosObj);
						if(tmp_cov.isEmpty()) {
							if(hasMoreLine[i]) {
								if( (line=cov_br[i].readLine())!=null ) {
									s = line.split("\\s+");
									a = Long.parseLong(s[1]);
									b = Long.parseLong(s[2]);
									c = Integer.parseInt(s[3]);
									for(long j=a; j<b; j++)
										tmp_cov.put(new RefPosObj(s[0], j), c);
								} else {
									hasMoreLine[i] = false;
								}
							}
						}
					}
				}
				if(cov>maxCov) 
					out_cov.write(min_RefPosObj.chr_id+"\t"+min_RefPosObj.position+"\t"+min_RefPosObj.position+"\n");
				min_RefPosObj = null;
				for(int i=0; i<cov_n; i++) {
					tmp_cov = cov_pool.get(i);
					if(tmp_cov.isEmpty()) continue;
					tmp_RefPosObj = tmp_cov.firstKey();
					if(min_RefPosObj==null || 
							tmp_RefPosObj.compareTo(min_RefPosObj)<0)
						min_RefPosObj = tmp_RefPosObj;
				}
				if(min_RefPosObj==null) break;
			}
			out_cov.close();
			for(int i=0; i<cov_n; i++) cov_br[i].close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
