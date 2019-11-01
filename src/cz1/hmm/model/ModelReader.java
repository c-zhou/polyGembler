package cz1.hmm.model;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.apache.commons.lang3.ArrayUtils;

import cz1.util.Utils;

public class ModelReader { 
	private final ZipFile in;
	private InputStream is;
	private BufferedReader br;

	public ModelReader(String in) {
		// TODO Auto-generated method stub
		try {
			this.in = new ZipFile(in);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			throw new RuntimeException(e.getMessage());
		}
	}

	private void setEntryReader(String target) {
		// TODO Auto-generated method stub
		try {
			String tfile = null;
			switch(target.toLowerCase()) {
			case "haplotype":
				tfile = "haplotype.txt";
				break;
			case "genotype":
				tfile = "genotype.txt";
				break;
			case "dosage":
				tfile = "dosage.txt";
				break;
			case "runinfo":
				tfile = "runinfo.txt";
				break;
			case "emission":
				tfile = "emission.txt";
				break;
			case "transition":
				tfile = "transition.txt";
				break;
			case "snp":
				tfile = "snp.txt";
				break;
			}
			if(tfile==null) throw new RuntimeException("!!!");
			ZipEntry entry;
			if( (entry=in.getEntry(tfile))==null) {
				throw new RuntimeException("!!!");
			}
			is = in.getInputStream(entry);
			br = Utils.getBufferedReader(is);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public void close() {
		// TODO Auto-generated method stub
		try {
			in.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void closeReader() {
		// TODO Auto-generated method stub
		try {
			br.close();
			is.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public double[][] getEmissionProbs(int ploidy) {
		// TODO Auto-generated method stub
		try {
			setEntryReader("emission");
			String line;
			List<String> lines = new ArrayList<>();
			while((line=br.readLine())!=null) {
				lines.add(line);
			}
			closeReader();
			
			int M = lines.size(), P = ploidy*2;
			double[][] emiss = new double[M][P];
			Pattern p = Pattern.compile("\\{(.*?)\\}");
			Matcher m;
			String[] prob_str;
			for(int i=0; i<M; i++) {
				line = lines.get(i);
				m = p.matcher(line);
				for(int j=0; j<P; j++) {
					m.find();
					prob_str = m.group(1).split(",|;");
					emiss[i][j] = Double.parseDouble(prob_str[1]);
				}
			}
			return emiss;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		throw new RuntimeException("!!!");
	}

	public int getMarkerNo() {
		// TODO Auto-generated method stub
		try {
			setEntryReader("runinfo");
			String line;
			while(!(line=br.readLine()).startsWith("##marker")) {}
			int no = Integer.parseInt(line.split("\\s+")[1]);
			closeReader();
			return no;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		throw new RuntimeException("!!!");
	}

	public int getSampleNo() {
		// TODO Auto-generated method stub
		try {
			setEntryReader("runinfo");
			String line;
			while(!(line=br.readLine()).startsWith("##sample")) {}
			int no = Integer.parseInt(line.split("\\s+")[1]);
			closeReader();
			return no;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		throw new RuntimeException("!!!");
	}

	public int getPloidy() {
		// TODO Auto-generated method stub
		try {
			setEntryReader("runinfo");
			String line;
			while(!(line=br.readLine()).startsWith("##ploidy")) {}
			int no = Integer.parseInt(line.split("\\s+")[1]);
			closeReader();
			return no;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		throw new RuntimeException("!!!");
	}

	public String[] getParents() {
		// TODO Auto-generated method stub
		try {
			setEntryReader("runinfo");
			String line;
			while(!(line=br.readLine()).startsWith("##parents")) {}
			String[] s = line.trim().split("\\s+");
			String[] parents = new String[2];
			for(int i=1; i<s.length; i++) parents[i-1] = s[i];
			closeReader();
			return parents;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		throw new RuntimeException("!!!");
	}

	public String[] getProgeny() {
		// TODO Auto-generated method stub
		try {
			setEntryReader("runinfo");
			String line;
			while(!(line=br.readLine()).startsWith("##progeny")) {}
			String[] s = line.trim().split("\\s+");
			String[] progeny = new String[s.length-1];
			System.arraycopy(s, 1, progeny, 0, progeny.length);
			closeReader();
			return progeny;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		throw new RuntimeException("!!!");
	}

	public long getRandomSeed() {
		// TODO Auto-generated method stub
		try {
			setEntryReader("runinfo");
			String line;
			while(!(line=br.readLine()).startsWith("##seed")) {}
			long seed = Long.parseLong(line.split("\\s+")[1]);
			closeReader();
			return seed;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		throw new RuntimeException("!!!");
	}

	public int getIterationNo() {
		// TODO Auto-generated method stub
		try {
			setEntryReader("runinfo");
			String line;
			while(!(line=br.readLine()).startsWith("##iteration")) {}
			int iter = Integer.parseInt(line.split("\\s+")[1]);
			closeReader();
			return iter;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		throw new RuntimeException("!!!");
	}

	public String[] getChrs() {
		// TODO Auto-generated method stub
		try {
			setEntryReader("runinfo");
			String line;
			while(!(line=br.readLine()).startsWith("##chrs")) {}
			String[] chrs = line.split("\\s+")[1].trim().split(",");
			closeReader();
			return chrs;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		throw new RuntimeException("!!!");
	}

	public boolean[] getChrsRev() {
		// TODO Auto-generated method stub
		try {
			setEntryReader("runinfo");
			String line;
			while(!(line=br.readLine()).startsWith("##chrs_rev")) {}
			String[] revs = line.split("\\s+")[1].trim().split(",");
			closeReader();
			boolean[] chrs_rev = new boolean[revs.length];
			for(int i=0; i<revs.length; i++) 
				chrs_rev[i] = Boolean.parseBoolean(revs[i]);
			return chrs_rev;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		throw new RuntimeException("!!!");
	}

	public int[] getModelLength() {
		// TODO Auto-generated method stub
		try {
			setEntryReader("runinfo");
			String line;
			while(!(line=br.readLine()).startsWith("##model_len")) {}
			String[] lens = line.split("\\s+")[1].trim().split(",");
			closeReader();
			int[] modelLength = new int[lens.length];
			for(int i=0; i<lens.length; i++) 
				modelLength[i] = Integer.parseInt(lens[i]);
			return modelLength;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		throw new RuntimeException("!!!");
	}


	public double[] getModelLoglik() {
		// TODO Auto-generated method stub
		try {
			setEntryReader("runinfo");
			String line;
			while(!(line=br.readLine()).startsWith("##model_ll")) {}
			String[] lls = line.split("\\s+")[1].trim().split(",");
			closeReader();
			double[] model_ll = new double[lls.length];
			for(int i=0; i<lls.length; i++) 
				model_ll[i] = Double.parseDouble(lls[i]);
			return model_ll;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		throw new RuntimeException("!!!");
	}
	
	public double getLoglik() {
		// TODO Auto-generated method stub
		try {
			setEntryReader("runinfo");
			String line;
			while(!(line=br.readLine()).startsWith("##loglik")) {}
			double ll = Double.parseDouble(line.split("\\s+")[1]);
			closeReader();
			return ll;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		throw new RuntimeException("!!!");
	}

	public int[] getSnpPosition(String scaff) {
		// TODO Auto-generated method stub
		try {
			setEntryReader("snp");
			List<Integer> pos = new ArrayList<>();
			String line;
			String[] s;
			while( (line=br.readLine())!=null ) {
				s = line.split("\\s+");
				if(s[2].equals(scaff)) 
					pos.add(Integer.parseInt(s[3]));
			}
			closeReader();
			if(pos.get(0)>pos.get(1)) Collections.reverse(pos);
			
			return ArrayUtils.toPrimitive(pos.toArray(new Integer[pos.size()]));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		throw new RuntimeException("!!!");		
	}

	public int[] getHapCounts() {
		// TODO Auto-generated method stub
		try {
			setEntryReader("haplotype");
			String line, states;
			String[] s;
			Map<Character, Integer> counts = new HashMap<>();
			while( (line=br.readLine())!=null ) {
				if(!line.startsWith("#")) continue;
				s = line.split("\\s+");
				states = s[s.length-1];
				if(states.startsWith("*")) continue;
				for(char h : states.toCharArray())
					counts.put(h, counts.getOrDefault(h, 0)+1);
			}
			closeReader();
			return ArrayUtils.toPrimitive(counts.values().toArray(new Integer[counts.size()]));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		throw new RuntimeException("!!!");
	}

	public Map<String, char[][]> getHaplotypeByPosition(int[] position, int ploidy) {
		// TODO Auto-generated method stub
		try {
			setEntryReader("haplotype");
			final Map<String, char[][]> haps = new HashMap<>();
			final int P = position.length;
			String line, states;
			String[] s;

			char[][] hap = new char[ploidy][P];
			int i = 0;
			while( (line=br.readLine())!=null ) {
				if(!line.startsWith("#")) continue;
				s = line.split("\\s+");
				states = s[s.length-1];
				for(int j=0; j<P; j++)
					hap[i][j] = states.charAt(position[j]);
				++i;
				if(i==ploidy) {
					haps.put(s[2].split(":")[0], hap);
					hap = new char[ploidy][P];
					i = 0;
				}
			}
			closeReader();
			return haps;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		throw new RuntimeException("!!!");
	}
	
	public Map<String, char[][]> getHaplotypeByPositionRange(int[] position, int ploidy) {
		// TODO Auto-generated method stub
		try {
			setEntryReader("haplotype");
			final Map<String, char[][]> haps = new HashMap<>();
			final int P = Math.abs(position[0]-position[1])+1;
			final boolean rev = position[0]>position[1];
			String line, states;
			String[] s;

			char[][] hap = new char[ploidy][P];
			int i = 0;
			while( (line=br.readLine())!=null ) {
				if(!line.startsWith("#")) continue;
				s = line.split("\\s+");
				states = s[s.length-1];
				if(rev) {
					for(int j=0; j<P; j++)
						hap[i][j] = states.charAt(position[0]-j);	
				} else {
					for(int j=0; j<P; j++)
						hap[i][j] = states.charAt(position[0]+j);
				}
				++i;
				if(i==ploidy) {
					haps.put(s[2].split(":")[0], hap);
					hap = new char[ploidy][P];
					i = 0;
				}
			}
			closeReader();
			return haps;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		throw new RuntimeException("!!!");
	}
}