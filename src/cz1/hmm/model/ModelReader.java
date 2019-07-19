package cz1.hmm.model;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
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
			String[] parents = new String[s.length-1];
			System.arraycopy(s, 1, parents, 0, parents.length);
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

	public int[] getDistance() {
		// TODO Auto-generated method stub
		try {
			setEntryReader("runinfo");
			String line;
			while(!(line=br.readLine()).startsWith("##distance")) {}
			String[] s = line.trim().split("\\s+");
			int[] distance = new int[s.length-1];
			for(int i=1; i<s.length; i++) 
				distance[i-1] = Integer.parseInt(s[i]);
			closeReader();
			return distance;
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

	public List<String> getSnpId() {
		// TODO Auto-generated method stub
		try {
			setEntryReader("snp");
			List<String> snps = new ArrayList<>();
			String line;
			while( (line=br.readLine())!=null )
				snps.add(line.split("\\s+")[2]);
			closeReader();
			return snps;
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
				for(char h : states.toCharArray()) {
					if(!counts.containsKey(h)) counts.put(h, 0);
					counts.put(h, counts.get(h)+1);
				}
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
}