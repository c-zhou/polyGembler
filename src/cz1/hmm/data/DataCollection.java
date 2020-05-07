package cz1.hmm.data;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.zip.ZipFile;

import org.apache.commons.lang3.ArrayUtils;

import cz1.hmm.tools.VCFtools;
import cz1.math.Combination;
import cz1.math.Algebra;
import cz1.util.Utils;

public class DataCollection {

	public static DataEntry[] readDataEntry(String zipFilePath,
			String[] contig, int ploidy) {
		DataEntry[] de = new DataEntry[contig.length];
		for(int i=0; i<de.length; i++) {
			double[] position = readPosition(zipFilePath, contig[i]);
			List<String[]> allele = readAllele(zipFilePath, contig[i]);
			List<List<int[]>> ad = readAlleleDepth(zipFilePath, 
					contig[i], allele);
			List<List<double[]>> gl = readGenotypeLikelihood(
					zipFilePath, contig[i], allele, ploidy);
			List<List<String[]>> gt = readGenotype(zipFilePath, contig[i],
					allele);
			List<String> sample = getSampleList(zipFilePath);
			de[i] = new DataEntry(contig[i], position, allele, ad, gl, gt, sample);
		}
		return de;
	}

	public static DataEntry[] readDataEntry(String zipFilePath,
			String[] contig, int[] startPos, int[] endPos, int ploidy) {
		
		DataEntry[] de = new DataEntry[contig.length];
		for(int i=0; i<de.length; i++) {
			double[] position = readPosition(zipFilePath, contig[i]);
			
			List<String[]> allele = readAllele(zipFilePath, contig[i]);
			List<List<int[]>> ad = readAlleleDepth(zipFilePath, 
					contig[i], allele);
			List<List<double[]>> gl = readGenotypeLikelihood(
					zipFilePath, contig[i], allele, ploidy);
			List<List<String[]>> gt = readGenotype(zipFilePath, contig[i],
					allele);
			List<String> sample = getSampleList(zipFilePath);
			
			int is=0, ie=0;
			for(int j=0; j!=position.length; j++) {
				if(position[j]<startPos[i]) is++; 
				if(position[j]<=endPos[i])  ie++;
			}
			
			double[] position2 = new double[ie-is];
			System.arraycopy(position, is, position2, 0, ie-is);
			
			de[i] = new DataEntry(contig[i], 
					position2, 
					allele.subList(is, ie), 
					ad.subList(is, ie), 
					gl.subList(is, ie), 
					gt.subList(is, ie), 
					sample);
		}
		return de;
	}
	
	public static void writeSOAPInputFile(String zipFilePath,
			int contigIndex, String output) {
		writeSOAPInputFile(zipFilePath, 
				getContigIDFromIndex(zipFilePath, contigIndex),
				output);
	}

	public static void writeSOAPInputFile(String zipFilePath,
			String contigId, String output) {
		writeSOAPInputFileMeta(zipFilePath, contigId, output, false);
	}

	public static void writeSOAPInputFileAll(String zipFilePath, String output) {
		File o = new File(output);
		if(o.exists() && o.isFile()) o.delete();
		List<String> contigs = getContigList(zipFilePath);
		for(String contig:contigs) writeSOAPInputFileMeta(zipFilePath, 
				contig, output, true);
	}

	public static void writeSOAPInputFileMeta(String zipFilePath,
			String contigId, String output, boolean append) {
		double[] position = readPosition(zipFilePath, contigId);
		List<String[]> allele = readAllele(zipFilePath, contigId);
		List<List<int[]>> ad = readAlleleDepth(zipFilePath, 
				contigId, allele);
		try {
			BufferedWriter bw = Utils.getBufferedWriter(output,append);
			for(int i=0; i<position.length; i++) {
				String line = ">"+contigId+"\t";
				line += position[i];
				line += "\tB";
				for(int j=2; j<allele.get(i).length; j++)
					line += "\tB"+j;
				bw.write(line+"\n");
				List<int[]> a = ad.get(i);
				for(int j=0; j<a.get(0).length; j++) {
					line = "";
					for(int k=0; k<a.size(); k++)
						line +=a.get(k)[j]+"\t";
					bw.write(line+"\n");
				}
			}
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
	}

	public static void reformatSOAPOutputFile(String template, String adFile,
			String SOAPoutput, String output) {

		Map<String, String> gpMap = new HashMap<String, String>();
		gpMap.put("0/0/0/0", "0,255,255,255,255");
		gpMap.put("0/0/0/1", "255,0,255,255,255");
		gpMap.put("0/0/1/1", "255,255,0,255,255");
		gpMap.put("0/1/1/1", "255,255,255,0,255");
		gpMap.put("1/1/1/1", "255,255,255,0,255");
		try {
			BufferedReader br, br2;
			br = Utils.getBufferedReader(template);
			String header = "";
			while(!header.startsWith("#CHROM")) header=br.readLine();
			br.close();
			br = Utils.getBufferedReader(adFile);
			br2 = Utils.getBufferedReader(SOAPoutput);
			BufferedWriter wr = Utils.getBufferedWriter(output);
			String line, _new_line;
			String[] s, _ad1, _ad2;
			int count = 0;
			while( (line=br2.readLine())!=null ) {
				if( !line.startsWith("#") ) {
					s = line.split("\\s+");
					_new_line = "";
					_new_line += s[0].replace(">", "");
					for(int i=1; i<8; i++) _new_line += "\t"+s[i];
					_new_line += "\t"+s[8]+":AD:PL";
					br.readLine();
					_ad1 = br.readLine().split("\\s+");
					_ad2 = br.readLine().split("\\s+");
					for(int i=0; i<_ad1.length; i++) {
						String g = s[i+9].split(":")[0];
						_new_line += "\t"+s[i+9]+":"+_ad1[i]+","+_ad2[i]+":"+gpMap.get(g);
					}
					wr.write(_new_line+"\n");
					count++;
					if(count%1000==0) System.err.println(count+" done.");
				} else {
					if(line.startsWith("##")) {
						wr.write(line+"\n");
					} else {
						wr.write(header+"\n");
					}
				}
			}
			br.close();
			br2.close();
			wr.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	private static List<List<String[]>> readGenotype(String zipFilePath,
			String contigId, List<String[]> allele) {
		// TODO Auto-generated method stub
		try {
			final ZipFile in = new ZipFile(zipFilePath);
			if(in.getEntry(contigId+"/GT")==null) {
				in.close();
				return null;
			}
			final InputStream is = in.getInputStream(
					in.getEntry(contigId+"/GT"));
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			List<String[]> entry;
			String[] genotype;
			String line;
			String[] s, s2;
			List<List<String[]>> gt = new ArrayList<List<String[]>>();
			int l=0;
			while( (line=br.readLine())!=null ) {
				String[] allele_l = allele.get(l++);
				s = line.split("\\s+");
				entry = new ArrayList<String[]>();
				for(int i=0; i<s.length; i++) {
					s2 = s[i].split("/");
					genotype = new String[s2.length];
					for(int j=0; j<s2.length; j++) 
						if(s2[j].equals("."))
							genotype[j] = ".";
						else
							genotype[j] = allele_l[Integer.parseInt(s2[j])];
					entry.add(genotype);
				}
				gt.add(entry);
			}
			in.close();
			return gt;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.err.println("Invalid contig index. Program halted.");
		System.exit(1);
		return null;
	}

	private static List<List<double[]>> readGenotypeLikelihood(
			String zipFilePath, String contigId, List<String[]> allele, int ploidy) {
		// TODO Auto-generated method stub
		try {
			final ZipFile in = new ZipFile(zipFilePath);
			if(in.getEntry(contigId+"/PL")==null) {
				in.close();
				return null;
			}
			final InputStream is = in.getInputStream(
					in.getEntry(contigId+"/PL"));
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			List<double[]> entry;
			double[] ll;
			String line;
			String[] s, s2;
			List<List<double[]>> gl = new ArrayList<List<double[]>>();
			int l = 0;
			while( (line=br.readLine())!=null ) {
				int n = allele.get(l++).length;
				s = line.split("\\s+");
				entry = new ArrayList<double[]>();
				for(int i=0; i<s.length; i++) {

					ll = new double[Combination.nmultichoosek(n, ploidy)];
					boolean miss = true;
					if(!s[i].equals(".")) {
						s2 = s[i].split(",");
						
						for(int j=0; j<s2.length; j++) {
							ll[j] = Math.pow(10,-
									Double.parseDouble(s2[j])/10);
							miss = miss&&ll[j]==1.0;
						}
					}

					if(miss) 
						Arrays.fill(ll, 0);
					else
						ll = Algebra.normalize(ll);

					entry.add(ll);
				}
				gl.add(entry);
			}
			in.close();
			return gl;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.err.println("Invalid contig index. Program halted.");
		System.exit(1);
		return null;
	}

	private static List<List<int[]>> readAlleleDepth(String zipFilePath,
			String contigId, List<String[]> allele) {
		// TODO Auto-generated method stub
		try {
			final ZipFile in = new ZipFile(zipFilePath);
			if(in.getEntry(contigId+"/AD")==null) {
				in.close();
				return null;
			}
			final InputStream is = in.getInputStream(
					in.getEntry(contigId+"/AD"));
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			List<int[]> entry;
			int[] depth;
			String line;
			String[] s, s2;
			List<List<int[]>> ad = new ArrayList<List<int[]>>();
			int l=0;
			while( (line=br.readLine())!=null ) {
				int n = allele.get(l++).length;
				s = line.split("\\s+");
				entry = new ArrayList<int[]>();
				for(int i=0; i<s.length; i++) {
					depth = new int[n];
					if(s[i].equals("."))
						for(int j=0; j<depth.length; j++)
							depth[j] = 0;
					else {
						s2 = s[i].split(",");
						// should be s2.length!!!
						for(int j=0; j<depth.length; j++)
							depth[j] = Integer.parseInt(s2[j]);
					}
					entry.add(depth);
				}
				ad.add(entry);
			}
			in.close();
			return ad;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.err.println("Invalid contig index. Program halted.");
		System.exit(1);
		return null;
	}

	private static List<String[]> readAllele(String zipFilePath, 
			String contigId) {
		// TODO Auto-generated method stub
		try {
			final ZipFile in = new ZipFile(zipFilePath);
			final InputStream is = in.getInputStream(
					in.getEntry(contigId+"/allele"));
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			List<String> entry;
			String line;
			String[] s;
			List<String[]> allele = new ArrayList<String[]>();
			while( (line=br.readLine())!=null ) {
				s = line.split("\\s+|,");
				entry = new ArrayList<String>();
				for(int i=0; i<s.length; i++) {
					entry.add(s[i]);
				}
				allele.add(entry.toArray(new String[entry.size()]));
			}
			in.close();
			return allele;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.err.println("Invalid contig index. Program halted.");
		System.exit(1);
		return null;
	}

	private static double[] readPosition(String zipFilePath, String contigId) {
		// TODO Auto-generated method stub
		try {
			final ZipFile in = new ZipFile(zipFilePath);

			System.out.println(zipFilePath);

			final InputStream is = in.getInputStream(
					in.getEntry(contigId+"/position"));
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			List<Double> pos = new ArrayList<Double>();
			String line;
			while( (line=br.readLine())!=null )
				pos.add(Double.parseDouble(line));
			in.close();
			return ArrayUtils.toPrimitive(pos.toArray(new Double[pos.size()]));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.err.println("Invalid contig index. Program halted.");
		System.exit(1);
		return null;
	}

	public static DataEntry[] readDataEntry(String zipFilePath,
			int[] contigIndex, int ploidy) {
		String[] contig = new String[contigIndex.length];
		for(int i=0; i<contig.length; i++)
			contig[i] = getContigIDFromIndex(zipFilePath,
					contigIndex[i]);
		return readDataEntry(zipFilePath, contig, ploidy);
	}

	public static String getContigIDFromIndex(String zipFilePath,
			int contigIndex) {
		// TODO Auto-generated method stub
		String contig = null;
		try {
			final ZipFile in = new ZipFile(zipFilePath);
			final InputStream is = in.getInputStream(
					in.getEntry("contig"));
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			while(contigIndex-- > 0) contig = br.readLine();
			br.close();
			in.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return contig.split("\\s+")[0];
	}

	public static int getSampleIndexFromId(String zipFilePath,
			String sampleId) {
		// TODO Auto-generated method stub
		int index = 0;
		try {
			final ZipFile in = new ZipFile(zipFilePath);
			final InputStream is = in.getInputStream(
					in.getEntry("samples"));
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			String sample;
			while( (sample=br.readLine())!=null ) 
				if( sample.equals(sampleId) ) { 
					br.close();
					in.close();
					return index;
				} else
					index++;
			br.close();
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return -1;
	}

	public static Map<String, Integer> getContigSizeMap(
			String zipFilePath) {
		try {
			final ZipFile in = new ZipFile(zipFilePath);
			final InputStream is = in.getInputStream(
					in.getEntry("contig"));
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			String line;
			Map<String, Integer> map = new HashMap<String, Integer>();
			while( (line=br.readLine())!=null ) {
				String[] s = line.split("\\s+");
				map.put(s[0], Integer.parseInt(s[1]));
			}
			in.close();
			return map;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.err.println("File corrupted. Program halted.");
		System.exit(1);
		return null;
	}

	public static List<String> getContigList(
			String zipFilePath) {
		try {
			final ZipFile in = new ZipFile(zipFilePath);
			final InputStream is = in.getInputStream(
					in.getEntry("contig"));
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			String line;
			List<String> contigs = new ArrayList<String>();
			while( (line=br.readLine())!=null ) {
				String[] s = line.split("\\s+");
				contigs.add(s[0]);
			}
			in.close();
			return contigs;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.err.println("File corrupted. Program halted.");
		System.exit(1);
		return null;
	}

	public static Map<Integer, String> getContigIndexMap(
			String zipFilePath) {
		try {
			final ZipFile in = new ZipFile(zipFilePath);
			final InputStream is = in.getInputStream(
					in.getEntry("contig"));
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			String line;
			Map<Integer, String> map = new HashMap<Integer, String>();
			int i=0;
			while( (line=br.readLine())!=null ) {
				String[] s = line.split("\\s+");
				map.put(++i,s[0]);
			}
			in.close();
			return map;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.err.println("File corrupted. Program halted.");
		System.exit(1);
		return null;
	}

	public static List<String> getSampleList(String zipFilePath) {
		try {
			final ZipFile in = new ZipFile(zipFilePath);
			final InputStream is = in.getInputStream(
					in.getEntry("samples"));
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			String line;
			List<String> samples = new ArrayList<String>();
			while( (line=br.readLine())!=null ) samples.add(line);
			in.close();
			return samples;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.err.println("File corrupted. Program halted.");
		System.exit(1);
		return null;
	}

	public static int[] getParentIndex(String zipFilePath, String[] parents) {
		List<String> samples = getSampleList(zipFilePath);
		List<Integer> list = new ArrayList<Integer>();
		for(String parent : parents) list.add(samples.indexOf(parent));
		return ArrayUtils.toPrimitive(
				list.toArray(new Integer[list.size()]));
	}

	public static int getSampleNumber(String zipFilePath) {
		return getSampleList(zipFilePath).size();
	}

	public static Map<String, Integer> readScaff(String zipFilePath) {
		// TODO Auto-generated method stub
		final Map<String, Integer> scaffs = new HashMap<String, Integer>();
		try {
			final ZipFile in = new ZipFile(zipFilePath);

			final InputStream is = in.getInputStream(in.getEntry("contig"));
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			String line;
			String[] s;
			while( (line=br.readLine())!=null ) {
				s = line.split("\\s+");
				scaffs.put(s[0], Integer.parseInt(s[1]));
			}
			in.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return scaffs;
	}

	public static Map<String, Integer> readScaff(String zipFilePath, 
			Set<String> targets) {
		// TODO Auto-generated method stub
		final Map<String, Integer> scaffs = new HashMap<String, Integer>();
		try {
			final ZipFile in = new ZipFile(zipFilePath);

			final InputStream is = in.getInputStream(in.getEntry("contig"));
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			String line;
			String[] s;
			while( (line=br.readLine())!=null ) {
				s = line.split("\\s+");
				if(targets.contains(s[0]))
					scaffs.put(s[0], Integer.parseInt(s[1]));
			}
			in.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return scaffs;
	}

	public static DataEntry[] readDataEntry(String in_dp_file, int ploidy) {
		// TODO Auto-generated method stub
		DataEntry[] de = new DataEntry[1];

		try {		
			BufferedReader br = Utils.getBufferedReader(in_dp_file);

			int N = Integer.parseInt(br.readLine());
			int M = Integer.parseInt(br.readLine());
			double[] position = new double[M];
			List<String> sample = new ArrayList<String>();
			String[] s = br.readLine().split("\\s+");
			for(int i=0; i!=N; i++) sample.add(s[i]);
			List<String[]> allele = new ArrayList<String[]>();
			List<List<int[]>> ad = new ArrayList<List<int[]>>();
			List<List<double[]>> gl = new ArrayList<List<double[]>>();
			br.readLine();
			br.readLine();
			String[] s0;
			for(int i=0; i!=M; i++) {
				s = br.readLine().split("\\s+");
				position[i] = Double.parseDouble(s[1]);
				String[] allele_i = s[2].split(",");
				allele.add(allele_i);
				List<int[]> ad_i = new ArrayList<int[]>();
				List<double[]> gl_i = new ArrayList<double[]>();
				for(int j=0; j!=N; j++) {
					s0 = s[j+3].split(",");
					int[] ad_ij = new int[2];
					if(s[j+3].equals(".")) {
						ad_ij[0] = 0;
						ad_ij[1] = 0;
					} else {
						ad_ij[0] = Integer.parseInt(s0[0]);
						ad_ij[1] = Integer.parseInt(s0[1]);
					}
					double[] pl_ij  = VCFtools.fit(ad_ij, ploidy);
					ad_i.add(ad_ij);
					gl_i.add(VCFtools.PL2GL(pl_ij));
				}
				ad.add(ad_i);
				gl.add(gl_i);
			}
			br.close();
			de[0] = new DataEntry(in_dp_file.split("\\.")[0]+"_"+N+"_"+M,
					position, allele, ad, gl, null, sample);
		} catch (IOException e) {
			e.printStackTrace();
			throw new RuntimeException("!!!");
		}
		return de;	
	}
}
