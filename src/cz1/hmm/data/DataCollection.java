package cz1.hmm.data;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import org.apache.commons.lang3.ArrayUtils;

import cz1.util.Algebra;
import cz1.util.Combination;
import cz1.util.Constants;
import cz1.util.IO;

public class DataCollection {
	
	public static void zip(String wd, String experimentId, 
			String vcfFilePath) {
		String zipFilePath = wd+Constants.file_sep+experimentId+".zip";
		List<String> info = Arrays.asList(Constants.
				_vcf_format_str.split(":"));
		int _idx_gt = info.indexOf("GT"),
				_idx_pl = info.indexOf("PL"),
				_idx_ad = info.indexOf("AD");
		try {
			
			ZipOutputStream out = new ZipOutputStream(new 
					FileOutputStream(zipFilePath));
			
			BufferedReader br = IO.getBufferedReader(vcfFilePath);
			String line;
			while( (line=br.readLine()).startsWith("##") ) {}
			String[] s = line.split("\\s+");
			List<String> sList = Arrays.asList(s);
			int _idx_start = sList.indexOf("FORMAT")+1,
					_idx_pos = sList.indexOf("POS"),
					_idx_ref = sList.indexOf("REF"),
					_idx_alt = sList.indexOf("ALT");
			out.putNextEntry(new ZipEntry("samples"));
			for(int i=_idx_start; i<s.length; i++) out.write((s[i]+
					Constants.line_sep).getBytes());
			List<Contig> contigs = new ArrayList<Contig>();
			line = br.readLine();
			while(line != null) {
				int i=1;
				String contig = line.split("\\s+")[0];
				List<String[]> snps = new ArrayList<String[]>();
				snps.add(line.split("\\s+"));
				while( (line=br.readLine())!=null && 
						line.split("\\s+")[0].equals(contig)) {
					i++;
					snps.add(line.split("\\s+"));
				}
				contigs.add(new Contig(contig,i));
				out.putNextEntry(new ZipEntry(contig+
						Constants.file_sep+"position"));
				for(int j=0; j<snps.size(); j++)
					out.write((snps.get(j)[_idx_pos]+Constants.line_sep).
							getBytes());
				out.putNextEntry(new ZipEntry(contig+
						Constants.file_sep+"allele"));
				for(int j=0; j<snps.size(); j++)
					out.write((snps.get(j)[_idx_ref]+
							"\t"+snps.get(j)[_idx_alt]+
							Constants.line_sep).
							getBytes());
				out.putNextEntry(new ZipEntry(contig+
						Constants.file_sep+"GT"));
				for(int j=0; j<snps.size(); j++) {
					for(int k=_idx_start; k<snps.get(j).length; k++)
						out.write((snps.get(j)[k].split(":")[_idx_gt]+"\t")
								.getBytes());
					out.write("\n".getBytes());
				}
				out.putNextEntry(new ZipEntry(contig+
						Constants.file_sep+"AD"));
				for(int j=0; j<snps.size(); j++) {
					for(int k=_idx_start; k<snps.get(j).length; k++)
						try {
							out.write((snps.get(j)[k].split(":")[_idx_ad]+"\t")
									.getBytes());
						} catch (ArrayIndexOutOfBoundsException e) {
							out.write(".\t".getBytes());
						}
					out.write("\n".getBytes());
				}
				out.putNextEntry(new ZipEntry(contig+
						Constants.file_sep+"PL"));
				for(int j=0; j<snps.size(); j++) {
					for(int k=_idx_start; k<snps.get(j).length; k++)
						try {
							out.write((snps.get(j)[k].split(":")[_idx_pl]+"\t")
									.getBytes());
						} catch (ArrayIndexOutOfBoundsException e) {
							out.write(".\t".getBytes());
						}
					out.write("\n".getBytes());
				}
			}
			br.close();
			
			Collections.sort(contigs);
			Collections.reverse(contigs);
			
			out.putNextEntry(new ZipEntry("contig"));
			for(Contig contig : contigs)
				out.write((contig.id+"\t"+contig.markers+
						Constants.line_sep).getBytes());
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static class Contig implements Comparable<Contig> {
		private String id;
		private int markers;
		
		public Contig(String id,
				int markers) {
			this.id = id;
			this.markers = markers;
		}
		
		@Override
		public int compareTo(Contig contig) {
			// TODO Auto-generated method stub
			return this.markers-contig.markers;
		}
	}
	
	public static DataEntry[] readDataEntry(String zipFilePath,
			String[] contig) {
		DataEntry[] de = new DataEntry[contig.length];
		for(int i=0; i<de.length; i++) {
			double[] position = readPosition(zipFilePath, contig[i]);
			List<String[]> allele = readAllele(zipFilePath, contig[i]);
			List<List<int[]>> ad = readAlleleDepth(zipFilePath, 
					contig[i], allele);
			List<List<double[]>> pl = readPhredScaledLikelihood(
					zipFilePath, contig[i], allele);
			List<List<char[]>> gt = readGenotype(zipFilePath, contig[i],
					allele);
			List<String> sample = getSampleList(zipFilePath);
			de[i] = new DataEntry(contig[i], position, allele, ad, pl, gt, sample);
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
			BufferedWriter bw = IO.getBufferedWriter(output,append);
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
			br = IO.getBufferedReader(template);
			String header = "";
			while(!header.startsWith("#CHROM")) header=br.readLine();
			br.close();
			br = IO.getBufferedReader(adFile);
			br2 = IO.getBufferedReader(SOAPoutput);
			BufferedWriter wr = IO.getBufferedWriter(output);
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
	
	private static List<List<char[]>> readGenotype(String zipFilePath,
			String contigId, List<String[]> allele) {
		// TODO Auto-generated method stub
		try {
			final ZipFile in = new ZipFile(zipFilePath);
			final InputStream is = in.getInputStream(
					in.getEntry(contigId+Constants.file_sep+"GT"));
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			List<char[]> entry;
			char[] genotype;
			String line;
			String[] s, s2;
			List<List<char[]>> gt = new ArrayList<List<char[]>>();
			int l=0;
			while( (line=br.readLine())!=null ) {
				Map<Integer, Character> gtMap = getUniversalAlleleMap(
						allele.get(l++).length);
				s = line.split("\\s+");
				entry = new ArrayList<char[]>();
				for(int i=0; i<s.length; i++) {
                    genotype = new char[Constants._ploidy_H];
                    s2 = s[i].split("/");
                    for(int j=0; j<s2.length; j++) 
                        if(s2[j].equals("."))
                            genotype[j] = 'N';
                        else
                            genotype[j] = gtMap.get(Integer.parseInt(s2[j]));
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

	private static Map<Integer, Character> getUniversalAlleleMap(int n) {
		// TODO Auto-generated method stub
		Map<Integer, Character> alleleMap = 
				new HashMap<Integer, Character>();
		for(int i=0; i<n; i++) {
			alleleMap.put(i, (char)('A'+i));
		}
		return alleleMap;
	}

	private static List<List<double[]>> readPhredScaledLikelihood(
			String zipFilePath, String contigId, List<String[]> allele) {
		// TODO Auto-generated method stub
		try {
			final ZipFile in = new ZipFile(zipFilePath);
			final InputStream is = in.getInputStream(
					in.getEntry(contigId+Constants.file_sep+"PL"));
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			List<double[]> entry;
			double[] ll;
			String line;
			String[] s, s2;
			List<List<double[]>> pl = new ArrayList<List<double[]>>();
			int l = 0;
			while( (line=br.readLine())!=null ) {
				int n = allele.get(l++).length;
				s = line.split("\\s+");
				entry = new ArrayList<double[]>();
				for(int i=0; i<s.length; i++) {
                    
					ll = new double[Combination.nmultichoosek(n, 
                            Constants._ploidy_H)];
                    boolean miss = true;
                    if(!s[i].equals(".")) {
                        s2 = s[i].split(",");
                        /**
                        for(int j=0; j<s2.length; j++) 
                            ll[j] = Math.pow(10,-
                                    Double.parseDouble(s2[j])/10);
                        **/
                        for(int j=0; j<s2.length; j++) {
                            ll[j] = Math.pow(10,-
                                    Double.parseDouble(s2[j])/10);
                            miss = miss&&ll[j]==1.0;
                        }
                    }
                    
                    if(miss) 
                    	Arrays.fill(ll, -1);
                    else
                    	ll = Algebra.normalize(ll);
                    
                    entry.add(ll);
                    
                    //int a = Algebra.maxIndex(ll);
                    //Arrays.fill(ll, 0);
                    //ll[a] = 1;
                    //entry.add(ll);
				}
				pl.add(entry);
			}
			in.close();
			return pl;
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
			final InputStream is = in.getInputStream(
					in.getEntry(contigId+Constants.file_sep+"AD"));
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
					in.getEntry(contigId+Constants.file_sep+"allele"));
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
            //System.out.println(contigId+Constants.file_sep+"position");

			final InputStream is = in.getInputStream(
					in.getEntry(contigId+Constants.file_sep+"position"));
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
			int[] contigIndex) {
		String[] contig = new String[contigIndex.length];
		for(int i=0; i<contig.length; i++)
			contig[i] = getContigIDFromIndex(zipFilePath,
					contigIndex[i]);
		return readDataEntry(zipFilePath, contig);
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
	
	public static int[] getParentIndex(String zipFilePath) {
		List<String> samples = getSampleList(zipFilePath);
		String[] parents = Constants._founder_haps.split(":");
		List<Integer> list = new ArrayList<Integer>();
		for(String parent : parents) list.add(samples.indexOf(parent));
	    return ArrayUtils.toPrimitive(
                list.toArray(new Integer[list.size()]));
    }
	
	public static int getSampleNumber(String zipFilePath) {
		return getSampleList(zipFilePath).size();
	}
}
