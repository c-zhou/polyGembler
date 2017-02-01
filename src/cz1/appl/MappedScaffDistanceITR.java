package cz1.appl;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.lang3.StringUtils;

import cz1.util.IO;

public class MappedScaffDistanceITR {

	public static void main(String[] args) throws IOException, InterruptedException {
		CommandLineParser parser = new PosixParser();

		// create the Options
		Options options = new Options();
		options.addOption( "f", "rf-file", true, "recombination frequency file." );
		options.addOption( "o", "output", true, "output file.");
		options.addOption( "m", "map-file", true, "scaffolds-chromosome colinear mapping file.");
		options.addOption( "g", "gm-file", true, "genetic linkage map file.");
		options.addOption( "s", "snp-file", true, "snp file.");
		String output=null, rfFile=null, clFile=null, gmFile=null, snpFile=null;
		try {
			// parse the command line arguments
			CommandLine line = parser.parse( options, args );
			if( line.hasOption("f") ) {
				rfFile = line.getOptionValue('f');
			} else throw new RuntimeException("!!!");
			if(line.hasOption("o")) {
				output = line.getOptionValue("o");
			} else throw new RuntimeException("!!!");
			if(line.hasOption("m")) {
				clFile = line.getOptionValue("m");
			} else throw new RuntimeException("!!!");
			if(line.hasOption("g")) {
				gmFile = line.getOptionValue("g");
			} else throw new RuntimeException("!!!");
			if(line.hasOption("s")) {
				snpFile = line.getOptionValue("s");
			} else throw new RuntimeException("!!!");
		}
		catch( ParseException exp ) {
			System.out.println( "Unexpected exception:" + exp.getMessage() );
		}

		MappedScaffDistanceITR.calculate(rfFile, clFile, gmFile, 
				snpFile, output);
	}

	private static class Scaffold implements Comparable<Scaffold> {
		private String id;
		private final double start;
		private final double end;
		
		public Scaffold(String id, double start, double end) {
			this.id = id;
			this.start = start;
			this.end = end;
		}
		
		@Override
		public int compareTo(Scaffold scaff) {
			// TODO Auto-generated method stub
			return scaff.id.compareTo(this.id);
		}
		
		public void setId(String id) {
			this.id = id;
		}
	}
	
	private static void calculate(String rfFile,
			String clFile, String gmFile, 
			String snpFile, 
			String output) {
		// TODO Auto-generated method stub
		Map<String, double[]> snp_position = new HashMap<String, double[]>();
		Map<String, Scaffold> map_info = new HashMap<String, Scaffold>();
		try {
			String line;
			String[] s;
			
			BufferedReader br = IO.getBufferedReader(snpFile);
			line = br.readLine();
			
			while( line!=null ) {
				s = line.split("_");
				double start = Double.parseDouble(s[s.length-1]),
						end = start;
				String[] ss = new String[s.length-1];
				System.arraycopy(s, 0, ss, 0, ss.length);
				String scf = StringUtils.join(ss,'_');
				
				while( (line=br.readLine())!=null ){
					s = line.split("_");
					ss = new String[s.length-1];
					System.arraycopy(s, 0, ss, 0, ss.length);
					if(scf.equals(StringUtils.join(ss,'_'))) 
						end = Double.parseDouble(s[s.length-1]);
					else
						break;
				}
				
				snp_position.put(scf.replace("scaffold", ""), new double[]{start, end});
			}
			br.close();
			
			br = IO.getBufferedReader(clFile);
			line = br.readLine();
			while( line!=null ) {
				if(line.startsWith("#")) {
					String ref = line.split("\\s+")[1];
					while( (line=br.readLine())!=null 
							&& !line.startsWith("#") ) {
						s = line.split("\\s+");
						if(!snp_position.containsKey(s[0]))
							continue;
						
						double[] np = snp_position.get(s[0]);
						double[] rp = new double[]{
								Double.parseDouble(s[1]), 
								Double.parseDouble(s[2])};
						double[] sp = new double[]{
								Double.parseDouble(s[3]), 
								Double.parseDouble(s[4])};
						if(np[0]>np[1] && sp[0]<sp[1] ||
								np[0]<np[1] && sp[0]>sp[1]) {
							double t = np[0];
							np[0] = np[1];
							np[1] = t;
						}
						
						double frac = Math.abs((rp[0]-rp[1])/(sp[0]-sp[1]));
						double start = Double.NEGATIVE_INFINITY,
								end = Double.NEGATIVE_INFINITY;
						if(sp[0]<sp[1]) {
							start = rp[0]+(np[0]-sp[0])*frac;
							end = rp[1]+(np[1]-sp[1])*frac;
						} else {
							start = rp[0]-(np[0]-sp[0])*frac;
							end = rp[1]-(np[1]-sp[1])*frac;
						}
						
						Scaffold scaff = new Scaffold(ref, 
								start, 
								end);
						map_info.put(s[0], scaff);
					}
 				}
			}
			br.close();
			
			br = IO.getBufferedReader(rfFile);
			String a, b;
			Map<String, String> rfs = new HashMap<String, String>();
			while( (line=br.readLine())!=null ) {
				s = line.split("\\s+");
				a = s[6].replace("scaffold", "");
				b = s[7].replace("scaffold", "");
				rfs.put(a+"_"+b, s[1]);
				rfs.put(b+"_"+a, s[1]);
			}
			br.close();

			
			Map<String, Set<Scaffold>> buckets = new HashMap<String, Set<Scaffold>>();
			
			br = IO.getBufferedReader(gmFile);
			while( !(line=br.readLine() ).startsWith("$order") ) {} ;
			while( (line=br.readLine() ).length()>0 ) {
				s = line.split("\\s+")[1].split("-");
				for(int i=0; i<s.length; i++) {
					s[i] = s[i].replace("scaffold", "");
					if(!map_info.containsKey(s[i])) continue;
					Scaffold si = map_info.get(s[i]);
					String id = si.id;
					if(!buckets.containsKey(id)) 
						buckets.put(id, new HashSet<Scaffold>());
					si.setId(s[i]);
					buckets.get(id).add(si);
				}
			}
			br.close();

			BufferedWriter bw = IO.getBufferedWriter(output+".ds0");
			for(String key : buckets.keySet()) {
				Set<Scaffold> bucket = buckets.get(key);
				
				for(Scaffold si : bucket) 
					for(Scaffold sj : bucket)
						if(rfs.containsKey(si.id+"_"+sj.id))
							bw.write(si.id+"\t"+si.id+"\t"+sj.id+"\t"+
									si.start+"\t"+si.end+"\t"+sj.start+"\t"+
									sj.end+"\t"+distance(si,sj)+"\t"+
									rfs.get(si.id+"_"+sj.id)+"\n");
			}
			bw.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return;
	}

	private static double distance(Scaffold x, Scaffold y) {
		// TODO Auto-generated method stub
		return Math.min(Math.abs(x.start-y.start), 
				Math.min(Math.abs(x.start-y.end),
						Math.min(Math.abs(x.end-y.start),
								Math.abs(x.end-y.end))));
	}
}
