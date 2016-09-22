package cz1.test;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.commons.io.input.ReversedLinesFileReader;
import org.apache.commons.lang3.ArrayUtils;

public class PairwiseCNVHapRFTest {

	public static void main(String[] args) 
			throws IOException, InterruptedException {
		CommandLineParser parser = new PosixParser();

		// create the Options
		Options options = new Options();
		options.addOption( "e", "experiment", true, "experiment name." );
		options.addOption( "w", "workspace", true, "directory contains input files." );
		options.addOption( "o", "output", true, "output file prefix." );
		options.addOption( "t", "thread", true, "number threads.");
		String experiment=null, workspace=null, output=null;
		int t=1;
		try {
			// parse the command line arguments
			CommandLine line = parser.parse( options, args );
			if( line.hasOption("e") ) {
				experiment = line.getOptionValue('e');
			}
			if(line.hasOption("w")) {
				workspace = line.getOptionValue("w");
			}
			if(line.hasOption("o")) {
				output = line.getOptionValue("o");
			}
			if(line.hasOption("t")) {
				t = Integer.parseInt(line.getOptionValue("t"));
			}
		}
		catch( ParseException exp ) {
			System.out.println( "Unexpected exception:" + exp.getMessage() );
		}

		PairwiseCNVHapRFTest pwRF = new PairwiseCNVHapRFTest(
				workspace,experiment,t);
		pwRF.calcRFsForAll2(output);

	}
	
	public PairwiseCNVHapRFTest(String workspace,
			String experiment,
			int threads) {
		this.wd = workspace;
		this.experiment = experiment;
		THREADS = threads;
	}

	private class rfCalculator implements Runnable {

		private final String str_c;

		public rfCalculator(String str_c) {
			this.str_c = str_c;
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub
			List<String> all = fileMap.get(str_c);
			double[] ll = new double[all.size()];
			double[] rf = new double[all.size()];
			
			Arrays.fill(ll, 0);
			Arrays.fill(rf, -1);
			
			for(int i=0; i<all.size(); i++) {

				File in = new File(wd+"/"+all.get(i)+"/stderr_true");
				if(in.exists()) {
					ReversedLinesFileReader rbr;
					try {
						rbr = new ReversedLinesFileReader(in);
						if(!rbr.readLine().startsWith("sampling from HMM")) {
							rbr.close();
							continue;
						}
						int k=0;
						String str = null;
						while( k++<6 && (str=rbr.readLine())!=null ) {}
						if(str==null) {
							rbr.close();
							continue;
						}
						String[] s0 = str.split("\\s+");
						double r = Double.parseDouble(s0[1]);
						rbr.readLine();
						if( (str=rbr.readLine()) == null) {
							rbr.close();
							continue;
						}
						s0 = str.split(" ");
						ll[i] = Double.parseDouble(s0[3]);
						rf[i] = r;
						rbr.close();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
			}
			
			ll = removeZERO(ll);
			rf = removeNEG(rf);
			int[] m = median(ll);
			
			if(m==null) return;
				
			String[] s = this.str_c.split(":");
			try {
				rfWriter.write(s[0]+
						"\t"+s[1]+
						"\t"+(rf[m[0]]+rf[m[1]])/2+
						"\t"+(ll[m[0]]+ll[m[1]])/2+
						"\t"+cat(rf, ",")+
						"\t"+cat(ll, ",")+
						"\n");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		private String cat(double[] ds, String sep) {
			if(ds.length==0) return "";
			StringBuilder sb = new StringBuilder(""+ds[0]);
			for(int i=1; i<ds.length; i++)
				sb.append(sep+ds[i]);
			return sb.toString();
		}
		
		private int[] median(double[] ds) {
			// TODO Auto-generated method stub
			double[] ds0 = removeZERO(ds);
			if(ds0==null) return null;
			Arrays.sort(ds0);
			int n = ds0.length;
			if (n % 2 == 0)
				return new int[]{n/2,n/2-1};
			else
				return new int[]{n/2,n/2};
		}

		private double[] removeNEG(double[] ds) {
			List<Double> ds0 = new ArrayList<Double>();
			for(double d : ds)
				if(d>=0)
					ds0.add(d);
			if(ds0.isEmpty()) return null;
			return ArrayUtils.toPrimitive(
					ds0.toArray(new Double[ds0.size()]));

		}
		
		private double[] removeZERO(double[] ds) {
			List<Double> ds0 = new ArrayList<Double>();
			for(double d : ds)
				if(d!=0)
					ds0.add(d);
			if(ds0.isEmpty()) return null;
			return ArrayUtils.toPrimitive(
					ds0.toArray(new Double[ds0.size()]));

		}
	}

	private static int THREADS = 1;
	private String wd;
	private static BufferedWriter rfWriter;
	private static Map<String, List<String>> fileMap;
	private static ExecutorService executor = 
			Executors.newFixedThreadPool(THREADS);
	private void reset() {
		executor = Executors.newFixedThreadPool(THREADS);
	}
	private String experiment;
	
	private void calcRFsForAll2( 
			String out) throws IOException, InterruptedException {

		File folder = new File(wd);
		File[] listFiles = folder.listFiles();
		fileMap = new ConcurrentHashMap<String, List<String>>();
		List<String> list;
		String[] s;
		for(File file:listFiles) {
			String name = file.getName();
			if(file.isDirectory() &&
					name.startsWith(experiment)) {
				name = name.replace(experiment,"experiment");
				s = name.split("\\.")[1].split("_");
				String str_c = s[0].compareTo(s[2])>0 ?
						s[0]+":"+s[2] : 
							s[2]+":"+s[0];
						if(fileMap.get(str_c)==null) {
							list = new ArrayList<String>();
							list.add(file.getName());
							fileMap.put(str_c, list);
						} else{
							fileMap.get(str_c).add(file.getName());
						}
			}
		}

		rfWriter = new BufferedWriter(new FileWriter(out));
		rfWriter.write("#scaff1\tscaff2\trf\tll\trfAll\tllAll\n");
		reset();
		for(String str_c : fileMap.keySet()) {
			executor.submit(new rfCalculator(str_c));
		}
		executor.shutdown();
		executor.awaitTermination(365, TimeUnit.DAYS);
		rfWriter.close();
	}
}
