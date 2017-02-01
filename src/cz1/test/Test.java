package cz1.test;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.commons.math3.stat.StatUtils;

import cz1.util.Constants;
import cz1.util.IO;

public class Test {

	public static void main(String[] args) throws IOException {

		Random random = new Random();
		Map<Long, Integer> map = new HashMap<Long, Integer>();
		for(int i=0; i<2400; i++)
			map.put(random.nextLong(), random.nextInt());
		long elapse = 0;
		List<Long> keys = new ArrayList<Long>(map.keySet());
		int nn = keys.size();
		for(int i=0; i<164157600; i++) {
			long x = keys.get(random.nextInt(nn));
			long tic = System.nanoTime();
			map.get(x);
			elapse += System.nanoTime()-tic;
		}
		System.out.println(elapse);
		
		System.exit(1);
		long tic;
		double time1 = 0, time2 = 0;
		double x;
		double x1 = random.nextDouble(),
				x2 = random.nextDouble();
		for(int i=0; i<10000000; i++) {

			tic = System.nanoTime();
			x = x1*x2;
			time2 += System.nanoTime()-tic;
			
			tic = System.nanoTime();
			x = x1+x2;
			time1 += System.nanoTime()-tic;
			
		}
		System.out.println(time1+ "\t" + time2);
		
		System.exit(1);
		System.out.println("Itr_sc0000111.1".replaceAll("^Itr_sc0{0,6}", "").
				replaceAll(".1$", ""));
		System.exit(1);
		
		System.out.println(Math.pow(2, 31)-1);
		System.out.println(Integer.MAX_VALUE);
		
		int n = 4;
		for(int i=0; i<n; i++) 
			for(int j=i; j<n; j++)
				System.out.println(j+" "+i);
		System.out.println((long) 10032030000000000.99931);
		
		System.out.println(Double.MAX_VALUE);
		System.out.println(Double.MIN_VALUE);
		System.out.println(Constants.logThreshMax);

		System.out.println(StatUtils.variance(new double[]{2,4,6}));
		
		if(false) {

			String[][] a = new String[10][5];
			if(a[9][1] == null) System.out.println("NULL");

			RandomGenerator rg = new Well19937c(1234);
			BetaDistribution beta = new BetaDistribution(rg,1,1);
			for(int i=0; i<100; i++)
				System.out.println(Math.floor(beta.sample()*100));


			String ge3 = "C:\\Users\\chenxi.zhou\\Desktop\\putty\\itr.bwa2m.b30.nnj.ge3.ds",
					pw = "C:\\Users\\chenxi.zhou\\Desktop\\putty\\rf_for_r_itr_pairwise.txt";

			int[] sep1 = new int[]{2, 33};
			int[] all1 = new int[] {642,609,362,1383,190,783,
					276,420,579,29,402,489,470,
					562,354,229,41,810,157,909,
					166,104,2377,17,467,139,64,
					93,523,135,100,83,440,77,6,
					206,335,103,180,744,2407,31,
					45,10,145,482,2001,932,175,
					307,595,528,663,247,2665};
			Map<String, Integer> all1Map = new HashMap<String, Integer>();
			for(int i=0; i<all1.length; i++)
				all1Map.put(""+all1[i], i);

			BufferedReader br = IO.getBufferedReader(ge3);
			Map<String, Double> dMap = new HashMap<String, Double>();
			String line, str_c;
			String[] s;
			while( (line=br.readLine())!=null ) {
				s = line.split("\\s+");
				str_c = s[1].compareTo(s[2])>0 ? s[2]+":"+s[1] : s[1]+":"+s[2];
				int scaff = Integer.parseInt(s[0]);
				if(scaff==1) {
					int contig1 = all1Map.get(s[1]),
							contig2 = all1Map.get(s[2]);
					if( (contig1<sep1[0]&&contig2<sep1[0]) || 
							(contig1>=sep1[0]&&contig1<sep1[1]&&
							contig2>=sep1[0]&&contig2<sep1[1]) ||
							(contig1>=sep1[1]&&contig2>=sep1[1]) ) {
						dMap.put(str_c, Double.parseDouble(s[6])-
								Double.parseDouble(s[5])+1);
					}
				} else {
					dMap.put(str_c, Double.parseDouble(s[6])-
							Double.parseDouble(s[5])+1);
				}
			}

			br.close();

			br = IO.getBufferedReader(pw);
			while( (line=br.readLine()) !=null ){
				if(line.startsWith("#"))
					continue;
				s = line.split("\\s+");
				str_c = s[0].compareTo(s[1])>0 ? s[1]+":"+s[0] : s[0]+":"+s[1];
				if(dMap.get(str_c)!=null)
					System.out.println(str_c+"\t"+dMap.get(str_c)+"\t"+s[2]);
			}
			br.close();



			System.out.println("Itrk_sc0000001.1".replaceAll("^Itrk{0,1}_sc0{0,9}", ""));
			System.out.println(50*Math.log(1/(1-2*0.31)));

			if(args.length==3) {
				br = new BufferedReader(
						new FileReader(args[0]));
				Map<String, Integer> size = new HashMap<String, Integer>();
				while( (line=br.readLine())!=null ) {
					if(line.startsWith("#"))
						continue;
					s = line.split("\\s+");
					size.put(s[0], Integer.parseInt(s[1]));
				}
				br.close();
				br = new BufferedReader(new FileReader(args[1]));
				BufferedWriter bw = new BufferedWriter(new FileWriter(args[2]));
				while( (line=br.readLine())!=null ) {
					s = line.split("\\s+");
					String id = "Trifida_";
					for(int i=0; i<5-s[0].length(); i++)
						id += "0";
					id += s[0];
					bw.write(id+"\t"+s[1]+"\t"+size.get(id)+"\t"+s[0]+"\n");
				}
				br.close();
				bw.close();
			}

			if(args.length==4) {
				br = new BufferedReader(
						new FileReader(args[0]));
				Map<String, Integer> size = new HashMap<String, Integer>();
				while( (line=br.readLine())!=null ) {
					if(line.startsWith("#"))
						continue;
					s = line.split("\\s+");
					size.put(s[0].replace("chr", ""), Integer.parseInt(s[1]));
				}
				br.close();
				br = new BufferedReader(new FileReader(args[1]));
				Map<Integer, String> ids = new HashMap<Integer, String>();
				int k=0;
				while( (line=br.readLine())!=null ) {
					if(line.startsWith("#"))
						continue;
					ids.put(k++, line);
				}
				br.close();
				br = new BufferedReader(new FileReader(args[2]));
				BufferedWriter bw = new BufferedWriter(new FileWriter(args[3]));
				while( (line=br.readLine())!=null ) {
					s = line.split("\\s+");
					bw.write(ids.get(Integer.parseInt(s[0]))+"\t"+
							s[1]+"\t"+size.get(s[0])+"\t"+s[0]+"\n");
				}
				br.close();
				bw.close();
			}
		}

	}
}
