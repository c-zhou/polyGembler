package cz1.tenx.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

public class FragmentUtils extends Executor {

	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+" -f/--fragment                Input fragment file.\n"
						+" -v/--vcf                     Input VCF file.\n"
						+" -r/--range                   Data range (chr/chr:start-end).\n"
						+" -c/--convert                 Type to convert.\n"
						+" -o/--out                     Output file.\n\n"
				);
	}

	private static enum Type {hapcompass};
	private String fragFile = null;
	private String vcfFile = null;
	private String outFile  = null;
	private Type convt = null;
	private String rangeChr;
	private int rangeLowerBound;
	private int rangeUpperBound;
	
	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		// create the command line parser

		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}

		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-f", "--fragment", true);
			myArgsEngine.add("-v", "--vcf", true);
			myArgsEngine.add("-r", "--range", true);
			myArgsEngine.add("-c", "--convert", true);
			myArgsEngine.add("-o", "--out", true);
			myArgsEngine.parse(args);
		}

		if(myArgsEngine.getBoolean("-f")) {
			this.fragFile = myArgsEngine.getString("-f");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your input fragment data file.");
		}

		if(myArgsEngine.getBoolean("-v")) {
			this.vcfFile = myArgsEngine.getString("-v");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your input VCF file.");
		}
		
		if(myArgsEngine.getBoolean("-r")) {
			this.setDataRange(myArgsEngine.getString("-r"));
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your data range.");
		}

		if(myArgsEngine.getBoolean("-c")) {
			switch(myArgsEngine.getString("-c").toLowerCase()) {
			case "hapcompass":
				this.convt = Type.hapcompass;
				break;
			default:
				throw new RuntimeException("!!!");
			}
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify the file type.");
		}

		if(myArgsEngine.getBoolean("-o")) {
			this.outFile = myArgsEngine.getString("-o");
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify your output file name.");
		}
	}

	private int var_start = Integer.MAX_VALUE, var_end = Integer.MIN_VALUE;
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		switch(this.convt) {
		case hapcompass:
			this.hapcompass();
			break;
		default:
			throw new IllegalArgumentException("!!!");
		}
	}
	
	private void hapcompass() {	
		// TODO Auto-generated method stub
		try {
			BufferedReader br1 = Utils.getBufferedReader(vcfFile);
			String line;
			String[] s;
			int position, index = 0;
			while( (line=br1.readLine())!=null ){
				if(line.startsWith("#")) continue;
				++index; // update variant index
				s = line.split("\\s+");
				if(!s[0].equals(rangeChr)) continue;
				position = Integer.parseInt(s[1]);
				if(position<rangeLowerBound) continue;
				if(position>rangeUpperBound) break;
				if(var_start>index) var_start = index;
				if(var_end  <index) var_end   = index;
			}
			br1.close();
		
			BufferedWriter bw  = Utils.getBufferedWriter(outFile);
			BufferedReader br2 = Utils.getBufferedReader(fragFile);
			int n, k, starti;
			String qual;
			String chr;
			while( (line=br2.readLine())!=null ) {
				s = line.split("\\s+");
				if(!s[1].split(":")[0].equals(rangeChr)) continue;
				n = s.length;
				if(Integer.parseInt(s[n-3]+s[n-2].length())<var_start) 
					continue;
				if(Integer.parseInt(s[5])>var_end) break;
				qual = s[n-1];
				chr  = s[1].split(":")[0];
				bw.write("1\t"+s[1]+"\t"+chr);
				k=0;
				for(int i=5; i<n-2; i+=2) {
					starti = Integer.parseInt(s[i])-1;
					for(char c : s[i+1].toCharArray()) {
						if(starti>=var_start&&starti<=var_end) {
							bw.write(" "+starti);
							bw.write(" "+(c-'0'));
							bw.write(" "+qual.charAt(k));
						}
						++starti;
						++k;
					}
				}
				bw.write("\n");
			}
			br2.close();
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void setDataRange(String dat_range) {
		// TODO Auto-generated method stub
		String[] s = dat_range.split(":");
		rangeChr = s[0];
		if(s.length==2) {
			s = s[1].split("-");
			rangeLowerBound = Integer.parseInt(s[0]);
			rangeUpperBound = Integer.parseInt(s[1]);
			if(rangeLowerBound>rangeUpperBound) 
				throw new RuntimeException("invalid data range!!!");
		} else {
			rangeLowerBound = Integer.MAX_VALUE;
			rangeUpperBound = Integer.MIN_VALUE;
			try {
				BufferedReader br = Utils.getBufferedReader(this.vcfFile);
				String line;
				int i;
				while( (line=br.readLine())!=null ){
					if(line.startsWith("#")) continue;
					s = line.split("\\s+");
					if(!s[0].equals(rangeChr)) continue;
					if( (i=Integer.parseInt(s[1]))<rangeLowerBound )
						rangeLowerBound = i;
					if( (i=Integer.parseInt(s[1]))>rangeUpperBound )
						rangeUpperBound = i;
				}
				br.close();		
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
}
