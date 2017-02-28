package cz1.hmm.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import cz1.util.ArgsEngine;
import cz1.util.Executor;

public class PseudoMoleculeConstructor extends Executor {
	
	private String assembly_file = null;
	private String out_file = null;
	private String mct_file = null;
	private double genome_size;
	private final Map<String, Scaffold> scaffolds = new HashMap<String, Scaffold>();
	private final Map<String, List<Scaffold>> pseudo_molecule = 
			new HashMap<String, List<Scaffold>>();
	
	public PseudoMoleculeConstructor(String assembly_file, 
			String mct_file,
			String out_file) {
		this.assembly_file = assembly_file;
		this.mct_file = mct_file;
		this.out_file = out_file;
		this.assemblyReader();
		this.geneticMapReader();
	}
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " Common:\n"
						+ "		-i/--input-mct				Input genetic linkage map file.\n"
						+ "		-a/--input-assembly			Input assembly fasta file.\n"
						+ "		-o/--prefix					Output pseudomolecule file.\n"
						);
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
			myArgsEngine.add("-i", "--input-mct", true);
			myArgsEngine.add("-a", "--input-assembly", true);
			myArgsEngine.add("-o", "--prefix", true);
		}
		
		if(myArgsEngine.getBoolean("-i")) {
			mct_file = myArgsEngine.getString("-i");
			this.geneticMapReader();
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your input genetic linkage map file.");
		}
		
		if(myArgsEngine.getBoolean("-a")) {
			assembly_file = myArgsEngine.getString("-a");
			this.assemblyReader();
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your assembly fasta file.");
		}

		if(myArgsEngine.getBoolean("-o")) {
			out_file = myArgsEngine.getString("-o");
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify your output file name.");
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		
		/**
		double bp_all=0, cm_all=0;
		for(String scf : scaffolds.keySet()) {
			scaffold = scaffolds.get(scf);
			bp_all += scaffold.size;
			cm_all += scaffold.gL;
		}
		double units = bp_all/cm_all;
		myLogger.info(units);
		List<Double> bp_cm = new ArrayList<Double>();
		for(String scf : scaffolds.keySet()) {
			scaffold = scaffolds.get(scf);
			if(scaffold.gL>0)
				bp_cm.add(scaffold.size/scaffold.gL);
		}
		Double[] unit_all = new Double[bp_cm.size()];
		bp_cm.toArray(unit_all);
		double sum = 0;
		for(int i=0; i<unit_all.length; i++) sum+=unit_all[i];
		double units2 = sum/unit_all.length;
		double rms = 0;
		for(int i=0; i<unit_all.length; i++) 
			rms+=Math.pow(unit_all[i]-units2,2);
		double unitsD = Math.sqrt(rms/(unit_all.length-1));
		myLogger.info(units2+" "+unitsD);
		**/
		
		double units = (genome_size-pseudoPhysicalSize())/pseudoGapSize();
		myLogger.info(units);
		Scaffold scaffold;
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(out_file));
			for(String lgid : pseudo_molecule.keySet()) {
				StringBuilder oos =  new StringBuilder();
				List<Scaffold> pm = pseudo_molecule.get(lgid);
				for(int i=0; i<pm.size(); i++) {
					scaffold = pm.get(i);
					int N = (int) (scaffold.gD*units);
					for(int j=0; j<N; j++) oos.append('N');
					if(scaffold.reverse)
						oos.append(scaffold.reverseSeq());
					else oos.append(scaffold.sequence);
				}
				bw.write(">"+lgid+"\n");
				bw.write(wrap(oos.toString(),50));
			}
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private class Scaffold {
		private final String name;
		private final String sequence;
		private final int size;
		private boolean reverse = false;
		private double gL = 0;
		private double gD = 0;
		
		public Scaffold(String name, String sequence) {
			this.name = name;
			this.sequence = sequence;
			this.size = this.sequence.length();
		}
		
		public void reverse() {
			this.reverse = !this.reverse;
		}
		
		public void setGL(double gL) {
			this.gL = gL;
		}
		
		public void setGD(double gD) {
			this.gD = gD;
		}
		
		public String reverseSeq() {
			return new StringBuilder(this.sequence).
					reverse().toString();
		}
		
		@Override
		public int hashCode() {
			return this.name.hashCode();
		}
		
		@Override
		public boolean equals(Object obj) {
			if(this==obj) return true;
			if(obj==null) return false;
			if(getClass()!=obj.getClass())
				return false;
			Scaffold scaff = (Scaffold) obj;
			return this.name.equals(scaff.name);
		}
	}
	
	public String toFasta(String name,
			String sequence, 
			int linewidth) {
		return ">"+name+
				"\n"+wrap(sequence, linewidth);
	}

	private String wrap(String sequence, 
			int linewidth) {
		// TODO Auto-generated method stub
		StringBuilder sb = new StringBuilder();
		int size = sequence.length();
		int n = size/linewidth, 
				residual=size%linewidth;
		if(n>0) sb.append(sequence.subSequence(0, 
				linewidth));
		for(int i=1; i<n; i++) {
			sb.append("\n");
			sb.append(sequence.subSequence(i*linewidth, 
					(i+1)*linewidth));
		}
		if(residual>0) {
			if(n>0) sb.append("\n");
			sb.append(sequence.subSequence(
					n*linewidth, size));
		}
		return sb.toString();
	}
	
	
	public void assemblyReader() {
		try {
			BufferedReader br = new BufferedReader(new FileReader(assembly_file));
			String line=br.readLine();
			String name;
			
			while( line!=null ) {
				if(line.startsWith(">")) {
					name = line.split("\\s+")[0].substring(1).
							replace('|', '_').replace('.', '_');
					StringBuilder sequence = new StringBuilder();
					while( (line=br.readLine())!=null &&
							!line.startsWith(">"))
						sequence.append(line);
					scaffolds.put(name, new Scaffold(name, 
							sequence.toString()));
					genome_size += sequence.length();
				} else line=br.readLine();
			}
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void geneticMapReader() {
		try {
			BufferedReader br = new BufferedReader(new FileReader(mct_file));
			String line = br.readLine();
			int lg = 0;
			String[] s;
			Scaffold scaffold;
			String name, scaffold_name;
			double start, end;
			while( line!=null ) {
				if(line.startsWith("group")) {
					List<Scaffold> scaff_list = new ArrayList<Scaffold>();
					double cm = 0;
					name = "LG"+(++lg);
					while( (line=br.readLine())!=null && 
							!line.startsWith("group") &&
							line.length()!=0) {
						s = line.split("\\(|\\)\\s+");
						scaffold_name = s[0];
						scaffold = scaffolds.get(scaffold_name);
						if(s[1].equals("-")) scaffold.reverse();
						start = Double.parseDouble(s[2]);
						s = br.readLine().split("\\(|\\)\\s+");
						end = Double.parseDouble(s[2]);
						scaffold.setGL(end-start);
						scaffold.setGD(start-cm);
						cm = end;
						scaff_list.add(scaffold);
					}
					pseudo_molecule.put(name, scaff_list);
					if(line.length()!=0)
						line = br.readLine();
				} else line=br.readLine();
			}
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public double pseudoPhysicalSize() {
		Scaffold scaffold;
		double bp_all = 0.0;
		for(String lgid : pseudo_molecule.keySet()) {
			List<Scaffold> pm = pseudo_molecule.get(lgid);
			for(int i=0; i<pm.size(); i++) {
				scaffold = pm.get(i);
				bp_all += scaffold.size;
			}
		}
		return bp_all;
	}
	
	public double pseudoGapSize() {
		Scaffold scaffold;
		double cm_all = 0.0;
		for(String lgid : pseudo_molecule.keySet()) {
			List<Scaffold> pm = pseudo_molecule.get(lgid);
			for(int i=0; i<pm.size(); i++) {
				scaffold = pm.get(i);
				cm_all += scaffold.gD;
			}
		}
		return cm_all;
	}
	
	
	public double coverage() {
		return pseudoPhysicalSize()/genomeSize();
	}
	
	public double genomeSize() {
		return this.genome_size;
	}

	public void setAssemblyError(Map<String, int[][]> errs) {
		// TODO Auto-generated method stub
		for(String errScaff : errs.keySet()) {
			int[][] err = errs.get(errScaff);
			String dnaSEQ = this.scaffolds.get(errScaff).sequence;
			int[][] chunk = findBPS(dnaSEQ, err);
			for(int i=0; i<chunk.length; i++)
				this.scaffolds.put(errScaff+"_"+(i+1),
						new Scaffold(errScaff+"_"+(i+1),
								dnaSEQ.substring(chunk[i][0], 
										chunk[i][1])));
		}
	}

	private int[][] findBPS(String errScaff, int[][] err) {
		// TODO Auto-generated method stub
		int[][] chunk = new int[err.length+1][2];
		for(int i=0; i<err.length; i++) {
			String seq = errScaff.substring(err[i][0], 
					err[i][1]);
			String[] gap = seq.split("[CAGT]+");
			if(gap.length==0) {
				// no gap found, chunk in between is discarded
				chunk[i][1] = err[i][0];
				chunk[i+1][0] = err[i][1];
			} else {
				// gaps found, select the longest one as break point
				String longest_gap = "";
				for(int j=0; j<gap.length; j++)
					if(longest_gap.length()<gap[j].length())
						longest_gap = gap[j];
				int star = seq.indexOf(longest_gap);
				chunk[i][1] = err[i][0]+star;
				chunk[i+1][0] = err[i][0]+longest_gap.length();
			}
		}
		chunk[err.length][1] = errScaff.length();
		return chunk;
	}
}



