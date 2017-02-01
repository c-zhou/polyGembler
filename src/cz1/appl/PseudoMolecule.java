package cz1.appl;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class PseudoMolecule {
	
	private final String name;
	private final String sequence;
	private final int size;
	
	public static void main(String[] args) {
		System.out.println(Integer.MAX_VALUE);
		construct(args[0], args[1], Integer.parseInt(args[2]), args[3]);
	}
	
	private static class Scaffold {
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
	
	public static void construct(String mct, String scaff, 
			int estimated_genome_size, String out) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(scaff));
			String line=br.readLine();
			String name;
			Map<String, Scaffold> scaffolds = new HashMap<String, Scaffold>();
			while( line!=null ) {
				if(line.startsWith(">")) {
					System.out.println(line);
					name = line.split("\\s+")[0].substring(1).
							replaceAll("Itr_sc0{0,9}", "").
							replaceAll("\\.1$", "");
					StringBuilder sequence = new StringBuilder();
					while( (line=br.readLine())!=null &&
							!line.startsWith(">"))
						sequence.append(line);
					scaffolds.put(name, new Scaffold(name, 
							sequence.toString()));
				} else line=br.readLine();
			}
			br.close();
			
			br = new BufferedReader(new FileReader(mct));
			line = br.readLine();
			int lg = 0;
			String[] s;
			Scaffold scaffold;
			String scaffold_name;
			double start, end;
			Map<String, List<Scaffold>> pseudo_molecule = 
					new HashMap<String, List<Scaffold>>();
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
			
			//calculate average length in base pairs for 1cM
			double bp_all=0, cm_all=0;
			for(String lgid : pseudo_molecule.keySet()) {
				List<Scaffold> pm = pseudo_molecule.get(lgid);
				for(int i=0; i<pm.size(); i++) {
					scaffold = pm.get(i);
					bp_all += scaffold.size;
					cm_all += scaffold.gD;
				}
			}
			double units = (estimated_genome_size-bp_all)/cm_all;
			System.out.println(bp_all+"\t"+units);
			
			/**
			double bp_all=0, cm_all=0;
			for(String scf : scaffolds.keySet()) {
				scaffold = scaffolds.get(scf);
				bp_all += scaffold.size;
				cm_all += scaffold.gL;
			}
			double units = bp_all/cm_all;
			System.out.println(units);
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
			System.out.println(units2+" "+unitsD);
			**/
			
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
				BufferedWriter bw = new BufferedWriter(new FileWriter(out+"/"+lgid+".fa"));
				bw.write(">"+lgid+"\n");
				bw.write(wrap(oos.toString(),50));
				bw.close();
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public PseudoMolecule(String name, String sequence) {
		this.name = name;
		this.sequence = sequence;
		this.size = this.sequence.length();
	}
	
	public String toFasta() {
		return ">"+this.name+
				"\n"+wrap(this.sequence, 50);
	}
	
	public String toFasta(int linewidth) {
		return ">"+this.name+
				"\n"+wrap(this.sequence, linewidth);
	}
	
	public static String toFasta(String name,
			String sequence, 
			int linewidth) {
		return ">"+name+
				"\n"+wrap(sequence, linewidth);
	}

	private static String wrap(String sequence, 
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
}
