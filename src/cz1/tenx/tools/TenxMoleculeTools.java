package cz1.tenx.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

public class TenxMoleculeTools extends Executor {

	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--in             Input file.\n"
						+ " -t/--threads        Threads.\n"
						+ " -o/--out            Output file.\n\n");	
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-i", "--in", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-o", "--out", true);
			myArgsEngine.parse(args);
		}
		if (myArgsEngine.getBoolean("-i")) {
			this.file_in = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the input file.");
		}
		
		if (myArgsEngine.getBoolean("-t")) {
			this.THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}
		
		if (myArgsEngine.getBoolean("-o")) {
			this.file_out = myArgsEngine.getString("-o");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the output file.");
		}
	}

	private String file_in = null;
	private String file_out = null;
	private final static Object lock = new Object();
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		Map<String, Molecule[]> mols = loadMolecule();
		BufferedWriter bw_out = Utils.getBufferedWriter(file_out);
		this.initial_thread_pool();
		for(String key : mols.keySet()) 
			this.executor.submit(new Runnable() {
				private Molecule[] mols;
				
				@Override
				public void run() {
					// TODO Auto-generated method stub
					try {
						// sort by alignment positions
						Arrays.sort(this.mols, new Comparator<Molecule>() {
							@Override
							public int compare(Molecule m1, Molecule m2) {
								// TODO Auto-generated method stub
								return  m1.chr_sn1.equals(m2.chr_sn1) ?
									  ( m1.chr_sn2.equals(m2.chr_sn2) ?
									  (	m1.s1==m2.s1  ? 
									  ( m1.e1==m2.e1  ? 
									  ( m1.s2==m2.s2  ? 
								      ( m1.e2-m2.e2 ) : 
								        m1.s2-m2.s2 ) : 
								        m1.e1-m2.e1 ) : 
								        m1.s1-m2.s1 ) :
								        m1.chr_sn2.compareTo(m2.chr_sn2)) :
								        m1.chr_sn1.compareTo(m2.chr_sn1)  ;
							}
						});
						
						// merge
						final int sz = mols.length;
						for(int i=0; i<sz; i++) {
							for(int j=i+1; j<sz; j++) {
								if( merge(i, j) ) {
									--i;
									break;
								}
							}
						}
						
						// write out
						synchronized(lock) {
							for(Molecule mol : this.mols) {
								if(mol!=null) 
									bw_out.write(mol.chr_sn1+":"+mol.s1+"-"+mol.e1+"\t"+mol.chr_sn2+":"+mol.s2+"-"+mol.e2+"\n");
							}
						}
					} catch (Exception e) {
						Thread t = Thread.currentThread();
						t.getUncaughtExceptionHandler().uncaughtException(t, e);
						e.printStackTrace();
						executor.shutdown();
						System.exit(1);
					}
				}

				private boolean merge(final int i, final int j) {
					// TODO Auto-generated method stub
					if(mols[i]==null||mols[j]==null||!overlap(mols[i],mols[j])) 
						return false;
					mols[i].s1 = Math.min(mols[i].s1, mols[j].s1);
					mols[i].e1 = Math.max(mols[i].e1, mols[j].e1);
					mols[i].s2 = Math.min(mols[i].s2, mols[j].s2);
					mols[i].e2 = Math.max(mols[i].e2, mols[j].e2);
					mols[j] = null;
					return true;
				}

				private boolean overlap(Molecule molecule, Molecule molecule2) {
					// TODO Auto-generated method stub
					System.out.println(molecule.chr_sn1.equals(molecule2.chr_sn1));
					System.out.println( molecule.chr_sn2.equals(molecule2.chr_sn2));
					System.out.println( overlap(molecule.s1, molecule.e1, molecule2.s1, molecule2.e1));
					System.out.println(overlap(molecule.s2, molecule.e2, molecule2.s2, molecule2.e2));
					return molecule.chr_sn1.equals(molecule2.chr_sn1) &&
						   molecule.chr_sn2.equals(molecule2.chr_sn2) &&
						   overlap(molecule.s1, molecule.e1, molecule2.s1, molecule2.e1) && 
						   overlap(molecule.s2, molecule.e2, molecule2.s2, molecule2.e2);
				}

				private boolean overlap(int s1, int e1, int s2, int e2) {
					// TODO Auto-generated method stub
					return ((double)s1-e2)*(s2-e1)>=0;
				}

				public Runnable init(Molecule[] molecules) {
					// TODO Auto-generated method stub
					this.mols = molecules;
					return this;
				}
				
			}.init(mols.get(key)));
		try {
			this.waitFor();
			bw_out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private Map<String, Molecule[]> loadMolecule() {
		// TODO Auto-generated method stub
		Map<String, List<Molecule>> mols = new HashMap<String, List<Molecule>>();
		try {
			BufferedReader br = Utils.getBufferedReader(file_in);
			String line;
			String[] s, s1, s2;
			while( (line = br.readLine())!=null ) {
				s = line.split("\\s+");
				s1 = s[0].split(":|-");
				s2 = s[1].split(":|-");
				if(!mols.containsKey(s1[0])) mols.put(s1[0], new ArrayList<Molecule>());
				mols.get(s1[0]).add(new Molecule(s1[0],Integer.parseInt(s1[1]), Integer.parseInt(s1[2]),
							          s2[0],Integer.parseInt(s2[1]), Integer.parseInt(s2[2])) );
			}
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		Map<String, Molecule[]> mols_arr = new HashMap<String, Molecule[]>();
		for(String key : mols.keySet())
			mols_arr.put(key, mols.get(key).toArray(new Molecule[mols.get(key).size()]));
		return mols_arr;
	}

	private final class Molecule {
		private String chr_sn1, chr_sn2;
		private int s1, s2, e1, e2;
		
		public Molecule(String chr_sn1,
				int s1,
				int e1,
				String chr_sn2,
				int s2,
				int e2) {
			this.chr_sn1 = chr_sn1;
			this.s1 = s1;
			this.e1 = e1;
			this.chr_sn2 = chr_sn2;
			this.s2 = s2;
			this.e2 = e2;
		}
	}
	
	public static void main(String[] args) {
		TenxMoleculeTools tenx = new TenxMoleculeTools();
		tenx.setParameters(new String[] {
				"-i", "c://users//chenxi.zhou//desktop//putty//440166_tanzania_0.9_collapsed4.txt",
				"-o", "c://users//chenxi.zhou//desktop//tt.txt"		
		});
		tenx.run();
	}
}

















