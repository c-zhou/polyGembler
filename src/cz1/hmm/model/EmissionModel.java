package cz1.hmm.model;

import java.io.BufferedOutputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.log4j.Logger;

import cz1.hmm.data.DataEntry;
import cz1.math.Combination;
import cz1.math.SaddlePointExpansion;
import cz1.util.Constants;
import cz1.util.Constants.Field;
import cz1.util.Utils;

public abstract class EmissionModel {
	private final static Logger myLogger = Logger.getLogger(EmissionModel.class);
	
	protected final static double mu_A_e = 10;
	protected final static double mu_A_m = 0.1;
	
	// at least 3 markers or 30% markers to keep a sample
	// at least 30 f1 progeny to run the program
	protected final double min_mf = 0.3;
	protected final double min_mn = 3;
	protected final int min_f1 = 30;
	
	protected static int iteration = 0;
	
	protected Field field;
	protected DataEntry de;

	protected int M; // #markers
	protected int N; // #individuals
	protected int H; // #founder haplotypes
	protected int K; // #compound hidden states
	protected int Nf1; // #f1 progeny
	protected String[] samples;
	protected String[] parents;
	protected boolean[] fi1ter;
	protected int[] parents_i;
	protected int[] progeny_i;
	protected double[] distance;
	protected double[] weights;
	
	protected List<Integer[]> sspace; // state space for each sample
	
	protected StateUnit state;
	protected ObUnit[][] obs;
	protected EmissionUnit[] emission;
	protected PathUnit[] pas;
	
	protected final List<Integer> conjs = new ArrayList<>(); // conjunctive position
	
	protected boolean logspace;
	protected String[] chrs;
	protected boolean[] chrs_rev;
	
	abstract double loglik();
	abstract double loglik(int fromMIndex, int toMIndex);
	abstract double findPath();
	abstract void write(String output, 
			String experiment, 
			String contig);
	
	public EmissionModel(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse,
			Field field,
			int ploidy,
			String[] parents,
			boolean logspace) {
		this.field = field;
		this.de = this.catDE(de, seperation, reverse);
		this.logspace = logspace;
		this.H = ploidy;
		this.parents = parents;
		this.initialise();
	}
	
	public EmissionModel(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse,
			Field field,
			int ploidy,
			String[] parents) {
		this(de, seperation, reverse, field, ploidy, parents, true);
	}
	
	public EmissionModel() {
		// TODO Auto-generated constructor stub
		super();
	}
	
	protected String[] true_id;
	protected int[] true_pos;
	protected int[] Ms;
	
	protected DataEntry catDE(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse) {
		// TODO Auto-generated method stub
		int M = 0;
		chrs = new String[de.length];
		chrs_rev = new boolean[de.length];
		System.arraycopy(reverse, 0, chrs_rev, 0, reverse.length);
		Ms = new int[de.length];
		for(int i=0; i<de.length; i++) {
			M += de[i].modelLength();
			chrs[i] = de[i].getId();
			Ms[i] = de[i].modelLength();
		}
		true_id = new String[M];
		true_pos = new int[M];
		for(int i=0; i<de.length; i++)
			if(reverse[i]) de[i].reverse();
		int m = 0;
		int[] pos = de[0].getIntegerPosition();
		Arrays.fill(true_id, 0, pos.length, de[0].getId());
		System.arraycopy(pos, 0, true_pos, m, pos.length);
		m += pos.length;
		for(int i=1; i<de.length; i++) {
			conjs.add(de[0].modelLength());
			de[0].addAll(de[i], seperation[i-1]);
			pos = de[i].getIntegerPosition();
			Arrays.fill(true_id, m, m+pos.length, de[i].getId());
			System.arraycopy(pos, 0, true_pos, m, pos.length);
			m += pos.length;
		}
		return de[0];
	}

	protected void initialise() {
		// TODO Auto-generated method stub
		this.samples = de.getSample();
		this.M = this.de.getAllele().size();
		this.N = this.samples.length;
		this.state = new StateUnit(H);
		this.K = state.hsc.length;
		
		List<Integer> pars = new ArrayList<Integer>();
		List<Integer> pros = new ArrayList<Integer>();
		for(int i=0; i<samples.length; i++) {
			boolean isprogeny = true;
			for(int j=0; j<parents.length; j++) {
				if(samples[i].equals(parents[j])) {
					pars.add(i);
					isprogeny = false;
					break;
				}
			}
			if(isprogeny) pros.add(i);
		}
		Nf1 = N-pars.size();
		this.parents_i = new int[2];
		Arrays.fill(parents_i, -1);
		if(pars.size()>2) throw new RuntimeException("More than TWO parental sample detected!!!");
		for(int i=0; i<pars.size(); i++)
			parents_i[i] = pars.get(i);
		if(pros.isEmpty()) throw new RuntimeException("No progeny sample provided!!!");
		this.progeny_i = ArrayUtils.toPrimitive(pros.toArray(new Integer[pros.size()]));
		this.fi1ter = new boolean[N];
		double[] position = de.getPosition();
		this.distance = new double[this.M-1];
		for(int i=0; i<distance.length; i++)
			distance[i] = Math.abs(position[i+1]-position[i]);
		this.sspace = new ArrayList<>(N);
		for(int i=0; i<N; i++) sspace.add(null);
		if(parents_i[0]!=-1) sspace.set(parents_i[0], new Integer[]{0});
		if(parents_i[1]!=-1) sspace.set(parents_i[1], new Integer[]{1});
		Integer[] progeny_s = new Integer[K-2];
		for(int i=0; i<K-2; i++) progeny_s[i] = i+2;
		for(int i : progeny_i) sspace.set(i, progeny_s);
		this.weights = this.field==Field.GT ? new double[] {Nf1/2.0, 1.0} : new double[]{1.0, 1.0};
		this.makeObUnits();
		this.makePathUnits();
		this.makeEmissionUnits();
	}

	private final int allele_d = 10; 
	private void makeObUnits() {
		// TODO Auto-generated method stub
		this.obs = new ObUnit[N][M];
		int miss_cnt = 0, miss_f1 = 0;
		double miss_max = M-Math.max(min_mn, M*min_mf);
		switch(this.field) {
		case AD:
			List<List<int[]>> ad = this.de.getAlleleDepth();
			if(ad==null) throw new RuntimeException("AD feild not available!!! Try GT (-G/--genotype) options.");
			int[] dp;
			for(int i=0; i<N; i++) {
				miss_cnt = 0;
				for(int j=0; j<M; j++) {
					dp = ad.get(j).get(i);
					obs[i][j] = new ObUnit(dp[0]+dp[1], dp[0], K);
					if(dp[0]+dp[1]==0) ++miss_cnt;
				}
				if(miss_cnt>miss_max) {
					fi1ter[i] = true;
					obs[i] = null;
					if(i!=parents_i[0]&&i!=parents_i[1]) ++miss_f1;
				}
			}
			break;
		case GT:
			List<List<String[]>> gt = this.de.getGenotype();
			if(gt==null) throw new RuntimeException("GT field not available!!! Try AD (-D/--allele-depth) option.");
			List<String[]> allele = this.de.getAllele();
			int acnt, bcnt;
			String[] a, g;
			for(int i=0; i<N; i++) {
				miss_cnt = 0;
				for(int j=0; j<M; j++) {
					a = allele.get(j);
					g = gt.get(j).get(i);
					acnt = 0;
					bcnt = 0;
					for(int k=0; k<H; k++) {
						acnt += (g[k].equals(a[0]) ? 1 : 0);
						bcnt += (g[k].equals(a[1]) ? 1 : 0);
					}
					acnt *= allele_d;
					bcnt *= allele_d;
					obs[i][j] = new ObUnit(acnt+bcnt, acnt, K);
					if(acnt+bcnt==0) ++miss_cnt;
				}
				if(miss_cnt>miss_max) {
					fi1ter[i] = true;
					obs[i] = null;
					if(i!=parents_i[0]&&i!=parents_i[1]) ++miss_f1;
				}
			}
			break;
		default:
			throw new RuntimeException("!!!");
		}
		
		int f1 = Nf1-miss_f1;
		if(f1<min_f1) {
			myLogger.info("f1 individuals left: "+f1+" (<"+min_f1+")");
			myLogger.info("program exit with status 0.");
			System.exit(0);
		} else {
			myLogger.info("f1 individuals left: "+f1);
		}
	}

	private void makeEmissionUnits() {
		// TODO Auto-generated method stub
		this.emission = new EmissionUnit[M];
		for(int i=0; i<M; i++) 
			emission[i] = new EmissionUnit(de.getAllele().get(i),
					H*2,
					K,
					bfrac(i));
	}
	
	private void makePathUnits() {
		// TODO Auto-generated method stub
		this.pas = new PathUnit[N];
		for(int i=0; i<N; i++) 
			pas[i] = fi1ter[i] ? null:new PathUnit();
	}
	
	protected void switchNumericalSpace() {
		// TODO Auto-generated method stub
		switchNumericalSpace(!logspace);
		logspace = !logspace;
	}
	
	protected void switchNumericalSpace(boolean logspace) {
		// TODO Auto-generated method stub
		if(logspace==this.logspace)
			return;
		for(int i=0; i<N; i++) {
			if(fi1ter[i]) continue;
			for(int j=0; j<M; j++) 
				obs[i][j].switchNumericSpace();
		}
		this.logspace = logspace;
	}
	
	private double bfrac(final int i) {
		// TODO Auto-generated method stub
		double baf = 0.5, acnt = 0, bcnt = 0;
		switch(this.field) {
		case AD:
			if(this.de.getAlleleDepth()==null)
				throw new RuntimeException("AD feild not available!!! Try GT (-G/--genotype) options.");
			List<int[]> ad = this.de.getAlleleDepth().get(i);
			for(int j=0; j<N; j++) {
				if(fi1ter[j]) continue;
				int[] aa = ad.get(j);
				acnt += aa[0];
				bcnt += aa[1];
			}
			break;
		case GT:
			if(this.de.getGenotype()==null)
				throw new RuntimeException("GT field not available!!! Try AD (-D/--allele-depth) option.");
			List<String[]> gt = this.de.getGenotype().get(i);
			String[] allele = this.de.getAllele().get(i);
			for(int j=0; j<N; j++) {
				if(fi1ter[j]) continue;
				String[] g = gt.get(j);
				for(int k=0; k<H; k++) {
					acnt += (g[k].equals(allele[0]) ? 1 : 0);
					bcnt += (g[k].equals(allele[1]) ? 1 : 0);
				}
			}
			break;
		default:
			throw new RuntimeException("!!!");
		}
		if(acnt!=0&&bcnt!=0) baf = bcnt/(acnt+bcnt);
		return baf;
	}

	protected void refresh() {
		// refresh emission probabilities for obs
		double[] emissc;
		for(int i=0; i<M; i++) {
			emissc = emission[i].emissc;
			for(int j=0; j<N; j++) {
				if(fi1ter[j]) continue;
				obs[j][i].updateEmiss(emissc);
			}
		}
	}
	
	protected static void clear(double[][][] matrix) {
		// TODO Auto-generated method stub
		int r = matrix.length, 
				c = matrix[0].length;
		for(int i=0; i<r; i++)
			for(int j=0; j<c; j++)
				Arrays.fill(matrix[i][j], 0);

	}

	protected static void clear(double[][] matrix) {
		// TODO Auto-generated method stub
		int r = matrix.length;
		for(int i=0; i<r; i++)
			Arrays.fill(matrix[i], 0);
	}

	protected static void clear(int[][] matrix) {
		// TODO Auto-generated method stub
		int r = matrix.length;
		for(int i=0; i<r; i++)
			Arrays.fill(matrix[i], 0);
	}
	
	protected static void clear(double[] array) {
		// TODO Auto-generated method stub
		Arrays.fill(array, 0);
	}
	
	public int iteration() {
		// TODO Auto-generated method stub
		return iteration;
	}
	
	protected class ObUnit {
		private final int cov;	// depth of coverage
		private final int aa;	// A-allele depth
		protected final double[] emiss; // place holder for emission probs
		protected final double[] logscale;
		
		public ObUnit(int cov, int aa, int k) {
			// TODO Auto-generated constructor stub
			this.cov = cov;
			this.aa = aa;
			this.emiss = new double[k];
			this.logscale = new double[3];
		}

		public double getLogScale(final int i) {
			// TODO Auto-generated method stub
			if(i==parents_i[0])
				return logscale[0];
			else if(i==parents_i[1])
				return logscale[1];
			else
				return logscale[2];
		}

		public void updateEmiss(double[] emissA) {
			if(emiss.length!=emissA.length) 
				throw new RuntimeException("!!!");
			if(cov==0) {
				Arrays.fill(emiss, 1.0/emiss.length);
				return;
			}
			for(int i=0; i<emiss.length; i++)
				emiss[i] = SaddlePointExpansion.logBinomialProbability(aa, cov, emissA[i]);
			
			if(!logspace) switchToNormalSpace();
		}

		public void switchNumericSpace() {
			// TODO Auto-generated method stub
			if(logspace)
				switchToNormalSpace();
			else
				switchToLogSpace();
		}

		private void switchToNormalSpace() {
			// TODO Auto-generated method stub
			// if not in log space
			logscale[0] = emiss[0];
			emiss[0] = 1.0;
			logscale[1] = emiss[1];
			emiss[1] = 1.0;
			logscale[2] = StatUtils.max(emiss, 2, K-2);
			for(int i=2; i<K; i++) {
				emiss[i] = Math.exp(emiss[i]-logscale[2]);
			}
		}

		private void switchToLogSpace() {
			// TODO Auto-generated method stub
			emiss[0] = logscale[0];
			emiss[1] = logscale[1];
			for(int i=2; i<K; i++) {
				emiss[i] = Math.log(emiss[i])+logscale[2];
			}
			Arrays.fill(logscale, 0);
		}

		public int getCov() {
			return this.cov;
		}

		public int getAa() {
			return this.aa;
		}

		public double[] getEmiss() {
			return this.emiss;
		}
	}

	protected class StateUnit {
		protected final char[] hs; // hidden states
		protected final int[][] hsc; // compound hidden states
		protected final String[] hsc_str; // compound hidden states str
		
		public StateUnit(final int h) {
			hs = makeHs(h*2);
			hsc = makeHsc();
			hsc_str = makeHscStr();
		}
		
		private char[] makeHs(int h) {
			// TODO Auto-generated method stub
			char[] hs = new char[h];
			for(int k=0; k<h; k++) {
				hs[k] = k<10 ? (char)('1'+k) : (char)('a'+k-10);
			}
			return hs;
		}

		private int[][] makeHsc() {
			// TODO Auto-generated method stub
			int mid1 = this.hs.length/2;
			int[] p1 = new int[mid1], p2 = new int[mid1];
			for(int i=0; i<mid1; i++) {
				p1[i] = i;
				p2[i] = mid1+i;
			}

			List<List<Integer>> com1 = Combination.combination(p1, mid1/2),
					com2 = Combination.combination(p2, mid1/2);
			int k = com1.size();
			int k2 = k*k+2;
			int[][] hsc = new int[k2][mid1];
			hsc[0] = p1;
			hsc[1] = p2;
			for(int i=0; i<k; i++) {
				for(int j=0; j<k; j++) {
					List<Integer> com = new ArrayList<Integer>(com1.get(i));
					com.addAll(com2.get(j));
					hsc[i*k+j+2] = ArrayUtils.toPrimitive(com.toArray(new Integer[mid1]));
				}
			}
			return hsc;
		}

		private String[] makeHscStr() {
			// TODO Auto-generated method stub
			int K = this.hsc.length;
			int h = this.hs.length/2;
			String[] hsc_str = new String[K];
			char[] st = new char[h];
			for(int i=0; i<K; i++) {
				for(int j=0; j<h; j++) {
					st[j] = hs[hsc[i][j]];
				}
				hsc_str[i] = StringUtils.join(st, '_');
			}
			return hsc_str;
		}
		
		protected char[] getHs() {
			return this.hs;
		}
		
		protected int[][] getHsc() {
			return this.hsc;
		}
		
		protected String[] getHscStr() {
			return this.hsc_str;
		}
	}

	protected class EmissionUnit {
		protected final String[] allele;
		protected final double[] emiss;
		protected final double[] emissc;
		protected final double bfrac;
		protected final double[][] count;
		protected final double[][][] cnts_prior; // priors for counting
		
		public EmissionUnit(final String[] allele,
				int H,
				int K,
				final double bfrac) {
			this.allele = allele;
			this.emiss = new double[H];
			this.emissc = new double[K];
			this.bfrac = bfrac;
			this.count = new double[H][allele.length];
			this.cnts_prior = new double[K][H/2][2];
			this.prior();
			this.updatec();
		}

		protected void prior() {
			// TODO Auto-generated method stub
			BetaDistribution beta = new BetaDistribution(Constants.rg, 
					(1-bfrac)*mu_A_e, bfrac*mu_A_e);
			for(int i=0; i<emiss.length; i++) {
				emiss[i] = beta.sample();
				if(emiss[i]==0) emiss[i] = 0.001;
				if(emiss[i]==1) emiss[i] = 0.999;
			}
		}

		public void addCount(int s, double acnt, double bcnt) {
			// TODO Auto-generated method stub
			int[] hs = state.hsc[s];
			double[][] prior = cnts_prior[s];
			for(int i=0; i<H; i++) {
				count[hs[i]][0] += acnt*prior[i][0];
				count[hs[i]][1] += bcnt*prior[i][1];
			}
		}
		
		protected void update() {
			// TODO Auto-generated method stub
			for(int i=0; i<emiss.length; i++) {
				emiss[i] = count[i][0]/
				(count[i][0]+count[i][1]);
			}
			this.updatec();
		}

		protected void updatec() {
			// TODO Auto-generated method stub
			int[][] hsc = state.hsc;
			for(int i=0; i<K; i++) {
				double p = 0, p1;
				for(int j=0; j<H; j++)
					p += emiss[hsc[i][j]];
				p1 = H-p;
				for(int j=0; j<H; j++) {
					cnts_prior[i][j][0] = emiss[hsc[i][j]]/p;
					cnts_prior[i][j][1] = (1-emiss[hsc[i][j]])/p1;
				}
				emissc[i] = p/H;
			}
		}
		
		protected String[] getAllele() {
			return this.allele;
		}
		
		protected double[] getEmiss() {
			return this.emiss;
		}
		
		protected double[] getEmissc() {
			return this.emissc;
		}
		
		protected double getBfrac() {
			return this.bfrac;
		}
		
		protected double[][] getCount() {
			return this.count;
		}
		
		protected void pseudo() {
			for(int i=0; i<count.length; i++) {
				count[i][0] = (1-bfrac)*mu_A_m;
				count[i][1] = bfrac*mu_A_m;
			}
		}
	}
	
	protected class PathUnit {
		protected String[] path_str;
		protected int[] path;
		
		public PathUnit() {
			// TODO Auto-generated constructor stub
			this.path_str = new String[M];
			this.path = new int[M];
		}
	}
	
	protected class ModelWriter {
		protected ZipOutputStream out = null;
		
		public ModelWriter(String file) {
			try {
				out = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream(file), 65536));
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		public void writeHaplotype() {
			try {
				out.putNextEntry(new ZipEntry("haplotype.txt"));
				out.write((""+loglik()+"\n").getBytes());
				out.write((""+M+"\n").getBytes());
				StringBuilder miss_str = new StringBuilder();
				for(int s=0; s<M; s++) miss_str.append("*");
				miss_str.append("\n");
				byte[] miss_byt = miss_str.toString().getBytes();
				for(int i=0; i<N; i++) {
					if(pas[i]==null) {
						for(int k=0; k<H; k++) {
							out.write(("# id "+samples[i]+":"+(k+1)+"\t\t\t").getBytes());
							out.write(miss_byt);
						}
					} else {
						String[] path = pas[i].path_str;
						List<String[]> path_s = new ArrayList<String[]>();
						for(int k=0; k<path.length; k++)
							path_s.add(path[k].split("_"));
						for(int k=0; k<H; k++) {
							out.write(("# id "+samples[i]+":"+(k+1)+"\t\t\t").getBytes());
							for(int s=0; s<path_s.size(); s++)
								out.write(path_s.get(s)[k].getBytes());
							out.write("\n".getBytes());
						}
					}
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		public void writeDosage() {
			try {
				out.putNextEntry(new ZipEntry("dosage.txt"));
				StringBuilder os = new StringBuilder();
				os.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
				for(String s : samples) {
					os.append("\t");
					os.append(s);
				}
				os.append("\n");
				out.write(os.toString().getBytes());
				
				List<String[]> allele = de.getAllele();
				double[] position = de.getPosition();
				String[] markers_id = de.getMarker();
				String id = de.getId();
				int[] states;
				double[] emiss;
				for(int i=0; i<M; i++) {
					os.setLength(0);
					os.append(id);
					os.append("\t");
					os.append((int) position[i]);
					os.append("\t");
					os.append(markers_id[i]);
					os.append("\t");
					os.append(allele.get(i)[0]);
					os.append("\t");
					os.append(allele.get(i)[1]);
					os.append("\t.\t.\t.\tGT");
					
					emiss = emission[i].getEmiss();
					for(int j=0; j<N; j++) {
						int dosa;
						if(fi1ter[j]) {
							dosa = -1;
						} else {
							states = state.getHsc()[pas[j].path[i]];
							dosa = 0;
							for(int k=0; k<H; k++) 
								if(emiss[states[k]]>=0.5)
									++dosa;
						}
						os.append("\t");
						os.append(dosa);
					}
					os.append("\n");
					out.write(os.toString().getBytes());
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		public void writeGenotype() {
			try {
				out.putNextEntry(new ZipEntry("genotype.txt"));
				StringBuilder os = new StringBuilder();
				os.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
				for(String s : samples) {
					os.append("\t");
					os.append(s);
				}
				os.append("\n");
				out.write(os.toString().getBytes());
				
				List<String[]> allele = de.getAllele();
				double[] position = de.getPosition();
				String[] markers_id = de.getMarker();
				String id = de.getId();
				int[] states;
				double[] emiss;
				for(int i=0; i<M; i++) {
					os.setLength(0);
					os.append(id);
					os.append("\t");
					os.append((int) position[i]);
					os.append("\t");
					os.append(markers_id[i]);
					os.append("\t");
					os.append(allele.get(i)[0]);
					os.append("\t");
					os.append(allele.get(i)[1]);
					os.append("\t.\t.\t.\tGT");
					
					emiss = emission[i].getEmiss();
					for(int j=0; j<N; j++) {
						os.append("\t");
						if(fi1ter[j]) {
							os.append(".");
							for(int k=1; k<H; k++) {
								os.append("|");
								os.append(".");
							}
						} else {
							states = state.getHsc()[pas[j].path[i]];
							os.append(emiss[states[0]]<0.5?1:0);
							for(int k=1; k<H; k++) {
								os.append("|");
								os.append(emiss[states[k]]<0.5?1:0);
							}
						}
					}
					os.append("\n");
					out.write(os.toString().getBytes());
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		public void writeEmissionModel() {
			try {
				out.putNextEntry(new ZipEntry("emission.txt"));
				for(int i=0; i<M; i++) {
					double[] emiss = emission[i].getEmiss();
					String[] allele = de.getAllele().get(i);
					out.write((de.getId()+"_"+de.getPosition()[i]+"\t\t\t").getBytes());
					for(int j=0; j<emiss.length; j++) {
						out.write((state.getHs()[j]+"-> {").getBytes());
						out.write((allele[0]+","+emiss[j]+";").getBytes());
						out.write((allele[1]+","+(1-emiss[j])+";").getBytes());
						out.write("} ".getBytes());
					}
					out.write("\n".getBytes());
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		public void writeSNP() {
			try {
				out.putNextEntry(new ZipEntry("snp.txt"));

				StringBuilder os = new StringBuilder();
				List<String[]> allele = de.getAllele();
				double[] position = de.getPosition();
				String[] markers_id = de.getMarker();
				String id = de.getId();
				for(int i=0; i<position.length; i++) {
					os.setLength(0);
					os.append(id);
					os.append("\t");
					os.append((int) position[i]);
					os.append("\t");
					os.append(true_id[i]);
					os.append("\t");
					os.append(true_pos[i]);
					os.append("\t");
					os.append(markers_id[i]);
					os.append("\t");
					os.append(allele.get(i)[0]);
					os.append("\t");
					os.append(allele.get(i)[1]);
					os.append("\n");
					out.write(os.toString().getBytes());
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		public void writeRunInfo() {
			try {
				out.putNextEntry(new ZipEntry("runinfo.txt"));
				out.write(("##marker: "+M+"\n").getBytes());
				out.write(("##sample: "+N+"\n").getBytes());
				out.write(("##ploidy: "+H+"\n").getBytes());
				out.write((("##parents: "+(parents[0]==null?"":parents[0])+
						" "+(parents[1]==null?"":parents[1])).trim()+"\n").getBytes());
				out.write(("##progeny:").getBytes());
				for(int i : progeny_i) out.write((" "+samples[i]).getBytes());
				out.write(("\n").getBytes());
				out.write(("##chrs: "+Utils.paste(chrs, ",")+"\n").getBytes());
				out.write(("##chrs_rev: "+Utils.paste(chrs_rev, ",")+"\n").getBytes());
				out.write(("##model_len: "+Utils.paste(Ms, ",")+"\n").getBytes());
				double[] model_ll = new double[Ms.length];
				int k = 0;
				for(int i=0; i<Ms.length; i++) {
					model_ll[i] = loglik(k, k+Ms[i]);
					k += Ms[i];
				}
				out.write(("##model_ll: "+Utils.paste(model_ll, ",")+"\n").getBytes());
				out.write(("##seed: "+Constants.seed+"\n").getBytes());
				out.write(("##iteration: "+iteration+"\n").getBytes());
				out.write(("##loglik: "+loglik()+"\n").getBytes());
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		public void close() {
			// TODO Auto-generated method stub
			try {
				out.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}
}
