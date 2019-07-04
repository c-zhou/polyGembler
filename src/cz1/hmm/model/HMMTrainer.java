package cz1.hmm.model;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.log4j.Logger;

import cz1.hmm.data.DataEntry;
import cz1.util.Constants;
import cz1.util.Constants.Field;

public class HMMTrainer extends HiddenMarkovModel {
	protected final static Logger myLogger = Logger.getLogger(HMMTrainer.class);
	private FBUnit[] forward, backward;
	private static int iteration = 0;
	
	public HMMTrainer(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse, 
			Field field) {
		// TODO Auto-generated constructor stub
		super(de, seperation, reverse, field);
		this.makeNaiveTrainer();
	}

	public HMMTrainer(DataEntry[] de,
			double[] seperation,
			boolean[] reverse,
			Field field, 
			String hmm_file) {
		// TODO Auto-generated constructor stub
		super(de, seperation, reverse, field);
		this.makeTrainerFromFile(hmm_file);
	}

	private void makeNaiveTrainer() {
		// TODO Auto-generated method stub
		this.forward = new FBUnit[N];
		for(int i=0; i<N; i++) 
			this.forward[i] = new FBUnit(false);
		this.backward = new FBUnit[N];
		for(int i=0; i<N; i++) 
			this.backward[i] = new FBUnit(true);
		return;
	}

	private void makeTrainerFromFile(String hmm_file) {
		// TODO Auto-generated method stub
		throw new RuntimeException("!!!");
	}
	
	@Override
	public void train() {
		// TODO Auto-generated method stub
		++iteration;
		refresh();
		//findViterbiPath();
		forward();
		backward();
		check();
		baumWelch();
	}

	private void baumWelch() {
		// TODO Auto-generated method stub
		
		FBUnit fw1, bw1;
		ObUnit ob1;
		double exp_c, exp, count;
		Integer[] ss;
		
		TransitionUnit t1;
		for(int i=0; i<M-1; i++) {
			// transitions
			t1 = transition[i];
			t1.pseudo();
			for(int j=0;j<N; j++) {
				ss = sspace.get(j);
				fw1 = forward[j];
				bw1 = backward[j];
				ob1 = obs[j][i+1];
				exp_c = fw1.logscale[i]+
						bw1.logscale[i+1]+
						ob1.getLogScale(j)-
						fw1.probability;
				
				if(exp_c>Constants.MAX_EXP_DOUBLE) { 
					for(int a : ss) {
						for(int b : ss) { 
							count = Math.exp(Math.log(
									fw1.probsMat[i][a]*
									t1.trans(a, b)*
									ob1.emiss[b]*
									bw1.probsMat[i+1][b])+
									exp_c);
							t1.addCount(a, b, count);
						}
					}
				} else {
					exp = Math.exp(exp_c);
					for(int a : ss) {
						for(int b : ss) { 
							count = fw1.probsMat[i][a]*
									t1.trans(a, b)*
									ob1.emiss[b]*
									bw1.probsMat[i+1][b]*
									exp;
							t1.addCount(a, b, count);
						}
					}
				}
			}
			//t1.update();
		}

		int acnt, bcnt;
		EmissionUnit e1;
		for(int i=0; i<M; i++) {
			e1 = emission[i];
			e1.pseudo();
			
			for(int j=0;j<N; j++) {
				ss = sspace.get(j);
				fw1 = forward[j];
				bw1 = backward[j];
				ob1 = obs[j][i];
				acnt = ob1.getAa();
				bcnt = ob1.getCov()-acnt;
				exp_c = fw1.logscale[i]+
						bw1.logscale[i]-
						fw1.probability;
				
				if(exp_c>Constants.MAX_EXP_DOUBLE) {
					for(int a : ss) {
						count = Math.exp(Math.log(
								fw1.probsMat[i][a]*
								bw1.probsMat[i][a])+
								exp_c);
						e1.addCount(a, count*acnt, count*bcnt);
					}
				} else {
					exp = Math.exp(exp_c);
					for(int a : ss) {
						count = fw1.probsMat[i][a]*
								bw1.probsMat[i][a]*
								exp;
						e1.addCount(a, count*acnt, count*bcnt);
					}
				}
			}
			e1.update();
		}
	}

	private void backward() {
		// TODO Auto-generated method stub
		for(int i=0; i<N; i++) {
			Integer[] ss = sspace.get(i);
			double[][] probsMat = backward[i].probsMat;
			double[] logscale = backward[i].logscale;
			ObUnit[] ob = obs[i];
			
			double[] emiss;
			TransitionUnit t;
			
			for(int k : ss) probsMat[M-1][k] = 1.0;
			logscale[M-1] = 0;
			double tmp; 
			
			for(int j=M-2; j>=0; j--) {	
				emiss = ob[j+1].emiss;
				t = transition[j];
				for(int k : ss) {
					tmp = 0;
					for(int z : ss) 
						tmp += t.trans(k, z)*emiss[z]*probsMat[j+1][z];
					probsMat[j][k] = tmp;
				}
				
				logscale[j] = ob[j+1].getLogScale(i);
				backward[i].scale(j);
			}
			
			double pi = 1.0/ss.length;
			double p = 0.0;
			emiss = ob[0].emiss;
			for(int z : ss)
				p += pi*emiss[z]*probsMat[0][z];
			backward[i].probability(p, ob[0].getLogScale(i));
		}
		return;
	}

	private void forward() {
		// TODO Auto-generated method stub
		for(int i=0; i<N; i++) {
			boolean isparent = i==parents_i[0]||i==parents_i[1];
			
			Integer[] ss = sspace.get(i);
			double pi = 1.0/ss.length;
			
			double[][] probsMat = forward[i].probsMat;
			double[] logscale = forward[i].logscale;
			ObUnit[] ob = obs[i];
			
			double[] emiss = ob[0].emiss;
			TransitionUnit t;
			
			for(int k : ss) probsMat[0][k] = pi*emiss[k];
			logscale[0] = ob[0].getLogScale(i);
			double tmp; 
			
			for(int j=1; j<M; j++) {
				
				emiss = ob[j].emiss;
				t = transition[j-1];
				
				for(int k : ss) {
					tmp = 0;
					for(int z : ss)
						tmp += probsMat[j-1][z]
								*t.trans(z, k);
					probsMat[j][k] = emiss[k]*tmp;
					
					if(isparent && probsMat[j][k]<Constants.threshMin) {
						
					}
				}
				
				logscale[j] = ob[j].getLogScale(i);
				forward[i].scale(j);
			}
			forward[i].probability(StatUtils.sum(probsMat[M-1]));
		}
		return;
	}

	private void check() {
		// TODO Auto-generated method stub
		if(iteration==0) return;
		for(int i=0; i<this.forward.length; i++) {
			double r = Math.abs(forward[i].probability-backward[i].probability);
			if(r>1e-6) 
				throw new RuntimeException("Different likelihood by forward and backward algorithm!");
		}
	}

	public double loglik1() {
		// TODO Auto-generated method stub
		if(iteration==0)
			return Double.NEGATIVE_INFINITY;
		else {
			double probability = 0;
			for(FBUnit bw : this.backward) probability += bw.probability;
			return probability;
		}
	}

	@Override
	public double loglik() {
		// TODO Auto-generated method stub
		if(iteration==0)
			return Double.NEGATIVE_INFINITY;
		else {
			double probability = 0;
			for(FBUnit fw : this.forward) probability += fw.probability;
			return probability;
		}
	}

	@Override
	public void write(String output, String experiment, String contig) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void write(String output, String experiment, String scaff, double loglik_diff) {
		// TODO Auto-generated method stub
		// TODO Auto-generated method stub
		
		this.findViterbiPath();

		String root = experiment+"."+
				scaff+"."+
				System.nanoTime()+"_1_1_all_0_200mb";

		try {
			ZipOutputStream out = new ZipOutputStream(new BufferedOutputStream(new 
					FileOutputStream(output+"/"+root+".zip"), 65536));

			out.putNextEntry(new ZipEntry("phasedStates/"+experiment+".txt"));
			out.write((""+this.loglik()+"\n").getBytes());
			out.write((""+M+"\n").getBytes());
			for(int i=0; i<N; i++) {
				String[] path = this.vbs[i].path_str;
				List<String[]> path_s = new ArrayList<String[]>();
				for(int k=0; k<path.length; k++)
					path_s.add(path[k].split("_"));
				for(int k=0; k<Constants._ploidy_H; k++) {
					out.write(("# id "+this.samples[i]+":"+(k+1)+"\t\t\t").getBytes());
					for(int s=0; s<path_s.size(); s++)
						out.write(path_s.get(s)[k].getBytes());
					out.write("\n".getBytes());
				}
			}
			
			out.putNextEntry(new ZipEntry("phased_genotypes"+experiment+".txt"));
			for(int i=0; i<N; i++) {
				int[] path = this.vbs[i].path;
				for(int j=0; j<M; j++) {
					int[] states = this.state.getHsc()[path[j]];
					double[] emiss = this.emission[j].getEmiss();
					int dosa = 0;
					for(int k=0; k<H; k++)
						if(emiss[states[k]]<0.5)
							++dosa;
					out.write((""+dosa+"\t").getBytes());
				}
				out.write("\n".getBytes());
			}

			out.putNextEntry(new ZipEntry("results_hmm/emissionModel.txt"));
			for(int i=0; i<M; i++) {
				double[] emiss = this.emission[i].getEmiss();
				String[] allele = this.de.getAllele().get(i);
				out.write((this.de.getId()+"_"+this.de.getPosition()[i]+"\t\t\t").getBytes());
				for(int j=0; j<emiss.length; j++) {
					out.write((this.state.getHs()[j]+"-> {").getBytes());
					out.write((allele[0]+","+emiss[j]+";").getBytes());
					out.write((allele[1]+","+(1-emiss[j])+";").getBytes());
					out.write("} ".getBytes());
				}
				out.write("\n".getBytes());
			}

			out.putNextEntry(new ZipEntry("results_hmm/transitionModel.txt"));
			for(int i=0; i<M-1; i++)
				out.write((this.de.getId()+"_"+this.de.getPosition()[i]+"\t\t\t"+this.transition[i].getJump()+"\n").getBytes());
			
			out.putNextEntry(new ZipEntry("stderr_true"));

			out.write("cz1.model.HiddenMarkovModel:\n".getBytes());
			out.write("cz1.model.HiidenMarkovModel$EM:\n".getBytes());
			out.write(("log prob is "+this.loglik()+" at "+iteration+"\n").getBytes());
			out.write(("random seed is "+Constants.seed).getBytes());

			out.putNextEntry(new ZipEntry("snp_"+experiment+".txt"));

			StringBuilder os = new StringBuilder();
			List<String[]> allele = de.getAllele();
			double[] position = de.getPosition();
			String[] markers_id = de.getMarker();
			String id = de.getId();
			for(int i=0; i<position.length; i++) {
				os.setLength(0);
				os.append(id);
				os.append("\t");
				int p = (int) position[i];
				os.append(p);
				os.append("\t");
				os.append(p);
				os.append("\t");
				os.append(markers_id[i]);
				os.append("\t");
				for(int j=0; j<allele.get(i).length; j++) {
					os.append(allele.get(i)[j]);
					os.append("\t");
				}
				os.append(experiment);
				os.append("\t");
				os.append("false\n");
				out.write(os.toString().getBytes());
			}			
			out.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
	}

	@Override
	public void print(boolean details) {
		// TODO Auto-generated method stub
		
	}
}
