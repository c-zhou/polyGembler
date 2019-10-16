package cz1.hmm.model;

import java.util.Arrays;

import org.apache.log4j.Logger;

import cz1.hmm.data.DataEntry;
import cz1.math.Algebra;
import cz1.util.Constants.Field;

public class ModelTrainer extends EmissionModel implements ForwardBackwardTrainer {
	
	private final static Logger myLogger = Logger.getLogger(ModelTrainer.class);
	
	private FBUnit[] forward, backward;
	
	public ModelTrainer(DataEntry[] de, 
			double[] seperation, 
			boolean[] reverse, 
			Field field,
			int ploidy,
			String[] parents) {
		// TODO Auto-generated constructor stub
		super(de, seperation, reverse, field, ploidy, parents, true);
		this.makeNaiveTrainer();
	}

	@Override
	public void train() {
		// TODO Auto-generated method stub
		++iteration;
		
		refresh();
		forward();
		backward();
		check();	
		em();
	}

	@Override
	public void makeNaiveTrainer() {
		// TODO Auto-generated method stub
		this.forward = new FBUnit[N];
		for(int i=0; i<N; i++) 
			this.forward[i] = new FBUnit(false);
		this.backward = new FBUnit[N];
		for(int i=0; i<N; i++) 
			this.backward[i] = new FBUnit(true);
		return;
	}


	@Override
	public double findPath() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	protected double loglik(int fromIndex, int toIndex) {
		// TODO Auto-generated method stub
		if(fromIndex<0||fromIndex>=toIndex||toIndex>M)
			throw new RuntimeException("!!!");
		
		int m = toIndex - fromIndex;
		double probability = 0;
		for(int i=0; i<N; i++) {

			Integer[] ss = sspace.get(i);
			double pi = Math.log(1.0/ss.length);

			double[][] probsMat = new double[m][K];
			for(int j=0; j<m; j++)
				Arrays.fill(probsMat[j], Double.NEGATIVE_INFINITY);
				
			ObUnit[] ob = obs[i];

			double[] emiss = ob[fromIndex].emiss;
			
			for(int k : ss) probsMat[0][k] = pi+emiss[k];

			for(int j=1; j<m; j++) {

				emiss = ob[fromIndex+j].emiss;
				for(int k : ss)
					probsMat[j][k] = emiss[k]+probsMat[j-1][k];	
			}
			probability += Algebra.sumExps(probsMat[m-1]);
		}
		return probability;	
	}
	
	@Override
	public void forward() {
		// TODO Auto-generated method stub
		for(int i=0; i<N; i++) {

			Integer[] ss = sspace.get(i);
			double pi = Math.log(1.0/ss.length);

			double[][] probsMat = forward[i].probsMat;
			ObUnit[] ob = obs[i];

			double[] emiss = ob[0].emiss;
			
			for(int k : ss) probsMat[0][k] = pi+emiss[k];

			for(int j=1; j<M; j++) {

				emiss = ob[j].emiss;
				for(int k : ss)
					probsMat[j][k] = emiss[k]+probsMat[j-1][k];	
			}
			forward[i].probability(Algebra.sumExps(probsMat[M-1]));
		}
		return;
	}

	@Override
	public void backward() {
		// TODO Auto-generated method stub
		for(int i=0; i<N; i++) {
			Integer[] ss = sspace.get(i);
			double[][] probsMat = backward[i].probsMat;
			ObUnit[] ob = obs[i];
			
			double[] emiss;
			for(int k : ss) probsMat[M-1][k] = 0;
			
			for(int j=M-2; j>=0; j--) {	
				emiss = ob[j+1].emiss;
				for(int k : ss) 
					probsMat[j][k] = emiss[k]+probsMat[j+1][k];
			}
			
			double pi = Math.log(1.0/ss.length);
			emiss = ob[0].emiss;
			double[] logs = new double[ss.length];
			int j=0;
			for(int k : ss) logs[j++] = pi+emiss[k]+probsMat[0][k];
			backward[i].probability(Algebra.sumExps(logs));
		}
		return;
	}

	@Override
	public void em() {
		// TODO Auto-generated method stub
		FBUnit fw1, bw1;
		ObUnit ob1;
		double count, coeff;
		Integer[] ss;
		
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
				coeff = weights[j==parents_i[0]||j==parents_i[1]?0:1];
				
				for(int a : ss) {
					count = coeff*
							Math.exp(fw1.probsMat[i][a]+
							bw1.probsMat[i][a]-
							fw1.probability);
					e1.addCount(a, count*acnt, count*bcnt);
				}
			}
			e1.update();
		}
	}


	@Override
	public double loglik() {
		// TODO Auto-generated method stub
		if(iteration==0)
			return Double.NEGATIVE_INFINITY;
		else {
			double probability = 0;
			for(int i=0; i<N; i++)
				probability += weights[i==parents_i[0]||i==parents_i[1]?0:1]*this.forward[i].probability;
			return probability;
		}
	}

	@Override
	public void check() {
		// TODO Auto-generated method stub
		if(iteration==0) return;
		for(int i=0; i<this.forward.length; i++) {
			double r = Math.abs(forward[i].probability-backward[i].probability);
			if(r>1e-6) 
				throw new RuntimeException("Different likelihood by forward and backward algorithm: forward, "+
						forward[i].probability+"; backward, "+backward[i].probability);
		}
	}

	@Override
	public void write(String output, String experiment, String contig) {
		// TODO Auto-generated method stub
		
	}
	
	protected class FBUnit { /** forward/backward unit */
		protected double[][] probsMat;
		protected double probability;
		protected final boolean backward;

		public FBUnit(boolean backward) {
			this.backward = backward;
			this.probability = 0;
			this.probsMat = new double[M][K];
			for(int i=0; i<M; i++)
				Arrays.fill(probsMat[i], Double.NEGATIVE_INFINITY);
		}

		public void probability(double p) {
			// TODO Auto-generated method stub
			this.probability = p;
		}

		public double probability() {
			// TODO Auto-generated method stub
			return this.probability;
		}	
	}
}
