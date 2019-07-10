package cz1.hmm.model;

public interface ForwardBackwardTrainer {
	
	void train();
	void makeNaiveTrainer();
	double findPath();
	void forward();
	void backward();
	void em();
	void check();
	double loglik();
	void write(String output, 
			String experiment, 
			String contig);

}
