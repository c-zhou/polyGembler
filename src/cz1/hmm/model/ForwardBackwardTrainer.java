package cz1.hmm.model;

public interface ForwardBackwardTrainer {
	
	void train();
	void makeNaiveTrainer();
	void forward();
	void backward();
	void em();
	void check();
}
