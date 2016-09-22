package cz1.model;

public abstract class HiddenMarkovModel {

	public abstract void train();
	public abstract double loglik();
	public abstract void write(String output, 
			String experiment, 
			String contig);
	public abstract void print(boolean details);
	public abstract void print();
}
