package cz1.test;

import cz1.util.Executor;

public class ThreadPool extends Executor {

	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		
	}
	
	public static void main(String[] args) {
		ThreadPool pool = new ThreadPool();
		myLogger.info(pool.executor==null ||
				pool.executor.isShutdown());
	}
}
