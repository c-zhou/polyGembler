package cz1.util;

import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.RejectedExecutionHandler;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;

import cz1.util.ArgsEngine;

public abstract class Executor {


	protected final static Logger myLogger = 
			Logger.getLogger(Executor.class);
	static {
		BasicConfigurator.configure();
	}
	protected ArgsEngine myArgsEngine = null;

	public abstract void printUsage();

	public abstract void setParameters(String[] args);
	
	public abstract void run();
	
	protected int THREADS = 1;
	protected static ExecutorService executor;
	protected BlockingQueue<Runnable> tasks = null;

	protected void initial_thread_pool() {
		tasks = new ArrayBlockingQueue<Runnable>(THREADS);
		executor = new ThreadPoolExecutor(THREADS, 
				THREADS, 
				1, 
				TimeUnit.SECONDS, 
				tasks, 
				new RejectedExecutionHandler(){
			@Override
			public void rejectedExecution(Runnable task,
					ThreadPoolExecutor arg1) {
				// TODO Auto-generated method stub
				try {
					tasks.put(task);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		});
	}
}
