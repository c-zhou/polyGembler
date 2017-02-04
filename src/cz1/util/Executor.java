package cz1.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
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

	protected final static int mb = 1024*1024;
	protected final static Runtime instance = Runtime.getRuntime();
	
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
	
	protected void require(String tool) {
		String command = "command -v "+tool+
				" >/dev/null 2>&1 && { echo \"true\"; } || { echo \"false\"; }";
		try {
			Process check = bash(command);
			BufferedReader in =
					new BufferedReader(new InputStreamReader(check.getInputStream()));
			String line = in.readLine();
			in.close();
			check.waitFor();
			if(line.equals("false")) 
				throw new RuntimeException(
						"\n\nTool "+tool+" is required but not available on this machine.\n"
						+ "In you have already installed it, you will need to add it\n"
						+ "to system path or use \"module load\" to load the tool.\n\n ");
		} catch (IOException | InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
	}
	
	protected Process bash(String command) {
		String[] runnable = new String[]{"bash","-c",command};
		try {
			return instance.exec(runnable);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			throw new RuntimeException("Runtime exception!!!");
		}
	}
	
	protected void bash(final String[] commands) {
		initial_thread_pool();
		
		for(int i=0; i<commands.length; i++) {
			executor.submit(new Runnable() {
				private int i;
				
				@Override
				public void run() {
					// TODO Auto-generated method stub
					String[] runnable = new String[]{"bash","-c",commands[i]};
					try {
						try {
							instance.exec(runnable).waitFor();
						} catch (InterruptedException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					} catch (IOException e) {
						// TODO Auto-generated catch block
						throw new RuntimeException("Runtime exception!!!");
					}
				}
				public Runnable init(final int i) {
			        this.i = i;
			        return(this);
			    }
			}.init(i));
		}
		
		try {
			executor.shutdown();
			executor.awaitTermination(365, TimeUnit.DAYS);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	protected static double maxMemory() {
		return instance.maxMemory() / mb;
	}
	
	protected static double totalMemory() {
		return instance.totalMemory() / mb;
	}
	
	protected static double freeMemory() {
		return instance.freeMemory() / mb;
	}
	
	protected static double usedMemory() {
		return totalMemory()-freeMemory();
	}
	
	protected static void usage() {
		myLogger.info("Max Memory: "+maxMemory());
		myLogger.info("Total Memory: "+totalMemory());
		myLogger.info("Free Memory: "+freeMemory());
		myLogger.info("Used Memory: "+usedMemory());
	}
}
