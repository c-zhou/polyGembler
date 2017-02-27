package cz1.test;

import java.util.Arrays;

import cz1.util.Executor;

public class TaskScheduler extends Executor {

	public static void main(String[] args) {
		TaskScheduler taskScheduler = new TaskScheduler();
		taskScheduler.THREADS = 32;
		taskScheduler.run();
	}
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		
	}
	
	final private int task_round = 5;
	final private int task_type = 3;
	final private int task_repeat = 4;
	
	// the first row of task_table is dummy tasks, will all set to zeros 
	final private int[][] task_table = new int[task_round+1][task_type];
	final private int[] task_progress = new int[task_type];
	
	private class Task implements Runnable {
		private final int task_round;
		private final int task_type;
		private final int task_repeat;
		
		public Task(int task_type,
				int task_round,
				int task_repeat) {
			this.task_type = task_type;
			this.task_round = task_round;
			this.task_repeat = task_repeat;
		}
		
		@Override
		public void run() {
			// TODO Auto-generated method stub
			try {
				Thread.sleep(10);
				
				myLogger.info(task_type+" "
						+ task_round+" "
						+ task_repeat);
				
				synchronized (task_table) {
					task_table[task_round][task_type]--;
					task_table.notify();
				}
			} catch (Exception e) {
				// TODO Auto-generated catch block
				Thread t = Thread.currentThread();
				t.getUncaughtExceptionHandler().uncaughtException(t, e);
				e.printStackTrace();
				executor.shutdown();
				System.exit(1);
			}
		}
	}
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		for (int i=1; i<task_table.length; i++)
		    Arrays.fill(task_table[i], task_repeat);
		this.initial_thread_pool();
		
		while(true) {
			this.scheduleTask();
			if(this.task_done()) break;
		}
		this.waitFor();
	}

	private void scheduleTask() {
		// TODO Auto-generated method stub
		for(int i=0; i<task_type; i++)
			if(task_progress[i]<=task_round &&
					task_table[task_progress[i]][i]==0) {
				task_progress[i]++;
				if(task_progress[i]<=task_round)
					for(int j=0; j<task_repeat; j++)
						executor.submit(new Task(i, 
								task_progress[i], j));
			}
	}

	private boolean task_done() {
		// TODO Auto-generated method stub
		for(int i=0; i<task_type; i++)
			if(task_progress[i]<=task_round)
				return false;
		return true;
	}

}
