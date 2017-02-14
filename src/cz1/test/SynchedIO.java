package cz1.test;

import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import cz1.util.Executor;
import cz1.util.Utils;

public class SynchedIO extends Executor {

	public static void main(String[] args) {
		SynchedIO synchedIO = new SynchedIO();
		synchedIO.THREADS = 4;
		synchedIO.run();
	}
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		
	}

	private final List<Integer> synchedBlock = 
			Collections.synchronizedList(new LinkedList<Integer>());
	
	public void write(final int i, final int[] Q, final int t) 
			throws InterruptedException {
		synchronized (synchedBlock) {
			while (synchedBlock.get(0)!=i)
				synchedBlock.wait();
			//IO.print("Thread "+t+" - ");
			Utils.print(Q);
			synchedBlock.remove(0);
			synchedBlock.notifyAll();
		}
	}
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		
		this.initial_thread_pool();
		for(int i=0; i<10000; i++) {
			
			int[] Q = new int[1];
			Arrays.fill(Q, i);
			synchronized (synchedBlock) {synchedBlock.add(i);}
			
			executor.submit(new Runnable() {
        		private int[] Q;
        		private int i;
        		
				@Override
				public void run() {
					// TODO Auto-generated method stub
					try {
						write(this.i, this.Q, 
								(int)Thread.currentThread().getId()%THREADS);
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				}
				
				public Runnable init(int[] Q, int i) {
			        this.Q = Q;
			        this.i = i;
			        return(this);
			    }
        	}.init(Q, i));
		}
		
		this.waitFor();
	}

}
