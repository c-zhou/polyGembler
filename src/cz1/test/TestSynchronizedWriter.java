package cz1.test;

import java.lang.Runnable;
import java.util.ArrayList;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.TimeUnit;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import cz1.model.SynchronizedWriter;

public class TestSynchronizedWriter {

	private static BufferedWriter bw;
	private static SynchronizedWriter synw;
	private final static int NUM_CORES = 
			Runtime.getRuntime().availableProcessors();

	public static void main(String[] args) throws IOException, InterruptedException {
		String[] a = new String[1];
		System.out.println(a[0]==null);
		
		synw = new SynchronizedWriter(new BufferedWriter(
				new FileWriter(
						"C:\\Users\\chenxi.zhou\\Desktop\\test.txt1")),
						new LinkedBlockingQueue<StringBuilder>());
		//bw = new BufferedWriter(new FileWriter(
		//	"C:\\Users\\chenxi.zhou\\Desktop\\test.txt1"));
		ExecutorService executor = Executors.newFixedThreadPool(NUM_CORES);
		System.out.println(NUM_CORES);
		executor.submit(synw);
		for (int i = 0; i < 50; i++) 
			executor.submit(r);
		synw.close();
		executor.shutdown();
		executor.awaitTermination(1, TimeUnit.SECONDS);
		//bw.close();
	}

	private static Runnable r = new Runnable() {

		@Override
		public void run() {
			// TODO Auto-generated method stub
			try {
				String s = "";
				for(int i=0; i<5000; i++) 
					s += "A";
				synw.write(Thread.currentThread().getId()+": "+s+"\n");
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	};
}
