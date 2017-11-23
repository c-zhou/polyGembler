package cz1.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import java.security.CodeSource;
import java.security.ProtectionDomain;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.RejectedExecutionHandler;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;

import cz1.util.ArgsEngine;

public abstract class Executor {

	protected final static int mb = 1024*1024;
	protected final static Runtime instance = Runtime.getRuntime();
	
	protected final static Logger myLogger = 
			Logger.getLogger(Executor.class);
	
	protected ArgsEngine myArgsEngine = null;

	public abstract void printUsage();

	public abstract void setParameters(String[] args);
	
	public abstract void run();
	
	protected int THREADS = 1;
	protected ExecutorService executor;
	protected BlockingQueue<Runnable> tasks = null;

	protected void initial_thread_pool(int thread_num) {
		if(this.tasks!=null&&tasks.size()>0)
			throw new RuntimeException("Thread pool initialisation "
					+ "exception!!! Task queue is not emputy!!!");
		if(this.executor!=null&&!executor.isShutdown())
			throw new RuntimeException("Thread pool initialisation "
					+ "exception!!! Executor is on the fly!!!");
		
		tasks = new ArrayBlockingQueue<Runnable>(thread_num);
		executor = new ThreadPoolExecutor(thread_num, 
				thread_num, 
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
	
	protected void initial_thread_pool() {
		this.initial_thread_pool(this.THREADS);
	}
	
	protected void waitFor() {
		// TODO Auto-generated method stub
		try {
			executor.shutdown();
			executor.awaitTermination(365, TimeUnit.DAYS);
		} catch (InterruptedException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
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
			this.consumeErr(check);
			if(line.equals("false")) 
				throw new RuntimeException(
						"\n\nTool "+tool+" is required but not available on this machine.\n"
						+ "In you have already installed it, you will need to add it\n"
						+ "to system path or use \"module load\" to load the tool.\n\n ");
		} catch (IOException e1) {
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
		this.initial_thread_pool();
		
		for(int i=0; i<commands.length; i++) {
			executor.submit(new Runnable() {
				private int i;
				
				@Override
				public void run() {
					// TODO Auto-generated method stub
					try {
							Process process = bash(commands[i]);
							consume(process);
					} catch (Exception e) {
						// TODO Auto-generated catch block
						Thread t = Thread.currentThread();
						t.getUncaughtExceptionHandler().uncaughtException(t, e);
						e.printStackTrace();
						executor.shutdown();
						System.exit(1);
					}
				}
				
				public Runnable init(final int i) {
			        this.i = i;
			        return(this);
			    }
			}.init(i));
		}
		
		this.waitFor();
	}

	private class BashStream implements Runnable {
		private final InputStream in;
		Thread thread;
		
		public BashStream(InputStream in) {
			this.in = in;
		}
		
		public void start () {
	        thread = new Thread (this);
	        thread.start();
	    }
		
		@Override
		public void run() {
			// TODO Auto-generated method stub
			try {
				BufferedReader reader = new BufferedReader(new InputStreamReader(in));
				String line;
				while( (line=reader.readLine())!=null )
					myLogger.info(line);
				in.close();
				reader.close();
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
	
	/**
	protected void consumeOnTheFly(Process process) {
		// TODO Auto-generated method stub
		try {
			BashStream stdout = new BashStream(process.getInputStream());
			BashStream stderr = new BashStream(process.getErrorStream());
			stdout.start();
			stderr.start();
			process.waitFor();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	**/
	
	protected void consume(Process process) {
		// TODO Auto-generated method stub
		try {
			BashStream stdout = new BashStream(process.getInputStream());
			BashStream stderr = new BashStream(process.getErrorStream());
			stdout.start();
			stderr.start();
			process.waitFor();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	protected void consumeErr(Process process) {
		// TODO Auto-generated method stub
		try {
			BashStream stderr = new BashStream(process.getErrorStream());
			stderr.start();
			process.waitFor();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	protected void consumeOut(Process process) {
		// TODO Auto-generated method stub
		try {
			BashStream stdout = new BashStream(process.getInputStream());
			stdout.start();
			process.waitFor();
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
	
	protected static URI uri() {
		try {
			return Executor.class.
					getProtectionDomain().
					getCodeSource().
					getLocation().
					toURI();
		} catch (URISyntaxException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
	
	protected static String makeExecutable(
			String executable,
			String out_prefix) {
        try {
        	final ZipFile zipFile = new ZipFile(new File(uri()));

            final File tempFile = new File(out_prefix+"/"+new File(executable).getName());
            //tempFile.deleteOnExit();
            final ZipEntry entry = zipFile.getEntry(executable);
            
            if(entry == null) {
            	zipFile.close();
            	throw new RuntimeException("Executable file "+executable+
            			" not found in "+zipFile.getName()+"!!!");
            }
            
            final InputStream zipStream  = zipFile.getInputStream(entry);
            final byte[] buf = new byte[65536];
            int i = 0;
            final OutputStream fileStream = new FileOutputStream(tempFile);
            while((i = zipStream.read(buf)) != -1) fileStream.write(buf, 0, i);
            
            fileStream.close();
            zipStream.close();
            zipFile.close();
            
            return tempFile.getAbsolutePath();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        
        return null;
	}
	
	protected static String makeExecutable(
			String executable,
			String out_prefix,
			boolean parent_only) {
		final String path =  makeExecutable(executable, out_prefix);
		return parent_only ? new File(path).getParent() : path;
	}
}
