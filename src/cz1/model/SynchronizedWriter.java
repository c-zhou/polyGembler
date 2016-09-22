package cz1.model;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.concurrent.BlockingQueue;

public class SynchronizedWriter implements Runnable {

	private final BufferedWriter writer;
	private final BlockingQueue<StringBuilder> queue;
	private boolean exit = false;
	
	public SynchronizedWriter(BufferedWriter writer,
			BlockingQueue<StringBuilder> queue) {
		this.writer = writer;
		this.queue = queue;
	}
	
	public void write(String msg) {
		try {
			queue.put(
					new StringBuilder(msg));
        }
        catch (InterruptedException e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
	
	public void close() {
		this.exit = true;
	}
	
	private void flush(String data) {
        try {
            writer.write(data);
            writer.flush();         
        }
        catch(IOException ex) {
            System.out.println(
            		"IOException: EventLog.flush: "+ 
            ex.getMessage());
        }
    }
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		while(true) {
            try {
            	if(this.exit && queue.isEmpty()) {
                	System.out.println(this.exit);
            		writer.flush();
            		writer.close();
            		return; 
            	}
                StringBuilder builder = queue.take(); 
                flush(builder.toString());
            } catch(InterruptedException | IOException e) {
                e.printStackTrace();
                System.exit(0);
            }
        }
	}
}
