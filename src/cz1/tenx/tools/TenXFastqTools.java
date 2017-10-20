package cz1.tenx.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

public class TenXFastqTools extends Executor {

	private String in_fastq;
	private String out_fastq;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--input-fastq    Input directory containing tag sequence files.\n"
						+ " -t/--threads        Threads (default is 1).\n"
						+ " -o/--output-fastq   Output directory.\n\n");
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}

		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-i", "--input-fastq", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-o", "--output-fastq", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-i")) {
			in_fastq = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your FASTQ files.");
		}

		if (myArgsEngine.getBoolean("-t")) {
			THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}
		
		if (myArgsEngine.getBoolean("-o")) {
			out_fastq = myArgsEngine.getString("-o");
		}
	}
	
	BufferedWriter out_bw;
	private final List<Long> synchedBlock = 
			Collections.synchronizedList(new LinkedList<Long>());
	
	public void write(final long i, final String o) 
			throws InterruptedException, IOException {
		synchronized (synchedBlock) {
			while (synchedBlock.get(0)!=i)
				synchedBlock.wait();
			out_bw.write(o);
			synchedBlock.remove(0);
			synchedBlock.notifyAll();
		}
	}
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		try {
			BufferedReader in_br = Utils.getBufferedReader(in_fastq);
			out_bw = Utils.getGZIPBufferedWriter(out_fastq);
			this.initial_thread_pool();
			
			final int block = 10000;
			String[][] Qs = new String[block][4];
			int k = 0;
			long bS = 0;
			String temp = in_br.readLine();
			while( temp!=null ) {
				Qs[k][0] = temp;
				Qs[k][1] = in_br.readLine();
				Qs[k][2] = in_br.readLine();
				Qs[k][3] = in_br.readLine();
				k++;
				temp = in_br.readLine();
				
				if(k==block || temp==null) {

					synchronized (synchedBlock) {synchedBlock.add(bS);}

					executor.submit(new Runnable() {
						private String[][] fastq;
						private long block_i;
						@Override
						public void run() {
							// TODO Auto-generated method stub
							try {
								final long start_i = block_i*block;
								final StringBuilder os_out = new StringBuilder();
								String[] s;
								for(int i=0; i<fastq.length; i++) {
									if(fastq[i][0]==null) break;

									s = fastq[i][0].split("\\s+");
									if(s.length==1) {
										os_out.append(fastq[i][0]);
									} else {
										os_out.append("@");
										os_out.append(s[1]);
										os_out.append(":");
										os_out.append(s[0].substring(1));
									}
									os_out.append("\n");
									os_out.append(fastq[i][1]);
									os_out.append("\n");
									os_out.append(fastq[i][2]);
									os_out.append("\n");
									os_out.append(fastq[i][3]);
									os_out.append("\n");
								}
								write(block_i, os_out.toString());
								if(start_i%1000000==0)
									myLogger.info(start_i+" reads processed");
							} catch (Exception e) {
								Thread t = Thread.currentThread();
								t.getUncaughtExceptionHandler().uncaughtException(t, e);
								e.printStackTrace();
								executor.shutdown();
								System.exit(1);
							}
						}

						public Runnable init(final String[][] Qs, final long bS) {
							this.fastq = Qs;
							this.block_i = bS;
							return(this);
						}
					}.init(Qs, bS));
					
					bS++;
					k = 0;
					Qs = new String[block][4];
				}
			}
			
			in_br.close();
			this.waitFor();
			out_bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
