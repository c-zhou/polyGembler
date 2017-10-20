package cz1.test;

import java.io.BufferedReader;
import java.io.IOException;

import cz1.util.Executor;
import cz1.util.Utils;

public class ThreadPool2 extends Executor {

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
	
	public static void main (String[] args) throws IOException {
		BufferedReader br = Utils.getBufferedReader("C:\\Users\\chenxi.zhou\\Desktop\\genetic mapping writing-up\\metafile\\comparison-new\\missing&error.txt");
		String line;
		String[][] rr = new String[24][30]; 
		while( (line=br.readLine()).startsWith("#") ) {}
		for(int i=0; i<30; i++) {
			for(int j=0; j<16; j++) {
				rr[j][i] = line.replace("Missing:", "").trim()+"/"
						+br.readLine().replace("Error:", "").trim();
				line = br.readLine();
			}
		}
		for(int i=0; i<30; i++) {
			for(int j=16; j<22; j++) {
				rr[j][i] = line.replace("Missing:", "").trim()+"/"
						+br.readLine().replace("Error:", "").trim();
				line = br.readLine();
			}
		}
		for(int i=0; i<30; i++) {
			for(int j=22; j<24; j++) {
				rr[j][i] = line.replace("Missing:", "").trim()+"/"
						+br.readLine().replace("Error:", "").trim();
				line = br.readLine();
			}
		}
		br.close();
		for(String[] r : rr) Utils.print(r);
	}
}
