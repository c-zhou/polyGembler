package cz1.util;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.BufferedReader;
import java.io.FilenameFilter;
import java.io.PrintStream;
import java.util.zip.GZIPInputStream;
import java.io.OutputStreamWriter;
import java.util.zip.GZIPOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.Calendar;
import java.text.SimpleDateFormat;

import org.apache.commons.lang3.ArrayUtils;

public class IO {
	public static BufferedReader getBufferedReader(String path) throws IOException {
		BufferedReader br = null;
		if(path.endsWith(".gz") || path.endsWith(".zip")){
			br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new 
					FileInputStream(path))), 65536);
		}else{
			br = new BufferedReader(new FileReader(path), 65536);
		}
		return br;
	}

	public static BufferedWriter getBufferedWriter(String path) throws IOException {
		return new BufferedWriter(new FileWriter(new File(path)));
	}
	
	public static BufferedWriter getBufferedWriter(String path, boolean append) 
			throws IOException {
		return new BufferedWriter(new FileWriter(new File(path), append));
	}

	public static BufferedWriter getGZIPBufferedWriter(String path) throws IOException {
		return new BufferedWriter(new OutputStreamWriter(new
				GZIPOutputStream(new FileOutputStream(path))));
	}

	public static String getSystemTime(){
		return new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").
				format(Calendar.getInstance().getTime());
	}

	public static void println() {
		System.out.println();
	}

	public static void println(char character) {
		System.out.println(character);
	}

	public static void print(char character) {
		System.out.print(character);
	}
	
	public static void print(String message) {
		System.out.print(message);
	}


	public static void println(String message) {
		print(message+"\n");
	}

	public static void print(ArrayList<Integer> array) {
		for(int i=0; i<array.size(); i++)
			print(array.get(i)+"\t");
		println();
	}

	public static void print(int[][] matrix) {
		for(int i=0; i<matrix.length; i++)
			print(matrix[i]);
	}

	public static void print(int a) {
		print(a+"");
	}

	public static void println(int a) {
		println(a+"");
	}

	public static void print(int[] array) {
		if(array==null) {
			IO.println("null");
			return;
		}
		for(int i=0; i<array.length; i++)
			print(array[i]+"\t");
		println();
	}

	public static void print(double a) {
		print(a+"\t");
	}

	public static void println(double a) {
		println(a+"");
	}

	public static void print(double[] array) {
		for(int i=0; i<array.length; i++)
			print(array[i]+"\t");
		println();
	}

	public static void print(Double[] array) {
		print(ArrayUtils.toPrimitive(array));
	}

	public static void print(double[][] matrix) {
		if(matrix!=null)
			for(int i=0; i<matrix.length; i++)
				print(matrix[i]);
	}

	public static void print(long[] array) {
		for(int i=0; i<array.length; i++)
			print(array[i]+"\t");
		println();
	}

	public static void print(Set<Character> set) {
		for(Character c : set) print(c);
		println();
	}

	public static void print(List<Character> list) {
		for(Character c : list) print(c);
		println();
	}

	public static void print(Character[] array) {
		for(Character c : array) print(c);
		println();
	}

	public static void print(Map<Character, Map<Character, Double>> map) {
		// TODO Auto-generated method stub
		Map<Character, Double> map0;
		for(Character c0 : map.keySet()) {
			print(c0+" - ");
			map0 = map.get(c0);
			for(Character c1 : map0.keySet())
				print(map0.get(c1)+"\t");
			println();
		}
	}

	public static String readFastaFileByID(String filePath, String id) {
		String sequence = null;
		try {
			BufferedReader br = getBufferedReader(filePath);
			String line;
			while( (line = br.readLine()) != null ) {
				if(line.startsWith(">") && line.equals(">"+id)) {
					StringBuilder sb = new StringBuilder();
					while( (line=br.readLine()) != null && !line.startsWith(">") ) 
						sb.append(br.readLine());
					sequence = sb.toString();
					return sequence;
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
		return sequence;
	}

	public static void print(char[] chars) {
		// TODO Auto-generated method stub
		for(char c : chars) IO.print(c+" ");
		IO.println();
	}
	
	public static void makeOutputDir(String dir) {
		// TODO Auto-generated method stub
		File out = new File(dir);
		if(!out.exists() || out.exists()&&!out.isDirectory()) {
			out.mkdir();
		}
	}

	public static void print(String[] strs) {
		// TODO Auto-generated method stub
		for(String str : strs) IO.print(str+" ");
		IO.println();
	}
}
