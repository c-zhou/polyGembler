package cz1.test;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.lang3.StringUtils;

public class TarGzioTest {

	public static void main(String[] args) {
		
		String targz = "C:\\Users\\chenxi.zhou\\Desktop\\"
				+ "genetic mapping writing-up\\metafile\\11.20.tar.gz";
		
		GZIPInputStream gzipInputStream=null;
		try {
			gzipInputStream = new GZIPInputStream( new FileInputStream(targz));
			TarArchiveInputStream is =  new TarArchiveInputStream(gzipInputStream);
			TarArchiveEntry entryx = null;

			int i=0;
			while((entryx = is.getNextTarEntry()) != null) {
			    if(entryx.isDirectory() && 
			    		entryx.getName().contains("tetraSIM10") &&
			    		StringUtils.countMatches(entryx.getName(), '/')==2) {
			    	System.out.println(++i+"\t"+entryx.getName());
			    }
			}

			entryx = new TarArchiveEntry("11.20/"
					+ "tetraSIM10.CHROM3_48148987_58833212.0423122615193808823_1_all_0_200mb/"
					+ "phasedStates/tetraSIM1.txt");
			if(is.canReadEntryData(entryx))
				System.out.println("entry exists");
			else 
				System.out.println("entry doesnot exist");
			
			OutputStream out = new FileOutputStream(new File("temp.txt"));

			out.close();
			is.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
}
