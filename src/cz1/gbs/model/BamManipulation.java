package cz1.gbs.model;

import java.io.File;
import java.io.FileReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.BufferedReader;
import java.util.zip.GZIPInputStream;
import java.io.OutputStreamWriter;
import java.util.zip.GZIPOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.FileInputStream;
import java.io.Reader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Calendar;
import java.text.SimpleDateFormat;

public class BamManipulation {

    public final static String NLS = System.getProperty("line.separator");
    public final static String SEP = System.getProperty("file.separator");

	public static void main(String[] args) throws Exception {
		
        BamManipulation bp = new BamManipulation();
        bp.extractUnmappedFile(args[0], args[1], args[2], args[3]);
        bp.runBwaAln(args[2],args[4],args[5],Integer.parseInt(args[6]),args[7],args[3]);
    } 

	public void extractUnmappedFile (String fastq_file, String bam_file, 
            String unmapped_file, String samtoolsDir) throws Exception {
        String command = samtoolsDir+SEP+"samtools view -f4 "+bam_file+" | awk '{print $1}'";
	//System.out.println(command);
        String[] runnable = new String[]{"bash","-c",command};
        Process p = Runtime.getRuntime().exec(runnable);
        BufferedReader br0 =
            new BufferedReader(new InputStreamReader(p.getInputStream()));
        BufferedReader br1 = getBufferedReader(fastq_file);
        BufferedWriter wr = getBufferedWriter(unmapped_file);
        String identifier, fastq;
        while((identifier=br0.readLine()) != null) {
            while((fastq=br1.readLine()) !=null && 
                    !fastq.startsWith("@"+identifier)){}
            wr.write(fastq+NLS);
            for(int i=0; i<3; i++) wr.write(br1.readLine()+NLS);
        }
        wr.close();
        br0.close();
        br1.close();
    }

    public void runBwaAln (String fastq_file, String mapped_file, String reference, 
            int threads, String bwaDir, String samtoolsDir) throws Exception {
        String[] commands = new String[] {
            bwaDir+SEP+"bwa aln -t "+threads+" "+reference+" "+fastq_file+" > "+mapped_file+".bwa",
            bwaDir+SEP+"bwa samse "+reference+" "+mapped_file+".bwa "+fastq_file+" > "+mapped_file+".sam",
            samtoolsDir+SEP+"samtools view -S -b "+mapped_file+".sam > "+mapped_file+".bam"};
        String[] runnable;
        for(int i=0; i<commands.length; i++) {
            //System.out.println(commands[i]);

            runnable = new String[]{"bash","-c",commands[i]};
            Runtime.getRuntime().exec(runnable).waitFor();
        }
        new File(fastq_file).delete();
        new File(mapped_file+".bwa").delete();
        new File(mapped_file+".sam").delete();
    }       

	public BufferedWriter getGZIPBufferedWriter(String path) throws IOException {
        return new BufferedWriter(new OutputStreamWriter(new
                    GZIPOutputStream(new FileOutputStream(path))));
    }
    
    public BufferedWriter getBufferedWriter(String path) throws IOException {
        return new BufferedWriter(new FileWriter(new File(path)));
    }

	public BufferedReader getBufferedReader(String path) throws IOException {
        	BufferedReader br = null;

         	if(path.endsWith(".gz")){
             		InputStream inputFileStream = new FileInputStream(path);
             		InputStream gzipStream = new GZIPInputStream(inputFileStream);
             		Reader decoder  = new InputStreamReader(gzipStream);
             		br = new BufferedReader(decoder, 65536);
         	}else{
             		br = new BufferedReader(new FileReader(path), 65536);
         	}
         	return br;
        }

    public static String getSystemTime(){
        return new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").
            format(Calendar.getInstance().getTime());
    }
}


