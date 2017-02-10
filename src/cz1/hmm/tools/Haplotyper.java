package cz1.hmm.tools;

import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.awt.print.PageFormat;
import java.awt.print.PrinterException;
import java.awt.print.PrinterJob;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;

import javax.imageio.ImageIO;
import javax.swing.JFrame;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

import com.itextpdf.awt.PdfGraphics2D;
import com.itextpdf.text.Document;
import com.itextpdf.text.DocumentException;
import com.itextpdf.text.Rectangle;
import com.itextpdf.text.pdf.PdfContentByte;
import com.itextpdf.text.pdf.PdfTemplate;
import com.itextpdf.text.pdf.PdfWriter;

import cz1.hmm.data.DataCollection;
import cz1.hmm.data.DataEntry;
import cz1.hmm.model.HiddenMarkovModel;
import cz1.hmm.model.HiddenMarkovModelBWT;
import cz1.hmm.model.HiddenMarkovModelVBT;
import cz1.hmm.swing.HMMFrame;
import cz1.hmm.swing.HMMPanel;
import cz1.hmm.swing.Printer;
import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Constants.Field;
import cz1.util.Executor;

public class Haplotyper extends Executor {

	private static HMMFrame hmmf = null;
	private static HMMPanel hmmp = null;
	private static String experiment = null;
	private static String workspace = null;
	private static String output = null;
	private static String[] contig = null;
	private static double[] seperation = null;
	private static boolean trainExp = false;
	private static boolean unifR = false;
	private static boolean vbt = false;
	private static boolean[] reverse = null;
	private static int max_iter = 30;
	private static int ploidy = 2;
	private static Field field = Field.PL;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		// create the command line parser

		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}

		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("e", "--experiment", true);
			myArgsEngine.add("w", "--workspace", true);
			myArgsEngine.add("o", "--output", true);
			myArgsEngine.add("c", "--contig", true);
			myArgsEngine.add("i", "--max-iteration", true);
			myArgsEngine.add("p", "--ploidy", true);
			myArgsEngine.add("P", "--parents", true);
			myArgsEngine.add("s", "--seed", true);
			myArgsEngine.add("S", "--initial-seperation", true);
			myArgsEngine.add("R", "--direction", true);
			myArgsEngine.add("G", "--genotype", false);
			myArgsEngine.add("D", "--allele-depth", false);
			myArgsEngine.add("L", "--genotype-likelihood", false);
			myArgsEngine.add("u", "--unif-rand", false);
			myArgsEngine.add("v", "--segmental-kmeans", false);
			myArgsEngine.parse(args);
		}

		if(myArgsEngine.getBoolean("-e")) {
			experiment = myArgsEngine.getString("-e");
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify your VCF file.");
		}

		if(myArgsEngine.getBoolean("-w")) {
			workspace = myArgsEngine.getString("-w");
		}

		if(myArgsEngine.getBoolean("-o")) {
			output = myArgsEngine.getString("-o");
		} else {
			output = workspace;
		}
		
		if(myArgsEngine.getBoolean("-c")) {
			contig = myArgsEngine.getString("-c").split(":");
		} else
			throw new RuntimeException("!!!");
		
		if(myArgsEngine.getBoolean("-S")) {
			String[] ss = myArgsEngine.getString("-S").split(":");
			if(ss.length<contig.length-1)
				throw new RuntimeException("!!!");
			seperation = new double[contig.length-1];
			for(int i=0; i<seperation.length; i++)
				seperation[i] = Double.parseDouble(ss[i]);
		} else {
			seperation = new double[contig.length-1];
			for(int i=0; i<seperation.length; i++)
				seperation[i] = Math.max(Math.round(
						Constants.rand.nextDouble()*
						Constants._max_initial_seperation),1);
		}
		
		if(myArgsEngine.getBoolean("-R")) {
			String[] dd = myArgsEngine.getString("-R").split(":");
			if(dd.length<contig.length)
				throw new RuntimeException("!!!");
			reverse = new boolean[contig.length];
			for(int i=0; i<reverse.length; i++)
				reverse[i] = Boolean.parseBoolean(dd[i]);
		} else {
			reverse = new boolean[contig.length];
			Arrays.fill(reverse, false);
		}
		
		if(myArgsEngine.getBoolean("-i")) {
			max_iter = Integer.parseInt(myArgsEngine.getString("-i"));
		}
		
		if(myArgsEngine.getBoolean("-t")) {
			trainExp = true;
		}
		
		if(myArgsEngine.getBoolean("-u")) {
			unifR = true;
		}
		
		if(myArgsEngine.getBoolean("-v")) {
			vbt = true;
			throw new RuntimeException("Viterbi training not supported yet!!!");
		}
		
		if(myArgsEngine.getBoolean("-p")) {
			ploidy = Integer.parseInt(myArgsEngine.getString("p"));
			Constants._ploidy_H = ploidy;
			Constants._haplotype_z = ploidy*2;
		}
		
		if(myArgsEngine.getBoolean("-P")) {
			Constants._founder_haps = myArgsEngine.getString("-P");
		} else
			throw new RuntimeException("!!!");
		if(myArgsEngine.getBoolean("s")) {
			Constants.seed = Long.parseLong(myArgsEngine.getString("s"));
			Constants.setRandomGenerator();
		}
		boolean isRF = Constants.isRF(seperation);
		if(isRF && unifR)
			for(int i=0; i<seperation.length; i++) 
				seperation[i] = Constants.haldane(seperation[i]*
						Constants.rand.nextDouble());
		else if(isRF)
			for(int i=0; i<seperation.length; i++) 
				seperation[i] = Constants.haldane(seperation[i]);
		
		int i = 0;
		if(myArgsEngine.getBoolean("-G")) {
			field = Field.GT;
			i++;
		}
		
		if(myArgsEngine.getBoolean("-D")) {
			field = Field.AD;
			i++;
		}
		
		if(myArgsEngine.getBoolean("-L")) {
			field = Field.PL;
			i++;
		}
		if(i>1) throw new RuntimeException("Options -G/--genotype, "
				+ "-D/--allele-depth, and -L/--genotype-likelihood "
				+ "are exclusive!!!");
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub

		System.out.println(Constants.seed);

		DataEntry[] de = DataCollection.readDataEntry(
				workspace+
				Constants.file_sep+
				experiment+".zip", 
				contig);

		//de.print();
		//System.exit(1);

		/*
		de.get(5);
		 */
		//for(int j=0; j<de.length; j++)
		//	for(int i=24; i<de[j].getSample().length; i++)
		//		de[j].remove(i);
		/**/
		//de.remove(new int[]{6,7,8,9,10,12,13,14});

		//if(!vbt)
		//	if(Constants._ploidy_H<6)
		//		hmm = new HiddenMarkovModelDupSLHEX(de, seperation, reverse, trainExp);
		//	else
		//		hmm = new HiddenMarkovModelDupSHEXA(de, seperation, reverse, trainExp);
		//else
		final HiddenMarkovModel hmm;
		if(vbt) {
			hmm = new HiddenMarkovModelBWT(de, seperation, reverse, trainExp, field);
		} else {
			hmm = new HiddenMarkovModelVBT(de, seperation, reverse, trainExp, field);
		}
		
		if(Constants.plot()){
			Runnable run = new Runnable(){
				public void run(){
					hmmf = new HMMFrame();
					hmmf.clearTabs();
					if(Constants.showHMM) 
						hmmp = hmmf.addHMMTab(hmm, hmm.de(), new File(workspace));
				}
			};
			Thread th = new Thread(run);
			th.run();
		}

		if(Constants.plot()){
			hmmf.pack();
			hmmf.setVisible(true);
		}

		double ll, ll0 = hmm.loglik();

		for(int i=0; i<max_iter; i++) {
			hmm.train();
			hmmp.update();
			ll = hmm.loglik();
			myLogger.info("----------loglik "+ll);
			if(ll<ll0) {
				throw new RuntimeException("!!!");
			}
			if( ll0!=Double.NEGATIVE_INFINITY && 
					Math.abs((ll-ll0)/ll0)< 1e-4)
				break;
			ll0 = ll;
		}

		if(Constants.printPlots()){
			try {
				float width = hmmf.jframe.getSize().width,
						height = hmmf.jframe.getSize().height;
				Document document = new Document(new Rectangle(width, height));
				PdfWriter writer = PdfWriter.getInstance(document, new FileOutputStream("C:\\Users\\chenxi.zhou\\Desktop\\hmmf.pdf"));
				document.open();
				PdfContentByte canvas = writer.getDirectContent();
				PdfTemplate template = canvas.createTemplate(width, height);
				Graphics2D g2d = new PdfGraphics2D(template, width, height);
				hmmf.jframe.paint(g2d);
				g2d.dispose();
				canvas.addTemplate(template, 0, 0);
				document.close();
			} catch (FileNotFoundException | DocumentException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		String contig_str = contig[0];
		for(int i=1; i<contig.length; i++) {
			if(contig_str.length()+contig[i].length()+8
					<=Constants.MAX_FILE_ID_LENGTH)
				contig_str += Constants.scaff_collapsed_str+contig[i];
			else {
				contig_str += Constants.scaff_collapsed_str+"etc"+contig.length;
				break;
			}
		}
		hmm.write(output, experiment.
				replace(".", "").
				replace("_", ""), contig_str);

	}
}
