package cz1.appl;

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

import cz1.data.DataCollection;
import cz1.data.DataEntry;
import cz1.model.HiddenMarkovModelBWT;
import cz1.model.HiddenMarkovModel;
import cz1.model.HiddenMarkovModelVBT;
import cz1.swing.HMMFrame;
import cz1.swing.HMMPanel;
import cz1.swing.Printer;
import cz1.util.Constants;

public class PolyMap {

	private static HMMFrame hmmf = null;
	private static HMMPanel hmmp = null;
	private static String experiment = null;
	private static String workspace = null;
	private static String output = null;
	private static String[] contig = null;
	private static double[] seperation = null;
	private static boolean dupS = false;
	private static boolean trainExp = false;
	private static boolean unifR = false;
	private static boolean vbt = false;
	private static boolean[] reverse = null;
	private static int max_iter = 30;
	private static int ploidy = 2;
	
	public static void main(String[] args) {

		// create the command line parser
		CommandLineParser parser = new PosixParser();
		
		// create the Options
		Options options = new Options();
		options.addOption( "e", "experiment", true, "experiment name." );
		options.addOption( "w", "workspace", true, "directory contains input files." );
		options.addOption( "o", "output", true, "output directory name." );
		options.addOption( "c", "contig", true, "contig id.");
		options.addOption( "i", "max-iteration", true, "maximum iterations.");
		options.addOption( "d", "duplicate-hs", false, "duplicate hidden states.");
		options.addOption( "p", "ploidy", true, "ploidy of genome.");
		options.addOption( "P", "parents", true, "parental haplotypes.");
		options.addOption( "s", "seed", true, "random seed.");
		options.addOption( "S", "initial-seperation", true, "initial seperation.");
		options.addOption( "R", "direction", true, "direction.");
		options.addOption( "t", "train-exp", false, "train the exp parameter for RF.");
		options.addOption( "u", "unif-rand", false, "sample initial RF from a uniform distribution.");
		options.addOption( "v", "segmental-kmeans", false, "using segmental k-means training.");
		
		try {
			// parse the command line arguments
			CommandLine line = parser.parse( options, args );
			if( line.hasOption("e") ) {
				experiment = line.getOptionValue('e');
			} else
				throw new RuntimeException("!!!");
			if(line.hasOption("w")) {
				workspace = line.getOptionValue("w");
			} else
				throw new RuntimeException("!!!");
			if(line.hasOption("o")) {
				output = line.getOptionValue("o");
			} else {
				output = workspace;
			}
			if(line.hasOption("c")) {
				contig = line.getOptionValue("c").split(":");
			} else
				throw new RuntimeException("!!!"); 
			if(line.hasOption("S")) {
				String[] ss = line.getOptionValue("S").split(":");
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
			if(line.hasOption("R")) {
				String[] dd = line.getOptionValue("R").split(":");
				if(dd.length<contig.length)
					throw new RuntimeException("!!!");
				reverse = new boolean[contig.length];
				for(int i=0; i<reverse.length; i++)
					reverse[i] = Boolean.parseBoolean(dd[i]);
			} else {
				reverse = new boolean[contig.length];
				Arrays.fill(reverse, false);
			}
			if(line.hasOption("i")) {
				max_iter = Integer.parseInt(line.getOptionValue("i"));
			}
			if(line.hasOption("d")) {
				dupS = true;
			}
			if(line.hasOption("t")) {
				trainExp = true;
			}
			if(line.hasOption("u")) {
				unifR = true;
			}
			if(line.hasOption("v")) {
				vbt = true;
			}
			if(line.hasOption("p")) {
				ploidy = Integer.parseInt(line.getOptionValue("p"));
				Constants._ploidy_H = ploidy;
				Constants._haplotype_z = ploidy*2;
			}
			if(line.hasOption("P")) {
				Constants._founder_haps = line.getOptionValue("P");
			} else
				throw new RuntimeException("!!!");
			if(line.hasOption("s")) {
				Constants.seed = Long.parseLong(line.getOptionValue("s"));
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
		} catch( ParseException exp ) {
			System.out.println( "Unexpected exception:" + exp.getMessage() );
		}

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
		final HiddenMarkovModelBWT hmm = 
				new HiddenMarkovModelBWT(de, seperation, reverse, trainExp);

		hmm.print(true);
		
		if(Constants.plot()>=1){
			//SwingUtilities.in
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
		
		if(Constants.plot()>=1){
			hmmf.pack();
			hmmf.setVisible(true);
		}
		
		double ll, ll0 = hmm.loglik();
		
		for(int i=0; i<max_iter; i++) {
			//System.out.println(i);
			hmm.train();
			hmmp.update();
			
			ll = hmm.loglik();
			
			System.err.println("----------loglik "+ll);
			//hmm.print(true);
			
			if(ll<ll0) {
				throw new RuntimeException("!!!");
			}
			
			if( ll0!=Double.NEGATIVE_INFINITY && 
					Math.abs((ll-ll0)/ll0)< 1e-4)
				break;
			ll0 = ll;
			//hmm.print(true);
		}
		hmm.print(true);
		
		/**
		PrinterJob pjob = PrinterJob.getPrinterJob();
		PageFormat preformat = pjob.defaultPage();
		preformat.setOrientation(PageFormat.LANDSCAPE);
		PageFormat postformat = pjob.pageDialog(preformat);
		//If user does not hit cancel then print.
		if (preformat != postformat) {
		    //Set print component
		    pjob.setPrintable(new Printer(hmmf.jframe), postformat);
		    if (pjob.printDialog()) {
		        try {
					pjob.print();
				} catch (PrinterException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
		    }
		}
		**/
		if(true){
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
