package cz1.test;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.geom.Ellipse2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.WindowConstants;

import org.apache.commons.lang3.StringEscapeUtils;
import org.apache.commons.math.stat.StatUtils;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.AxisSpace;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.imagemap.ImageMapUtilities;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.Range;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import com.itextpdf.awt.PdfGraphics2D;
import com.itextpdf.text.Document;
import com.itextpdf.text.DocumentException;
import com.itextpdf.text.Rectangle;
import com.itextpdf.text.pdf.PdfContentByte;
import com.itextpdf.text.pdf.PdfTemplate;
import com.itextpdf.text.pdf.PdfWriter;

import cz1.util.Executor;
import cz1.util.Utils;

public class JfreeChart extends Executor {

	private JFrame jframe;
	private Color background_color = Color.WHITE,
			line_color = Color.BLACK,
			font_color = Color.BLACK;
	
	public JfreeChart() {
		jframe = new JFrame();
		jframe.getContentPane().setLayout(new BorderLayout());
		jframe.setSize(new Dimension(1100,400));
		jframe.setMinimumSize(new Dimension(1100,400));
		jframe.setPreferredSize(new Dimension(1100,400));
		jframe.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
		jframe.setVisible(false);
	}
	
	public JfreeChart(Dimension dim) {
		jframe = new JFrame();
		jframe.getContentPane().setLayout(new BorderLayout());
		jframe.setSize(dim);
		jframe.setMinimumSize(dim);
		jframe.setPreferredSize(dim);
		jframe.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
		jframe.setVisible(false);
	}
	
	public JfreeChart(Dimension dim, GridLayout gridLayout) {
		// TODO Auto-generated constructor stub
		jframe = new JFrame();
		jframe.getContentPane().setLayout(gridLayout);
		jframe.setSize(dim);
		jframe.setMinimumSize(dim);
		jframe.setPreferredSize(dim);
		jframe.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
		jframe.setVisible(false);
	}

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
	
	private String[] files_in;

	final Font font3 = new Font("Dialog", Font.BOLD, 14); 
	final Font font4 = new Font("Dialog", Font.PLAIN, 14); 

	final static String formatStr2= "%1.3f";
	private void plotLineChart(String file_in,
			Dimension dim,
			Color color,
			String position,
			String title) {
		JPanel collinearChart = new JPanel(); 
		collinearChart.setSize(dim);
		collinearChart.setMinimumSize(dim);
		collinearChart.setPreferredSize(dim);
		collinearChart.setBackground(color);
		
		final XYSeriesCollection lineSeries = new XYSeriesCollection();
		int n1 = 0;
		double x_max = Double.NEGATIVE_INFINITY, 
				y_max = Double.NEGATIVE_INFINITY;
		try {
			BufferedReader br = Utils.getBufferedReader(file_in);
			String line = br.readLine();
			String[] s;
			int series = 0;
			while( (line=br.readLine())!=null &&
					line.length()>0) {
				s = line.split("\\s+");
				series = Integer.parseInt(s[2]);
				if(lineSeries.getSeriesCount()<series) {
					lineSeries.addSeries(new XYSeries(series));
					n1 = 0;
				}
				double x = Double.parseDouble(s[1])/1000,
						y = Double.parseDouble(s[0]);
				lineSeries.getSeries(series-1).add(x,y);
				if(x>x_max) x_max = x;
				if(y>y_max) y_max = y;
				n1++;
			}
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		final int n = n1;
		final int m = lineSeries.getSeriesCount();
		// create error bars
		final XYSeriesCollection errorBarSeries2 = new XYSeriesCollection();
		for(int i=0; i<n; i++) {
			double[] f = new double[m];
			for(int j=0; j<m; j++)
				f[j] = lineSeries.getSeries(j).getY(i).doubleValue();
			double mean = StatUtils.mean(f),
					sd = Math.sqrt(StatUtils.variance(f));
			double x = lineSeries.getSeries(0).getX(i).doubleValue();
			
			errorBarSeries2.addSeries(new XYSeries(i));
			errorBarSeries2.getSeries(i).add(x, mean-sd);
			errorBarSeries2.getSeries(i).add(x, mean);
			errorBarSeries2.getSeries(i).add(x, mean+sd);
		}
		
		final XYSeriesCollection upperBoundSeries = new XYSeriesCollection();
		final XYSeriesCollection upperQ3Series = new XYSeriesCollection();
		final XYSeriesCollection boxSeriesE = new XYSeriesCollection();
		final XYSeriesCollection boxSeriesN = new XYSeriesCollection();
		final XYSeriesCollection boxSeriesW = new XYSeriesCollection();
		final XYSeriesCollection boxSeriesS = new XYSeriesCollection();
		final XYSeriesCollection q2Series = new XYSeriesCollection();
		final XYSeriesCollection lowerQ2Series = new XYSeriesCollection();
		final XYSeriesCollection lowerBoundSeries = new XYSeriesCollection();		
		final XYSeriesCollection meanDotSeries = new XYSeriesCollection();
		final XYSeriesCollection fillingSeries = new XYSeriesCollection();
		
		for(int i=0; i<n; i++) {
			double[] f = new double[m];
			for(int j=0; j<m; j++)
				f[j] = lineSeries.getSeries(j).getY(i).doubleValue();
			double q1 = StatUtils.percentile(f, 25),
					q2 = StatUtils.percentile(f, 50),
					q3 = StatUtils.percentile(f, 75);
			double IQR = q3-q1;
			double[] boundary = new double[]{q1-1.5*IQR, q3+1.5*IQR};
			final Set<Integer> outliers = new HashSet<Integer>();
			for(int k=0; k<m; k++) 
				if(f[k]<boundary[0] || f[k]>boundary[1])
					outliers.add(k);
			double lower_bound = Double.POSITIVE_INFINITY, 
					upper_bound = Double.NEGATIVE_INFINITY;
			for(int k=0; k<m; k++) {
				if(outliers.contains(k)) continue;
				if(f[k]<lower_bound) lower_bound = f[k];
				if(f[k]>upper_bound) upper_bound = f[k];
			}
			double x = lineSeries.getSeries(0).getX(i).doubleValue();
			double mean = StatUtils.mean(f);
			
			if(mean>=thres) {
				int w = upperBoundSeries.getSeriesCount();
				upperBoundSeries.addSeries(new XYSeries(i));
				upperBoundSeries.getSeries(w).add(x-half_width, upper_bound);
				upperBoundSeries.getSeries(w).add(x+half_width, upper_bound);
				upperQ3Series.addSeries(new XYSeries(i));
				upperQ3Series.getSeries(w).add(x, upper_bound);
				upperQ3Series.getSeries(w).add(x, q3);
				boxSeriesE.addSeries(new XYSeries(i));
				boxSeriesE.getSeries(w).add(x+half_width, q1);
				boxSeriesE.getSeries(w).add(x+half_width, q3);
				boxSeriesN.addSeries(new XYSeries(i));
				boxSeriesN.getSeries(w).add(x+half_width, q3);
				boxSeriesN.getSeries(w).add(x-half_width, q3);
				boxSeriesW.addSeries(new XYSeries(i));
				boxSeriesW.getSeries(w).add(x-half_width, q3);
				boxSeriesW.getSeries(w).add(x-half_width, q1);
				boxSeriesS.addSeries(new XYSeries(i));
				boxSeriesS.getSeries(w).add(x-half_width, q1);
				boxSeriesS.getSeries(w).add(x+half_width, q1);
				q2Series.addSeries(new XYSeries(i));
				q2Series.getSeries(w).add(x-half_width, q2);
				q2Series.getSeries(w).add(x+half_width, q2);
				lowerQ2Series.addSeries(new XYSeries(i));
				lowerQ2Series.getSeries(w).add(x, q1);
				lowerQ2Series.getSeries(w).add(x, lower_bound);
				lowerBoundSeries.addSeries(new XYSeries(i));
				lowerBoundSeries.getSeries(w).add(x-half_width, lower_bound);
				lowerBoundSeries.getSeries(w).add(x+half_width, lower_bound);
				
				double x_inter = x-half_width, x_end = x+half_width;
				while( x_inter<x_end) {
					XYSeries tmp = new XYSeries("shading");
					tmp.add(x_inter, q1);
					tmp.add(x_inter, q3);
					fillingSeries.addSeries(tmp);
					x_inter += 0.01;
				}
			}
			meanDotSeries.addSeries(new XYSeries(i));
			meanDotSeries.getSeries(i).add(x, mean);
		}
		
		final JFreeChart chart = ChartFactory.createXYLineChart(
				null,
				null, // domain axis label
				null, // range axis label
				null, // data
				PlotOrientation.VERTICAL, 
				false, // include legend
				true, // tooltips?
				false // URL generator? Not required...
				);
		final NumberAxis domainAxis = new NumberAxis("SNP Physical Position (Kb) ["+n+" SNPs]");
		final NumberAxis rangeAxis = new NumberAxis("Estimated RF");
		domainAxis.setTickUnit(new NumberTickUnit(50));
		domainAxis.setRange(-x_max/20/dim.getWidth()*dim.getHeight(), 
				x_max*(1+0.1/dim.getWidth()*dim.getHeight()));
		rangeAxis.setRange(-y_max/20, y_max*1.1);
		domainAxis.setLabelFont(font3);
		rangeAxis.setLabelFont(font3);
		domainAxis.setTickLabelFont(font4);		
		rangeAxis.setTickLabelFont(font4);
		
		XYPlot plot = chart.getXYPlot();
		plot.setDomainAxis(domainAxis);
		plot.setRangeAxis(rangeAxis);
		
		for(int i=0; i<meanDotSeries.getSeriesCount(); i++) {
			double mean = meanDotSeries.getSeries(i).getY(0).doubleValue();
			if(mean<thres) continue;
			double x = meanDotSeries.getSeries(i).getX(0).doubleValue();
			XYTextAnnotation anno = new XYTextAnnotation(
					String.format(formatStr2, mean)+"\u00B1"+
							String.format(formatStr2,
									mean-errorBarSeries2.getSeries(i).getY(0).doubleValue()),
					x>x_max/2?x-30:x+30, mean);
			anno.setFont(font3);
			plot.addAnnotation((XYAnnotation) anno);
		}
		final ChartPanel cp =   new ChartPanel(
				chart,
				collinearChart.getWidth(),
				collinearChart.getHeight(),
				collinearChart.getWidth(),
				collinearChart.getHeight(),
				collinearChart.getWidth(),
				collinearChart.getHeight(),
				ChartPanel.DEFAULT_BUFFER_USED,
				true, 
				true,
				true,  
				true,
				true
				);
		

		plot.setDataset(0, meanDotSeries);
		plot.setRenderer(0, getRenderer1());
		
		plot.setDataset(1, upperQ3Series);
		plot.setRenderer(1, getRenderer3());
		plot.setDataset(2, lowerQ2Series);
		plot.setRenderer(2, getRenderer3());
		
		plot.setDataset(3, boxSeriesE);
		plot.setRenderer(3, getRenderer2());
		plot.setDataset(4, boxSeriesN);
		plot.setRenderer(4, getRenderer2());
		plot.setDataset(5, boxSeriesW);
		plot.setRenderer(5, getRenderer2());
		plot.setDataset(6, boxSeriesS);
		plot.setRenderer(6, getRenderer2());
		plot.setDataset(7, q2Series);
		plot.setRenderer(7, getRenderer2());
		plot.setDataset(8, upperBoundSeries);
		plot.setRenderer(8, getRenderer2());
		plot.setDataset(9, lowerBoundSeries);
		plot.setRenderer(9, getRenderer2());
		
		plot.setDataset(10, fillingSeries);
		plot.setRenderer(10, getRenderer4());
		
		plot.setDataset(99, lineSeries);
		plot.setRenderer(99, getRenderer0());
		
		/**
		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer(){
			Stroke soild = new BasicStroke(2.0f);
			Stroke dashed =  new BasicStroke(1.0f,BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, new float[] {10.0f}, 0.0f);
			@Override
			public Stroke getItemStroke(int row, int column) {
				if (row>=N[0]&&row<N[1]){
					return new BasicStroke(stroke*2);
				}  else if(row>=N[1]&&row<N[2]) {
					return new BasicStroke(stroke*2);
				} else {
					return new BasicStroke(stroke/2);
				}
			}
		};
		
		
		//XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		for(int i=0; i<seriesCollection.getSeriesCount(); i++) {
			if(i>=N[3]) {
				renderer.setSeriesPaint(i, Color.gray);
			} else if(i>=N[2]) {
				renderer.setSeriesPaint(i, Color.decode(hex_str[0]));
			} else if(i>=N[1]) {
				renderer.setSeriesPaint(i, Color.decode(hex_str[0]));
			} else if(i>=N[0]) {
				renderer.setSeriesPaint(i, Color.black);
			} else {
				renderer.setSeriesPaint(i, Color.decode(hex_str[4]));
				renderer.setSeriesShapesVisible(i,true);
				renderer.setSeriesShapesFilled(i,true);
				renderer.setSeriesShape(i, shape);
			}
				
		}
		**/
		
		chart.setBackgroundPaint(Color.WHITE);
		chart.setBorderPaint(Color.WHITE);
		chart.setTitle(title);
		cp.setPreferredSize(dim);
		cp.setBorder(null);
		collinearChart.add(cp);
		jframe.add(collinearChart,position);
	}
	

	final Shape shape = new Ellipse2D.Double(-4.0,-4.0, 8.0, 8.0);
	final String[] hex_str = new String[]{"#F8766D","#7CAE00","#00BFC4","#C77CFF",
			"#B55D60", "#857AAA", "#5F9E6E", "#5975A4"};
	final double half_width = 3;
	final double thres = 0.05;
	final int stroke = 2;
	
	private XYLineAndShapeRenderer getRenderer0() {
		// TODO Auto-generated method stub
		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		renderer.setShapesVisible(false); 
		renderer.setShapesFilled(false);
		renderer.setPaint(Color.gray);
		return renderer;
	}
	
	private XYLineAndShapeRenderer getRenderer1() {
		// TODO Auto-generated method stub
		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		renderer.setPaint(Color.decode(hex_str[4]));
		renderer.setShapesVisible(true);
		renderer.setShapesFilled(true);
		renderer.setShape(shape);
		return renderer;
	}
	
	private XYLineAndShapeRenderer getRenderer2() {
		// TODO Auto-generated method stub
		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		renderer.setStroke(new BasicStroke(stroke));
		renderer.setPaint(Color.DARK_GRAY);
		renderer.setShapesVisible(false);
		renderer.setShapesFilled(false);
		return renderer;
	}
	
	
	private XYLineAndShapeRenderer getRenderer3() {
		// TODO Auto-generated method stub
		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		renderer.setStroke(new BasicStroke(stroke));
		renderer.setPaint(Color.DARK_GRAY);
		renderer.setShapesVisible(false);
		renderer.setShapesFilled(false);
		return renderer;
	}

	private XYItemRenderer getRenderer4() {
		// TODO Auto-generated method stub
		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		renderer.setStroke(new BasicStroke(stroke/4));
		Color c = Color.decode(hex_str[0]);
		renderer.setPaint(c);
		renderer.setShapesVisible(false);
		renderer.setShapesFilled(false);
		return renderer;
	}
	
	private XYLineAndShapeRenderer getRenderer5() {
		// TODO Auto-generated method stub
		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		renderer.setStroke(new BasicStroke(3));
		renderer.setPaint(Color.DARK_GRAY);
		renderer.setShapesVisible(false);
		renderer.setShapesFilled(false);
		return renderer;
	}
	
	private void plotLineChart2(String file_in,
			Dimension dim,
			Color color,
			String position) {
		JPanel collinearChart = new JPanel(); 
		collinearChart.setSize(dim);
		collinearChart.setMinimumSize(dim);
		collinearChart.setPreferredSize(dim);
		collinearChart.setBackground(color);
		
		final XYSeriesCollection lineSeries = new XYSeriesCollection();
		int n1 = 0;
		try {
			BufferedReader br = Utils.getBufferedReader(file_in);
			String line = br.readLine();
			String[] s;
			int series = 0;
			while( (line=br.readLine())!=null &&
					line.length()>0) {
				s = line.split("\\s+");
				series = Integer.parseInt(s[2]);
				if(lineSeries.getSeriesCount()<series) {
					lineSeries.addSeries(new XYSeries(series));
					n1 = 0;
				}
				lineSeries.getSeries(series-1).add(
						Double.parseDouble(s[1])/1000, 
						Double.parseDouble(s[0]));
				n1++;
			}
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		final int n = n1;
		final int m = lineSeries.getSeriesCount();
		// create error bars
		final XYSeriesCollection errorBarSeries2 = new XYSeriesCollection();
		for(int i=0; i<n; i++) {
			double[] f = new double[m];
			for(int j=0; j<m; j++)
				f[j] = lineSeries.getSeries(j).getY(i).doubleValue();
			double mean = StatUtils.mean(f),
					sd = Math.sqrt(StatUtils.variance(f));
			double x = lineSeries.getSeries(0).getX(i).doubleValue();
			
			errorBarSeries2.addSeries(new XYSeries(i));
			errorBarSeries2.getSeries(i).add(x, mean-sd);
			errorBarSeries2.getSeries(i).add(x, mean);
			errorBarSeries2.getSeries(i).add(x, mean+sd);
		}
		
		final XYSeriesCollection errorBarSeries = new XYSeriesCollection();
		
		final XYSeriesCollection upperBoundSeries = new XYSeriesCollection();
		final XYSeriesCollection upperQ3Series = new XYSeriesCollection();
		final XYSeriesCollection boxSeries = new XYSeriesCollection();
		final XYSeriesCollection q2Series = new XYSeriesCollection();
		final XYSeriesCollection lowerQ2Series = new XYSeriesCollection();
		final XYSeriesCollection lowerBoundSeries = new XYSeriesCollection();		
		final XYSeriesCollection boxplotSeries = new XYSeriesCollection();	
		final XYSeriesCollection meanDotSeries = new XYSeriesCollection();
		
		final double half_width = 1;
		final double thres = 0.05;
		
		final int stroke = 2;
		
		for(int i=0; i<n; i++) {
			double[] f = new double[m];
			for(int j=0; j<m; j++)
				f[j] = lineSeries.getSeries(j).getY(i).doubleValue();
			double q1 = StatUtils.percentile(f, 25),
					q2 = StatUtils.percentile(f, 50),
					q3 = StatUtils.percentile(f, 75);
			double IQR = q3-q1;
			double[] boundary = new double[]{q1-1.5*IQR, q3+1.5*IQR};
			final Set<Integer> outliers = new HashSet<Integer>();
			for(int k=0; k<m; k++) 
				if(f[k]<boundary[0] || f[k]>boundary[1])
					outliers.add(k);
			double lower_bound = Double.POSITIVE_INFINITY, 
					upper_bound = Double.NEGATIVE_INFINITY;
			for(int k=0; k<m; k++) {
				if(outliers.contains(k)) continue;
				if(f[k]<lower_bound) lower_bound = f[k];
				if(f[k]>upper_bound) upper_bound = f[k];
			}
			double x = lineSeries.getSeries(0).getX(i).doubleValue();
			double mean = StatUtils.mean(f);
			
			if(mean>=thres) {
				int w = upperBoundSeries.getSeriesCount();
				upperBoundSeries.addSeries(new XYSeries(i));
				upperBoundSeries.getSeries(w).add(x-half_width, upper_bound);
				upperBoundSeries.getSeries(w).add(x+half_width, upper_bound);
				upperQ3Series.addSeries(new XYSeries(i));
				upperQ3Series.getSeries(w).add(x, upper_bound);
				upperQ3Series.getSeries(w).add(x, q3);
				boxSeries.addSeries(new XYSeries(i));



				boxSeries.getSeries(w).add(x-half_width, q1);
		
				boxSeries.getSeries(w).add(x-half_width, q3);
				boxSeries.getSeries(w).add(x+half_width, q3);
				boxSeries.getSeries(w).add(x+half_width, q1);	
				
				q2Series.addSeries(new XYSeries(i));
				q2Series.getSeries(w).add(x-half_width, q2);
				q2Series.getSeries(w).add(x+half_width, q2);
				lowerQ2Series.addSeries(new XYSeries(i));
				lowerQ2Series.getSeries(w).add(x, q1);
				lowerQ2Series.getSeries(w).add(x, lower_bound);
				lowerBoundSeries.addSeries(new XYSeries(i));
				lowerBoundSeries.getSeries(w).add(x-half_width, lower_bound);
				lowerBoundSeries.getSeries(w).add(x+half_width, lower_bound);
			}
			meanDotSeries.addSeries(new XYSeries(i));
			meanDotSeries.getSeries(i).add(x, mean);
		}
		
		final XYSeriesCollection seriesCollection = new XYSeriesCollection();
		final int[] N = new int[10];
		
		N[0]+=addSeries(seriesCollection, meanDotSeries);
		N[1] = N[0];
		N[1]+=addSeries(seriesCollection, boxplotSeries);
		N[2] = N[1];
		N[2]+=addSeries(seriesCollection, upperBoundSeries);
		N[2]+=addSeries(seriesCollection, boxSeries);
		N[2]+=addSeries(seriesCollection, q2Series);
		N[2]+=addSeries(seriesCollection, lowerBoundSeries);
		N[3] = N[2];
		N[3]+=addSeries(seriesCollection, upperQ3Series);
		N[3]+=addSeries(seriesCollection, lowerQ2Series);
		addSeries(seriesCollection, lineSeries);
		
		final JFreeChart chart = ChartFactory.createXYLineChart(
				null,
				null, // domain axis label
				null, // range axis label
				seriesCollection, // data
				PlotOrientation.VERTICAL, 
				false, // include legend
				true, // tooltips?
				false // URL generator? Not required...
				);
		final NumberAxis domainAxis = new NumberAxis("SNP Physical Position (Kb) ["+n+" SNPs]");
		final NumberAxis rangeAxis = new NumberAxis("Estimated RF");
		domainAxis.setTickUnit(new NumberTickUnit(50));
		rangeAxis.setRange(-0.02, 0.48);
		chart.getXYPlot().setDomainAxis(domainAxis);
		chart.getXYPlot().setRangeAxis(rangeAxis);
		
		final ChartPanel cp =   new ChartPanel(
				chart,
				collinearChart.getWidth(),
				collinearChart.getHeight(),
				collinearChart.getWidth(),
				collinearChart.getHeight(),
				collinearChart.getWidth(),
				collinearChart.getHeight(),
				ChartPanel.DEFAULT_BUFFER_USED,
				true, 
				true,
				true,  
				true,
				true
				);
		
		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer(){
			Stroke soild = new BasicStroke(2.0f);
			Stroke dashed =  new BasicStroke(1.0f,BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, new float[] {10.0f}, 0.0f);
			@Override
			public Stroke getItemStroke(int row, int column) {
				if (row>=N[0]&&row<N[1]){
					return new BasicStroke(stroke*2);
				}  else if(row>=N[1]&&row<N[2]) {
					return new BasicStroke(stroke*2);
				} else {
					return new BasicStroke(stroke/2);
				}
			}
		};
		
		final Shape shape = new Ellipse2D.Double(-3.0,-3.0, 6.0, 6.0);
		final String[] hex_str = new String[]{"#F8766D","#7CAE00","#00BFC4","#C77CFF",
				"#B55D60", "#857AAA", "#5F9E6E", "#5975A4"};
		//XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		for(int i=0; i<seriesCollection.getSeriesCount(); i++) {
			if(i>=N[3]) {
				renderer.setSeriesPaint(i, Color.gray);
			} else if(i>=N[2]) {
				renderer.setSeriesPaint(i, Color.decode(hex_str[0]));
			} else if(i>=N[1]) {
				renderer.setSeriesPaint(i, Color.decode(hex_str[0]));
			} else if(i>=N[0]) {
				renderer.setSeriesPaint(i, Color.black);
			} else {
				renderer.setSeriesPaint(i, Color.decode(hex_str[4]));
				renderer.setSeriesShapesVisible(i,true);
				renderer.setSeriesShapesFilled(i,true);
				renderer.setSeriesShape(i, shape);
			}
				
		}
		renderer.setBaseShapesVisible(false);
		renderer.setBaseShapesFilled(false);
		chart.getXYPlot().setRenderer(renderer);
		//chart.getXYPlot().setDatasetRenderingOrder(DatasetRenderingOrder.FORWARD);
		
		chart.setBackgroundPaint(Color.WHITE);
		chart.setBorderPaint(Color.WHITE);
		cp.setPreferredSize(dim);
		cp.setBorder(null);
		collinearChart.add(cp);
		jframe.add(collinearChart,position);
	}
	
	private int addSeries(XYSeriesCollection seriesCollection,
			XYSeriesCollection extraSeriesCollection) {
		// TODO Auto-generated method stub
		for(int i=0; i<extraSeriesCollection.getSeriesCount(); i++)
			seriesCollection.addSeries(extraSeriesCollection.getSeries(i));
		return extraSeriesCollection.getSeriesCount();
	}

	private void plotCollinearChart(final String file_in,
			Dimension dim,
			Color color,
			String position,
			String title) {
		JPanel collinearChart = new JPanel(); 
		collinearChart.setSize(dim);
		collinearChart.setMinimumSize(dim);
		collinearChart.setPreferredSize(dim);
		collinearChart.setBackground(color);
		
		final XYSeriesCollection blocks = new XYSeriesCollection();
		try {
			BufferedReader br = Utils.getBufferedReader(file_in);
			String line = br.readLine();
			String[] s;
			int series = 0;
			while( (line=br.readLine())!=null &&
					line.length()>0) {
				blocks.addSeries(new XYSeries(series));
				s = line.split("\\s+");
				blocks.getSeries(series).add(
						Double.parseDouble(s[0])/1000000, 
						Double.parseDouble(s[1])/1000);
				line = br.readLine();
				s = line.split("\\s+");
				blocks.getSeries(series).add(
						Double.parseDouble(s[0])/1000000, 
						Double.parseDouble(s[1])/1000);
				br.readLine();
				series++;
			}
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		final JFreeChart chart = ChartFactory.createXYLineChart(
				null,
				null, // domain axis label
				null, // range axis label
				blocks, // data
				PlotOrientation.VERTICAL, 
				false, // include legend
				true, // tooltips?
				false // URL generator? Not required...
				);
		
		String[] s = new File(file_in).getName().split("\\.");
		
		final NumberAxis domainAxis = new NumberAxis("NSP306v2 scaffold "
				+ s[1].replaceAll("^0{0,5}", "")+"(Mb)");
		final NumberAxis rangeAxis = new NumberAxis("ITR_r1.0 scaffold "
				+ s[5].replaceAll("^sc0{0,6}", "").replaceAll(".1$", "")+"(Kb)");
		rangeAxis.setRange(-10, 340);
		chart.getXYPlot().setDomainAxis(domainAxis);
		chart.getXYPlot().setRangeAxis(rangeAxis);
		domainAxis.setLabelFont(font3);
		rangeAxis.setLabelFont(font3);
		domainAxis.setTickLabelFont(font4);		
		rangeAxis.setTickLabelFont(font4);
		
		final ChartPanel cp =   new ChartPanel(
				chart,
				collinearChart.getWidth(),
				collinearChart.getHeight(),
				collinearChart.getWidth(),
				collinearChart.getHeight(),
				collinearChart.getWidth(),
				collinearChart.getHeight(),
				ChartPanel.DEFAULT_BUFFER_USED,
				true, 
				true,
				true,  
				true,
				true
				);
		chart.getXYPlot().setRenderer(getRenderer5());
		chart.setBackgroundPaint(Color.WHITE);
		chart.setBorderPaint(Color.WHITE);
		chart.setTitle(title);
		cp.setPreferredSize(dim);
		cp.setBorder(null);
		collinearChart.add(cp);
		jframe.add(collinearChart,position);
	}
	
	private void show() {
		this.jframe.pack();
		this.jframe.setVisible(true);
	}
	
	public static void main(String[] args) {
		JfreeChart jf = new JfreeChart(new Dimension(1170,600));
		jf.plotLineChart("C:\\Users\\chenxi.zhou\\Desktop\\"
				+ "genetic mapping writing-up\\meta_final\\itr_misass\\"
				+ "15.txt"
				, new Dimension(630,250),Color.white,BorderLayout.WEST,
				"(a)");
		jf.plotCollinearChart("C:\\Users\\chenxi.zhou\\Desktop\\"
				+ "genetic mapping writing-up\\meta_final\\itr_misass\\"
				+ "Trifida.00004.fa.vs.Itr.sc000015.fa.chain.rdotplot"
				, new Dimension(270,250),Color.white,BorderLayout.CENTER,
				"(b)");
		jf.plotCollinearChart("C:\\Users\\chenxi.zhou\\Desktop\\"
				+ "genetic mapping writing-up\\meta_final\\itr_misass\\"
				+ "Trifida.40002.fa.vs.Itr.sc000015.fa.chain.rdotplot"
				, new Dimension(270,250),Color.white,BorderLayout.EAST,
				"(c)");
		jf.show();
		jf.print("C:\\Users\\chenxi.zhou\\Desktop\\test.pdf");
	}
	
	public static void main2(String[] args) {
		JfreeChart jf = new JfreeChart(new Dimension(600, 2000), 
				new GridLayout(6,1));
		jf.plotLineChart("C:\\Users\\chenxi.zhou\\Desktop\\"
				+ "genetic mapping writing-up\\meta_final\\itr_misass\\"
				+ "91.txt"
				, new Dimension(600,220),Color.white,BorderLayout.WEST,
				"(a) Itr_sc000091.1");
		jf.plotLineChart("C:\\Users\\chenxi.zhou\\Desktop\\"
				+ "genetic mapping writing-up\\meta_final\\itr_misass\\"
				+ "96.txt"
				, new Dimension(600,220),Color.white,BorderLayout.WEST,
				"(b) Itr_sc000096.1");
		jf.plotLineChart("C:\\Users\\chenxi.zhou\\Desktop\\"
				+ "genetic mapping writing-up\\meta_final\\itr_misass\\"
				+ "141.txt"
				, new Dimension(600,220),Color.white,BorderLayout.WEST,
				"(c) Itr_sc000141.1");
		jf.plotLineChart("C:\\Users\\chenxi.zhou\\Desktop\\"
				+ "genetic mapping writing-up\\meta_final\\itr_misass\\"
				+ "246.txt"
				, new Dimension(600,220),Color.white,BorderLayout.WEST,
				"(d) Itr_sc000246.1");
		jf.plotLineChart("C:\\Users\\chenxi.zhou\\Desktop\\"
				+ "genetic mapping writing-up\\meta_final\\itr_misass\\"
				+ "314.txt"
				, new Dimension(600,220),Color.white,BorderLayout.WEST,
				"(e) Itr_sc000314.1");
		jf.plotLineChart("C:\\Users\\chenxi.zhou\\Desktop\\"
				+ "genetic mapping writing-up\\meta_final\\itr_misass\\"
				+ "380.txt"
				, new Dimension(600,220),Color.white,BorderLayout.WEST,
				"(f) Itr_sc000380.1");
		jf.show();
		jf.print("C:\\Users\\chenxi.zhou\\Desktop\\test.pdf");
	}
	
	
	public void print(String plot_pdf) {
		try {
			float width = jframe.getSize().width,
					height = jframe.getSize().height;
			Document document = new Document(new Rectangle(width, height));
			PdfWriter writer = PdfWriter.getInstance(document, new FileOutputStream(plot_pdf));
			document.open();
			PdfContentByte canvas = writer.getDirectContent();
			PdfTemplate template = canvas.createTemplate(width, height);
			Graphics2D g2d = new PdfGraphics2D(template, width, height);
			jframe.paint(g2d);
			g2d.dispose();
			canvas.addTemplate(template, 0, 0);
			document.close();
		} catch (FileNotFoundException | DocumentException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}



