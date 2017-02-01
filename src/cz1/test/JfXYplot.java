package cz1.test;

import java.awt.Color;

import javax.swing.JFrame;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYDifferenceRenderer;
import org.jfree.data.function.Function2D;
import org.jfree.data.function.NormalDistributionFunction2D;
import org.jfree.data.general.DatasetUtilities;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

public class JfXYplot {
	double mean = 0.0, sd = 1.0;

	public JfXYplot(JFrame frame) {
		XYDataset dataset = initDataset();
		JFreeChart chart = ChartFactory.createXYLineChart(
				"Normal Distribution",
				"X",
				"PDF",
				dataset,
				PlotOrientation.VERTICAL,
				false,
				false,
				false
				);

		XYPlot plot = chart.getXYPlot();
		NumberAxis domain = (NumberAxis) plot.getDomainAxis();
		domain.setAutoRangeStickyZero(false); //Fixes the margin issue with 0
		domain.setTickUnit(new NumberTickUnit(sd)); //Spacing on X-axis should be standard deviation + mean   

		XYSeriesCollection dataset1 = (XYSeriesCollection) dataset;

		XYSeries area1 = new XYSeries("area1");
		area1.add(-4, 0);
		area1.add(-1, 0);
		area1.add(-1, 0.25);
		area1.add(1, 0.25);
		area1.add(1, 0);
		area1.add(4, 0);
		dataset1.addSeries(area1);
		XYDifferenceRenderer renderer1 = new XYDifferenceRenderer();
		Color none = new Color(0, 0, 0, 0);
		renderer1.setNegativePaint(none);//hide the area where the values for the second series are higher
		renderer1.setSeriesPaint(1, none);//hide the outline of the "difference area"

		plot.setRenderer(0, renderer1);

		final ChartPanel chartPanel = new ChartPanel(chart);
		frame.getContentPane().add(chartPanel);
	}

	private XYDataset initDataset() {
		double minX = mean - (4 * sd), maxX = mean + (4 * sd);  //Minimum and Maximum values on X-axis (4 deviations)
		Function2D normal = new NormalDistributionFunction2D(mean, sd);
		XYDataset dataset = DatasetUtilities.sampleFunction2D(normal, minX, maxX, 100, "Normal");
		return dataset;
	}

	public static void main(String[] args) {
		JFrame frame = new JFrame("Normal Distribution Demo");
		new JfXYplot(frame);
		frame.pack();
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.setVisible(true);
	}
}