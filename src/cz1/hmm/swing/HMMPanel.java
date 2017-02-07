package cz1.hmm.swing;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Stroke;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;
import java.util.logging.Logger;

import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JSplitPane;
import javax.swing.SwingUtilities;

import org.freehep.graphics2d.VectorGraphics;
import org.freehep.graphicsio.ImageConstants;
import org.freehep.graphicsio.ImageGraphics2D;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartRenderingInfo;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYAnnotation;
import org.jfree.chart.annotations.XYPolygonAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.AxisLocation;
import org.jfree.chart.axis.AxisSpace;
import org.jfree.chart.axis.LogAxis;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.PlotRenderingInfo;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.Range;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.Layer;
import org.jfree.ui.RectangleEdge;

import cz1.hmm.model.HiddenMarkovModel;
import cz1.hmm.model.HiddenMarkovModel.TP;
import cz1.util.Algebra;
import cz1.util.Constants;
import cz1.util.IO;

/*@Author Lachlan Coin*/
public class HMMPanel extends JPanel implements PropertyChangeListener{
	final List<Integer> location;
	final List<Character> major;
	final List<Character> minor;
	Color  background_color = Color.WHITE,
			line_color = Color.BLACK,
			font_color = Color.BLACK;
	// int x_start =0;
	//  int x_end;



	int mult = 25;
	//final File hmmF;

	/* public void write(){
        try{
        for(int i=0; i<this.getComponentCount(); i++){
            InnerPanel p = (InnerPanel) this.getComponent(i);
            this.writeToZipFile(p,hmmF,  p.getName());
        }
        }catch(Exception exc){
            exc.printStackTrace();
        }
    }
    public void writeToZipFile(Component charts, File dir, String id) throws Exception{
        File out = new File(dir, (id)+".png");
        ImageGraphics2D g = new ImageGraphics2D(out,charts, ImageConstants.PNG); 
        g.startExport();
        charts.print(g);
        g.endExport();
    }*/
	final double noIndiv;
	Font small_font,
	large_font,
	small_italic_font,
	large_italic_font;
	FontMetrics fm_small,
	fm_large,
	fm_small_italic,
	fm_large_italic;

	Logger logger = Logger.global.getAnonymousLogger();
	// int[] y_loc;
	//  int domain_thickness = 30;
	// double ratio;

	//   double[][] mat;
	//double[][] hittingProb;
	public void update(){
		if(plot==2|| Constants.plot()>=2){
			for(int i=0; i<this.mat.length; i++)
				this.hmm.ep(i).getCounts(this.mat[i]);;
				this.ratePanel.removeAll();
				this.ratePanel.add(this.getChart());
				SwingUtilities.invokeLater(new Runnable(){
					public void run(){
						updateUI();
					}
				});
		}

	}

	InnerPanel ip;
	LocPanel locP;
	public HMMPanel(Dimension dim, HiddenMarkovModel hmm, List<Integer> location, List<Character> major, List<Character> minor, int noIndiv, File outdir){

		{

			ca = ColorAdapter.get(hmm);

		}
		this.maxLineWidth = maxLineWidth * (4.0/(double)hmm.hs());
		this.ratePanel = new JPanel();

		//for(int k=0; k<hmm.hs()-1;k++){
		//	this.rates.addSeries(new XYSeries(""+k));
		//}
		this.location = location;
		this.major = major;
		this.minor = minor;
		this.noIndiv = noIndiv;
		this.chartBF = IndividualPlot.makeOrDelete(outdir, "hmmF");
		mat = new double[hmm.noSnps()][hmm.hs()];
		//hmm.getCounts(hmm.noSnps, mat);
		for(int i=0; i<mat.length; i++) 
			hmm.ep(i).getCounts(mat[i]);
		/* for(int i=0; i<mat.length; i++){
        for(int j=1; j<mat[i].length; j++){
            mat[i][j] = transform(mat[i][j]); 
        }
    }*/
		this.numFounders = hmm.hs()-1;
		this.noSnps = hmm.noSnps();
		this.hmm = hmm;
		this.ratePanel.add(this.getChart());
		this.dim = dim;
		this.setMinimumSize(dim);
		this.setPreferredSize(dim);
		this.setSize(dim);
		this.minloc = location.get(0);
		this.reg_length = location.get(noSnps-1)-minloc;
		// this.hittingProb = hmm.getHittingProb(hmm.noSnps);
		this.setLayout(new BorderLayout());
		
		double[] x_coor = new double[this.noSnps];
		for(int i=0; i<this.noSnps; i++) x_coor[i] = this.getX(i, false)+shift;
		this.locP = new LocPanel(location,
				dim.width, 
				50,
				new double[]{x_coor[0], this.dim.width-x_coor[this.noSnps-1]},
				x_coor,
				"hmm");

		/**
		JSplitPane sjp2 = new JSplitPane(JSplitPane.VERTICAL_SPLIT,locP, ratePanel);
		sjp2.setDividerLocation(0.2);
		JSplitPane sjp =   new JSplitPane1(JSplitPane.VERTICAL_SPLIT,  
				ip = new InnerPanel(), sjp2
				);
		sjp.setDividerLocation(0.6); 
		this.add(sjp);
		*/
		this.setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
		this.add(ip=new InnerPanel());
		this.add(locP);
		this.add(ratePanel);
		ratePanel.setBackground(Color.white);

		/* if(stToP.length==0){
        this.addTab("no phen", new InnerPanel());

    }
    else{
        for(int i=0; i<stToP.length; i++){
            this.addTab(stToP[i], new InnerPanel(stToP[i]));
        }

    }*/
		//  emStSp =Emiss.getSpaceForNoCopies(Constants.backgroundCount());
		// ((EmissionState)hmm.getState(1)).getEmissionStateSpace();
		//  List<Integer> cn = emStSp.copyNumber;
		/* for(int i=0; i<cn.size(); i++){
        int[] geno = emStSp.getGenoForCopyNo(cn.get(i));
        for(int k=0; k<geno.length; k++){
            ca1.getColor(geno[k]+"");
        }
    }*/
	}


	Dimension dim;
	HiddenMarkovModel hmm;

	// final int index=0;

	//  final double textOffset;

	public double getX(int i, boolean useLoc){
		if(useLoc) return offset_x +(((double)location.get(i)-minloc)/this.reg_length) *(this.getWidth()-(offset_x+offset_x)) ;
		else return offset_x +((((double)i)/(double)this.noSnps)) *(this.getWidth()-(offset_x+offset_x)) ;
	}

	public double getY(int i){
		return ((double)i) *height*3/4 +offset_y;
	}



	public double getStartX(){
		return 5.0;
	}

	public double getStartY(){
		return (double)  (this.getHeight())/2.0;
	}

	public static Font font16 = new Font("SansSerif", Font.PLAIN, Constants.hmmFontSize()+2);
	public static Font font16b = new Font("SansSerif", Font.BOLD, Constants.hmmFontSize()+2); 
	double wid;
	double wid1;
	double height;
	double shape;
	double max_shape;


	double  prev_loc =0;

	final int minloc;
	final double reg_length;
	double maxLineWidth = 35;
	static int len = 300;
	static int width = 1000;
	final int noSnps;
	final int numFounders;
	final double offset_x = mult;
	// final double offset_xR = mult;
	double offset_y = 50;

	static Color LIGHTGREEN = new Color(0,100,0);
	static Color DARKGREEN  = new Color(189, 183, 107);


	public  ColorAdapter ca;


	// static ColorAdapter ca1 = new ColorAdapter(new Color[] {Color.BLUE, Color.RED, Color.RED, Color.GREEN, Color.GREEN, Color.GREEN});

	final JPanel ratePanel;
	XYSeriesCollection rates;
	// XYSeriesCollection ratesSeries = new XYSeriesCollection();
	
	public Range updateRates(){
		rates = new XYSeriesCollection();
		rates.addSeries(new XYSeries("probs"));
		for(int k=0; k<rates.getSeriesCount(); k++){
			rates.getSeries(k).clear();
		}
		double upper = Double.NEGATIVE_INFINITY,
				lower = Double.POSITIVE_INFINITY;
		for(int i=1; i<hmm.noSnps(); i++){
			double[][] probs = hmm.tp(i).probs();
			double r = 0;
			for(int j=1; j<probs.length; j++)
				r += 1-probs[j][j];
			r /= (probs.length-1);
			if(r>upper) upper = r;
			if(r<lower) lower = r;
			this.rates.getSeries(0).add(i+.5, r);
		}
		lower = Math.floor(lower*10)/10-offset_yR;
		upper = Math.ceil(upper*10)/10+offset_yR;
		
		//shading objects
		//rates.addSeries(new XYSeries("shade1"));
		//rates.getSeries(1).add(5,0);
		//rates.getSeries(1).add(6,.3);
		XYSeries xy = rates.getSeries(0);
		for(int i=0; i<xy.getItemCount(); i++) {
			if(xy.getY(i).doubleValue()<shading_thresh)
				continue;
			double x = xy.getX(i).doubleValue();
			double x_start = Math.floor(x),
					x_end = Math.ceil(x),
					y_start = lower+offset_yR,
					y_end = upper-offset_yR;
			if(shading_slope==0) {
				double y_inter = y_start+shading_step;
				while( y_inter<y_end) {
					XYSeries tmp = new XYSeries(
							"shading "+(rates.getSeriesCount()+1));
					tmp.add(x_start, y_inter);
					tmp.add(x_end, y_inter);
					rates.addSeries(tmp);
					y_inter += shading_step;
				}
			} else if(Double.isInfinite(shading_slope)) {
				double x_inter = x_start+shading_step;
				while( x_inter<x_end) {
					XYSeries tmp = new XYSeries(
							"shading "+(rates.getSeriesCount()+1));
					tmp.add(x_inter, y_start);
					tmp.add(x_inter, y_end);
					rates.addSeries(tmp);
					x_inter += shading_step;
				}
			} else {
				double x_start2 = Double.NaN,
						x_end2 = Double.NaN;
				if(shading_slope>0) {
					x_start2 = x_start-y_end/shading_slope;
					x_end2 = x_end-y_start/shading_slope;
				} else {
					x_start2 = x_start-y_start/shading_slope;
					x_end2 = x_end-y_end/shading_slope;
				}
				double x_inter = x_start2+shading_step;
				while( x_inter<x_end2 ) {
					XYSeries tmp = new XYSeries(
							"shading "+(rates.getSeriesCount()+1));
					double y0 = shading_slope*(x_start-x_inter);
					double y1 = shading_slope*(x_end-x_inter);
					double x0 = y_start/shading_slope+x_inter;
					double x1 = y_end/shading_slope+x_inter;
					if(y0>=y_start && y0<=y_end) 
						tmp.add(x_start, y0);
					if(y1>=y_start && y1<=y_end) 
						tmp.add(x_end, y1);
					if(x0>=x_start && x0<=x_end) 
						tmp.add(x0, y_start);
					if(x1>=x_start && x1<=x_end) 
						tmp.add(x1, y_end);
					rates.addSeries(tmp);
					x_inter += shading_step;
				}
			}
		}
		return new Range(lower, upper);
	}

	private final double offset_yR = 0.03;
	private final double shading_step = 0.05;
	private final double shading_thresh = 0.1;
	private final double shading_slope = -0.3;
	
	private final int shift = 50;
	private final int shift0 = 75;
	static String formatStr= "%1.3e";
	static String formatStr2= "%1.4f";
	
	public ChartPanel getChart(){
		Range range =  this.updateRates();
		final JFreeChart chart = ChartFactory.createXYLineChart(
				null,
				null, // domain axis label
				null, // range axis label
				rates, // data
				PlotOrientation.VERTICAL, 
				false, // include legend
				true, // tooltips?
				false // URL generator? Not required...
				);
		
        final NumberAxis domainAxis = new NumberAxis("SNP relative physical position");
        //final LogarithmicAxis rangeAxis = new LogarithmicAxis("Trans-out probability");
        final NumberAxis rangeAxis = new NumberAxis("Trans-out probability");
        AxisSpace space = new AxisSpace();
        space.setRight(22); //reserved space on the left side of the plot
        space.setLeft(67);
        chart.getXYPlot().setFixedRangeAxisSpace(space);
        //double fixedDemension = HMMPanel.this.getX(0, false)+shift;
        //rangeAxis.setFixedDimension(fixedDemension);
        domainAxis.setTickUnit(new NumberTickUnit(1));
        //rangeAxis.setLog10TickLabelsFlag(true);
        rangeAxis.setRange(range);
        domainAxis.setRange(new Range(1-offset_yR, noSnps+offset_yR));
        chart.getXYPlot().setDomainAxis(domainAxis);
        chart.getXYPlot().setRangeAxis(rangeAxis);
		
		final XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		for(int i=0; i<rates.getSeriesCount(); i++){
			Color c;
			if(i==0) c=ca.getColor(i);
			else c=Color.lightGray;
			renderer.setSeriesPaint(i, c);
			if(i==0)
				renderer.setSeriesShapesVisible(i, true);
			else
				renderer.setSeriesShapesVisible(i, false);
		}
		/**
		XYPolygonAnnotation a = new XYPolygonAnnotation(new double[] {5.0,
				   0.0, 5.0, 0.5, 6.0, 0.5, 6.0, 0.0}, null, null,Color.black) {
            @Override
            public void draw(Graphics2D g2, XYPlot plot, Rectangle2D dataArea, ValueAxis domainAxis, ValueAxis rangeAxis, int rendererIndex, PlotRenderingInfo info) {
                Graphics2D g22 = (Graphics2D) g2.create();
                g22.setXORMode(new Color(0xff0000));
                super.draw(g22, plot, dataArea, domainAxis, rangeAxis, rendererIndex, info); 
            }
        };
		renderer.addAnnotation(a, Layer.BACKGROUND);
		*/
		chart.getXYPlot().setRenderer(renderer);
		double anno_offset_x = 0.5;
		
		XYSeries d = rates.getSeries(0);
		for(int i=0; i<d.getItemCount(); i++){
			double x = d.getX(i).doubleValue();
			double y = d.getY(i).doubleValue();
			if(y>0.1)
				chart.getXYPlot().addAnnotation(
					(XYAnnotation) new XYTextAnnotation(String.format(formatStr2,y),
							x>noSnps/2?(x-anno_offset_x):(x+anno_offset_x), y));
			else
				chart.getXYPlot().addAnnotation(
						(XYAnnotation) new XYTextAnnotation(String.format(formatStr,y),
								x, y+.05));
		}
		
		
		//  chart.getXYPlot().getRangeAxis().get
		//chart.getXYPlot().setRangeAxis(new LogarithmicAxis("rates"));
		final ChartPanel cp =   new ChartPanel(
				chart,
				this.getWidth(), //width
				150, //height
				this.getWidth(), //mindrawWidth
				150, //mindrawHeight
				this.getWidth(), //maxDrawWith
				150,//maxDrawHeight
				ChartPanel.DEFAULT_BUFFER_USED,
				true,  // properties
				true,  // save
				true,  // print
				true,  // zoom
				true   // tooltips		

				);
		{
			Range r =  chart.getXYPlot().getDomainAxis().getRange();
			chart.getXYPlot().getDomainAxis().setAutoRange(false);
			chart.getXYPlot().getDomainAxis().setRange(1,r.getUpperBound());
		}
		{
			Range r =  chart.getXYPlot().getRangeAxis().getRange();
			if(r.getLowerBound() < range.getUpperBound()){
				chart.getXYPlot().getRangeAxis().setAutoRange(false);
				chart.getXYPlot().getRangeAxis().setRange(r.getLowerBound(), range.getUpperBound());
			}
		}
		chart.setBackgroundPaint(Color.WHITE);
		chart.setBorderPaint(Color.WHITE);
		cp.setBorder(null);
		
		/**
		RectangleEdge edge = chart.getXYPlot().getDomainAxisEdge();
		Rectangle2D edge2 = cp.getScreenDataArea();
		double x = RectangleEdge.coordinate(edge2, edge);
		ValueAxis edge3 = chart.getXYPlot().getDomainAxisForDataset(0);
		double xx = domainAxis.valueToJava2D(3, edge2, edge);
		**/
	    return cp;
	}

	class InnerPanel extends JPanel{
		VectorGraphics vg;
		final int phenIndexToPaint;// =0;
		final double[] angle;
		final double[] angle1;

		public InnerPanel(){
			this.phenIndexToPaint = -1;
			this.angle = null;
			this.angle1 = null;
			this.setPreferredSize(new Dimension(dim.width,250));
			this.setMinimumSize(this.getPreferredSize());
			this.setSize(this.getPreferredSize());
		}

		public void paint0(){
			/*   vg.setColor(Color.BLACK);
      drawOval(getStartX(), getStartY(),max_shape, 0, 360);
       if( hmm!=null){
          // double[] trans  = new double[numFounders];
           for(int k=0; k<numFounders; k++){
            //   double shape_k = max_shape* mat[1][k+1];
               double val  = transform2(hmm.getTransitionScore(0, k+1, 0));
               vg.setColor(Color.black);//getLineColor(val,  Color.BLACK));
               vg.setLineWidth(maxLineWidth *val);

//               vg.drawLine(getStartX(),getStartY(), HMMPanel.this.getX(0, false),  HMMPanel.this.getY(k));
           }
       }*/
		}
		void drawOval(double xcen, double ycen, double shape, double startAngle, double endAngle){
			vg.drawArc(xcen - shape/2.0, ycen - shape/2.0, shape, shape, startAngle, endAngle);
		}
		void fillOval(double xcen, double ycen, double shape, double startAngle, double endAngle){
			vg.fillArc(xcen - shape/2.0, ycen - shape/2.0, shape, shape, startAngle, endAngle);
		}
		public void write(int i){
			double x_start = HMMPanel.this.getX(i, false)+shift;
			// double textOffset =  ((((double)0.5)/(double)noSnps)) *(this.getWidth()-(offset_x+offset_x)) ;
			vg.setFont(font16);
			for(int j=0; j< numFounders;j++){
				double y_start = HMMPanel.this.getY(j);
				vg.setColor(Color.BLACK);
				//vg.setColor(Color.WHITE);
				double[] probs= hmm.ep(i).probs()[j+1];
				double rem =1.0;
				for(int k=0; k<probs.length; k++){

					{
						double p =  Math.round((k==probs.length-1 ? rem :probs[k])*100.0)/100.0;
						rem -=p;
						String s = p+"";
						if(p>0){
							vg.drawString(
									//emStSp.getGenotypeString(emStSp.get(k))+
									s,
									x_start-(s.length()-.5)*Constants.hmmFontSize/4, 
									y_start+k*Constants.hmmFontSize+(k==0?-2.5:2.5));
						}
						// sb.append(emStSp.get(k)+":"+probs[k]+"\n");
					}
				}
			}
		}

		public void write0(){
			double x_start = shift0/2;
			// double textOffset =  ((((double)0.5)/(double)noSnps)) *(this.getWidth()-(offset_x+offset_x)) ;
			vg.setFont(font16b);
			for(int j=0; j< numFounders; j++){
				double y_start = HMMPanel.this.getY(j);
				vg.setColor(Color.BLACK);
				
				//String[] genoList = hmm.ob(0).allele();
				String[] genoList = new String[]{"A","B"};
				for(int k=0; k<genoList.length; k++){
					vg.drawString(
							genoList[k]
									,
									x_start, y_start+k*Constants.hmmFontSize+(k==0?-2.5:2.5));
				}

			}
			vg.setFont(font16);
		}

		public void paint( int i){
			double x_start = HMMPanel.this.getX(i, false)+shift;
			vg.setColor(Color.BLACK);
			//  vg.setLineWidth(0.5);
			double x_loc = HMMPanel.this.getX(i, true);
			//  vg.drawLine(x_loc, y_0, x_loc, y_1);
			// vg.drawLine(x_loc, y_1, x_start, HMMPanel.this.getY(0));
			/*  if(i==0 || x_loc> 100+prev_loc){
           String str = HMMPanel.this.location.get(i)+"";
           vg.setColor(Color.BLACK);
           vg.drawString(str, x_loc, 10);
           prev_loc = x_loc;

       }*/
			double x_start1 = i<noSnps-1 ? HMMPanel.this.getX(i+1,false) : HMMPanel.this.getX(i,false);
			x_start1 += shift;
			for(int j=0; j< numFounders; j++){
				if(Constants.CHECK){
					double sum = Algebra.sum(mat[i]);
					try{
						if(Math.abs(1.0-sum)>0.001) throw new RuntimeException("!! "+sum);
					}catch(Exception exc){
						exc.printStackTrace();
						System.exit(0);
					}
				}


				double y_start = HMMPanel.this.getY(j);
				Color c = ca.getColor(j+"");
				//   vg.setColor(c);
				double shape_j =max_shape* transform2( mat[i][j+1]);
				//  vg.setLineWidth(line_width);
				//  drawOval(x_start, y_start, shape_j);

				int max;
				Color c1 ;

				double[] probs =  hmm.ep(i).probs()[j+1];
				max = Algebra.maxIndex(probs);
				int[] comp = new int[]{0, 1};
				double frac = mod(probs, max, comp);
				//System.err.println(emStSp.get(max)+" "+comp.length+" "+frac);
				c1 = HMMPanel.this.modify(c, frac);
				vg.setColor(c1);
				//vg.setColor(Color.red);

				fillOval(x_start, y_start, shape_j, 0, 360);
				
				// vg.fillOval(x_start, y_start, shape_j,shape_j);
				// vg.setColor(getColor());

				// vg.setLineWidth(HMMPanel.this.maxLineWidth);
				//  double[][] mat1 = mat;
				double mult  = Constants.plotFlux() ? mat[i][j+1] : 0.2;
				if(i<noSnps-1 && hmm!=null){
					for(int k=0; k<numFounders; k++){
						//  hmm.getHittingProb(hmm.noSnps);
						double pr = hmm.tp(i+1).probs()[j+1][k+1];
						//double val = transform2(pr*(mult));
						double val = pr;
						vg.setColor(c);//getLineColor(val, c));
						//vg.setColor(Color.red);
						double wid = HMMPanel.this.maxLineWidth*val;

						//  vg.setStroke(getStroke(Math.pow(pr,1.0)));

						if(wid>0) {
							vg.setLineWidth(wid);
							vg.drawLine(x_start,HMMPanel.this.getY(j), x_start1,  HMMPanel.this.getY(k));
						}
					}
				}


			}
		}
		float[] dash = new float[2];
		private Stroke getStroke(double d) {
			dash[0] = Math.max(0,(float) (d*10.0));
			dash[1] = Math.max(0,(float) ((1-d)*10.0));
			return new BasicStroke(1.0f, BasicStroke.CAP_SQUARE, BasicStroke.JOIN_MITER, 10.0f, dash, 0.0f);
		}
		public void paint( Graphics g ) {
			//  double x_0 = HMMPanel.this.getX(0, true);
			//   double x_1 = HMMPanel.this.getX(noSnps-1, true);
			//offset_y = ((double)getHeight())*0.1;

			g.setColor(HMMPanel.this.background_color);
			g.fillRect(0,0, getWidth(), getHeight());
			
			vg = VectorGraphics.create(g);
			//    vg.drawRect(x_0, y_0, x_1 - x_0, y_1-y_0);//(arg0, arg1, arg2, arg3)
			wid =  ((double)getWidth() - (2*offset_x) ) / ((double) reg_length);
			wid1 =  ((double)getWidth() - (2*offset_x) ) / ((double) noSnps);
			height  = ((double)getHeight()-offset_y ) / ((double) numFounders-1.0);
			shape =Math.min(wid1, height);
			max_shape = Math.min(shape*5, Math.max(wid1, height));

			vg.setFont(font16);
			paint0();
			for(int i=0; i<noSnps; i++){
				paint(i);
				write(i);
			}
			write0();
			for(int i=0; i<noSnps; i++){
				write(i);
			} 


		}
	}
	static final double log2 = Math.log(2);
	double[][] mat;

	public double mod( double[] probs, int max, int[] comp ){
		double sum = 0;
		for(int i=0; i<comp.length; i++){
			sum+= -probs[comp[i]]*(Math.log(probs[comp[i]])/log2);
		}
		double info = -probs[max]*(Math.log(probs[max])/log2);
		return  transform1((sum-info)/sum);
	}
	public static Color modify(Color c, double frac){

		try{
			return new Color(c.getRed(), c.getGreen(), c.getBlue(), (int) Math.floor(c.getAlpha()*(frac)));
		}catch(Exception exc){
			exc.printStackTrace();
			System.err.println(frac);
			return Color.BLACK;
			// System.exit(0);
		}
		//  return null;
	}


	// final EmissionStateSpace emStSp;


	Character uncertain = '-';



	//   double y_0 = 10;
	//  double y_1 = 20;



	private double transform(double d) {
		//System.err.println(" t "+d);
		if(Constants.showAll(1))return d;
		if(Constants.plotFlux && d<1e-4 || ! Constants.plotFlux && d <0.05) return 0;
		else return d;
		//       return Math.pow(d,2);//Math.pow(10, Math.log(d)/Math.log(2));
	}

	private double transform2(double d) {
		//System.err.println(" t "+d);
		// if(d<0.0001) return 0;
		// else return d;
		return Math.pow(d, Constants.hmmBubblePow());//Math.pow(d,2);//Math.pow(10, Math.log(d)/Math.log(2));
	}
	private double transform1(double d) {
		return Math.pow(10, Math.log10(d)/Math.log(1.5));
	}
	private Color getLineColor(double d, Color c) {
		try{
			return new Color(c.getRed(),  c.getGreen(), c.getBlue(),(int)Math.round(d*255));
		}catch(Exception exc){
			System.err.println(d);
			exc.printStackTrace();
			System.exit(0);
			//           return null;
		}
		return null;
	}


	int plot = 0;
	public void setToPlot(int i){
		plot= i;
	}
	final File chartBF;


	public synchronized void propertyChange(PropertyChangeEvent arg0) {
		String nme = arg0.getPropertyName();
		/*if(nme.equals("emiss")){
        if(plot<2 && Constants.plot()<2) return;

	   Object[] obj = (Object[]) arg0.getNewValue();

	   StateDistribution dist = (StateDistribution) obj[0];
	   Integer l = (Integer) obj[1];
	   if(l==0) clear();
	   Integer i = (Integer)obj[2];
	   this.addInformation(dist,i);
   }*/    
		if(nme.equals("setToPlot")){
			int level = (Integer) arg0.getNewValue();
			setToPlot(level);
			return;
		}
		if(nme.equals("done") ){
			// if(plot<2 && Constants.plot()<2) return;
			//this.update();
			try{

				ip.setMinimumSize(dim);


				ip.setSize(dim);
				ratePanel.setSize(dim);
				int lp_height = locP.getSize().height;
				if(Constants.printPlots() && plot==2){
					{

						System.err.println("printing hmm");
						// int maxSNPS =Constants.maxSNPS();
						//  int n = (int) Math.ceil((double) this.noSnps/ maxSNPS);
						//  for(int i=0; i<n; i++){
						//  for(int j=0; j<this.getTabCount(); j++){
						Properties p = new Properties();
						p.setProperty("PageSize", "A4");
						boolean inclLocP = false;
						//int len1 = maxSNPS;//Math.min(maxSNPS, location.size()  - i*maxSNPS);
						Dimension dim1 = new Dimension(dim.width,dim.height+this.ratePanel.getHeight()+(inclLocP? lp_height : 0));	
						AffineTransform at = new AffineTransform();
						ImageGraphics2D g = new ImageGraphics2D(new File(this.chartBF, "hmm.png")
						,dim1, ImageConstants.PNG) ;
						g.setProperties(p);
						g.startExport();
						at.setToTranslation(0, 0);
						g.setTransform(at);
						this.ip.print(g);
						at.setToTranslation(0,dim.height );
						g.setTransform(at);
						this.ratePanel.print(g);
						if(inclLocP){ 
							at.setToTranslation(0,dim.height +ratePanel.getHeight());
							g.setTransform(at);

							this.locP.print(g);
						}
						g.endExport();
						System.err.println("done printing hmm");
						// }
						//  }
					}
					if(false)	{
						Properties p = new Properties();
						p.setProperty("PageSize", "A4");
						ChartPanel cp = (ChartPanel) this.ratePanel.getComponent(0);
						cp.setSize(dim);
						ImageGraphics2D g = new ImageGraphics2D(new File(this.chartBF, "rates.png"), cp,
								ImageConstants.PNG);

						g.setProperties(p);
						g.startExport();
						cp.print(g);
						g.endExport();

					}

				}
			}catch(Exception exc){
				exc.printStackTrace();
			}
		}


		/* if(nme.equals("init")){

    }
    else if(nme.equals("emiss")){

    }

    else 
    else if(nme.equals("expec_i")){

    }
    else if(nme.equals("expectation1")){

      //  this.updateUI();
    }
    else if(nme.equals("dist_maximisation")){

    }
    else*/ 
		if(nme.equals("hmm_maximisation")){
			update();
		}
		// else throw new RuntimeException("!! "+nme);

	}

	private void clear() {
		for(int i=0; i<this.mat.length; i++){
			Arrays.fill(mat[i], 0.0);
		}

	}


}