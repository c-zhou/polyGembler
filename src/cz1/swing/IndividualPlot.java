package cz1.swing;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;

import javax.swing.JTabbedPane;


public class IndividualPlot extends JTabbedPane implements PropertyChangeListener{

	@Override
	public void propertyChange(PropertyChangeEvent arg0) {
		// TODO Auto-generated method stub
		
	}

	public static File makeOrDelete(File outdir, String st){
		File fi = new File(outdir, st);
		if(!fi.exists()) fi.mkdir();
		else{
			File[] f = fi.listFiles();
			for(int i=0; i<f.length; i++){
				f[i].delete();
			}
		}
		return fi;
	}
	//AbstractVectorGraphicsIO[] gB_; //= new AbstractVectorGraphicsIO[];
	// AbstractVectorGraphicsIO[] gR_;// = new HashMap<Integer, AbstractVectorGraphicsIO>();



}
