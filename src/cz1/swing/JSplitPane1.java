package cz1.swing;

import java.awt.Graphics;

import javax.swing.JComponent;
import javax.swing.JSplitPane;

public class JSplitPane1 extends JSplitPane{
	
	boolean hasProportionalLocation=false;
	double proportionalLocation;
	boolean isPainted;
	public JSplitPane1(int split, JComponent splitPane, JComponent splitPane2) {
		super(split,splitPane,splitPane2);
	}

	public void setDividerLocation(double proportionalLocation) {
        if (!isPainted) {       
            hasProportionalLocation = true;
            this.proportionalLocation = proportionalLocation;
        }
        else
            super.setDividerLocation(proportionalLocation);
    }

    public void paint(Graphics g) {
        if (!isPainted) {       
            if (hasProportionalLocation)
                super.setDividerLocation(proportionalLocation);
            isPainted = true;
        }
        super.paint(g);
    } 
	
}
