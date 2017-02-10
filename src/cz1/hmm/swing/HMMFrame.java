package cz1.hmm.swing;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.SwingUtilities;

import org.apache.commons.lang3.ArrayUtils;

import cz1.hmm.data.DataEntry;
import cz1.hmm.model.HiddenMarkovModel;
import cz1.util.Constants;

public class HMMFrame {
	public JFrame jframe;
	Headless jframeH;
	List<HMMPanel> hmm_tab = new ArrayList<HMMPanel>();
	JTabbedPane tabs = new JTabbedPane();
	List<String> tab_names = new ArrayList<String>();
	public String[] getPoss(){
		return tab_names.toArray(new String[0]);
	}
	public void setTab(String name){
		tabs.setSelectedIndex(tab_names.indexOf(name));
	}
	public String currentTab(){
		return tab_names.get(tabs.getSelectedIndex());
	}
	/*  public JButton jbutton_it = new JButton("it");
	public JButton jbutton_sw = new JButton("sw");
	public JButton jbutton_plot1 = new JButton("plot -1");
	public JButton jbutton_plot2 = new JButton("plot -2");
	public JButton jbutton_finish = new JButton("finish");*/
	JPanel buttons = new JPanel();
	public HMMFrame(){
		//super();
		if(Constants.plot()){
			jframe = new JFrame();
			/*    buttons.setLayout(new BoxLayout(buttons, BoxLayout.X_AXIS));
            buttons.add(jbutton_it);
            buttons.add(jbutton_sw);
            buttons.add(jbutton_plot1);
            buttons.add(jbutton_plot2);
            buttons.add(jbutton_finish);*/
			jframe.getContentPane().setLayout(new BorderLayout());
			jframe.add(BorderLayout.CENTER,tabs);
			/*     jframe.add(BorderLayout.BEFORE_FIRST_LINE, buttons);
            jbutton_it = new JButton("train");
            jbutton_sw = new JButton("expand");
            jbutton_plot1 = new JButton("plot -1");
            jbutton_plot2 = new JButton("plot -2");*/
			jframe.setSize(new Dimension(600,450));
			//jframe.setExtendedState(JFrame.MAXIMIZED_BOTH);
			jframe.setDefaultCloseOperation(jframe.DISPOSE_ON_CLOSE);
		}
		else{
			jframeH = new Headless(tabs);
		}
	}

	public HMMPanel addHMMTab(HiddenMarkovModel hmm, DataEntry de, File out){
		int[] d = de.getIntegerPosition();
		HMMPanel hmmp = new HMMPanel(
				jframe.getSize(),
				hmm, 
				Arrays.asList(ArrayUtils.toObject(d)), 
				de.getAlleleA(), 
				de.getAlleleB(), 
				0, 
				out);
		JScrollPane jscp = new JScrollPane(hmmp,  JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
		this.tabs.add("HMM",jscp);
		this.pack();
		this.hmm_tab.add(hmmp);
		//  if(Constants.plot()>=1) setVisible(true);
		return hmmp;
	}

	public void addTab(String name, JComponent ip){
		this.tabs.add(name,ip);
		tab_names.add(name);
		if(Constants.plot())
		{
			setVisible(true);
			this.pack();
		}
	}



	public void clearTabs() {
		tabs.removeAll();
		if(hmm_tab.size()>0){
			hmm_tab.get(0).removeAll();
			// hmm_tab.get(0).re
		}

		this.hmm_tab.clear();
	}
	public void pack() {
		if(jframe!=null) jframe.pack();
		else this.jframeH.pack();
	}
	public void setVisible(final boolean b){
		SwingUtilities.invokeLater(
				new Runnable()
				{
					public void run(){
						if(jframe!=null) jframe.setVisible(b);
						else jframeH.setVisible(b);
					}
				});

	}


}
