/**
 * 
 */
package cz1.swing;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.gicentre.utils.colour.CIELuv;

import cz1.model.HiddenMarkovModel;
import cz1.util.Constants;

public class ColorAdapter{


	final Color[] colors;

	public ColorAdapter(){
		this.transf = new HashMap<Integer, Integer>();
		char[] c = Constants.modify(0);
		if(c!=null){
			int cnt =2;
			for(int i=0; i<c.length; i++){

				if(c[i]=='0' || c[i]=='-' || c[i] =='N' || c[i]=='_'){
					transf.put(i,0);
				}
				else if(c[i]=='2' || c[i]=='3' || c[i]=='X' || c[i] =='Y' || c[i] =='Z') {
					transf.put(i,1);
				}
				else{
					transf.put(i, cnt);
					cnt++;
				}
			}
		}
		this.colors = this.defaultC.toArray(new Color[0]);
		// mod= Constants.modify0;
	}
	static int[] cols= new int[5];
	static List<Color> defaultC = new ArrayList<Color>();
	{
		for(int i=0; i<cols.length; i++){
			cols[i] =(int) Math.floor((double)i*(255.0/((double) cols.length)));
		}
		for(int i=cols.length-1;i>=0; i--){
			for(int j=0; j<cols.length; j++){
				for(int k=cols.length-1; k>=0; k--){
					if(cols[i]>230 && cols[j]>230 & cols[k]>230) continue;
					defaultC.add(new Color(cols[i], cols[j], cols[k]));
				}
			}
		}
		try{
			defaultC = randomise(defaultC);
		}catch(Exception exc){}
	}
	public ColorAdapter(Color[] array) {
		this.colors = array;
		//  mod= Constants.modify0;
	}

	public ColorAdapter(HiddenMarkovModel hmm) {
		switch(Constants.colorPalette) {
		case DEFAULT:
			colors = getDefaultColors(hmm);
			break;
		case GGPLOT2:
			colors = getGGPlot2Colors(hmm);
			break;
		default:
			throw new RuntimeException("!!!");
		}
	}
	
	public static Color[] colo = new Color[]{Color.red, Color.gray, Color.green, Color.blue};
	static{
		if(colo.length<Constants.maxCopies()+2){
			Color[] colo1 = new Color[Constants.maxCopies()+2];
			Arrays.fill(colo1, Color.blue);
			System.arraycopy(colo, 0,colo1, 0, colo.length);
			colo = colo1;
		}
	}
	private static Color[] getDefaultColors(HiddenMarkovModel hmm) {
		Color[] colors = new Color[hmm.hs()];
		colors[0] = Color.GRAY;
		//boolean firstNormal = true;
		List<Integer>[] normal = new List[colo.length];
		for(int i=0; i<colo.length; i++){
			normal[i] = new ArrayList<Integer>();
		}
		// int[] count = new int[colo.length];
		for(int i=1; i<colors.length; i++){
			Integer nocop = 1;
			if(nocop==null) nocop=1;
			colors[i] = colo[nocop];
			normal[nocop].add(i);
		}
		for(int k=0; k<normal.length; k++){
			if(normal[k].size()>1 ){
				if(k==1){
					float[] d = partition(normal[k].size());
					for(int i=0; i<normal[k].size(); i++){
						colors[normal[k].get(i)] =  Color.getHSBColor(d[i], (float)1.0, (float) Constants.state_brightness);
					}
				}
				else if(k==0){
					for(int i=0; i<normal[k].size(); i++){
						colors[normal[k].get(i)] =  Color.red;
					}

				}
				else if(k==2){
					for(int i=0; i<normal[k].size(); i++){
						colors[normal[k].get(i)] =  Color.GREEN;
					}

				}

			}

		}
		return colors;

	}

	private Color[] getGGPlot2Colors(HiddenMarkovModel hmm) {
		// TODO Auto-generated method stub
		Color[] colors = new Color[hmm.hs()];
		colors[0] = Color.GRAY;
		String[] hex_str = null; 
		switch(hmm.hs()) {
		case 5:
			hex_str = new String[]{"#F8766D","#7CAE00","#00BFC4","#C77CFF"};
			break;
		case 9:
			hex_str = new String[]{"#F8766D","#CD9600","#7CAE00","#00BE67",
									"#00BFC4","#00A9FF","#C77CFF","#FF61CC"};
			break;
		case 13:
			hex_str = new String[]{"#F8766D","#DE8C00","#B79F00","#7CAE00",
									"#00BA38","#00C08B","#00BFC4","#00B4F0",
									"#619CFF","#C77CFF","#F564E3","#FF64B0"};
			break;
		default:
			throw new RuntimeException("!!!");
		}
		for(int i=1; i<colors.length; i++)
			colors[i] = Color.decode(hex_str[i-1]);
		return colors;
	}
	
	private static float[] partition(int n) {
		float[] res = new float[n];
		float jump = (float)(1.0/(double) (n));
		for(int i=0; i<res.length; i++){
			res[i] = (float) ((double)i)*jump;
		}
		return res;
	}
	//char[] mod;
	static Color[] col = new Color[] {Color.RED, Color.GREEN, Color.CYAN, Color.BLUE, Color.PINK, Color.YELLOW, Color.MAGENTA, Color.ORANGE};
	private List<Color> randomise(List<Color> defaultC2) {
		List<Color> ne = new ArrayList<Color>();
		for(int i=0; i<col.length; i++){
			ne.add(col[i]);
		}
		while(defaultC2.size()>0){
			int i = Constants.nextInt(defaultC2.size());
			ne.add(defaultC2.get(i));
			defaultC2.remove(i);
		}
		return ne;
	}
	Map<String, Color> m = new HashMap<String, Color>();

	Map<Integer, Integer> transf  = new HashMap<Integer, Integer>();
	public Color getColor(int sh1){
		if(sh1<0) return Color.black;
		/*  Integer sh = transf.get(sh1);
     if(sh==null){
    	 Color c = getColor(sh1+"");
    	 transf.put(sh1,sh);
     }*/
		return getColor(sh1+"");
		/*     if(mod[sh]=='0') return getColor("0");
        else if (mod[sh]=='2') return getColor("1");
        else return getColor(sh+"");*/
	}

	public Color getColor(String st){
		//int ind = st1.indexOf('_');
		//	String st = ind>0 ? st1.substring(0,ind) : st1;
		try{
			String[] str = st.split(":");
			Color[] color = new Color[str.length];
			if(str.length==1){
				color[0] =m.get(str[0]);
				if(color[0]==null ){
					if(m.size()>=colors.length){
						m.put(str[0],color[0] =  getRandomColor());
					}
					else{
						m.put(str[0], color[0] =colors[m.size()] );
					}
				}

			}
			else{
				for(int i=0; i<str.length; i++){
					color[i] = getColor(str[i]);
				}
			}
			return merge(color);
		}catch(Exception exc){
			System.err.println("warning - problem with "+exc.getMessage());
			return Color.BLACK;
		}
	}
	public static int rand(){
		return (int) Math.floor(Math.random()*255.0);
	}
	public static Color getRandomColor() {
		// TODO Auto-generated method stub
		return Color.getHSBColor((float)Math.random(), (float)Math.random(), (float) Math.random());
		//return new Color(rand(), rand(),rand());
	}



	/* public static Color getRandomColor1() {
		// TODO Auto-generated method stub
    	int[] d = new int[] {rand(), rand(),rand()};
    	int sum = Constants.su
		return new Color(d[0], d[1], d[2]);
	}*/

	public static Color merge(Color[] color2) {
		/*if(color2.length==2){

    		if(	(color2[0].equals(Color.RED) && color2[1].equals(Color.GREEN) ) || 
    			(color2[1].equals(Color.RED) && color2[0].equals(Color.GREEN) )){
    		return Color.orange;
    	}
    	}*/
		int[] rgb = new int[3];
		Arrays.fill(rgb,0);
		double sum=0;
		// float sum=(float) 0.0;
		//    Arrays.fill(hsb, 0);
		for(int i=0; i<color2.length; i++){

			if(color2[i].equals(Color.black)) sum += 0.4;
			else{
				sum+=1;
				rgb[0]+=color2[i].getRed();
				rgb[1]+=color2[i].getGreen();
				rgb[2]+=color2[i].getBlue();
			}

		}
		for(int i=0; i<rgb.length; i++){
			rgb[i] = (int) Math.floor((double)rgb[i]/(double) sum);//Math.min(255, rgb[i]);
		}
		Color c =  new Color(rgb[0], rgb[1], rgb[2]);
		if(c.equals(Color.yellow)) return Color.orange;
		return c;
	}

	Color[] color2;

	public Color get(double[] d ){
		if(color2==null) 
			color2 = new Color[d.length];
		for(int i=0; i<color2.length; i++){
			color2[i] = this.getColor(i+"");
		}
		return merge1(color2, d);
	}
	public static Color merge1(Color[] color2, double[] d) {
		int[] rgb = new int[3];
		Arrays.fill(rgb,0);
		double sum= 0.0;
		//    Arrays.fill(hsb, 0);
		for(int i=0; i<color2.length; i++){
			rgb[0]+=color2[i].getRed()*d[i];
			rgb[1]+=color2[i].getGreen()*d[i];
			rgb[2]+=color2[i].getBlue()*d[i];
			sum+=d[i];
		}
		for(int i=0; i<rgb.length; i++){
			rgb[i] = (int)Math.round((double)rgb[i]/sum);
		}
		return new Color(rgb[0], rgb[1], rgb[2]);
	}

	static  int no_states;
	static  ColorAdapter ca;
	public static ColorAdapter get(HiddenMarkovModel hmm) {
		// TODO Auto-generated method stub
		if(ca==null || no_states!=hmm.hs()) {
			ca = new ColorAdapter(hmm);
			Color c1 = ca.getColor("start");
			for(int i=0; i<hmm.hs(); i++){
				Color c2 =  ca.getColor(i);
			}
			// Color c3 = ca.getColor(1);
			//Color c4  = ca.getColor(2);
			// Color c5  = ca.getColor(3);
			no_states = hmm.hs();

		}
		return ca;
	}

	public Color getColor3(short s) {
		return this.colors[s];
	}


}