package cz1.gbs.core;

import java.util.ArrayList;
import java.util.Random;

public class Digest {

    private final long RANDOM_SEED = 112019890314L;
    private final Random random= new Random(RANDOM_SEED);
    
    private Enzyme enzyme;
    private double rate;
    private Digestion digestion;

    public Digest(Enzyme enzyme, String sequence, double rate) {
        this.enzyme = enzyme;
        this.rate = rate;
        this.digestion = digest(sequence);
    }

    public Digest(Enzyme enzyme, String sequence) {
        this(enzyme, sequence, 0.95);
        this.digestion = digest(sequence);
    }

	public Digestion digest(String sequence) {
        String rs[] = enzyme.getRecognizationSequence();
        int enzyme_length = rs[0].length();
        int overhang = enzyme.getOverhangLength();
        char[] buffer = new char[enzyme_length];
        int buffered = 0, r;
        char base;

        ArrayList<Integer> cut = new ArrayList<Integer>();
        ArrayList<Integer> recognization = new ArrayList<Integer>();
        cut.add(0); recognization.add(-1);
        for(int i=0; i<sequence.length(); i++) {
            base = sequence.charAt(i);
            buffer[buffered++] = Character.toUpperCase(base);
            if (buffered == enzyme_length) {
                r = recognize(String.valueOf(buffer),rs);
                if(r!=-1 && random.nextDouble()<rate) { 
                    cut.add(i-overhang+1);
                    recognization.add(r);            
                    buffered = 0;
                } else {
                    buffered--;
                    System.arraycopy(buffer, 1, buffer, 0, enzyme_length-1);
                }
            }
        }
        cut.add(sequence.length()); recognization.add(-1);
        return new Digestion(cut, recognization);
	}

    public void statistics (ArrayList<Integer> cut) {
        int[] stats = new int[15];
        int s;
        for(int i=1; i<cut.size(); i++) {
            s = cut.get(i)-cut.get(i-1);
            if(s<100) stats[0]++;
            else if(s<200) stats[1]++;
            else if(s<300) stats[2]++;
            else if(s<400) stats[3]++;
            else if(s<500) stats[4]++;
            else if(s<600) stats[5]++;
            else if(s<700) stats[6]++;
            else if(s<800) stats[7]++;
            else if(s<900) stats[8]++;
            else if(s<1000) stats[9]++;
            else if(s<10000) stats[10]++;
            else if(s<100000) stats[11]++;
            else if(s<1000000) stats[12]++;
            else if(s<10000000) stats[13]++;
            else stats[14]++;
        }

        System.out.println("fragments 0-100 nt: "+stats[0]);
        System.out.println("fragments 100-200 nt: "+stats[1]);
        System.out.println("fragments 200-300 nt: "+stats[2]);
        System.out.println("fragments 300-400 nt: "+stats[3]);
        System.out.println("fragments 400-500 nt: "+stats[4]);
        System.out.println("fragments 500-600 nt: "+stats[5]);
        System.out.println("fragments 600-700 nt: "+stats[6]);
        System.out.println("fragments 700-800 nt: "+stats[7]);
        System.out.println("fragments 800-900 nt: "+stats[8]);
        System.out.println("fragments 900-1kb nt: "+stats[9]);
        System.out.println("fragments 1kb-10kb nt: "+stats[10]);
        System.out.println("fragments 10kb-100kb nt: "+stats[11]);
        System.out.println("fragments 100kb-1Mb nt: "+stats[12]);
        System.out.println("fragments 1Mb-10Mb nt: "+stats[13]);
        System.out.println("fragments >10Mb nt: "+stats[14]);
    }

    private static int recognize(String hit, String[] rs) {
        for(int i=0; i<rs.length; i++) 
            if(hit.equals(rs[i])) return i;
        return -1;
    }

    public ArrayList<Integer> getCutsite() {
        return digestion.cut;
    }

    public ArrayList<Integer> getRecognization() {
        return digestion.recognization;
    }

    private class Digestion {
        ArrayList<Integer> cut;
        //to record which regconization sequence it is,
        //for enzyme that recognize only one sequence
        //all the elements are 0 (except the 1st and the last one, -1)
        ArrayList<Integer> recognization;

        public Digestion(ArrayList<Integer> cut, ArrayList<Integer> recognization) {
            this.cut = cut;
            this.recognization = recognization;
        }
    }
}
