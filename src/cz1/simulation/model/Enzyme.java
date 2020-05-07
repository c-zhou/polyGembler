package cz1.simulation.model;

public class Enzyme {

    private String name;
    private String[] recognizationSequence;
    private int overhangLength;
    private String[][] overhang;

    public Enzyme(String name) {
        switch(name.toLowerCase()) {
            case "apeki" :
                this.name = "ApeKI";
                this.recognizationSequence = new String[]{"GCAGC","GCTGC"};
                this.overhangLength = 4;
                break;
            case "psti" :
                this.name = "PstI";
                this.recognizationSequence = new String[]{"CTGCAG"};
                this.overhangLength = 1;
                break;
            case "ecot22i" :
                this.name = "EcoT22I";
                this.recognizationSequence = new String[]{"ATGCAT"};
                this.overhangLength = 1;
                break;
            case "mspi" :
                this.name = "MspI";
                this.recognizationSequence = new String[]{"CCGG"};
                this.overhangLength = 3;
                break;
            case "tsei" :
            	this.name = "TseI";
            	this.recognizationSequence = new String[]{"GCAGC", "GCTGC"};
            	this.overhangLength = 4;
            	break;
            default:
                throw new RuntimeException("Unrecognized enzyme. Program exit.");
        }
        int n = recognizationSequence.length;
        int l = recognizationSequence[0].length();
        overhang = new String[n][2];
        String r;
        for(int i=0; i<n; i++) {
            if(overhangLength*2>l) {
                r = recognizationSequence[i].substring(l-overhangLength,overhangLength);
                overhang[i] = new String[]{"",r};
            } else {
                r = recognizationSequence[i].substring(overhangLength,l-overhangLength);
                overhang[i] = new String[]{r,""};
            }
        }
    }

    public String getName() {
        return name;
    }

    public String[] getRecognizationSequence() {
        return recognizationSequence;
    }

    public int getOverhangLength() {
        return overhangLength;
    }

    public String[][] getOverhang() {
        return overhang;
    }
}
