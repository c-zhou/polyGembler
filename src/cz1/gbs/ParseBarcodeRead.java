package cz1.gbs;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import org.apache.commons.lang3.StringUtils;

/**
 * Takes a key file and then sets up the methods to decode a read from the sequencer.
 * The key file decribes how barcodes are related to their taxon.  Generally, a keyfile
 * with all flowcells is included, and then the flowcell and lane to be processed are
 * indicated in the constructor.
 *
 * @author Ed Buckler, Jeff Glaubitz, and James Harriman
 *
 */
public class ParseBarcodeRead {

	private static int chunkSize = BaseEncoder.chunkSize;
    protected int maximumMismatchInBarcodeAndOverhang = 0;
    protected static String[] initialCutSiteRemnant = null;
    protected static int readEndCutSiteRemnantLength;
    static String nullS = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    protected static String[] likelyReadEnd = null;
    protected static String theEnzyme = null;
    static int maxBarcodeLength = 10;
    private Barcode[] theBarcodes;
    private long[] quickBarcodeList;
    private HashMap<Long, Integer> quickMap;
    private static Field field;
	static {
		try {
			field = String.class.getDeclaredField("value");
		} catch (NoSuchFieldException | SecurityException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		field.setAccessible(true);
	}
	
    public int getMaximumMismatchInBarcodeAndOverhang(){
        return maximumMismatchInBarcodeAndOverhang;
    }

    public String[] getInitialCutSiteRemnant(){
        return initialCutSiteRemnant;
    }

    public int getReadEndCutSiteRemnantLength(){
        return readEndCutSiteRemnantLength;
    }

    public String[] getLikelyReadEnd(){
        return likelyReadEnd;
    }

    public String getTheEnzyme(){
        return theEnzyme;
    }

    public int getMaxBarcodeLength(){
        return maxBarcodeLength;
    }

    public Barcode[] getTheBarcodes(){
        return theBarcodes;
    }

    public long[] getQuickBarcodeList(){
        return quickBarcodeList;
    }

    public HashMap<Long, Integer> getQuickMap(){
        return quickMap;
    }

    /**
     * Create the barcode parsing object
     * @param keyFile file location for the keyfile
     * @param enzyme name of the enzyme
     * @param flowcell name of the flowcell to be processed
     * @param lane name of the lane to be processed
     */
    public ParseBarcodeRead(String keyFile, String enzyme, String flowcell, String lane) {
        if (enzyme != null) {
            chooseEnzyme(enzyme);
        } else {
            chooseEnzyme(getKeyFileEnzyme(keyFile));
        }

        int totalBarcodes = setupBarcodeFiles(new File(keyFile), flowcell, lane);
        System.out.println("Total barcodes found in lane:" + totalBarcodes);
    }

    public void setMaxDivergence(int maximumMismatch){
        this.maximumMismatchInBarcodeAndOverhang = maximumMismatch;
    }

    public int getMaxDivergence(){
        return this.maximumMismatchInBarcodeAndOverhang;
    }

    /**
     * Determines which cut sites to look for, and sets them, based on the
     * enzyme used to generate the GBS library. For two-enzyme GBS both enzymes
     * MUST be specified and separated by a dash "-". e.g. PstI-MspI, SbfI-MspI
     *
     * @param enzyme The name of the enzyme (case insensitive)
     */
    //TODO these should all be private static final globals, then just use this set which one is active.
    public static void chooseEnzyme(String enzyme) {
        // Check for case-insensitive (?i) match to a known enzyme
        // The common adapter is: [readEndCutSiteRemnant]AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG  
        if (enzyme.matches("(?i)apek[i1]")) {
            theEnzyme = "ApeKI";
            initialCutSiteRemnant = new String[]{"CAGC", "CTGC"};
            likelyReadEnd = new String[]{"GCAGC", "GCTGC", "GCAGAGAT", "GCTGAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 4;
        } else if (enzyme.matches("(?i)pst[i1]")) {
            theEnzyme = "PstI";
            initialCutSiteRemnant = new String[]{"TGCAG"};
            //likelyReadEnd = new String[]{"CTGCAG", "CTGCAAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            likelyReadEnd = new String[]{"CTGCAAGATCGGAAGA"};
            readEndCutSiteRemnantLength = 5;
        } else if (enzyme.matches("CviAII")) {
        	theEnzyme = "CviAII";
        	initialCutSiteRemnant = new String[]{"ATG"};
        	likelyReadEnd = new String[]{"CTGCAAGATCGGAAGA"};
        	readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("TseI")) {
        	theEnzyme = "TseI";
        	initialCutSiteRemnant = new String[]{"CAGC","CTGC"};
        	likelyReadEnd = new String[]{"CTGCAAGATCGGAAGA"};
        	readEndCutSiteRemnantLength = 4;
        } else if (enzyme.matches("(?i)ecot22[i1]")) {
            theEnzyme = "EcoT22I";
            initialCutSiteRemnant = new String[]{"TGCAT"};
            likelyReadEnd = new String[]{"ATGCAT", "ATGCAAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        } else if (enzyme.matches("(?i)pas[i1]")) {
            theEnzyme = "PasI";
            initialCutSiteRemnant = new String[]{"CAGGG", "CTGGG"};
            likelyReadEnd = new String[]{"CCCAGGG", "CCCTGGG", "CCCTGAGAT", "CCCAGAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        } else if (enzyme.matches("(?i)hpaii|(?i)hpa2")) {
            theEnzyme = "HpaII";
            initialCutSiteRemnant = new String[]{"CGG"};
            likelyReadEnd = new String[]{"CCGG", "CCGAGATCGG"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)msp[i1]")) {
            theEnzyme = "MspI";  // MspI and HpaII are isoschizomers (same recognition seq and overhang)
            initialCutSiteRemnant = new String[]{"CGG"};
            likelyReadEnd = new String[]{"CCGG", "CCGAGATCGG"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)pst[i1]-apek[i1]")) {
            theEnzyme = "PstI-ApeKI";
            initialCutSiteRemnant = new String[]{"TGCAG"};
            likelyReadEnd = new String[]{"GCAGC", "GCTGC", "CTGCAG", "GCAGAGAT", "GCTGAGAT"}; // look for ApeKI site, PstI site, or common adapter for ApeKI
            readEndCutSiteRemnantLength = 4;
        } else if (enzyme.matches("(?i)pst[i1]-ecot22[i1]")) {
            theEnzyme = "PstI-EcoT22I";
            initialCutSiteRemnant = new String[]{"TGCAG", "TGCAT"};
            likelyReadEnd = new String[]{"ATGCAT", "CTGCAG", "CTGCAAGAT", "ATGCAAGAT"}; // look for EcoT22I site, PstI site, or common adapter for PstI/EcoT22I
            readEndCutSiteRemnantLength = 5;
        } else if (enzyme.matches("(?i)pst[i1]-msp[i1]")) {
            theEnzyme = "PstI-MspI";
            initialCutSiteRemnant = new String[]{"TGCAG"};
            likelyReadEnd = new String[]{"CCGG", "CTGCAG", "CCGAGATC"}; // look for MspI site, PstI site, or common adapter for MspI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)pst[i1]-msp[i1]-GDFcustom")) {
            theEnzyme = "PstI-MspI-GDFcustom";
            initialCutSiteRemnant = new String[]{"TGCAG"};
            // changed from  CCGAGAT to CCGCTCAGG, as IGD/GDF used a custom Y adapter for MspI
            likelyReadEnd = new String[]{"CCGG", "CTGCAG", "CCGCTCAGG"}; // look for MspI site, PstI site, or GDF custom common adapter for MspI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)pst[i1]-taq[i1]")) {
            theEnzyme = "PstI-TaqI";
            initialCutSiteRemnant = new String[]{"TGCAG"};
            likelyReadEnd = new String[]{"TCGA", "CTGCAG", "TCGAGATC"}; // look for TaqI site, PstI site, or common adapter for TaqI
            readEndCutSiteRemnantLength = 3;
        } else if(enzyme.matches("(?i)PaeR7[i1]-Hha[i1]")) {
            theEnzyme = "PaeR7I-HhaI";
            // Requested by Ye, Songqing, use same Y adapter as Polland paper  -QS
            initialCutSiteRemnant=new String[]{"TCGAG"};
            likelyReadEnd = new String[]{"GCGC", "CTCGAG", "GCGAGATC"}; // look for HhaI site, PaeR7I site, or common adapter for HhaI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)sbf[i1]-msp[i1]")) {
            theEnzyme = "SbfI-MspI";
            initialCutSiteRemnant = new String[]{"TGCAGG"};
            likelyReadEnd = new String[]{"CCGG", "CCTGCAGG", "CCGAGATC"}; // look for MspI site, SbfI site, or common adapter for MspI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)asis[i1]-msp[i1]")) {
            theEnzyme = "AsiSI-MspI";
            initialCutSiteRemnant = new String[]{"ATCGC"};
            likelyReadEnd = new String[]{"CCGG", "GCGATCGC", "CCGAGATC"}; // look for MspI site, AsiSI site, or common adapter for MspI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)bsshii-msp[i1]|(?i)bssh2-msp[i1]")) {
            theEnzyme = "BssHII-MspI";
            initialCutSiteRemnant = new String[]{"CGCGC"};
            likelyReadEnd = new String[]{"CCGG", "GCGCGC", "CCGAGATC"}; // look for MspI site, BssHII site, or common adapter for MspI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)fse[i1]-msp[i1]")) {
            theEnzyme = "FseI-MspI";
            initialCutSiteRemnant = new String[]{"CCGGCC"};
            likelyReadEnd = new String[]{"CCGG", "GGCCGGCC", "CCGAGATC"}; // look for MspI site, FseI site, or common adapter for MspI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)sal[i1]-msp[i1]")) {
            theEnzyme = "SalI-MspI";
            initialCutSiteRemnant = new String[]{"TCGAC"};
            likelyReadEnd = new String[]{"CCGG", "GTCGAC", "CCGAGATC"}; // look for MspI site, SalI site, or common adapter for MspI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)ecor[i1]-msp[i1]")) {
            theEnzyme = "EcoRI-MspI";   //  G^AATTC  C^CGG
            initialCutSiteRemnant = new String[]{"AATTC"};
            likelyReadEnd = new String[]{"CCGG", "GAATTC", "CCGAGATC"}; // look for MspI site, EcoRI site, or Poland et al. 2012 Y-adapter for MspI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)hindiii-msp[i1]|(?i)hind3-msp[i1]")) {
            theEnzyme = "HindIII-MspI"; // A^AGCTT   C^CGG
            initialCutSiteRemnant = new String[]{"AGCTT"};
            likelyReadEnd = new String[]{"CCGG", "AAGCTT", "CCGAGATC"}; // look for MspI site, HindIII site, or Poland et al. 2012 Y-adapter for MspI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)sexa[i1]-sau3a[i1]")) {
            theEnzyme = "SexAI-Sau3AI";  // A^CCWGGT   ^GATC (not blunt)
            initialCutSiteRemnant = new String[]{"CCAGGT", "CCTGGT"};
            likelyReadEnd = new String[]{"GATC", "ACCAGGT", "ACCTGGT", "GATCAGATC"}; // look for SexAI site, Sau3AI site, or Poland et al. 2012 Y-adapter for Sau3AI
            readEndCutSiteRemnantLength = 4;
        } else if (enzyme.matches("(?i)bamh[i1l]-mluc[i1]")) {
            theEnzyme = "BamHI-MluCI";  // G^GATCC   ^AATT (not blunt)
            initialCutSiteRemnant = new String[]{"GATCC"};
            likelyReadEnd = new String[]{"AATT", "GGATCC", "AATTAGATC"}; // look for MluCI site, BamHI site, or Poland et al. 2012 Y-adapter for MluCI
            readEndCutSiteRemnantLength = 4;
        } else if (enzyme.matches("(?i)psti-msei|(?i)pst1-mse1")) {
            theEnzyme = "PstI-MseI"; // CTGCA^G   T^TAA
            initialCutSiteRemnant = new String[]{"TGCAG"};
            likelyReadEnd = new String[]{"TTAA", "CTGCAG", "TTAAGATC"}; // look for MseI site, PstI site, or Poland et al. 2012 Y-adapter for MseI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)avaii-msei|(?i)ava2-mse1")) {
            theEnzyme = "AvaII-MseI"; // G^GWCC   T^TAA  W=AorT
            initialCutSiteRemnant = new String[]{"GACC", "GTCC"};
            likelyReadEnd = new String[]{"TTAA", "GGACC", "GGTCC", "TTAAGATC"}; // look for MseI site, AvaII site, or Poland et al. 2012 Y-adapter for MseI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)ecori-msei|(?i)ecor1-mse1")) {
            theEnzyme = "EcoRI-MseI"; // G^AATTC   T^TAA
            initialCutSiteRemnant = new String[]{"AATTC"};
            likelyReadEnd = new String[]{"TTAA", "GAATTC", "TTAAGATC"}; // look for MseI site, EcoRI site, or Poland et al. 2012 Y-adapter for MseI
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)ecori-avaii|(?i)ecor1-ava2")) {
            theEnzyme = "EcoRI-AvaII"; // G^AATTC   G^GWCC
            initialCutSiteRemnant = new String[]{"AATTC"};
            likelyReadEnd = new String[]{"GGACC", "GGTCC", "GAATTC", "GGACAGATC", "GGTCAGATC"}; // look for AvaII site, EcoRI site, or Poland et al. 2012 Y-adapter for AvaII
            readEndCutSiteRemnantLength = 4;
        } else if (enzyme.matches("(?i)ecori-hinfi|(?i)ecor1-hinf1")) {
            theEnzyme = "EcoRI-HinfI"; // G^AATTC   G^ANTC
            initialCutSiteRemnant = new String[]{"AATTC"};
            likelyReadEnd = new String[]{"GAATC", "GACTC", "GAGTC", "GATTC", "GAATTC", "GAATAGATC", "GACTAGATC", "GAGTAGATC", "GATTAGATC"}; // look for HinfI site, EcoRI site, or Poland et al. 2012 Y-adapter for HinfI
            readEndCutSiteRemnantLength = 4;
        } else if (enzyme.matches("(?i)bbvci-mspi|(?i)bbvc1-msp1")) {
            theEnzyme = "BbvCI-MspI"; // CCTCAGC (-5/-2)   C^CGG
            initialCutSiteRemnant = new String[]{"TCAGC"};
            likelyReadEnd = new String[]{"CCGG", "CCTCAGC", "CCGAGATC"}; // look for MspI site, BbvCI site, or Poland et al. 2012 Y-adapter for MspI
            readEndCutSiteRemnantLength = 3;
        } else if(enzyme.matches("(?i)apo[i1]")){
            theEnzyme = "ApoI";
            initialCutSiteRemnant=new String[]{"AATTC","AATTT"};
            likelyReadEnd = new String[]{"AAATTC","AAATTT","GAATTC","GAATTT","AAATTAGAT","GAATTAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        } else if (enzyme.matches("(?i)BamH[i1l]")) {
            theEnzyme = "BamHI";
            initialCutSiteRemnant = new String[]{"GATCC"};
            likelyReadEnd = new String[]{"GGATCC", "GGATCAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            // likelyReadEnd = new String[]{"GGATCC", "AGATCGGAA", "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG"}; // <-- corrected from this by Jeff Glaubitz on 2012/09/12
            readEndCutSiteRemnantLength = 5;
        } else if (enzyme.matches("(?i)mse[i1]")) {
            theEnzyme = "MseI";
            initialCutSiteRemnant = new String[]{"TAA"};
            likelyReadEnd = new String[]{"TTAA", "TTAAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 3;
        } else if(enzyme.matches("(?i)Sau3A[i1]")){
            theEnzyme = "Sau3AI";
            initialCutSiteRemnant=new String[]{"GATC"};
            likelyReadEnd = new String[]{"GATC","GATCAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 4;
        } else if(enzyme.matches("(?i)nde[i1]")){
            theEnzyme = "NdeI";
            initialCutSiteRemnant=new String[]{"TATG"};
            likelyReadEnd = new String[]{"CATATG","CATAAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 4;
        } else if(enzyme.matches("(?i)hinp1[i1]")){
            theEnzyme = "HinP1I";
            initialCutSiteRemnant=new String[]{"CGC"};
            likelyReadEnd = new String[]{"GCGC","GCGAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 3;
        } else if(enzyme.matches("(?i)sbf[i1]")){
            theEnzyme = "SbfI";
            initialCutSiteRemnant=new String[]{"TGCAGG"};
            likelyReadEnd = new String[]{"CCTGCAGG","CCTGCAAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 6;
        } else if (enzyme.matches("(?i)hindiii|(?i)hind3")) {
            theEnzyme = "HindIII"; // A^AGCTT
            initialCutSiteRemnant=new String[]{"AGCTT"};
            likelyReadEnd = new String[]{"AAGCTT","AAGCTAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        } else if(enzyme.matches("(?i)ecor[i1]")) {
            theEnzyme = "EcoRI";  // G^AATTC
            initialCutSiteRemnant= new String[]{"AATTC"};
            likelyReadEnd = new String[]{"GAATTC","GAATTAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        } else if(enzyme.matches("(?i)cviq[i1]")){  
            theEnzyme = "CviQI";  // CviQI and Csp6I are isoschizomers (same recognition seq and overhang)
            initialCutSiteRemnant=new String[]{"TAC"};
            likelyReadEnd = new String[]{"GTAC","GTAAGATCGG"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 3;
        } else if(enzyme.matches("(?i)csp6[i1]")){  
            theEnzyme = "Csp6I";  // Csp6I and CviQI are isoschizomers (same recognition seq and overhang)
            initialCutSiteRemnant=new String[]{"TAC"};
            likelyReadEnd = new String[]{"GTAC","GTAAGATCGG"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)nlaiii|(?i)nla3")) {
            theEnzyme = "NlaIII"; // CATG^
            initialCutSiteRemnant=new String[]{"CATG"};
            likelyReadEnd = new String[]{"CATG","CATGAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 4;
        } else if(enzyme.matches("(?i)sph[i1]")){
            theEnzyme = "SphI";  // GCATG^C
            initialCutSiteRemnant=new String[]{"CATGC"};
            likelyReadEnd = new String[]{"GCATGC","GCATGAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        } else if(enzyme.matches("(?i)nsp[i1]")){
            theEnzyme = "NspI";  // RCATG^Y
            initialCutSiteRemnant=new String[]{"CATGC","CATGT"};
            likelyReadEnd = new String[]{"ACATGT","GCATGC","ACATGAGAT","GCATGAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        } else if(enzyme.matches("(?i)kpn[i1]")){
            theEnzyme = "KpnI";  // GGTAC^C
            initialCutSiteRemnant=new String[]{"GTACC"};
            likelyReadEnd = new String[]{"GGTACC","GGTACAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        } else if(enzyme.matches("(?i)sty[i1]")){
            theEnzyme = "StyI";  // C^CWWGG
            initialCutSiteRemnant=new String[]{"CAAGG","CATGG","CTAGG","CTTGG"};
            likelyReadEnd = new String[]{"CCAAGG","CCATGG","CCTAGG","CCTTGG","CCAAGAGAT","CCATGAGAT","CCTAGAGAT","CCTTGAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        } else if(enzyme.matches("(?i)styi-msei|(?i)sty1-mse1")){
            theEnzyme = "StyI-MseI";  // C^CWWGG & T^TAA
            initialCutSiteRemnant=new String[]{"CAAGG","CATGG","CTAGG","CTTGG"};
            likelyReadEnd = new String[]{"TTAA","CCAAGG","CCATGG","CCTAGG","CCTTGG","TTAAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 3;
        } else if(enzyme.matches("(?i)fse[i1]")){
            theEnzyme = "FseI";  // GGCCGG^CC
            initialCutSiteRemnant=new String[]{"CCGGCC"};
            likelyReadEnd = new String[]{"GGCCGGCC","AGATCGGAAG"}; // full cut site (from partial digest or chimera) or Morishige et al (BMC Genomics, 2013) T adapter start
            readEndCutSiteRemnantLength = 0;  // assumes that common T adapter is far more likely than a second full cut site
        } else if(enzyme.matches("(?i)NgoMIV|(?i)NgoM4")){
            theEnzyme = "NgoMIV";  // G^CCGGC
            initialCutSiteRemnant=new String[]{"CCGGC"};
            likelyReadEnd = new String[]{"GCCGGC","AGATCGGAAG"}; // full cut site (from partial digest or chimera) or Morishige et al (BMC Genomics, 2013) T adapter start
            readEndCutSiteRemnantLength = 0;  // assumes that common T adapter is far more likely than a second full cut site
        } else if(enzyme.matches("(?i)msl[i1]")){
            theEnzyme = "MslI";  // CAYNN^NNRTG  -- has 32 different cut sites (assuming constrained to palindromic YNN -- 32^2 otherwise)
            initialCutSiteRemnant=new String[]{""};
            likelyReadEnd = new String[]{"AGATCGGA"}; // common adapter start only (too many possible cut sites!)
            readEndCutSiteRemnantLength = 0;  
        } else if (enzyme.matches("(?i)RBSTA")) {
            theEnzyme = "RBSTA";
            initialCutSiteRemnant = new String[]{"TA"};
            likelyReadEnd = new String[]{"TTAA", "GTAC", "CTAG", "TTAAGAT", "GTAAGAT", "CTAAGAT"}; // full cut site (from partial digest or chimera) of MseI, CVIQi, XspI or common adapter start
            readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("(?i)RBSCG")) {
            theEnzyme = "RBSCG";
            initialCutSiteRemnant = new String[]{"CG"};
            likelyReadEnd = new String[]{"CCGC", "TCGA", "GCGC", "CCGG", "ACGT", "CCGAGAT", "TCGAGAT", "GCGAGAT", "ACGAGAT"}; // full cut site (from partial digest or chimera) of AciI, TaqaI, HinpI, HpaII, HpyCH4IV or common adapter start
            readEndCutSiteRemnantLength = 3;
        } else if(enzyme.matches("(?i)ignore")){
            theEnzyme = "unspecified";  // can be used for new enzymes -- only looks for barcodes and common adapter starts
            initialCutSiteRemnant=new String[]{""};
            likelyReadEnd = new String[]{"AGATCGGA"}; // common adapter start only
            readEndCutSiteRemnantLength = 0;  
        } else {
            System.out.println("The software didn't recognize your restriction enzyme (-e option).\n"
                    +"Currently, only the following enzymes are recognized for single enzyme digests:\n"
                    +"  ApeKI"    +"\n"
                    +"  ApoI"     +"\n"
                    +"  BamHI"    +"\n"
                    +"  Csp6I"    +"\n"
                    +"  CviQI"    +"\n"
                    +"  EcoRI"    +"\n"
                    +"  EcoT22I"  +"\n"
                    +"  FseI"     +"\n"
                    +"  HindIII"  +"\n"
                    +"  HinP1I"   +"\n"
                    +"  HpaII"    +"\n"
                    +"  KpnI"     +"\n"
                    +"  MseI"     +"\n"
                    +"  MslI"     +"\n"
                    +"  MspI"     +"\n"
                    +"  NdeI"     +"\n"
                    +"  NgoMIV"   +"\n"
                    +"  NlaIII"   +"\n"
                    +"  NspI"     +"\n"
                    +"  PasI"     +"\n"
                    +"  PstI"     +"\n"
                    +"  Sau3AI"   +"\n"
                    +"  SbfI"     +"\n"
                    +"  SphI"     +"\n"
                    +"  StyI"     +"\n"
                    +"  RBSTA"    +"\n"
                    +"  RBSCG"    +"\n"
                    +"  ignore"    +"\n"
                    +"Or the following for two-enzyme digests:\n"
                    +"  AsiSI-MspI"   +"\n"
                    +"  AvaII-MseI"   +"\n"
                    +"  BamHI-MluCI"  +"\n"
                    +"  BbvCI-MspI"   +"\n"
                    +"  BssHII-MspI"  +"\n"
                    +"  EcoRI-AvaII"  +"\n"
                    +"  EcoRI-HinfI"  +"\n"
                    +"  EcoRI-MseI"   +"\n"
                    +"  EcoRI-MspI"   +"\n"
                    +"  FseI-MspI"    +"\n"
                    +"  HindIII-MspI" +"\n"
                    +"  HindIII-NlaIII" +"\n"
                    +"  PaeR7I-HhaI"  +"\n"
                    +"  PstI-ApeKI"   +"\n"
                    +"  PstI-EcoT22I" +"\n"
                    +"  PstI-MseI"    +"\n"
                    +"  PstI-MspI"    +"\n"
                    +"  PstI-MspI-GDFcustom"+"\n"
                    +"  PstI-TaqI"    +"\n"
                    +"  SalI-MspI"    +"\n"
                    +"  SbfI-MspI"    +"\n"
                    +"  SexAI-Sau3AI" +"\n"
                    +"  StyI-MseI"    +"\n"
                    +"  ignore"    +"\n"
            );
            System.out.println("For two-enzyme digest, enzyme names should be separated by a dash, e.g. -e PstI-MspI");
            System.out.println("\nIf your enzyme is not on the above list you can use \"-e ignore\". In this case\n"
                    +"   barcodes and common adapter start sequences will be recognized, but chimeric DNA\n"
                    +"   fragments (or partial digests) will not be trimmed.");
        }
        System.out.println("Enzyme: " + theEnzyme);
    }

    private String getKeyFileEnzyme(String keyFileName) {
        String result = null;
        try {
            BufferedReader br = new BufferedReader(new FileReader(keyFileName), 65536);

            String temp;
            int currLine = 0;
            while (((temp = br.readLine()) != null)) {
                String[] s = temp.split("\\t");  //split by whitespace
                if (currLine > 0) {
                    String enzymeName;
                    if (s.length < 9) {
                        enzymeName = "";
                    } else {
                        enzymeName = s[8];
                    }
                    if (!enzymeName.equals("")) {
                        result = enzymeName;
                        break;
                    }
                }
                currLine++;
            }
        } catch (Exception e) {
            System.out.println("Couldn't open key file to read Enzyme: " + e);
        }
        return result;
    }

    /**
     * Reads in an Illumina key file, creates a linear array of {@link Barcode} objects
     * representing the barcodes in the key file, then creates a hash map containing
     * indices from the linear array indexed by sequence.  The names of barcode objects
     * follow the pattern samplename:flowcell:lane:well, since sample names alone are not unique.
     *
     * @param keyFile Illumina key file.
     * @param flowcell Only barcodes from this flowcell will be added to the array.
     * @param lane Only barcodes from this lane will be added to the array.
     * @return Number of barcodes in the array.  
     */
    private int setupBarcodeFiles(File keyFile, String flowcell, String lane) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(keyFile), 65536);
            ArrayList<Barcode> theBarcodesArrayList = new ArrayList<Barcode>();
            String temp;
            int k = 0;
            while (((temp = br.readLine()) != null)) {
                String[] s = temp.split("\\t");  //split by whitespace
                Barcode theBC = null;
                if (s[0].equals(flowcell) && s[1].equals(lane)) {
                    String well = (s[6].length() < 2) ? (s[5] + '0' + s[6]) : s[5] + s[6];
                    if (s.length < 8 || s[7] == null || s[7].equals("")) {  // use the plate and well
                        theBC = new Barcode(s[2], initialCutSiteRemnant, s[3] + ":" + s[0] + ":" + s[1] + ":" + s[4] + ":" + well, flowcell, lane, k++);
                    } else {  // use the "libraryPlateWellID" or whatever is in column H of the key file, IF it is an integer
                        try {
                            int libPrepID = Integer.parseInt(s[7]);
                            theBC = new Barcode(s[2], initialCutSiteRemnant, s[3] + ":" + s[0] + ":" + s[1] + ":" + libPrepID, flowcell, lane, k++);
                        } catch (NumberFormatException nfe) {
                            theBC = new Barcode(s[2], initialCutSiteRemnant, s[3] + ":" + s[0] + ":" + s[1] + ":" + s[4] + ":" + well, flowcell, lane, k++);
                        }
                    }
                    theBarcodesArrayList.add(theBC);
                    System.out.println(theBC.barcodeS + " " + theBC.taxaName);
                    theBC = new Barcode("N"+theBC.barcodeS.substring(1), initialCutSiteRemnant, theBC.taxaName, theBC.flowcell, theBC.lane, k++);
                    theBarcodesArrayList.add(theBC);
                    System.out.println(theBC.barcodeS + " " + theBC.taxaName);
                }
            }
            theBarcodes = new Barcode[theBarcodesArrayList.size()];
            theBarcodesArrayList.toArray(theBarcodes);
            Arrays.sort(theBarcodes);
            int nBL = theBarcodes[0].barOverLong.length;
            
            quickBarcodeList = new long[theBarcodes.length * nBL];
            quickMap = new HashMap<Long, Integer>();
            for (int i = 0; i < theBarcodes.length; i++) {
                for (int j = 0; j < nBL; j++) {
                    quickBarcodeList[i * nBL + j] = theBarcodes[i].barOverLong[j];
                    quickMap.put(theBarcodes[i].barOverLong[j], i);
                }
            }
            Arrays.sort(quickBarcodeList);
        } catch (Exception e) {
            System.out.println("Error with setupBarcodeFiles: " + e);
        }

        return theBarcodes.length;
    }

    /**
     * Returns the best barcode match for a given sequence.  
     * @param queryS query sequence to be tested against all barcodes
     * @param maxDivergence maximum divergence to permit
     * @return best barcode match (null if no good match)
     */
    Barcode findBestBarcode(String queryS, int maxDivergence) {
        long query = BaseEncoder.getLongFromSeq(queryS.substring(0,chunkSize));
        
        //note because the barcodes are polyA after the sequence, they should always
        //sort ahead of the hit, this is the reason for the -(closestHit+2)
        
        
        int closestHit = Arrays.binarySearch(quickBarcodeList, query);
        /*      THIS IS THE NEW PIPELINE APPROACH THAT DOES NOT WORK
         if(closestHit>-2) return null; //hit or perfect
         if((query&quickBarcodeList[-(closestHit+2)])!=quickBarcodeList[-(closestHit+2)]) return null;
         int index =quickMap.get(quickBarcodeList[-(closestHit+2)]);
         //      System.out.println(theBarcodes[index].barcodeS);
         return theBarcodes[index];
         //note to see if it is a perfect match you can just bit AND
         */

        //  Below is the old pipeline approach, which works (at least for maxDivergence of 0)
        if (closestHit < -1) {  // should always be true, as the barcode+overhang is padded to 32 bases with polyA
            int index = quickMap.get(quickBarcodeList[-(closestHit + 2)]);
            if (theBarcodes[index].compareSequence(query, 1) == 0) {
                return theBarcodes[index];
            } else if (maxDivergence == 0) {
                return null;  // return null if not a perfect match
            }
        } else {
            return null;  // should never go to this line
        }
        int maxLength = 0, minDiv = maxDivergence + 1, countBest = 0;
        Barcode bestBC = null;
        for (Barcode bc : theBarcodes) {
            int div = bc.compareSequence(query, maxDivergence + 1);
            if (div <= minDiv) {
                if ((div < minDiv) || (bc.barOverLength > maxLength)) {
                    minDiv = div;
                    maxLength = bc.barOverLength;
                    bestBC = bc;
                    countBest = 1;
                } else {  //it is a tie, so return that not resolvable
                    bestBC = null;
                    countBest++;
                }
            }
        }
        return bestBC;
    }

    /**
     * The barcode libraries used for this study can include two types of
     * extraneous sequence at the end of reads. The first are chimeras created
     * with the free ends. These will recreate the restriction site. The second
     * are short regions (less than 64bp), so that will they will contain a
     * portion of site and the universal adapter. This finds the first of site
     * in likelyReadEnd, keeps the restriction site overhang and then sets
     * everything to polyA afterwards
     *
     * @param seq An unprocessed tag sequence.
     * @param maxLength The maximum number of bp in the processed sequence.
     * @return returnValue A ReadBarcodeResult object containing the unprocessed
     * tag, Cut site position, Processed tag, and Poly-A padded tag.
     */


    private static String commonAdapter = "AGATCGGAAGAGCGGTTC";
    
    public static int removeSeqAfterSecondCutSite(String seq) {
        //this looks for a second restriction site or the common adapter start, and then turns the remaining sequence to AAAA
        int cutSitePosition = 9999;

        String processedSequence = null;
        //Look for cut sites, starting at a point past the length of the initial cut site remnant that all reads begin with
        String match = null;
        for (String potentialCutSite : likelyReadEnd) {
            int p = seq.indexOf(potentialCutSite, 1);
            if ((p > 1) && (p < cutSitePosition)) {
                cutSitePosition = p;
                match = potentialCutSite;
            }
        }
        
        if (false && theEnzyme.equalsIgnoreCase("ApeKI") && cutSitePosition == 2
                && (match.equalsIgnoreCase("GCAGC") || match.equalsIgnoreCase("GCTGC"))) {  // overlapping ApeKI cut site: GCWGCWGC
            seq = seq.substring(3, seq.length());  // trim off the initial GCW from GCWGCWGC
            cutSitePosition = 9999;
            for (String potentialCutSite : likelyReadEnd) {
                int p = seq.indexOf(potentialCutSite, 1);
                if ((p > 1) && (p < cutSitePosition)) {
                    cutSitePosition = p;
                }
            }
        }
        
        if (cutSitePosition < 9999) {  // Cut site found
            //Trim tag to sequence up to & including the cut site
            cutSitePosition += readEndCutSiteRemnantLength;
        } else {
            if (seq.length() <= 0) {
                //If cut site is missing because there is no sequence
            	cutSitePosition = -1;
            } else {
                //If cut site is missing because it is beyond the end of the sequence (or not present at all)
                cutSitePosition = seq.length();
            }
        }

        int i = seq.substring(0, cutSitePosition).indexOf(commonAdapter);
        return i<0 ? cutSitePosition : i;
    }


	/**
     * Return a {@link ReadBarcodeResult} that captures the processed read and taxa
     * inferred by the barcode
     * @param seqS DNA sequence from the sequencer
     * @param qualS quality score string from the sequencer
     * @param fastq (fastq = true?; qseq=false?)
     * @param minQual minimum quality score
     * @return If barcode and cut site was found returns the result and
     * processed sequence, if the barcode and cut site were not found return
     * null
     */
    public ReadBarcodeResult parseReadIntoTagAndTaxa2(String seqS, String qualS, int minQual) {
    	if(seqS==null) return null;
    	
        Barcode bestBarcode = findBestBarcode(seqS, maximumMismatchInBarcodeAndOverhang);
        if (bestBarcode == null) {
            return null;  //overhang missing so skip
        }
        
        int cutSitePosition = removeSeqAfterSecondCutSite(seqS.substring(bestBarcode.barLength));
        String processedSequence = null;
        //if(cutSitePosition-bestBarcode.barLength<36) return null;
        if(cutSitePosition==-1) {
        	return null;
        } else {
        	if ((minQual > 0) && (qualS != null)) {
        		StringBuilder seqSB = new StringBuilder(seqS);
        		for(int i=bestBarcode.barLength; i<cutSitePosition; i++)
        			if(qualS.charAt(i)-33 < minQual) 
        				seqSB.setCharAt(i,'N');
        		processedSequence = seqSB.toString().substring(bestBarcode.barLength, cutSitePosition);
            } else {
            	processedSequence = seqS.substring(bestBarcode.barLength, 
            			bestBarcode.barLength+cutSitePosition);
            }
        }
        processedSequence = StringUtils.strip(processedSequence, "N");
        
        ReadBarcodeResult rbr = new ReadBarcodeResult(BaseEncoder.getBitSetFromSeq(
        		processedSequence), 
        		bestBarcode.getTaxaId()/2);
        return rbr;
    }
    
    public ReadBarcodeResult parseReadIntoTagAndTaxa(String seqS, String qualS, int minQual) {
    	if(seqS==null) return null;
    	
    	if(minQual!=0 && seqS!=null) {
    		try {
				final char[] seqC = (char[]) field.get(seqS),
						qualC = (char[]) field.get(qualS);
				int len = seqC.length;
				for(int i=0; i<len; i++)
					if(qualC[i]-33<minQual)
						seqC[i] = 'N';
				seqS = String.valueOf(seqC);
			} catch (IllegalArgumentException | IllegalAccessException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
    	}
    	
        Barcode bestBarcode = findBestBarcode(seqS, maximumMismatchInBarcodeAndOverhang);
        if (bestBarcode == null) {
            return null;  //overhang missing so skip
        }
        
        seqS = StringUtils.strip(seqS.substring(bestBarcode.barLength),"N");
        if(seqS.length()<16) return null;
        int cutSitePosition = Math.min(SWSearch(seqS, likelyReadEnd[0]), seqS.length());
        String processedSeqS = StringUtils.strip(seqS.substring(0, cutSitePosition),"N");
        if(processedSeqS.length()>=16)
        	return new ReadBarcodeResult(BaseEncoder.getBitSetFromSeq(
        			processedSeqS), 
        			bestBarcode.getTaxaId()/2);
        return null;
    }
    
   /**Returns the number of barcodes for the flowcell and lane*/
    public int getBarCodeCount() {
        return theBarcodes.length;
    }

    /**Returns the {@link Barcode} for the flowcell and lane*/
    public Barcode getTheBarcodes(int index) {
        return theBarcodes[index];
    }

    /**Returns the taxaNames for the flowcell and lane*/
    public String[] getTaxaNames() {
        String[] result = new String[getBarCodeCount()];
        for (int i = 0; i < result.length; i++) {
            result[i] = getTheBarcodes(i).getTaxaName();
        }
        return result;
    }

	public String[] getSortedTaxaNames() {
		// TODO Auto-generated method stub
		String[] result = new String[getBarCodeCount()];
		for (int i = 0; i < theBarcodes.length; i++) 
            result[theBarcodes[i].taxonId] = theBarcodes[i].taxaName;
		String[] taxa = new String[getBarCodeCount()/2];
		for(int i=0; i<taxa.length; i++)
			taxa[i] = result[2*i];
        return taxa;
	}
	
	public static int SWSearch(final String refS, final String queryS) {
		try {
			final char[] ref = (char[]) field.get(refS),
					query = (char[]) field.get(queryS);
			int nr = ref.length, nq = query.length;
			int[][] mD = new int[nq+1][nr+1];
			for(int i=1; i<=nq-1; i++) 
				for(int j=i, b=Math.min(nr-6+i, nr-1); j<=b; j++)
					mD[i][j] = mD[i-1][j-1] + 
							(query[i-1]==ref[j-1] ? 1 : 0);
			double score = 0.0, s;
			int ijs = -1, ijn = -1;
			int nonelead = 0, n;
			for(int i=6; i<=nq; i++) {
				mD[i][nr] = mD[i-1][nr-1]+(query[i-1]==ref[nr-1] ? 1 : 0);
				s = ((double)mD[i][nr])/i;
				if(s>=score) {
					score = s;
					ijs = nr-i+5;
				}
				n = mD[i][nr]-mD[5][nr-i+5];
				if(n>=nonelead) {
					nonelead = n;
					ijn = nr-i+5;
				}
			}
			for(int j=nr-1; j>=nq; j--) {
				mD[nq][j] = mD[nq-1][j-1]+(query[nq-1]==ref[j-1] ? 1 : 0);
				s = ((double)mD[nq][j])/nq;
				if(s>=score) {
					score = s;
					ijs = j-nq+5;
				}
				n = mD[nq][j] -mD[5][j-nq+5];
				if(n>=nonelead) {
					nonelead = n;
					ijn = j-nq+5;
				}
			}
			//matrix(ref, query, mD);
			//System.out.println(score);
			if(score>=.75) return ijs;
			else if(nonelead>8) {
				for(int i=0; i<Math.min(nr-ijn-11,7); i++)
					nonelead += ref[ijn+11+i]
							==commonAdapter.charAt(i+11) ? 1 : 0;
				if(nonelead>15) return ijn;
			}
		} catch (IllegalArgumentException | IllegalAccessException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return 9999;
	}
	
	private static void matrix(char[] ref, char[] query, int[][] mD) {
		System.out.print("\t\t");
		for(int i=0; i<ref.length; i++)
			System.out.print(ref[i]+"\t");
		System.out.println();
		for(int i=-1; i<query.length; i++) {
			if(i==-1) System.out.print("\t");
			else System.out.print(query[i]+"\t");
			for(int j=0; j<mD[i+1].length; j++)
				System.out.print(mD[i+1][j]+"\t");
			System.out.println();
		}
	}
}
