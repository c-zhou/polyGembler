package cz1.gbs.model;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.lang.reflect.Field;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import org.apache.commons.lang3.StringUtils;

import cz1.gbs.core.Barcode;
import cz1.gbs.core.BaseEncoder;
import cz1.gbs.core.ReadBarcodeResult;

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
            likelyReadEnd = new String[]{"CTGCAG", "CTGCAAGAT"}; // full cut site (from partial digest or chimera) or common adapter start
            readEndCutSiteRemnantLength = 5;
        } else if (enzyme.matches("CviAII")) {
        	theEnzyme = "CviAII";
        	initialCutSiteRemnant = new String[]{"ATG"};
        	likelyReadEnd = new String[]{"CATG", "CATAGATCG"};
        	readEndCutSiteRemnantLength = 3;
        } else if (enzyme.matches("TseI")) {
        	theEnzyme = "TseI";
        	initialCutSiteRemnant = new String[]{"CAGC","CTGC"};
        	likelyReadEnd = new String[]{"GCAGC", "GCTGC", "GCAGAGAT", "GCTGAGAT"};
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
        } else {
            System.out.println("The software didn't recognize your restriction enzyme (-e option).\n"
                    +"Currently, only the following enzymes are recognized for single enzyme digests:\n"
                    +"  ApeKI"    +"\n"
                    +"  CviAII"    +"\n"
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
                    +"  TseI"    +"\n"
                    +"  RBSTA"    +"\n"
                    +"  RBSCG"    +"\n"
            );
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
                    System.out.println(theBC.getBarcodeString() + " " + theBC.getTaxaName());
                }
            }
            br.close();
            
            theBarcodes = new Barcode[theBarcodesArrayList.size()];
            theBarcodesArrayList.toArray(theBarcodes);
            Arrays.sort(theBarcodes);
            int nBL = theBarcodes[0].getBarOverLong().length;
            
            quickBarcodeList = new long[theBarcodes.length * nBL];
            quickMap = new HashMap<Long, Integer>();
            for (int i = 0; i < theBarcodes.length; i++) {
                for (int j = 0; j < nBL; j++) {
                    quickBarcodeList[i * nBL + j] = theBarcodes[i].getBarOverLong()[j];
                    quickMap.put(theBarcodes[i].getBarOverLong()[j], i);
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
        int maxLength = 0, minDiv = maxDivergence + 1;
        Barcode bestBC = null;
        for (Barcode bc : theBarcodes) {
            int div = bc.compareSequence(query, maxDivergence + 1);
            if (div <= minDiv) {
                if ((div < minDiv) || (bc.getBarOverLength() > maxLength)) {
                    minDiv = div;
                    maxLength = bc.getBarOverLength();
                    bestBC = bc;
                } else {  //it is a tie, so return that not resolvable
                    bestBC = null;
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
        
        seqS = StringUtils.strip(seqS.substring(bestBarcode.getBarLength()),"N");
        if(seqS.length()<32) return null;
        int cutSitePosition = seqS.length();
        for(String end : likelyReadEnd) {
        	int w = seqS.indexOf(end);
        	if(w>=0 && w<cutSitePosition) cutSitePosition = w; 
        }
        String processedSeqS = StringUtils.strip(seqS.substring(0, cutSitePosition),"N");
        if(processedSeqS.length()>=32)
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
            result[theBarcodes[i].getTaxaId()] = theBarcodes[i].getTaxaName();
		String[] taxa = new String[getBarCodeCount()/2];
		for(int i=0; i<taxa.length; i++)
			taxa[i] = result[2*i];
        return taxa;
	}
}
