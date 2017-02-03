package cz1.gbs.core;

import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;

public class BaseEncoder {

	public static final int chunkSize = 21;
    public static final char[] bases = {'A', 'C', 'G', 'T', 'N'};
    public static final Map<Integer, Character> baseMap;
    static
    {
    	baseMap = new HashMap<Integer, Character>();
    	baseMap.put(1, 'A');
    	baseMap.put(6, 'C');
    	baseMap.put(4, 'G');
    	baseMap.put(3, 'T');
    	baseMap.put(0, 'N');
    	baseMap.put(7, 'N');
    }
    
    public static final Map<Character, Character> cBaseMap;
    static
    {
    	cBaseMap = new HashMap<Character, Character>();
    	cBaseMap.put('A', 'T');
    	cBaseMap.put('C', 'G');
    	cBaseMap.put('T', 'A');
    	cBaseMap.put('G', 'C');
    	cBaseMap.put('N', 'N');
    }
    
    private BaseEncoder() {
    }

    /**
     * Returns a BitSet for a sequence
     * @param seq
     * @return 3-bit encoded sequence
     *	1		001			A
     *	6		110			C
     *	4		100			G
     *	3		011			T
     *	5	 	111/000		N
     */
    public static BitSet getBitSetFromSeq(String seq) {
        int seqLength = seq.length();
        BitSet v = new BitSet(seqLength*3+1);
        v.set(seqLength*3);
        for (int i = 0; i < seqLength; i++) {
            switch (seq.charAt(i)) {
            	case 'A':
                case 'a':
                	v.set(i*3+2);
                    break;
                case 'C':
                case 'c':
                    v.set(i*3);
                    v.set(i*3+1);
                    break;
                case 'G':
                case 'g':
                	v.set(i*3);
                    break;
                case 'T':
                case 't':
                	v.set(i*3+1);
                	v.set(i*3+2);
                    break;
                case 'N':
                case 'n':
                case '.':
                	break;
                default:
                    return null;
            }
        }
        return v;
    }
        
    /**
     * Returns a long for a sequence (demultiplexing)
     * @param seq
     * @return 2-bit encode sequence
     */
    public static long getLongFromSeq(String seq) {
        int seqLength = seq.length();
        long v = 0;
        for (int i = 0; i < seqLength; i++) {
            switch (seq.charAt(i)) {
            	case 'N':
            	case 'n':
            		v=v<<3;
            		break;
                case 'A':
                case 'a':
                    v = (v << 3) + (byte) 1;
                    break;
                case 'C':
                case 'c':
                    v = (v << 3) + (byte) 2;
                    break;
                case 'G':
                case 'g':
                    v = (v << 3) + (byte) 3;
                    break;
                case 'T':
                case 't':
                    v = (v << 3) + (byte) 4;
                    break;
                default:
                    return -1;
            }
        }
        if (seqLength == chunkSize) {
            return v;
        }
        if (seqLength > chunkSize) {
            return -1;
        }
        v = (v << (3 * (chunkSize - seqLength))); //if shorter fill with NNNN
        return v;
    }

    /**
     * Returns the reverse complement of a sequence
     * @param seq DNA sequence
     * @return reverse complement sequence
     */
    public static String getReverseComplement(String seq) {
        StringBuilder sb = new StringBuilder(seq.length());
        for (int i = seq.length() - 1; i >= 0; i--) {
            sb.append(cBaseMap.get(seq.charAt(i)));
        }
        return sb.toString();
    }

    /**
     * Returns the reverse complement of a encoded sequence
     * @param seq 3-bit encoded sequence
     * @return 3-bit encoded reverse complement sequence
     */
    public static BitSet getReverseComplement(BitSet seq) {
        int len = seq.length();
    	BitSet rev = new BitSet(len);
        for (int i = 0; i < len-1; i++) 
        	if( !seq.get(len-i-2) )
        		rev.set(i);
        rev.set(len-1);
        return rev;
    }
    
     /**
     * Return a sequence represented by a 3-bit encoded BitSet
     * @param val 3-bit encoded sequence
     * @return sequence
     */
    public static String getSequenceFromBitSet(BitSet val) {
        StringBuilder seq = new StringBuilder( (val.length()-1)/3 );
        int mask;
        for (int i=0;i<val.length()-1; i+=3) {
            mask = 0;
        	for(int j=0; j<3; j++) 
        		mask = (mask << 1) + (val.get(i+j) ? 1 : 0);
        	seq.append(baseMap.get(mask));
        }
        return seq.toString();
    }
    
    /**
     * Return the divergence of two long-encoded sequence (demultiplexing)
     * @param val1 long encoded sequence
     * @param val2 long encoded sequence
     * @param lengthOfComp number of bases for comparison
     * @param maxDivergence stops when divergence exceeds this number
     * @return divergence
     */
    public static byte seqDifferencesForSubset(long val1, long val2, 
    		int lengthOfComp, int maxDivergence) {
        long mask = 7;
        byte cnt = 0;
        long diff = val1 ^ val2;
        diff = diff >> (3 * (chunkSize - lengthOfComp));  //shift to 5' end of sequence
        for (int x = 0; x < lengthOfComp && cnt < maxDivergence; x++) {
            if ((diff & mask) > 0) {
                cnt++;
            }
            diff = diff >> 3;
        }
        return cnt;
    }
    
    /**
     * Returns the position of the first low quality positions based on a quality string
     * @param quality fastq quality string
     * @param minQual minimum quality threshold
     * @return position of first low quality position (quality length is returned if no low
     * quality base is found).
     */
    public static int getFirstLowQualityPos(String quality, int minQual) {
        int qualInt = 0;
        for (int i = 0; i < quality.length(); i++) {
            qualInt = (int) quality.charAt(i) - 33;
            if (qualInt < minQual) {
                return i;
            }
        }
        return quality.length();
    }
}
