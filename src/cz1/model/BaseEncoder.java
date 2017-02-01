package cz1.model;

import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;

/**
 * Utility class for encoding tags into longs.
 * <p>
 * Sequencing reads are chunked into 32bp and recorded in a 64-bit long.  Only 
 * A (00), C (01), G (10), T (11) are encoded.  Any other character sets the entire long to -1.
 * Missing data at the end is padded with poly-A or (0).  This missing end, is tracked
 * by the tag length attribute.
 * <p>
 * Some of these methods should be transitioned to {@link net.maizegenetics.pal.alignment.NucleotideAlignmentConstants},
 * however, BaseEncoder only supports four states, while NucleotideAlignment includes gaps, insertions, and missing.
 * 
 * @author Ed Buckler
 */
public class BaseEncoder {

	public static final int chunkSize = 21;
	/** defines the number of bases fitting with a long */
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
     * Returns a long for a sequence in a String
     * @param seq
     * @return 2-bit encode sequence (-1 if an invalid sequence state is provided e.g. N)
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
     * Returns a long for a sequence in a String
     * @param seq
     * @return 3-bit encode sequence (-1 if an invalid sequence state is provided e.g. N)
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
     * Returns the reverse complement of a sequence already encoded in a 2-bit long.
     * <p>
     * Note: polyA is used represent unknown, but reverse complement will change it to polyT which does not mean the same
     * sometimes it is best to reverseComplement by text below
     * @param seq  2-bit encoded sequence
     * @param len  length of the sequence
     * @return  2-bit reverse complement
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
     * Returns a string based reverse complement.  Get around issues with the poly-A tailing in the 2-bit encoding approach.
     *
     * @param seq  DNA sequence
     * @return  reverse complement DNA sequence
     */
    public static String getReverseComplement(String seq) {
        StringBuilder sb = new StringBuilder(seq.length());
        for (int i = seq.length() - 1; i >= 0; i--) {
            sb.append(cBaseMap.get(seq.charAt(i)));
        }
        return sb.toString();
    }

     /**
     * Return a string representation of the 2-bit encoded long.
     * @param val 2-bit encoded sequence
     * @param len length of the sequence
     * @return DNA sequence as a string
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
    
    public static byte seqDifferencesForSubset(long seq1, long seq2, int lengthOfComp, int maxDivergence) {
        long mask = 7;
        byte cnt = 0;
        long diff = seq1 ^ seq2;
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
     * Returns the position of the first low quality positions based on a quality
     * fastq (?) string.
     * @param quality fastq quality string
     * @param minQual minimum quality threshold
     * @return position of first low quality position (quality length is returned is not low
     * quality base is found.
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
