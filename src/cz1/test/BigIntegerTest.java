package cz1.test;

import java.lang.reflect.Field;
import java.math.BigInteger;
import java.util.BitSet;
import java.util.Random;

import org.apache.commons.math3.stat.StatUtils;

public class BigIntegerTest {
	
	public static void main2(String[] args) {
		for(int i=0; i<8; i++)
			System.out.println((char) i);
		
		Random rand = new Random();
		for(int r=0; r<9; r++) {	
			int k=0;
			long start, end;
			
			boolean[] b3 = new boolean[576];
			boolean[] b4 = new boolean[576];
			
			for(int i=0; i<576; i++) {
				b3[rand.nextInt(576)] = true;
				b4[rand.nextInt(576)] = true;
			}
			
			int[] check = new int[300];
			for(int i=0; i<300; i++)
				check[i] = rand.nextInt(576);
			start = System.nanoTime();
			for(int i=0; i<10000000; i++) {
				k=0;
				for(int j=0; j<300; j++)
					if(b3[j]&&b4[j])
						k++;
			}
			end = System.nanoTime();
			System.out.println(k+"	"+(end-start));

			BitSet b5 = new BitSet(4000);
			BitSet b6 = new BitSet(4000);
			
			for(int i=0; i<576; i++) {
				b5.set(rand.nextInt(4000));
				b6.set(rand.nextInt(4000));
			}
			
			start = System.nanoTime();
			for(int i=0; i<10000000; i++) {
				b5.or(b6);
				k=b5.cardinality();
			}
			end = System.nanoTime();
			System.out.println(k+"	"+(end-start));

			BigInteger b1 = new BigInteger(576, rand);
			BigInteger b2 = new BigInteger(576, rand);
			start = System.nanoTime();
			for(int i=0; i<10000000; i++)
				k = b1.and(b2).bitCount();
			end = System.nanoTime();
			System.out.println(k+"	"+(end-start));
		
		}
	}
	

	public static void main(String[] args) {
		
		double[][] ds = new double[3][];
		ds[0] = null;
		ds[1] = new double[]{1,2,3};
		
		double[][] ds2 = new double[][]{ds[1]};
		System.out.println(StatUtils.sum(ds2[0],1,1));
		
		System.out.println("etc123".matches("^etc[0-9]{1,4}$"));
		System.out.println("etc123".replaceAll("^etc",""));
		System.out.println("etc123_12121_212121_21212".replaceAll("_[0-9]{1,}$",""));
		System.out.println("etc123_12121_212121_c333321212".replaceAll(".*[^\\d](\\d+).*", "$1"));
		
		int[] key4 = new int[]{1164, 1164, 1035, 1035};
		int shift_bits2 = 12;
		long key = ((((((long) key4[0]<<shift_bits2)+key4[1])<<shift_bits2)+key4[2])<<shift_bits2)+key4[3];
		System.out.println(key);
		System.out.println(new BigInteger("-1941916661"));
		System.exit(1);
		
		long start, end;
		
		char a = 1;
		System.out.println(a<<20);
		System.out.println((int) Math.ceil(Math.log(4+1)/Math.log(2)));
		Random rand = new Random();
		BitSet b5 = new BitSet1(999);
		BitSet b6 = new BitSet1(3);
		BitSet b7;
		
		for(int i=0; i<1000; i++) 
			b5.set(rand.nextInt(999));
		for(int i=0; i<3; i++)
			b6.set(i);
		
		for(int r=0; r<9; r++) {	
			int k=0;
			start = System.nanoTime();
			for(int i=0; i<100000; i++) {
				for(int j=0; j<333; j++) {
					b7 = b5.get(j*3, j*3+3);
					//b7.and(b6);
					//k=b7.cardinality();
				}
			}
			end = System.nanoTime();
			System.out.println("	"+b5.length()+"	"+(end-start));
		}
	}
	
	public static void main3(String[] args) {
		long start, end;
		
		Random rand = new Random();
		BigInteger b5 = new BigInteger(4003, rand);
		BigInteger b6 = new BigInteger("3");
		BigInteger b7;
		
		for(int r=0; r<9; r++) {	
			int k=0;
			b7 = b5;
			start = System.nanoTime();
			for(int i=0; i<1000; i++) {
				//b7 = b5;
				for(int j=0; j<1334; j++) {
					b7.and(b6);
					b7.shiftRight(3);
				}
			}
			end = System.nanoTime();
			System.out.println(k+"	"+b5.bitCount()+"	"+(end-start));
		}
	}
}

class BitSet1 extends BitSet {

    private long[] words;
    private static Field wordsField;

    static {
        try {
            wordsField = BitSet.class.getDeclaredField("words");
            wordsField.setAccessible(true);
        } catch (NoSuchFieldException e) {
            throw new IllegalStateException(e);
        }
    }
    
    public BitSet1(final int regLength) {
        super(regLength);
        try {
            words = (long[]) wordsField.get(this);
        } catch (IllegalAccessException e) {
            throw new IllegalStateException(e);
        }
    }
    
    public byte[] toByteArray(int bits){
    	byte[] bytes = new byte[(int) Math.
    	                        ceil(1.0*this.length()/bits)];
    	byte mask = (byte)(Math.pow(2,bits)-1);
    	int r = 64%bits, carry = 0, shift;
    	long w;
    	int k = 0;
    	byte residual = 0;
    	for(int i=0; i<words.length; i++) {
    		w = words[i];
    		if(carry>0) {
    			shift = bits-carry;
    			byte maskR = (byte)(Math.pow(2,shift)-1);
    			bytes[k++] = (byte)(((maskR&w)<<shift)+
    					residual);
    			w>>>=shift;
    		}
    		bytes[k++] = (byte)(mask&w);
    		
    	}
    	return null;
    }
}
