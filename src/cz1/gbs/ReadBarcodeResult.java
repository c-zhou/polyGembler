/*
 * ReadBarcodeResult
 */
package cz1.gbs;

import java.util.BitSet;

/**
 * Container class for returning the results of parsed barcoded sequencing read.
 * <p>
 * {@link #unprocessedSequence} is the original sequence with the barcode still attached
 * <p>
 * {@link #processedSequence} is the sequence with the barcode removed and cut down to the given {@link length}.
 * <p>
 * {@link #paddedSequence} is the {@link #processedSequence} padded with polyA to length if the 
 * {@link #processedSequence} was shorter than {@link length}.
 * 
 * @author Ed Buckler
 */
public class ReadBarcodeResult {
    /**Original sequence from sequencer*/
    public String unprocessedSequence = null;
    /**Processed sequence with barcode removed*/
    public String processedSequence = null;
    /**length of the processed sequence*/
    byte length;
    /**Sequence encoded in 2-bit long array*/
    public BitSet read;
    /**Taxon name implied by the barcode sequence*/
    String taxonName;
    /**Taxon index implied by the barcode sequence*/
    public int taxonId;
    
    //TODO this instantiation should also include the orginal unprocessedSequence, processedSequence, and paddedSequence - the the object encode it
    public ReadBarcodeResult(BitSet read, String taxon) {
        this.read = read;
        this.taxonName = taxon;
    }
    
    public ReadBarcodeResult(BitSet read, int taxon) {
        this.read = read;
        this.taxonId = taxon;
    }

    public ReadBarcodeResult(String sequence) {
        unprocessedSequence = sequence;
    }

    public byte getLength() {
        return length;
    }

    public BitSet getRead() {
        return read;
    }

    /**Return taxon name implied by the barcode sequence*/
    public String getTaxonName() {
        return taxonName;
    }
}
