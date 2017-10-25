package cz1.breeding.data;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.inference.TestUtils;

import cz1.hmm.data.DataEntry;
import cz1.util.Algebra;
import cz1.util.Combination;
import cz1.util.Constants;

public abstract class Segregation {
	
	protected int ploidy;
	protected double[][] geno_freq_conf = null;
	protected boolean[][] f1_allele_dosage = null;
	protected int[][] parental_allele_dosage = null;
    protected double[] geno_ad_freq = null;
    protected final Map<Integer, Integer> index = 
    		new HashMap<Integer, Integer>();
    
    public Segregation() {
    	this.ploidy = Constants._ploidy_H;
    	this.configure();
    }
    
    public Segregation(int ploidy) {
    	this.ploidy = ploidy;
    	this.configure();
    }
    
    public double[][] calcDosaConfig (DataEntry data, int[] founder_hap) {
    	int no = data.modelLength();
    	double[][] pval = new double[no][geno_freq_conf.length];
		List<String[]> alleles = data.getAllele();
		if(data.getGenotype()!=null) {
			List<List<String[]>> genotypes = data.getGenotype();
			for(int i=0; i!=no; i++) {
				long[] observation = new long[this.ploidy+1];
				List<String[]> genotype = genotypes.get(i);
				String a = alleles.get(i)[1];
				int j = -1;
				for(String[] geno : genotype) {
					j++;
					if(j==founder_hap[0]||j==founder_hap[1])
						continue;
					int iz = baCount(geno, a);
					if(iz!=-1) observation[iz]++;
				}
				
				double[] chisq_p = new double[geno_freq_conf.length];
                for(int z=0; z!=geno_freq_conf.length; z++)
                    chisq_p[z] = TestUtils.chiSquareTest(geno_freq_conf[z], observation);
                pval[i] = chisq_p;
            }
		} else if(data.getGenotypeLikelihood()!=null) {
			throw new RuntimeException("Not implemented yet!!!");
		} else if(data.getAlleleDepth()!=null) {
			throw new RuntimeException("Not implemented yet!!!");
		} else {
			throw new RuntimeException("Truncated data entry!!!");
		}
		return pval;
    }
    
    public double[] calcSegs(DataEntry data, int[] founder_hap) {
		// TODO Auto-generated method stub
		int no = data.modelLength();
		double[] segs = new double[no];
		List<String[]> alleles = data.getAllele();
		if(data.getGenotype()!=null) {
			List<List<String[]>> genotypes = data.getGenotype();
			for(int i=0; i!=no; i++) {
				long[] observation = new long[this.ploidy+1];
				List<String[]> genotype = genotypes.get(i);
				String a = alleles.get(i)[1];
				int j = -1;
				for(String[] geno : genotype) {
					j++;
					if(j==founder_hap[0]||j==founder_hap[1])
						continue;
					int iz = baCount(geno, a);
					if(iz!=-1) observation[iz]++;
				}
				
				double[] chisq_p = new double[geno_freq_conf.length];
                for(int z=0; z!=geno_freq_conf.length; z++)
                    chisq_p[z] = TestUtils.chiSquareTest(geno_freq_conf[z], observation);
                
                /***
                int iz = Algebra.maxIndex(chisq_p);
                if(chisq_p[iz]==0) throw new RuntimeException("!!!");
                int[] dosa = this.parental_allele_dosage[iz];
                segs[i] = Math.abs(ploidy/2-dosa[0])+Math.abs(ploidy/2-dosa[1]);
			    **/

                segs[i] = 1-StatUtils.max(chisq_p);
            }
		} else if(data.getGenotypeLikelihood()!=null) {
			throw new RuntimeException("Not implemented yet!!!");
		} else if(data.getAlleleDepth()!=null) {
			throw new RuntimeException("Not implemented yet!!!");
		} else {
			throw new RuntimeException("Truncated data entry!!!");
		}
		return segs;
	}
    
    private int baCount(String[] genotype, String allele) {
		// TODO Auto-generated method stub
		int ba = 0;
		if(genotype[0].equals(".")) 
			return -1;
		for(String g : genotype) 
			if(g.equals(allele)) ba++;
		return ba;
	}
    
	private void configure() {
		// TODO Auto-generated method stub
		this.configuration(this.ploidy);
		this.geno_ad_freq = new double[this.geno_freq_conf.length];
		for(int j=0; j!=this.geno_freq_conf.length; j++) {
			double af = 0;
			for(int z=0; z<=this.ploidy; z++) 
				af += (this.ploidy-z)*this.geno_freq_conf[j][z];
			this.geno_ad_freq[j] = af/this.ploidy;
		}
	}
	
	private void configuration(int h) {
		
		final List<int[]> tmp_parental_allele_dosage = new ArrayList<int[]>();
        int g = h+1;
        char[][] genotypes = new char[g][h];
        for(int i=0; i<g; i++) {
            Arrays.fill(genotypes[i], 0, h-i, 'A');
            Arrays.fill(genotypes[i], h-i, h, 'B');
        }
        final List<double[]> configs = new ArrayList<double[]>();
        final List<boolean[]> dosa = new ArrayList<boolean[]>();
        for(int i=0; i<g; i++) {
            int[] gameteA = this.gamete(genotypes[i]);
            for(int j=i; j<g; j++) {
                //if( i==0&&j==0 ||
                //  i==0&&j==(g-1)||
                //  i==(g-1)&&j==(g-1) )
                //  //two homozygotes
                //  continue;
                //else {
                int[] gameteB = this.gamete(genotypes[j]);
                double[] os = new double[g];
                for(int k=0; k<gameteA.length; k++)
                    for(int l=0; l<gameteB.length; l++)
                        os[gameteA[k]+gameteB[l]] += 1.0;

                addDosa(dosa, os);
                addConf(configs, 
                		normalize(Algebra.normalize(os),
                        Constants.seq_err), 
                		h-i, 
                		h-j, 
                		tmp_parental_allele_dosage);
                //}
            }
        }
        this.geno_freq_conf = new double[configs.size()][g];
        this.f1_allele_dosage = new boolean[dosa.size()][this.ploidy+1];
        for(int i=0; i<configs.size(); i++) this.geno_freq_conf[i] = configs.get(i);
        for(int i=0; i<dosa.size(); i++) this.f1_allele_dosage[i] = dosa.get(i);
        this.parental_allele_dosage = new int[tmp_parental_allele_dosage.size()][2];
        for(int i=0; i<this.parental_allele_dosage.length; i++) {
        	int[] pa = tmp_parental_allele_dosage.get(i);
        	this.parental_allele_dosage[i] = pa;
        	this.index.put( (pa[0]*this.ploidy+pa[1])*(pa[0]>pa[1]?1:-1), i);
        	this.index.put( (pa[1]*this.ploidy+pa[0])*(pa[1]>pa[0]?1:-1), i);
        }
        return;
    }
    
    private void addDosa(List<boolean[]> dosa,
    		double[] os) {
		// TODO Auto-generated method stub
		boolean[] ds = new boolean[os.length];
		for(int i=0; i<os.length; i++)
			if(os[i]>0) ds[i] = true;
		dosa.add(ds);
		return;
	}

	private double[] normalize(double[] norm_arr, double seqErr) {
        // TODO Auto-generated method stub
        int a = norm_arr.length;
        double[] arr_copy = new double[a];
        System.arraycopy(norm_arr, 0, arr_copy, 0, a);
        for(int i=0; i!=a; i++)
            arr_copy[i] = arr_copy[i]*(1-seqErr)+(1-norm_arr[i])*seqErr/(a-1);
        return arr_copy;
    }

    private void configuration() {
        this.configuration(this.ploidy);
    }

    private void addConf(List<double[]> configs,
            double[] norm, int p1, int p2,
            final List<int[]> parental_allele_dosage) {
        // TODO Auto-generated method stub
        boolean exist = false;
        for(int i=0; i<configs.size(); i++) {
            double[] a = configs.get(i);
            boolean y = true;
            for(int j=0; j<a.length; j++) {
                y = y && a[j]==norm[j];
            }
            exist = y;
        }
        if(!exist) {
            configs.add(norm);
            parental_allele_dosage.add(new int[]{p1, p2});
        }
    }

    private int[] gamete(char[] genotype) {
        // TODO Auto-generated method stub
        List<List<Character>> gametes =
                Combination.combination(
                        ArrayUtils.toObject(genotype),
                        genotype.length/2);
        int[] a = new int[gametes.size()];
        for(int j=0; j<a.length; j++)
            for(int k=0; k<gametes.get(j).size(); k++)
                if(gametes.get(j).get(k)=='A')
                    a[j] += 1.0;
        return a;
    }
    
    public int[][] parentDosa() {
		// TODO Auto-generated method stub
		return this.parental_allele_dosage;
	}

	public boolean[][] f1Dosa() {
		// TODO Auto-generated method stub
		return this.f1_allele_dosage;
	}
	
	public int getIndexFromParentalDosa(int[] parentDosa) {
		// TODO Auto-generated method stub
		return this.index.get((parentDosa[0]*this.ploidy+parentDosa[1])
				*(parentDosa[0]>parentDosa[1]?1:-1));
	}
	
	public int[] getParentalDosaByIndex(final int i) {
		// TODO Auto-generated method stub
		return this.parental_allele_dosage[i];
	}

	public boolean[] getF1DosaByIndex(int i) {
		// TODO Auto-generated method stub
		return this.f1_allele_dosage[i];
	}
}
