package cz1.util;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Random;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class Constants {
	public final static double eps = 1e-12;
	public static long seed = System.nanoTime();
	public static Random rand = new Random(seed);
	public static RandomGenerator rg = new Well19937c(seed);
	//public final static String file_sep = System.getProperty("file.separator");
	public final static String file_sep = "/";
	public final static String line_sep = System.getProperty("line.separator");
	public final static String user_dir = System.getProperty("user.dir");
	public final static String user_home = System.getProperty("user.home");
	public final static String user_name = System.getProperty("user.name");
	public final static int MAX_FILE_ID_LENGTH = 128;
	public final static NumberFormat formatter = new DecimalFormat("#0.000000");
	public final static double MIN_EXP_DOUBLE = Math.log(Double.MIN_VALUE);
	public final static double MAX_EXP_DOUBLE = Math.log(Double.MAX_VALUE);
	public final static double log2 = Math.log(2);
	/*** constants for Dirichlet and Beta distribution
	 * e for E-step
	 * m for M-step 
	 */
	public final static String _use_field_="PL";
	
	public final static double _mu_alpha_e=1;
	public final static double _mu_theta_e=1;
	public final static double _mu_J_e=1e5;
	public final static double[] _pseudo_ = new double[]{1E10, 1, 1, 1 ,1};
	/*** background probability of recombination between consecutive bases **/
	public final static double _con_base_r=1e-8;
	public final static double _soften_ = 0.05;
	public final static double _max_initial_seperation = 1e7; 
	/*** ploidy **/
	public static int _ploidy_H=2;
	/*** number of ancestral haplotypes - hidden states **/
	public static int _haplotype_z=4;
	public final static String _comment_syntax ="#";
	public final static String _vcf_format = "GT:AD:DP:GQ:PL";
	public static String _founder_haps = "P1:P2";
	//public final static String _founder_haps = "Tetra_TrifidaP1:Tetra_TrifidaP2";
	//public final static String _founder_haps = "Atlantic:B1829";
	//public final static String _founder_haps = "Trifida_D_P1:Trifida_D_P2";
	public final static double _beta_alpha=0.5;
	public final static double _beta_beta=0.5;
	public final static int _beta_mean=1000000;
	public final static double _max_af=0;
	//public final static double _max_af=0.1;
	
	public final static int _n_samples = 190;
	public final static double _outlier_thresh = 0.05;
	public final static double _drop_thres = 0.5;
	
	public final static double threshMax = 1e-100d;
	public final static double threshMin = 1e-200d;
	//public final static double threshMax = 1e-10d;
	//public final static double threshMin = 1e-20d;
	public final static double logThreshMax = Math.log(threshMax);
	public final static double logThreshMin = Math.log(threshMin);
	
	public final static String scaff_collapsed_str = "____";
	
	public static int iteration = 0;

	private final static Logger logger = 
			LogManager.getLogger(Constants.class.getName());
	
	@SuppressWarnings("unused")
	private static void check() {
		if(_haplotype_z<_ploidy_H)
			logger.error("Haplotypes should NOT "
					+ "be less than ploidy. Program halted.");
		return;
	}

	public static void setRandomGenerator() {
		// TODO Auto-generated method stub
		rand = new Random(seed);
		rg = new Well19937c(seed);
	}
}
