package cz1.ngs.assembly;

import cz1.util.Executor;

public class SNPeraser extends Executor {

	public static void main() {
		SNPeraser eraser = new SNPeraser();
		eraser.setParameters(new String[] {
				"-t", "3",
				"-b", "C:\\Users\\chenxi.zhou\\Desktop\\snp_error_correctrion\\bt2_alignment\\art_hiseq2000.bam",
				"-v", "C:\\Users\\chenxi.zhou\\Desktop\\snp_error_correctrion\\discoSNP\\unitig_k_31_c_auto_D_100_P_10_b_0_coherent.vcf",
				"-o", "C:\\Users\\chenxi.zhou\\Desktop\\snp_error_correctrion\\art_hiseq2000.fa"
		});
		eraser.run();
	}
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		
	}
}
