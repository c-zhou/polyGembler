package cz1.test;

import java.util.List;

import org.apache.commons.lang3.StringUtils;

import cz1.util.Combination;
import cz1.util.IO;
import cz1.util.Permutation;

public class PermutationTest {
	public static void main(String[] args) {
		List<List<String>> perms = Permutation.multiPermutation(new String[]{"A","B"},4);
		
		System.err.println(perms.size());
		for(List<String> perm : perms) IO.println(StringUtils.join(perm.toArray(),"_"));
		
		double a = Double.NEGATIVE_INFINITY-Double.NEGATIVE_INFINITY;
		System.out.println( Math.abs((a-a)/a)<1e-4 );
	}
}
