package cz1.appl;

import java.io.File;

import cz1.data.DataCollection;
import cz1.tools.PolyGembler;

public class Main {

        public static void main(String[] args) {
                switch(args[0].toLowerCase()) {
                case "datapreparation":
                        String vcf = args[1];
                        String zip = new File(vcf).getName().
                replaceAll(".vcf.gz$", "").
                replaceAll(".vcf$", "");
                        String out_dir = args[2];
                        DataCollection.zip(out_dir, zip, vcf);
                        break;
                case "haplotypephasing":
                        String[] args2 = new String[args.length-1];
                        System.arraycopy(args, 1, args2, 0, args2.length);
                        PolyGembler.main(args2);
                        break;
                default:
                        throw new RuntimeException("Undefined tool!!!");
                }
        }

}
