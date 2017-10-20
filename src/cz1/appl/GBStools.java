package cz1.appl;

import java.io.File;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;

import cz1.gbs.tools.FastqToTagSequence;
import cz1.gbs.tools.GBSpileup;
import cz1.gbs.tools.MergeBedGenCov;
import cz1.gbs.tools.MergeTagSequence;
import cz1.gbs.tools.SamFileFilter;
import cz1.gbs.tools.TagSequenceToFastq;
import cz1.hmm.data.DataCollection;
import cz1.hmm.tools.ZipDataCollection;
import cz1.hmm.tools.DataPreparation;
import cz1.hmm.tools.Gembler;
import cz1.hmm.tools.Haplotyper;
import cz1.simulation.tools.GBSSimulator;
import cz1.simulation.tools.PopulationSimulator;

public class GBStools {
	
	protected final static Logger myLogger = 
			Logger.getLogger(GBStools.class);
	
	public static void main(String[] args) {
	
		if(args.length<1) {
			printUsage();
			throw new RuntimeException("Undefined tool!!!");
		}
		String[] args2 = new String[args.length-1];
		System.arraycopy(args, 1, args2, 0, args2.length);
		switch(args[0].toLowerCase()) {
		case "fastq2tag":
			FastqToTagSequence fastq2TagSequence = new FastqToTagSequence();
			fastq2TagSequence.setParameters(args2);
			fastq2TagSequence.run();
			break;
		case "mergetag":
			MergeTagSequence mergeTagSequence = new MergeTagSequence();
			mergeTagSequence.setParameters(args2);
			mergeTagSequence.run();
			break;
		case "tag2fastq":
			TagSequenceToFastq tagSequence2Fastq = new TagSequenceToFastq();
			tagSequence2Fastq.setParameters(args2);
			tagSequence2Fastq.run();
			break;
		case "mergebedgencov":
			MergeBedGenCov mergeBedGenCov = new MergeBedGenCov();
			mergeBedGenCov.setParameters(args2);
			mergeBedGenCov.run();
			break;
		case "samfilefilter":
			SamFileFilter samFileFilter = new SamFileFilter();
			samFileFilter.setParameters(args2);
			samFileFilter.run();
			break;
		default:
			printUsage();
			throw new RuntimeException("Undefined tool!!!");
		}
	}
	
	private static void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " fastq2tag           Convert GBS reads to tag sequences. \n"
						+ " mergetag            Merge tag sequences. \n"
						+ " mergebedgencov      Merge bedtools genome coverage files. \n"
						+ " samfilefilter       Filter abnormally high coverage sites from sam/bam files. \n"
						+ " tag2fastq           Convert tag sequences to GBS reads. \n\n");
	}
}
