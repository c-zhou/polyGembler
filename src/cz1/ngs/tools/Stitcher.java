package cz1.ngs.tools;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.DiscreteDomain;
import com.google.common.collect.ImmutableRangeSet;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;

import cz1.ngs.model.AlignmentSegment;
import cz1.ngs.model.CompoundAlignmentSegment;
import cz1.ngs.model.Sequence;
import cz1.util.ArgsEngine;
import cz1.util.Constants;
import cz1.util.Executor;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class Stitcher extends Executor {

	private static enum Task {all, parse, stitch, zzz}
	private Task task_list = Task.zzz;

	private String align_file;
	private boolean debug  = false;
	private boolean ddebug = false;
	private String out_prefix = null;
	private String subject_file;
	private String query_file;
	private int minQual = 20;

	private int num_threads = Runtime.getRuntime().availableProcessors();

	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		switch(this.task_list) {
		case zzz:
			myLogger.info(
					"\n\nChoose a task from listed:\n"
							+ " parse                   Count links.\n"
							+ " stitch                  Anchor contigs to generate scaffolds. \n"
							+ " all                     Run parse followed by stitch.\n"
							+ "\n");
			break;

		case parse:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -s/--subject            The FASTA file contain subject/reference sequences. \n"
							+ " -a/--align              The alignment file. \n"
							+ " -t/--threads            Number of threads to use (default: all available cores).\n"
							+ " -d/--debug              Debugging mode will have extra information printed out.\n"
							+ " -dd/--debug-debug       Debugging mode will have more information printed out than -d mode.\n"
							+ " -o/--out-prefix         Prefix of the output files.\n"
							+ "\n");
			break;
		case stitch:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -s/--subject            The FASTA file contain subject/reference sequences. \n"
							+ " -q/--query              The FASTA file contain query/scaffold/contig sequences. \n"
							+ " -a/--align              The alignment file with query sequences aligned to subject sequences.\n"
							+ " -t/--threads            Number of threads to use (default: all available cores).\n"
							+ " -d/--debug              Debugging mode will have extra information printed out.\n"
							+ " -dd/--debug-debug       Debugging mode will have more information printed out than -d mode.\n"
							+ " -o/--out-prefix         Prefix of the output files.\n"
							+ "\n");
			break;
		case all:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -s/--subject            The FASTA file contain subject/reference sequences. \n"
							+ " -q/--query              The FASTA file contain query/scaffold/contig sequences. \n"
							+ " -a/--align              The alignment file with query sequences aligned to subject sequences.\n"
							+ " -t/--threads            Number of threads to use (default: all available cores).\n"
							+ " -d/--debug              Debugging mode will have extra information printed out.\n"
							+ " -dd/--debug-debug       Debugging mode will have more information printed out than -d mode.\n"
							+ " -o/--out-prefix         Prefix of the output files.\n"
							+ "\n");
			break;
		default:
			throw new RuntimeException("!!!");
		}
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}

		switch(args[0].toLowerCase()) {
		case "parse":
			this.task_list = Task.parse;
			break;
		case "stitch":
			this.task_list = Task.stitch;
			break;
		case "all":
			this.task_list = Task.all;
			break;
		default:
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");	
		}

		String[] args2 = new String[args.length-1];
		System.arraycopy(args, 1, args2, 0, args2.length);


		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-a", "--align", true);
			myArgsEngine.add("-s", "--subject", true);
			myArgsEngine.add("-q", "--query", true);
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-d", "--debug", false);
			myArgsEngine.add("-dd", "--debug-debug", false);
			myArgsEngine.add("-o", "--out-prefix", true);
			myArgsEngine.parse(args2);
		}

		switch(this.task_list) {
		case parse:
			break;
		case stitch:
		case all:
			if (myArgsEngine.getBoolean("-q")) {
				this.query_file = myArgsEngine.getString("-q");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the query file.");
			}
			break;
		default:
			throw new RuntimeException("!!!");
		}
		
		if (myArgsEngine.getBoolean("-s")) {
			this.subject_file = myArgsEngine.getString("-s");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the subject file.");
		}
		
		if(myArgsEngine.getBoolean("-a")) {
			this.align_file = myArgsEngine.getString("-a");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the alignment BAM file.");
		}
		
		if (myArgsEngine.getBoolean("-t")) {
			int t = Integer.parseInt(myArgsEngine.getString("-t"));
			if(t<this.num_threads) {
				this.num_threads = t;
			}
		}
		this.THREADS = this.num_threads;
		Constants.omp_threads = this.num_threads;
		
		if (myArgsEngine.getBoolean("-d")) {
			this.debug = true;
		}

		if (myArgsEngine.getBoolean("-dd")) {
			this.debug  = true;
			this.ddebug = true;
		}

		if (myArgsEngine.getBoolean("-o")) {
			this.out_prefix = myArgsEngine.getString("-o");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the prefix of output files.");
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		switch(this.task_list) {
		case parse:
			this.parse();
			break;
		case stitch:
			this.stitch();
			break;
		case all:
			this.parse();
			this.align_file = this.out_prefix+"_parsed.bam";
			this.stitch();
			break;
		default:
			printUsage();
			throw new IllegalArgumentException("\n\nUnkonwn task list.\n\n");	
		}
	}

	private final SamReaderFactory factory =
			SamReaderFactory.makeDefault()
			.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, 
					SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
			.validationStringency(ValidationStringency.SILENT);

	private class BAMBarcodeIterator {

		private final SamReader samReader;
		private final SAMRecordIterator iter;
		private SAMRecord samRecord = null;

		public BAMBarcodeIterator(String bam_file) {
			this.samReader = factory.open(new File(bam_file));
			this.iter = this.samReader.iterator();
			this.samRecord = iter.hasNext() ? iter.next() : null;
		}

		public boolean hasNext() {
			return samRecord != null;
		}

		public List<SAMRecord> next() {

			if(!this.hasNext()) throw new RuntimeException("!!!");

			List<SAMRecord> records = new ArrayList<SAMRecord>();
			String query = samRecord.getReadName();

			while( samRecord!=null && samRecord.getReadName().equals(query) ) {
				if(!samRecord.getReadUnmappedFlag() && 
						!samRecord.isSecondaryAlignment())
					records.add(samRecord);
				samRecord = iter.hasNext()?iter.next():null;
			}

			return records;
		}

		public void close() {
			try {
				this.iter.close();
				this.samReader.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}

	private final static double mxgr = 3.0;
	private final static double diff = 0.1;
	private final static int[] bs = new int[]{10000, 5000, 1000, 500};

	private void stitch() {
		// TODO Auto-generated method stub

		// final Map<String, Sequence> qry_seqs = Sequence.parseFastaFileWithRevCmpAsMap(query_file);
		final Map<String, Sequence>  sub_seqs = Sequence.parseFastaFileAsMap(subject_file);
		
		/***
		final Map<String, ImmutableRangeSet<Integer>> sub_gaps = new HashMap<>();

		for(Map.Entry<String, Sequence> entry : sub_seqs.entrySet()) {
			String seq_sn = entry.getKey();
			String seq_str = entry.getValue().seq_str();

			final RangeSet<Integer> range_set = TreeRangeSet.create();
			char c;
			for(int j=0; j<seq_str.length(); j++) {
				c = seq_str.charAt(j);
				if(c=='N'||c=='n')
					range_set.add(Range.closed(j, j).
							canonical(DiscreteDomain.integers()));
			}
			sub_gaps.put(seq_sn, ImmutableRangeSet.copyOf(range_set));
		}
		**/
		
		final Map<String, List<CompoundAlignmentSegment>> initPlace = new HashMap<>();

		// read alignment file and place the query sequences
		final BAMBarcodeIterator iter = new BAMBarcodeIterator(this.align_file);

		List<SAMRecord> query_recs;
		String qry_id, sub_id;
		final Map<String, List<SAMRecord>> blocks = new HashMap<>();
		List<SAMRecord> block;
		RangeSet<Integer> qryCov, subCov;
		int qstart, qend, sstart, send;
		Range<Integer> span;
		String bx;
		
		while(iter.hasNext()) {
			query_recs = iter.next();
			if(query_recs.isEmpty()) continue;
			
			blocks.clear();
			for(SAMRecord rec : query_recs) {
				bx = rec.getStringAttribute("BX");
				if(!blocks.containsKey(bx)) 
					blocks.put(bx, new ArrayList<>());
				blocks.get(bx).add(rec);
			}
			
			for(String key : blocks.keySet()) {
				block = blocks.get(key);
				qryCov = getQryCovs(block);
				subCov = getSubCovs(block);
				qry_id = block.get(0).getReadName();
				if(reverse(block)) qry_id+="'";
				sub_id = block.get(0).getReferenceName();

				if(!initPlace.containsKey(sub_id))
					initPlace.put(sub_id, new ArrayList<>());
				span = qryCov.span();
				qstart = span.lowerEndpoint();
				qend   = span.upperEndpoint();
				span = subCov.span();
				sstart = span.lowerEndpoint();
				send   = span.upperEndpoint();
				initPlace.get(sub_id).add(new CompoundAlignmentSegment(
						key,      // a universal key for the compount alignment segment
						qry_id,   // query (e.g., gene) sequence id
						sub_id,   // subject (e.g., reference genome) sequence id
						qstart,   // start of alignment in query
						qend,     // end of alignment in query
						sstart,   // start of alignment in subject
						send,     // end of alignment in subject
						subCov,
						qryCov,
						getCovFromCovs(subCov),
						getCovFromCovs(qryCov)));
			}

		}
		iter.close();

		// now process each subject
		this.initial_thread_pool();
		// for(final String sub : initPlace.keySet()) {
		for(final String sub : new String[]{"chr01_itf", "chr02_itf"}) {

			this.executor.submit(new Runnable() {
				private String sub_id;

				@Override
				public void run() {
					// TODO Auto-generated method stub
					try {
						final List<CompoundAlignmentSegment> init = initPlace.get(sub_id);
						final List<CompoundAlignmentSegment> asmb = new ArrayList<>();
						final RangeSet<Integer> molSubCov = TreeRangeSet.create();
						final int subL = sub_seqs.get(this.sub_id).seq_ln();
						final int[] subDepth = new int[subL];
						RangeSet<Integer> cov;
						int lb, ub;
						for(CompoundAlignmentSegment as : init) {
							cov = as.getSubCov();
							for(Range<Integer> r : cov.asRanges()) {
								lb = r.lowerEndpoint()-1;
								ub = r.upperEndpoint()-1;
								for(int i=lb; i<ub; i++)
									++subDepth[i];
							}
						}
						ImmutableRangeSet<Integer> subCov, subExt;

						for(int iter = 0; iter<bs.length; iter++) {
							// block size from large to small
							// extract segments with size not smaller than b
							final int b = bs[iter];
							myLogger.info("Iteration #"+iter+", block size "+b);
							final List<CompoundAlignmentSegment> segs = new ArrayList<>();
							for(int i=0; i<init.size(); i++) {
								CompoundAlignmentSegment seg = init.get(i);
								if(seg!=null&&seg.getQryLen()>=b) {
									segs.add(seg);
									init.set(i, null);
								}
							}
							myLogger.info("  #segments "+segs.size());

							// try to fill gaps on the molecule with these segments
							Collections.sort(segs, new AlignmentSegment.SubjectCoordinationComparator());
							final List<ImmutableRangeSet<Integer>> subExts = new ArrayList<>();
							final RangeSet<Integer> complMolSubCov = molSubCov.complement();
							for(CompoundAlignmentSegment seg : segs) {
								subCov = ImmutableRangeSet.copyOf(seg.getSubCov());
								subExt = subCov.intersection(complMolSubCov);
								subExts.add(subExt);
								myLogger.info(seg.qseqid()+"\t"+subExt.toString());
							}


							myLogger.info("");
						}

					} catch (Exception e) {
						// TODO Auto-generated catch block
						Thread t = Thread.currentThread();
						t.getUncaughtExceptionHandler().uncaughtException(t, e);
						e.printStackTrace();
						executor.shutdown();
						System.exit(1);
					}
				}

				public Runnable init(String sub_id) {
					// TODO Auto-generated method stub
					this.sub_id = sub_id;
					return this;
				}

			}.init(sub));
		}
		this.waitFor();
	}

	private void mergeBlocks(List<List<SAMRecord>> blocks, Map<String, ImmutableRangeSet<Integer>> sub_gaps) {
		// TODO Auto-generated method stub
		List<SAMRecord> prev, next;
		final int n = blocks.size();
		RangeSet<Integer> prevCov, nextCov, subGap;
		Range<Integer> prevSpan, nextSpan;
		int cov, gap;

		outerloop:
		for(int i=n-1; i>0; ) {
			next = blocks.get(i);
			nextCov = getSubCovs(next);
			nextSpan = nextCov.span();

			for(int j=i-1; j>=0; j--) {
				prev = blocks.get(j);
				if(prev.get(0).getReferenceIndex().intValue()!=
						next.get(0).getReferenceIndex().intValue())
					return;
				prevCov = getSubCovs(prev);
				prevSpan = prevCov.span();

				gap = nextSpan.lowerEndpoint()-prevSpan.upperEndpoint();
				if(gap<0) {
					for(int k=i; k>j; k--) {
						prev.addAll(blocks.get(k));
						blocks.remove(k);
					}
					sortBySubPos(prev);
					i = j;
					continue outerloop;
				} else {
					cov = Math.min(getCovFromCovs(prevCov), getCovFromCovs(nextCov));
					subGap = TreeRangeSet.create();
					subGap.add(Range.closed(prevSpan.upperEndpoint(), nextSpan.lowerEndpoint()));
					gap -= getCovFromCovs(sub_gaps.get(prev.get(0).getReferenceName()).intersection(subGap));
					if(gap<=cov*mxgr) {
						for(int k=i; k>j; k--) {
							prev.addAll(blocks.get(k));
							blocks.remove(k);
						}
						sortBySubPos(prev);
						i = j;
						continue outerloop;
					}
				}
			}
			return;
		}
	}
	
	private final Object lock = new Object();
	private long molCount = 0;
	private final int altBlockSize = 1000;

	private void parse() {
		// TODO Auto-generated method stub
				final Map<String, Sequence>  sub_seqs = Sequence.parseFastaFileAsMap(subject_file);
		final Map<String, ImmutableRangeSet<Integer>> sub_gaps = new HashMap<>();

		for(Map.Entry<String, Sequence> entry : sub_seqs.entrySet()) {
			String seq_sn = entry.getKey();
			String seq_str = entry.getValue().seq_str();

			final RangeSet<Integer> range_set = TreeRangeSet.create();
			char c;
			for(int j=0; j<seq_str.length(); j++) {
				c = seq_str.charAt(j);
				if(c=='N'||c=='n')
					range_set.add(Range.closed(j, j).
							canonical(DiscreteDomain.integers()));
			}
			sub_gaps.put(seq_sn, ImmutableRangeSet.copyOf(range_set));
		}

		final BAMBarcodeIterator iter = new BAMBarcodeIterator(this.align_file);

		final SAMFileWriter outputSam =  new SAMFileWriterFactory().
				makeSAMOrBAMWriter(iter.samReader.getFileHeader(), true, 
						new File(this.out_prefix+"_parsed.bam") );

		List<SAMRecord> query_recs;

		this.initial_thread_pool();
		while(iter.hasNext()) {
			query_recs = iter.next();
			if(query_recs.isEmpty()) continue;

			this.executor.submit(new Runnable() {
				private List<SAMRecord> query_recs;

				@Override
				public void run() {
					// TODO Auto-generated method stub

					try {
						final List<List<SAMRecord>> blocks = new ArrayList<>();
						List<SAMRecord> block;
						SAMRecord prev, next;
						RangeSet<Integer> qryCov, subCov, subGap;
						int cov, gap;

						sortBySubPos(query_recs);
						
						blocks.clear();
						next = query_recs.get(0);
						blocks.add(new ArrayList<>());
						blocks.get(blocks.size()-1).add(next);
						subCov = TreeRangeSet.create();

						for(int i=1; i<query_recs.size(); i++) {
							prev = next;
							subCov.add(Range.closed(prev.getAlignmentStart(), prev.getAlignmentEnd()));
							next =  query_recs.get(i);    ;

							if(prev.getReferenceIndex().intValue()==next.getReferenceIndex().intValue()) {
								gap = next.getAlignmentStart()-prev.getAlignmentEnd();
								if(gap>0) {
									cov = Math.min(getCovFromCovs(subCov), next.getAlignmentEnd()-next.getAlignmentStart());
									subGap = TreeRangeSet.create();
									subGap.add(Range.closed(prev.getAlignmentEnd(), next.getAlignmentStart()));
									gap -= getCovFromCovs(sub_gaps.get(prev.getReferenceName()).intersection(subGap));
									if(gap>cov*mxgr) {
										mergeBlocks(blocks, sub_gaps);
										blocks.add(new ArrayList<>());
										subCov.clear();
									}
								}
							} else {
								mergeBlocks(blocks, sub_gaps);
								blocks.add(new ArrayList<>());
								subCov.clear();
							}
							blocks.get(blocks.size()-1).add(next);
						}

						final double[] covs = new double[blocks.size()];
						final List<RangeSet<Integer>> qryCovs = new ArrayList<>();
						double ext = Double.MIN_VALUE;
						for(int i=0; i<blocks.size(); i++) {
							block = blocks.get(i);
							qryCov = getQryCovs(block);
							cov = getCovFromCovs(qryCov);
							covs[i] = cov;
							qryCovs.add(qryCov);
							if(cov>ext) {
								ext = cov;
							}
						}

						ext = Math.min(altBlockSize, ext*(1-diff));
						synchronized(lock) {
							for(int i=0; i<blocks.size(); i++) {
								if(covs[i]<ext) continue;
								block = blocks.get(i);
								++molCount;
								for(SAMRecord rec : block) {
									rec.setAttribute("BX", "b"+String.format("%08d", molCount));
									outputSam.addAlignment(rec);
								}
							}
						}

					} catch (Exception e) {
						// TODO Auto-generated catch block
						Thread t = Thread.currentThread();
						t.getUncaughtExceptionHandler().uncaughtException(t, e);
						e.printStackTrace();
						executor.shutdown();
						System.exit(1);
					}
				}

				public Runnable init(List<SAMRecord> query_recs) {
					// TODO Auto-generated method stub
					this.query_recs = query_recs;
					return this;
				}

			}.init(query_recs));
		}

		this.waitFor();
		iter.close();
		outputSam.close();
	}

	protected void sortBySubPos(List<SAMRecord> recs) {
		// TODO Auto-generated method stub
		Collections.sort(recs, new Comparator<SAMRecord>() {

			@Override
			public int compare(SAMRecord o1, SAMRecord o2) {
				// TODO Auto-generated method stub
				int r1 = o1.getReferenceIndex().intValue(),
						r2 = o2.getReferenceIndex();
				return r1==r2 ? o1.getAlignmentStart()-o2.getAlignmentStart() : r1-r2;
			}

		});

	}

	private boolean reverse(List<SAMRecord> sams) {
		// TODO Auto-generated method stub
		RangeSet<Integer> forward = TreeRangeSet.create();
		RangeSet<Integer> reverse = TreeRangeSet.create();
		for(SAMRecord sam : sams) {
			int offset = getHardClipAt5End(sam);
			int lower = sam.getReadPositionAtReferencePosition(sam.getAlignmentStart(), true);
			int upper = sam.getReadPositionAtReferencePosition(sam.getAlignmentEnd(),   true);
			if(sam.getReadNegativeStrandFlag()) {
				reverse.add(Range.closed(lower+offset, upper+offset).canonical(DiscreteDomain.integers()));
			} else {
				forward.add(Range.closed(lower+offset, upper+offset).canonical(DiscreteDomain.integers()));
			}
		}

		return getCovFromCovs(forward)<getCovFromCovs(reverse);
	}

	private int getHardClipAt5End(SAMRecord sam) {
		// TODO Auto-generated method stub
		CigarElement cigar = sam.getCigar().getCigarElement(0);
		int clip = 0;
		if(cigar.getOperator()==CigarOperator.HARD_CLIP)
			clip += cigar.getLength();
		return clip;
	}

	private int getCovFromCovs(RangeSet<Integer> rs) {
		// TODO Auto-generated method stub
		int cov = 0;
		for(Range<Integer> r : rs.asRanges())
			cov += r.upperEndpoint()-r.lowerEndpoint();
		return cov;
	}

	private int getSpanFromCovs(RangeSet<Integer> rs) {
		// TODO Auto-generated method stub
		Range<Integer> span = rs.span();
		return span.upperEndpoint()-span.lowerEndpoint();
	}

	private RangeSet<Integer> getQryCovs(List<SAMRecord> recs) {
		// TODO Auto-generated method stub
		RangeSet<Integer> covs = TreeRangeSet.create();
		for(SAMRecord rec : recs) {
			int offset = getHardClipAt5End(rec);
			int lower = rec.getReadPositionAtReferencePosition(rec.getAlignmentStart(), true);
			int upper = rec.getReadPositionAtReferencePosition(rec.getAlignmentEnd(),   true);
			covs.add(Range.closed(lower+offset, upper+offset).canonical(DiscreteDomain.integers()));
		}
		return covs;
	}

	private RangeSet<Integer> getSubCovs(List<SAMRecord> recs) {
		// TODO Auto-generated method stub
		RangeSet<Integer> covs = TreeRangeSet.create();
		for(SAMRecord rec : recs) 
			covs.add(Range.closed(rec.getAlignmentStart(), rec.getAlignmentEnd())
					.canonical(DiscreteDomain.integers()));
		return covs;
	}

	private int distance(RangeSet<Integer> span, Range<Integer> ext) {
		// TODO Auto-generated method stub

		Range<Integer> s = span.span();
		int ls = s.lowerEndpoint(), us = s.upperEndpoint(),
				le = ext.lowerEndpoint(), ue = ext.upperEndpoint();
		if(ue<=ls) 
			return ls-ue;
		else if(us<=le)
			return le-us;
		return 0;
	}

}
