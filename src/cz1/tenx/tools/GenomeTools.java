package cz1.tenx.tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMFileReader.ValidationStringency;

public class GenomeTools extends Executor {

	private String in_bam;
	private String out_bed;
	private int match_thres = 100;
	private int insert_thres = 10;
	private int delete_thres = 10;
	private int clip_thres = 10;
	private int min_size = 200;
	private boolean skip_primary = false;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--input-bam      Input directory containing tag sequence files.\n"
						+ " -0/--match          Match threshold (default is 100).\n"
						+ " -1/--insert         Insert threshold (default is 10).\n"
						+ " -2/--delete         Delete threshold (default is 10).\n"
						+ " -3/--clip           Clip threshold (default is 10).\n"
						+ " -l/--min-size       Min region size to output (default is 200). \n"
						+ " -z/--skip-primary   Skip the primary alignment (default false). \n"
						+ " -o/--output-bed     Output directory.\n\n");
	}

	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		if (args.length == 0) {
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
		}

		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-i", "--input-bam", true);
			myArgsEngine.add("-0", "--match", true);
			myArgsEngine.add("-1", "--insert", true);
			myArgsEngine.add("-2", "--delete", true);
			myArgsEngine.add("-3", "--clip", true);
			myArgsEngine.add("-l", "--min-size", true);
			myArgsEngine.add("-z", "--skip-primary", false);
			myArgsEngine.add("-o", "--output-bam", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-i")) {
			in_bam = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your FASTQ files.");
		}
		
		if (myArgsEngine.getBoolean("-o")) {
			out_bed = myArgsEngine.getString("-o");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your FASTQ files.");
		}
		
		if (myArgsEngine.getBoolean("-0")) {
			match_thres = Integer.parseInt(myArgsEngine.getString("-0"));
		}
		
		if (myArgsEngine.getBoolean("-1")) {
			insert_thres = Integer.parseInt(myArgsEngine.getString("-1"));
		}
		
		if (myArgsEngine.getBoolean("-2")) {
			delete_thres = Integer.parseInt(myArgsEngine.getString("-2"));
		}
		
		if (myArgsEngine.getBoolean("-3")) {
			clip_thres = Integer.parseInt(myArgsEngine.getString("-3"));
		}
		
		if (myArgsEngine.getBoolean("-l")) {
			min_size = Integer.parseInt(myArgsEngine.getString("-l"));
		}
		
		if (myArgsEngine.getBoolean("-z")) {
			skip_primary = true;
		}
	}
	
	@Override
	public void run() {
		// TODO Auto-generated method stub

		final SAMFileReader inputSam = new SAMFileReader(new File(in_bam));
		inputSam.setValidationStringency(ValidationStringency.SILENT);
		BufferedWriter bed_bw = Utils.getBufferedWriter(out_bed);
		
		
		SAMRecordIterator iter=inputSam.iterator();
		Cigar cigar;
		List<CigarElement> cigar_ele;
		CigarElement e;
		String bed_region;
		List<SAMRecord> records = new ArrayList<SAMRecord>();
		String[] ss;
		
		int seq_len, nc, nr, pos_cursor;
		List<List<Integer>> deletion_pos = new ArrayList<List<Integer>>();
		SAMRecord record = iter.next(), record_i;
		
		try {
			while(true) {
				records.clear();
				records.add(record);
				bed_region = record.getReadName();
				
				while( iter.hasNext() && 
						(record=iter.next()).getReadName().equals(bed_region) )
					records.add(record);
				nr = records.size();
				
				seq_len = records.get(0).getReadLength();
				boolean[][] house_keeper = new boolean[nr][seq_len];
				deletion_pos.clear();
				for(int i=0; i<nr; i++) deletion_pos.add(new ArrayList<Integer>());
				
				for(int i=0; i<nr; i++) {
					record_i = records.get(i);
					
					if( record_i.getReadUnmappedFlag() ||
							skip_primary && 
							!record_i.getSupplementaryAlignmentFlag() &&
							!record_i.getNotPrimaryAlignmentFlag() )
						continue;
					
					cigar = record_i.getCigar();
					cigar_ele = new ArrayList<CigarElement>(cigar.getCigarElements());
					nc = cigar_ele.size();
					pos_cursor = 0;
					if(record_i.getReadNegativeStrandFlag())
						Collections.reverse(cigar_ele);
					
					for(int j=0; j<nc; j++) {
						e = cigar_ele.get(j);
						switch( CigarOperator.enumToCharacter(e.getOperator()) ) {
						case 'H':
						case 'S':
							if(e.getLength()<=clip_thres)
								Arrays.fill(house_keeper[i], pos_cursor, 
										pos_cursor+e.getLength(), true);
							pos_cursor += e.getLength();
							break;
						case 'M':
							Arrays.fill(house_keeper[i], pos_cursor, 
									pos_cursor+e.getLength(), true);
							pos_cursor += e.getLength();
							break;
						case 'D':
							if(e.getLength()>delete_thres)
								deletion_pos.get(i).add(pos_cursor);
							break;
						case 'I':
							if(e.getLength()<=insert_thres)
								Arrays.fill(house_keeper[i], pos_cursor, 
										pos_cursor+e.getLength(), true);
							pos_cursor += e.getLength();
							break;
						default:
						}
					}
				}
				
				boolean[] out = new boolean[seq_len], house_keeper_i;
				List<Integer> deletion_pos_i;
				for(int i=0; i<nr; i++) {
					house_keeper_i = house_keeper[i];
					deletion_pos_i = deletion_pos.get(i);
					deletion_pos_i.add(0, 0);
					deletion_pos_i.add(seq_len);
					for(int j=0; j<deletion_pos_i.size()-1; j++)
						render(out, 
								house_keeper_i, 
								deletion_pos_i.get(j), 
								deletion_pos_i.get(j+1));
				}
				
				ss = bed_region.split(":|-");
				int piv = Integer.parseInt(ss[3]);
				int ps = 0, pe = 0;
				List<String> split_regions = new ArrayList<String>();
				for(int i=0; i!=seq_len; i++) {
					if(out[i]) {
						if(pe-ps>=min_size) {
							bed_bw.write(ss[2]+"\t"+(piv+ps)+"\t"+(piv+pe)+"\n");
							split_regions.add("::"+ss[2]+":"+(piv+ps)+"-"+(piv+pe));
						}
						ps = i+1;
						pe = i+1;
					} else {
						pe = i+1;
					}
				}
				if(pe-ps>=min_size) {
					bed_bw.write(ss[2]+"\t"+(piv+ps)+"\t"+(piv+pe)+"\n");
					split_regions.add("::"+ss[2]+":"+(piv+ps)+"-"+(piv+pe));
				}
				if(split_regions.size()==0) {
					myLogger.info(bed_region+" -> deleted.");
				} else if(!bed_region.equals(split_regions.get(0))) {
					String oo = split_regions.get(0);
					for(int i=1; i!=split_regions.size(); i++)
						oo += ", "+split_regions.get(i);
					myLogger.info(bed_region+" -> "+oo);
				}
				
				if(record.getReadName().equals(bed_region)) break;
			}
			

			bed_bw.close();
			iter.close();
			inputSam.close();
			
		} catch (IOException ex) {
			ex.printStackTrace();
			System.exit(1);
		}
	}

	private void render(final boolean[] out, final boolean[] house_keeper, final int a, final int b) {
		// TODO Auto-generated method stub
		int s = -1, e = a;
		for(int i=a; i!=b; i++) {
			if(house_keeper[i]) {
				if(s==-1) {
					s = a;
					e = a+1;
				} else e++;
			} else {
				if(s==0&&s<e || 
						e==house_keeper.length&&s<e || 
						s>=a&&(e-s)>match_thres) 
					Arrays.fill(out, s, e, true);
				s = i+1;
				e = i+1;
			}
		}
		if(s==0&&s<e ||
				e==house_keeper.length&&s<e || 
				s>=a&&(e-s)>match_thres) 
			Arrays.fill(out, s, e, true);
	}
	
	public static void main(String[] args) {
		
		System.out.println("a".equals(null)+"xxx");
		
		GenomeTools genomeTools = new GenomeTools();
		genomeTools.setParameters(new String[]{"-i","null","-o","null","-0","1","-l","1"});
		boolean[] out = new boolean[12];
		boolean[] house_keeper = new boolean[]{true, true, true, true, false, true, true, true, false, true, false, false};
		genomeTools.render(out, house_keeper, 0, 1);
		genomeTools.render(out, house_keeper, 1, 7);
		genomeTools.render(out, house_keeper, 7, 12);
		for(int i=0; i!=12; i++) System.out.println(out[i]);
		
		String[] ss = "::Chr01:2263357-2263844".split(":|-");
		int piv = Integer.parseInt(ss[3]);
		int ps = 0, pe = 0;
		for(int i=0; i!=12; i++) {
			if(out[i]) {
				if(pe-ps>=genomeTools.min_size)
					System.out.println(ss[2]+"\t"+(piv+ps)+"\t"+(piv+pe));
				ps = i+1;
				pe = i+1;
			} else {
				pe = i+1;
			}
		}
		if(pe-ps>=genomeTools.min_size) System.out.println(ss[2]+"\t"+(piv+ps)+"\t"+(piv+pe));
	}
}
