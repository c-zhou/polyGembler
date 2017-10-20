package cz1.tenx.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.zip.GZIPOutputStream;

import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;

import org.apache.commons.lang3.StringUtils;
import org.renjin.eval.Context;
import org.renjin.primitives.io.serialization.RDataWriter;
import org.renjin.primitives.matrix.DoubleMatrixBuilder;
import org.renjin.sexp.ListVector;

import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;

public class TenXMoleculeStats extends Executor {

	private static final boolean debug = false;
	
	private String in_bam;
	private String out_stats;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " -i/--input-bam      Input directory containing tag sequence files.\n"
						+ " -t/--threads        Threads (default is 1).\n"
						+ " -o/--output-stats   Output file.\n\n");
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
			myArgsEngine.add("-t", "--threads", true);
			myArgsEngine.add("-o", "--output-stats", true);
			myArgsEngine.parse(args);
		}

		if (myArgsEngine.getBoolean("-i")) {
			in_bam = myArgsEngine.getString("-i");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your FASTQ files.");
		}

		if (myArgsEngine.getBoolean("-t")) {
			THREADS = Integer.parseInt(myArgsEngine.getString("-t"));
		}
		
		if (myArgsEngine.getBoolean("-o")) {
			out_stats = myArgsEngine.getString("-o");
		}
	}
	
	private static final int max_gap = 20000;
	private static final int oob_thresh = 20000;
	private static BufferedWriter oos;
	private static SAMFileWriter oob;
	private static final Object lock = new Object(),
			lockoos = new Object();
	private static long counter = 0, last_counter = 0;
	
	@Override
	public void run() {
		// TODO Auto-generated method stub

		final SAMFileReader inputSam = new SAMFileReader(new File(in_bam));
		final SAMSequenceDictionary seqDic = inputSam.getFileHeader().getSequenceDictionary();
		inputSam.setValidationStringency(ValidationStringency.SILENT);
		/***
		try {
			oos = Utils.getGZIPBufferedWriter(out_stats);
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		**/
		oos = Utils.getBufferedWriter(out_stats);
		SAMFileHeader header = inputSam.getFileHeader();
		header.setAttribute("SO","unsorted");
		oob = new SAMFileWriterFactory().
				makeSAMOrBAMWriter(header, true, new File(out_stats+".bam"));
		
		List<SAMRecord> bc_records = new ArrayList<SAMRecord>();
		SAMRecordIterator iter = inputSam.iterator();
		String bc_str = null, bc_tmp;
		SAMRecord record = null;
		String read_name;
		
		this.initial_thread_pool();
		while(iter.hasNext()) {
			
			while(bc_str==null && iter.hasNext()) {
				record = iter.next();
				read_name = record.getReadName();
				if(!read_name.startsWith("BX:Z:"))
					continue;
				bc_str = read_name.substring(5, 23);
			}
			if(bc_str==null) break;
			bc_records.add(record);
			
			while(iter.hasNext()) {
				record = iter.next();
				read_name = record.getReadName();
				if(!read_name.startsWith("BX:Z:")) {
					bc_str = null;
					break;
				}
				bc_tmp = read_name.substring(5, 23);
				if(bc_str.equals(bc_tmp))
					bc_records.add(record);
				else {
					bc_str = bc_tmp;
					break;
				}
			}
			
			executor.submit(new Runnable() {
				private List<SAMRecord> bc_records;
				
				private final StringBuilder os = new StringBuilder();
				private String barcode;
				private int chr, p0, p1, read_count;
				private final List<Integer> gap = new ArrayList<Integer>();
				private final List<SAMRecord> os_record = new ArrayList<SAMRecord>();
				private int molecule_id = 0;
				
				@Override
				public void run() {
					// TODO Auto-generated method stub
					try {
						synchronized (lock) {
							counter += bc_records.size();
							if(counter-last_counter>10000) {
								myLogger.info("["+Thread.currentThread().getName()+"] "+counter+" records processed.");
								last_counter = counter;
							}
						}
						
						barcode = bc_records.get(0).getReadName().substring(5, 23);
						int n_bc = bc_records.size();
						SAMRecord record;
						for(int i=0; i!=n_bc;i++) {
							record = bc_records.get(i);
							if(record.getMappingQuality()<20 ||
									record.getDuplicateReadFlag() || 
									record.getNotPrimaryAlignmentFlag() ||
									record.getSupplementaryAlignmentFlag() ) {
								bc_records.remove(i);
								--i;
								--n_bc;
							}
						}
						if(n_bc==0) return;
						Collections.sort(bc_records, new Comparator<SAMRecord>() {
							@Override
							public int compare(SAMRecord r, SAMRecord r2) {
								// TODO Auto-generated method stub
								int f = r.getReferenceIndex() - r2.getReferenceIndex();
								return f==0 ? 
										(r.getAlignmentStart()-r2.getAlignmentStart()) : f;
							}
						});
						
						if(debug) {
							myLogger.info("["+Thread.currentThread().getName()+"] "+barcode+", "+bc_records.size());
							for(int i=0; i!=n_bc; i++) {
								record = bc_records.get(i);
								myLogger.info("    "+record.getReferenceName()+": "+record.getAlignmentStart());
							}
						}
						
						Iterator<SAMRecord> it = bc_records.iterator();
						record = it.next();
						this.initNewMolecule(record);
						int g;
						String read_name, mol_id;
						
						while( it.hasNext() ) {
							record = it.next();
							if(record.getAlignmentStart()-
									p1>max_gap || 
									chr!=record.getReferenceIndex()) {
								this.addNewMoleculeRecord();
								if(p1-p0+1>=oob_thresh) {
									mol_id = StringUtils.leftPad(""+molecule_id, 5, '0');
									synchronized (lockoos) {
										for(SAMRecord rc : os_record) {
											read_name = rc.getReadName();
											rc.setReadName(barcode+":"+mol_id+read_name.substring(23));
											rc.setAttribute("RX", barcode+"-"+mol_id);
											oob.addAlignment(rc);
										}
									}
									++molecule_id;
								}
								this.initNewMolecule(record);
							} else {
								g = record.getAlignmentStart()-p1;
								if(g>0) gap.add(g);
								p1 = Math.max(record.getAlignmentEnd(), p1);
								os_record.add(record);
								++read_count;
							}
						}
						this.addNewMoleculeRecord();
						oos.write(os.toString());
						if(p1-p0+1>=oob_thresh) {
							mol_id = StringUtils.leftPad(""+molecule_id, 5, '0');
							synchronized (lockoos) {
								for(SAMRecord rc : os_record) {
									read_name = rc.getReadName();
									rc.setReadName(barcode+":"+mol_id+read_name.substring(23));
									rc.setAttribute("RX", barcode+"-"+mol_id);
									oob.addAlignment(rc);
								}
							}
							++molecule_id;
						}
					} catch (Exception e) {
						Thread t = Thread.currentThread();
						t.getUncaughtExceptionHandler().uncaughtException(t, e);
						e.printStackTrace();
						executor.shutdown();
						System.exit(1);
					}
				}

				private void initNewMolecule(SAMRecord record) {
					// TODO Auto-generated method stub
					chr = record.getReferenceIndex();
					p0 = record.getAlignmentStart();
					p1 = record.getAlignmentEnd();
					read_count = 1;
					gap.clear();
					os_record.clear();
					os_record.add(record);
				}

				private void addNewMoleculeRecord() throws IOException {
					// TODO Auto-generated method stub
					int maxG, l, gS;
					double cvg;
					l = p1-p0+1;
					if(gap.isEmpty()) {
						maxG = 0;
						gS = 0;
						cvg = 1;
					} else {
						maxG = gap.get(0);
						gS = 0;
						for(Integer i : gap) {
							if(i>maxG) maxG=i;
							gS += i;
						}
						cvg = (l-gS)/(double) l;
					}
					//os.setLength(0);
					os.append(barcode);
					os.append("\t");
					os.append(seqDic.getSequence(chr).getSequenceName());
					os.append(":");
					os.append(p0);
					os.append("-");
					os.append(p1);
					os.append("\t");
					os.append(l);
					os.append("\t");
					os.append(read_count);
					os.append("\t");
					os.append(gS);
					os.append("\t");
					os.append(cvg);
					os.append("\t");
					os.append(maxG);
					os.append("\t");
					if(gap.isEmpty()) os.append("-");
					else {
						os.append(gap.get(0));
						int n = gap.size();
						for(int i=1; i!=n; i++) {
							os.append(",");
							os.append(gap.get(i));
						}
					}
					os.append("\n");
					//oos.write(os.toString());
				}

				public Runnable init(final List<SAMRecord> bc_records) {
					this.bc_records = bc_records;
					return(this);
				}
			}.init(bc_records));
			bc_records = new ArrayList<SAMRecord>();
		}
		iter.close();
		inputSam.close();
		this.waitFor();
		try {
			oos.close();
			oob.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		myLogger.info("Completed. "+counter+" bam records processed.");
	}
	
	private void r() {
		// TODO Auto-generated method stub
		ScriptEngineManager manager = new ScriptEngineManager();
		ScriptEngine engine = manager.getEngineByName("Renjin"); 
		if(engine == null) { 
			throw new RuntimeException("Renjin not found!!!"); 
		}
		try {

			Process p = Runtime.getRuntime().exec("wc -l " + out_stats);
			p.waitFor();
			BufferedReader pbr = new BufferedReader(new InputStreamReader(p.getInputStream()));
			String line;
			int mol_n = 0;
			while ((line = pbr.readLine()) != null) 
				mol_n = Integer.parseInt(line.split("\\s+")[0]);
			pbr.close();
			
			DoubleMatrixBuilder mol_len = new DoubleMatrixBuilder(1, mol_n);
			final Map<String, Integer> bc_stats = new HashMap<String, Integer>();
			BufferedReader br = Utils.getBufferedReader(out_stats);
			String[] s;
			int i = 0;
			while( (line=br.readLine())!=null ) {
				s = line.split("\\s+");
				if(bc_stats.containsKey(s[0])) 
					bc_stats.put(s[0], bc_stats.get(s[0])+1);
				else bc_stats.put(s[0], 1);
				mol_len.set(0, i++, Integer.parseInt(s[2]));
			}
			br.close();
			
			int bc_n = bc_stats.size();
			DoubleMatrixBuilder bc_nstats = new DoubleMatrixBuilder(1, bc_n);
			i = 0;
			for(Map.Entry<String, Integer> entry : bc_stats.entrySet()) 
				bc_nstats.set(0, i++, entry.getValue());
			
			Context context = Context.newTopLevelContext();
			FileOutputStream fos = new FileOutputStream(out_stats+".RData");
			GZIPOutputStream zos = new GZIPOutputStream(fos);
			RDataWriter writer = new RDataWriter(context, zos);
			
			ListVector.NamedBuilder Rdat = new ListVector.NamedBuilder();
			Rdat.add("mol_stats", mol_len.build());
			Rdat.add("bc_stats", bc_nstats.build());
			writer.save(Rdat.build());
			writer.close();
			
		} catch (IOException | InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private void hist() {
		// TODO Auto-generated method stub
		
		try {
			BufferedReader br = Utils.getBufferedReader(out_stats);
			String line;
			String[] s;
			final Map<Integer, Integer> mol_hist = new HashMap<Integer, Integer>();
			final Map<Integer, Long> mol_L = new HashMap<Integer, Long>();
			final Map<String, Integer> bc_hist = new HashMap<String, Integer>();
			mol_hist.put(100,0);
			mol_hist.put(200,0);
			mol_hist.put(500,0);
			mol_hist.put(1000,0);
			mol_hist.put(2000,0);
			mol_hist.put(5000,0);
			mol_hist.put(10000,0);
			mol_hist.put(20000,0);
			mol_hist.put(50000,0);
			mol_hist.put(100000,0);
			mol_hist.put(Integer.MAX_VALUE,0);
			mol_L.put(100,0L);
			mol_L.put(200,0L);
			mol_L.put(500,0L);
			mol_L.put(1000,0L);
			mol_L.put(2000,0L);
			mol_L.put(5000,0L);
			mol_L.put(10000,0L);
			mol_L.put(20000,0L);
			mol_L.put(50000,0L);
			mol_L.put(100000,0L);
			mol_L.put(Integer.MAX_VALUE,0L);
			int mol_len;
			counter = 0;
			while( (line=br.readLine())!=null ) {
				if(++counter%1000000==0) myLogger.info(counter+" items processed.");;
				s = line.split("\\s+");
				if(bc_hist.containsKey(s[0])) 
					bc_hist.put(s[0], bc_hist.get(s[0])+1);
				else bc_hist.put(s[0], 1);
				if(s[1].startsWith("Chr00")) continue;
				mol_len = Integer.parseInt(s[2]);
				if(mol_len<100) {
					mol_hist.put(100, mol_hist.get(100)+1);
					mol_L.put(100, mol_L.get(100)+mol_len);
				} else if(mol_len<200) {
					mol_hist.put(200, mol_hist.get(200)+1);
					mol_L.put(200, mol_L.get(200)+mol_len);
				} else if(mol_len<500) {
					mol_hist.put(500, mol_hist.get(500)+1);
					mol_L.put(500, mol_L.get(500)+mol_len);
				} else if(mol_len<1000) {
					mol_hist.put(1000, mol_hist.get(1000)+1);
					mol_L.put(1000, mol_L.get(1000)+mol_len);
				} else if(mol_len<2000) {
					mol_hist.put(2000, mol_hist.get(2000)+1);
					mol_L.put(2000, mol_L.get(2000)+mol_len);
				} else if(mol_len<5000) {
					mol_hist.put(5000, mol_hist.get(5000)+1);
					mol_L.put(5000, mol_L.get(5000)+mol_len);
				} else if(mol_len<10000) {
					mol_hist.put(10000, mol_hist.get(10000)+1);
					mol_L.put(10000, mol_L.get(10000)+mol_len);
				} else if(mol_len<20000) {
					mol_hist.put(20000, mol_hist.get(20000)+1);
					mol_L.put(20000, mol_L.get(20000)+mol_len);
				} else if(mol_len<50000) {
					mol_hist.put(50000, mol_hist.get(50000)+1);
					mol_L.put(50000, mol_L.get(50000)+mol_len);
				} else if(mol_len<100000) {
					mol_hist.put(100000, mol_hist.get(100000)+1);
					mol_L.put(100000, mol_L.get(100000)+mol_len);
				} else {
					mol_hist.put(Integer.MAX_VALUE, mol_hist.get(Integer.MAX_VALUE)+1);
					mol_L.put(Integer.MAX_VALUE, mol_L.get(Integer.MAX_VALUE)+mol_len);
				}
			}
			br.close();
			
			BufferedWriter mol_bw = Utils.getBufferedWriter(out_stats+".molL");
			mol_bw.write("<100bp\t"+mol_hist.get(100)+"\t"+mol_L.get(100)+"bp\n");
			mol_bw.write("<200bp\t"+mol_hist.get(200)+"\t"+mol_L.get(200)+"bp\n");
			mol_bw.write("<500bp\t"+mol_hist.get(500)+"\t"+mol_L.get(500)+"bp\n");
			mol_bw.write("<1Kbp\t"+mol_hist.get(1000)+"\t"+mol_L.get(1000)+"bp\n");
			mol_bw.write("<2Kbp\t"+mol_hist.get(2000)+"\t"+mol_L.get(2000)+"bp\n");
			mol_bw.write("<5Kbp\t"+mol_hist.get(5000)+"\t"+mol_L.get(5000)+"bp\n");
			mol_bw.write("<10Kbp\t"+mol_hist.get(10000)+"\t"+mol_L.get(10000)+"bp\n");
			mol_bw.write("<20Kbp\t"+mol_hist.get(20000)+"\t"+mol_L.get(20000)+"bp\n");
			mol_bw.write("<50Kbp\t"+mol_hist.get(50000)+"\t"+mol_L.get(50000)+"bp\n");
			mol_bw.write("<100Kbp\t"+mol_hist.get(100000)+"\t"+mol_L.get(100000)+"bp\n");
			mol_bw.write(">100Kbp\t"+mol_hist.get(Integer.MAX_VALUE)+"\t"+mol_L.get(Integer.MAX_VALUE)+"bp\n");
			mol_bw.close();
			
			BufferedWriter bc_bw = Utils.getBufferedWriter(out_stats+".bc");
			for(Map.Entry<String, Integer> entry : bc_hist.entrySet()) 
				bc_bw.write(entry.getKey()+"\t"+entry.getValue()+"\n");
			bc_bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		TenXMoleculeStats stats = new TenXMoleculeStats();
		stats.setParameters(new String[]{
				"-i", "C:\\Users\\chenxi.zhou\\Desktop\\10x_igv\\xx.bam",
				"-t", "3",
				"-o", "C:\\Users\\chenxi.zhou\\Desktop\\10x_igv\\out.xx.txt"
		});
		stats.run();
		//stats.r();
		stats.hist();
	}
}
