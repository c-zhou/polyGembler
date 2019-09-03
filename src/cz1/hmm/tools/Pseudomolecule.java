package cz1.hmm.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import cz1.ngs.model.Sequence;
import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

public class Pseudomolecule extends Executor {
	private String out_file = null;
	private int gap_size = 1000;
	private Map<String, Sequence> sequences;
	private Map<String, List<Sequence>> molecules;
	
	public Pseudomolecule(String map_file,
			String seq_file, 
			String out_file) {
		this.sequences = Sequence.parseFastaFileAsMap(seq_file);
		this.out_file = out_file;
		this.molecules = this.readGeneticMap(map_file);
	}
	
	public Pseudomolecule(String map_file,
			String seq_file, 
			String out_file,
			Map<String, int[][]> errs) {
		this.sequences = Sequence.parseFastaFileAsMap(seq_file);
		this.parseAssemblyErrors(errs);
		this.molecules = this.readGeneticMap(map_file);
		this.out_file = out_file;
	}
	
	public Pseudomolecule(String map_file,
			String seq_file,
			String err_file,
			String out_file) {
		this.sequences = Sequence.parseFastaFileAsMap(seq_file);
		this.parseAssemblyErrors(err_file);
		this.readGeneticMap(map_file);
		this.out_file = out_file;
	}

	public Pseudomolecule() {
		// TODO Auto-generated constructor stub
		super();
	}

	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		myLogger.info(
				"\n\nUsage is as follows:\n"
						+ " Common:\n"
						+ "     -i/--map                    Input genetic linkage map file.\n"
						+ "     -a/--assembly               Input assembly FASTA file.\n"
						+ "     -e/--error                  Assembly error file.\n"
						+ "     -n/--gap                    Gap size between sequences (default 1000). \n"
						+ "                                 The gaps will be filled with character 'n'.\n"
						+ "     -o/--prefix                 Output pseudomolecule file.\n"
						);
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
			myArgsEngine.add("-i", "--map", true);
			myArgsEngine.add("-a", "--assembly", true);
			myArgsEngine.add("-e", "--error", true);
			myArgsEngine.add("-n", "--gap", true);
			myArgsEngine.add("-o", "--prefix", true);
			myArgsEngine.parse(args);
		}
		
		if(myArgsEngine.getBoolean("-a")) {
			this.sequences = Sequence.parseFastaFileAsMap(myArgsEngine.getString("-a"));
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your assembly file in FASTA format.");
		}
		
		if(myArgsEngine.getBoolean("-e")) {
			this.parseAssemblyErrors(myArgsEngine.getString("-e"));
		}
		
		if(myArgsEngine.getBoolean("-i")) {
			this.molecules = this.readGeneticMap(myArgsEngine.getString("-i"));
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify your input genetic linkage map file.");
		}
		
		if(myArgsEngine.getBoolean("-n")) {
			this.gap_size = Integer.parseInt(myArgsEngine.getString("-n"));
		}
		
		if(myArgsEngine.getBoolean("-o")) {
			out_file = myArgsEngine.getString("-o");
		}  else {
			printUsage();
			throw new IllegalArgumentException("Please specify your output file name.");
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		try {
			BufferedWriter bw = Utils.getBufferedWriter(out_file);
			for(String id : molecules.keySet()) {
				StringBuilder oos =  new StringBuilder();
				List<Sequence> seqs = molecules.get(id);
				oos.append(seqs.get(0).seq_str());
				for(int i=1; i<seqs.size(); i++) {
					oos.append(Sequence.polyn(gap_size));
					oos.append(seqs.get(i).seq_str());
				}
				bw.write(Sequence.formatOutput(id, oos.toString()));
			}
			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public Map<String, List<Sequence>> readGeneticMap(String map_file) {
		Map<String, List<Sequence>> molecules = new HashMap<>();
		try {
			BufferedReader br = Utils.getBufferedReader(map_file);
			String line = br.readLine();
			int lg = 0;
			String[] s;
			Sequence sequence;
			String name, seqid;
			while( line!=null ) {
				if(line.startsWith("group")) {
					List<Sequence> seq_list = new ArrayList<>();
					name = "LG"+(++lg);
					while( (line=br.readLine())!=null && 
							!line.startsWith("group") &&
							line.length()!=0) {
						s = line.split("\\(|\\)\\s+");
						seqid = s[0];
						sequence = sequences.get(seqid);
						if(s[1].equals("-"))
							sequence = sequence.revCompSeq();
						seq_list.add(sequence);
						br.readLine();
					}
					molecules.put(name, seq_list);
					if(line!=null&&line.length()!=0)
						line = br.readLine();
				} else line=br.readLine();
			}
			br.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return molecules;
	}
	
	private Map<String, int[][]> parseAssemblyErrorFile(String err_file) {
		// TODO Auto-generated method stub
		Map<String, List<int[]>> errList = new HashMap<>();
		try {
			BufferedReader br_err = Utils.getBufferedReader(err_file);
			String line;
			String[] s;
			while( (line=br_err.readLine())!=null) {
				s = line.split("\\s+");
				errList.putIfAbsent(s[0], new ArrayList<>());
				errList.get(s[0]).add(new int[]{Integer.parseInt(s[1]), Integer.parseInt(s[2])});
			}
			br_err.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		Map<String, int[][]> errs = new HashMap<>();
		for(Map.Entry<String, List<int[]>> entry : errList.entrySet()) {
			List<int[]> err = entry.getValue();
			Collections.sort(err, new Comparator<int[]>() {

				@Override
				public int compare(int[] i0, int[] i1) {
					// TODO Auto-generated method stub
					return Integer.compare(i0[0], i1[0]);
				}
				
			});
		
			int[][] arr = new int[err.size()][];
			for(int i=0; i<arr.length; i++)
				arr[i] = err.get(i);
			errs.put(entry.getKey(), arr);
		}
		
		return errs;
	}
	
	private void parseAssemblyErrors(String err_file) {
		// TODO Auto-generated method stub
		this.parseAssemblyErrors(this.parseAssemblyErrorFile(err_file));
	}
	

	private void parseAssemblyErrors(Map<String, int[][]> errs) {
		// TODO Auto-generated method stub
		for(Map.Entry<String, int[][]> entry : errs.entrySet()) {
			String seqid = entry.getKey();
			if(seqid.equals("contig_516"))
				myLogger.info("");
			String seqstr = sequences.get(seqid).seq_str();
			int[][] bps = findBPS(seqstr, entry.getValue());
			for(int i=0; i<bps.length; i++) {
				String bid = seqid+"_"+(i+1);
				sequences.put(bid, new Sequence(bid, seqstr.substring(bps[i][0], bps[i][1])));
			}
		}
	}
	
	private int[][] findBPS(String seqstr, int[][] err) {
		// TODO Auto-generated method stub
		int[][] chunk = new int[err.length+1][2];
		for(int i=0; i<err.length; i++) {
			String seq = seqstr.substring(err[i][0], err[i][1]);
			String gap = seq.replaceAll("[^ATCGatcg]+", "").replaceAll("[ATCGatcg]+$","");
			if(gap.length()>0&&gap.replace('n', 'N').replaceAll("[N]+","").length()==0) {
				// exactly one gap found, choose as break point
				int star = Math.min(seq.indexOf('n'), seq.indexOf('N'));
				chunk[i][1] = err[i][0]+star;
				chunk[i+1][0] = err[i][0]+gap.length();
			} else {
				// no gap or more than one gap found, chunk in between is discarded
				chunk[i][1] = err[i][0];
				chunk[i+1][0] = err[i][1];
			}
		}
		chunk[err.length][1] = seqstr.length();
		return chunk;
	}
	
	public void setGapSize(int gap_size) {
		this.gap_size = gap_size;
	}
}



