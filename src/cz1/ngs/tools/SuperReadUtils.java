package cz1.ngs.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import cz1.ngs.model.Sequence;
import cz1.util.ArgsEngine;
import cz1.util.Executor;
import cz1.util.Utils;

public class SuperReadUtils extends Executor {

	private static enum Task {create, clean, overlap, zzz}
	private Task task_list = Task.zzz;
	
	private String out_prefix = null;
	private String utg_file = null;
	private String sr_file  = null;
	private String work_dir = null;
	private int _k_ = 0;
	
	@Override
	public void printUsage() {
		// TODO Auto-generated method stub
		switch(this.task_list) {
		case zzz:
			myLogger.info(
					"\n\nChoose a task from listed:\n"
							+ " create                  Create super reads sequence file and check. \n"
							+ " clean                   Clean redundant super reads. \n"
							+ " overlap                 Make overlap graph. \n"
							+ "\n");
			break;
		case create:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -wd/--work-dir          Working directory.\n"
							+ " -k/--kmer-size          K-mer size. \n"
							+ " -u/--unitig-file        Unitig file in FASTA format. \n"
							+ " -o/--out-prefix         Prefix of the output files. \n"
							+ "\n");
			break;
		case clean:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -s/--sequence-file      Sequence file in FASTA format. \n"
							+ " -o/--out-prefix         Prefix of the output files. \n"
							+ "\n");
			break;
		case overlap:
			myLogger.info(
					"\n\nUsage is as follows:\n"
							+ " -s/--sequence-file      Sequence file in FASTA format. \n"
							+ " -u/--unitig-file        Unitig file in FASTA format. \n"
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
		case "create":
			this.task_list = Task.create;
			break;
		case "clean":
			this.task_list = Task.clean;
			break;
		case "overlap":
			this.task_list = Task.overlap;
			break;
		default:
			printUsage();
			throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");	
		}

		String[] args2 = new String[args.length-1];
		System.arraycopy(args, 1, args2, 0, args2.length);
		
		if (myArgsEngine == null) {
			myArgsEngine = new ArgsEngine();
			myArgsEngine.add("-wd", "--work-dir", true);
			myArgsEngine.add("-s", "--sequence-file", true);
			myArgsEngine.add("-u", "--unitig-file", true);
			myArgsEngine.add("-k", "--kmer-size", true);
			myArgsEngine.add("-o", "--out-prefix", true);
			myArgsEngine.parse(args2);
		}
		
		if (myArgsEngine.getBoolean("-o")) {
			this.out_prefix = myArgsEngine.getString("-o");
		} else {
			printUsage();
			throw new IllegalArgumentException("Please specify the prefix of output files.");
		}
		
		switch(this.task_list) {
		case create:
			if (myArgsEngine.getBoolean("-wd")) {
				this.work_dir = myArgsEngine.getString("-wd");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the working directory.");
			}
			if (myArgsEngine.getBoolean("-k")) {
				this._k_ = Integer.parseInt(myArgsEngine.getString("-k"));
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the kmer size.");
			}
			if (myArgsEngine.getBoolean("-u")) {
				this.utg_file = myArgsEngine.getString("-u");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the unitig sequence file.");
			}
			break;
		case clean:
			if (myArgsEngine.getBoolean("-s")) {
				this.sr_file = myArgsEngine.getString("-s");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the super read sequence file.");
			}
			break;
		case overlap:
			if (myArgsEngine.getBoolean("-s")) {
				this.sr_file = myArgsEngine.getString("-s");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the super read sequence file.");
			}
			if (myArgsEngine.getBoolean("-u")) {
				this.utg_file = myArgsEngine.getString("-u");
			} else {
				printUsage();
				throw new IllegalArgumentException("Please specify the unitig sequence file.");
			}
			break;
		default:
			throw new RuntimeException("!!!");	
		}
	}

	@Override
	public void run() {
		// TODO Auto-generated method stub
		switch(this.task_list) {
		case zzz:
			myLogger.info("Task list is empty!!!");
			break;
		case create:
			this.create();
			break;
		case clean:
			this.clean();
			break;
		case overlap:
			this.overlap();
			break;
		default:
			throw new RuntimeException("!!!");
		}
		return;
	}

	private void create() {
		// TODO Auto-generated method stub
		final Map<String, Sequence> utgs = Sequence.parseFastaFileWithRevCmpAsMap(utg_file);
		// sr all sequence file    : $wd/superReadSequences.fasta.all
		// sr all name file        : $wd/superReadNames.txt
		// sr reduce sequence file : $wd/superReadSequences.fasta
		// sr reduce name file     : $wd/reduce.tmp.renamed
		
		// 1. check if seq names are correct
		myLogger.info("#1.check sequence...");
		final Map<String, String> names = new HashMap<>();
		try {
			BufferedReader br_alls = Utils.getBufferedReader(this.work_dir+"/superReadSequences.fasta.all");
			BufferedReader br_alln = Utils.getBufferedReader(this.work_dir+"/superReadNames.txt");
			String line, seqid, utgn;
			int seqn = 0;
			while( (line=br_alls.readLine())!=null ) {
				++seqn;
				if(seqn%100000==0) myLogger.info("#1.sequence check done: "+seqn);
				if(!line.startsWith(">")) throw new RuntimeException("!!!");
				seqid = line.trim().substring(1);
				utgn = br_alln.readLine();
				if(!br_alls.readLine().equals(getSequence(utgs,utgn)))
					throw new RuntimeException("!!!");
				names.put(seqid, utgn);
			}
			br_alls.close();
			br_alln.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// 1 check done
		myLogger.info("#1.sequence check done: ALL");
		
		// 2. check if reduced seq names are correct
		myLogger.info("#2.check sequence reduced...");
		List<String> finals = new ArrayList<>(); 
		try {
			BufferedReader br_redn = Utils.getBufferedReader(this.work_dir+"/reduce.tmp.renamed");
			String line, seqid;
			Set<String> redn = new HashSet<>();
			while( (line=br_redn.readLine())!=null )
				redn.add(line.split("\\s+")[0]);
			br_redn.close();
			Set<String> alln = new HashSet<>();
			alln.addAll(names.keySet());
			alln.removeAll(redn);
			
			BufferedReader br_reds = Utils.getBufferedReader(this.work_dir+"/superReadSequences.fasta");
			int seqn = 0;
			while( (line=br_reds.readLine())!=null ) {
				++seqn;
				if(seqn%100000==0) myLogger.info("#2.sequence reduced check done: "+seqn);
				if(!line.startsWith(">")) throw new RuntimeException("!!!");
				seqid = line.trim().substring(1);
				if(!br_reds.readLine().equals(getSequence(utgs,names.get(seqid))))
					throw new RuntimeException("!!!");
				finals.add(seqid);
			}
			br_reds.close();
			
			if(alln.size()!=finals.size()||!alln.containsAll(finals))
				throw new RuntimeException("!!!");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		myLogger.info("#2.sequence reduced check done: ALL");
		
		// 3. writting super reads sequence file
		myLogger.info("#3.writting sr sequences...");
		try {
			BufferedWriter bw_finals = Utils.getBufferedWriter(this.out_prefix+".fasta");
			int seqn = 0;
			for(String seqid : finals) {
				++seqn;
				if(seqn%100000==0) myLogger.info("#3.writting sr sequences done: "+seqn);
				bw_finals.write(">"+names.get(seqid)+"\n");
				//bw_finals.write(">"+seqid+" "+names.get(seqid)+"\n");
				bw_finals.write(getSequence(utgs,names.get(seqid))+"\n");
			}
			bw_finals.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		myLogger.info("#3.sr sequence writting done.");
	}

	private String getSequence(Map<String, Sequence> utgs, String utgn) {
		// TODO Auto-generated method stub
		String[] uns = utgn.split("_");
		String un;
		for(int i=0; i<uns.length; i++) {
			un = uns[i];
			if(un.endsWith("F")) 
				uns[i] = un.substring(0, un.length()-1);
			else if(un.endsWith("R"))
				uns[i] = un.substring(0, un.length()-1)+"'";
			else throw new RuntimeException("!!!");
		}
		StringBuilder seq = new StringBuilder(utgs.get(uns[0]).seq_str());
		for(int i=1; i<uns.length; i++)
			seq.append(utgs.get(uns[i]).seq_str().substring(_k_-1));
		return seq.toString();
	}

	final Object lock = new Object();
	
	private void overlap() {
		// TODO Auto-generated method stub
		Map<String, Integer> utgs = new HashMap<>();
		List<String> sr = new ArrayList<>();
		try {
			BufferedReader br_utg = Utils.getBufferedReader(this.utg_file);
			String line, utg;
			int index = 0;
			while( (line=br_utg.readLine())!=null ) {
				if(line.startsWith(">")) {
					utg = line.split("\\s+")[0].substring(1);
					utgs.put(utg+"F", index++);
					utgs.put(utg+"R", index++);
				}
			}
			br_utg.close();
			
			BufferedReader br_sr = Utils.getBufferedReader(this.sr_file);
			while( (line=br_sr.readLine())!=null ) {
				if(line.startsWith(">")) {
					sr.add(line.split("\\s+")[0].substring(1));
				}
			}
			br_sr.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		final int nseq = sr.size();
		final int[][] utgmat= new int[nseq][];
		String[] s;
		for(int i=0; i<nseq; i++) {
			s = sr.get(i).split("_");
			int[] u = new int[s.length];
			for(int j=0; j<u.length; j++) 
				u[j] = utgs.get(s[j]);
			utgmat[i] = u;
		}
		
		Map<Integer, Set<Integer>> utg_coords = new HashMap<>();
		int[] u;
		for(int i=0; i<nseq; i++) {
			u = utgmat[i];
			for(int j : u) {
				if(!utg_coords.containsKey(j))
					utg_coords.put(j, new HashSet<>());
				utg_coords.get(j).add(i);
			}
		}
		
		// get all potential pairs
		final Set<Long> pairs = new HashSet<>();
		long key;
		for(Set<Integer> v : utg_coords.values()) {
			for(int i : v) {
				for(int j : v) {
					if(i<j) {
						key =  i;
						key <<= 32;
						key += j;
						pairs.add(key);
					}
				}
			}
		}
		
		// process each pair
		myLogger.info("#sequence pairs to process: "+pairs.size());
		final long[] counter = new long[10];
		this.initial_thread_pool();
		for(long k : pairs) {
			this.executor.submit(new Runnable() {
				private long pair;
				
				@Override
				public void run() {
					// TODO Auto-generated method stub
					
					
					synchronized(lock) {
						++counter[0];
						if(counter[0]%1000000==0)
							myLogger.info("#sequence pairs processed: "+counter[0]);
					}
				}

				public Runnable init(long pair) {
					// TODO Auto-generated method stub
					this.pair = pair;
					return this;
				}

			}.init(k));
		}
		this.waitFor();
	}

	private void clean() {
		// TODO Auto-generated method stub
		final List<Sequence> seqs = Sequence.parseFastaFileAsList(sr_file);
		
		final List<String> seqids = new ArrayList<>();
		for(Sequence seq : seqs) {
			seqids.add(seq.seq_sn());
			seqids.add(rev_sr(seq.seq_sn()));
		}
		
		final int nseq = seqids.size()/2;
		final int[] nlens = new int[nseq];
		for(int i=0; i<nseq; i++) 
			nlens[i]  = seqids.get(i*2).length();
		
		Map<String, Set<Integer>> utg_coords = new HashMap<>();
		String[] utgs;
		String seqid;
		
		myLogger.info("#1.sequences to process: "+nseq);
		for(int i=0; i<nseq; i++) {
			if(i%10000==0) myLogger.info("#1.sequences processed: "+i);
			seqid = seqids.get(i*2);
			utgs = seqid.split("_");
			for(String utg : utgs) {
				if(!utg_coords.containsKey(utg)) 
					utg_coords.put(utg, new HashSet<>());
				utg_coords.get(utg).add(i);
			}
		}
		
		Set<Integer> sr = new HashSet<>();
		Set<Integer> reduns = new HashSet<>();
		
		myLogger.info("#2.sequences to process: "+nseq);
		outerloop:
			for(int i=0; i<nseq; i++) {
				if(i%10000==0) myLogger.info("#2.sequences processed: "+i+
						"; sequences contained: "+reduns.size());
				for(int j=0; j<2; j++) {
					seqid = seqids.get(i*2+j);
					utgs = seqid.split("_");
					sr.clear();
					for(String utg : utgs) {
						if(!utg_coords.containsKey(utg)) 
							continue;
						for(int u : utg_coords.get(utg)) {
							if(i!=u && 
									nlens[i]<=nlens[u] && 
									!reduns.contains(u))
								sr.add(u*2);
						}
					}
					for(int s : sr) {
						if(seqids.get(s).contains(seqid)) {
							reduns.add(i);
							myLogger.error("####"+(j==0?"":"reversed ")+"sequence "+
									seqid+" is contained in "+seqids.get(s));
							continue outerloop;
						}
					}
				}
			}
		
		myLogger.info("#sequences cleaned: "+reduns.size());
		
		try {
			BufferedWriter bw_fa = Utils.getBufferedWriter(this.out_prefix+".fa");
			for(int i=0; i<nseq; i++) {
				if(reduns.contains(i)) continue;
				bw_fa.write(seqs.get(i).formatOutput(Integer.MAX_VALUE));
			}
			bw_fa.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private String rev_sr(String seqid) {
		// TODO Auto-generated method stub
		String[] utgs = seqid.split("_");
		StringBuilder rev = new StringBuilder();
		for(int i=utgs.length-1; i>=0; i--) {
			rev.append(rev_utg(utgs[i]));
		}
		return rev.toString();
	}
	
	private String rev_utg(String utgid) {
		// TODO Auto-generated method stub
		int n = utgid.length();
		int c = utgid.charAt(n-1);
		if(c!='F'&&c!='R') throw new RuntimeException("!!!");
		return utgid.substring(0, n-1)+(c=='F'?'R':'F');
	}
}


