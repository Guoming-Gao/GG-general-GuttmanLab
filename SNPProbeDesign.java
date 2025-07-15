package guttmanlab.core.probes;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import guttmanlab.core.annotation.Annotation;
import guttmanlab.core.annotation.Gene;
import guttmanlab.core.annotation.SingleInterval;
import guttmanlab.core.annotation.Annotation.Strand;
import guttmanlab.core.annotation.io.BEDFileIO;
import guttmanlab.core.datastructures.IntervalTree;
import guttmanlab.core.datastructures.IntervalTree.Node;
import guttmanlab.core.math.Statistics;
import guttmanlab.core.sequence.FastaFileIOImpl;
import guttmanlab.core.sequence.Sequence;

public class SNPProbeDesign {
	
	int maxNumberPerIntron;
	String PP1_left="TTTCGTCCGCGAGTGACCAG";
	String PP1_right="CAACGTCCATGTCGGGATGC";
	String PP2_left="TCAGGGCACGAGGACATTCG";
	String PP2_right="TCCGGCAAGATTGCTCTCCC";
	List<String> probes;

	public SNPProbeDesign(File geneFile, File snpFile, int probeLength, String save, File genomeInputPath, int maxNum) throws NumberFormatException, IOException {
		this.maxNumberPerIntron=maxNum;
		Map<String, Sequence> genomeSequenceMap=FastaFileIOImpl.readFromFilesByName(genomeInputPath.listFiles());
		Map<String, IntervalTree<String>> snpTree=parseTree(snpFile);
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(geneFile.getAbsolutePath());
		Collection<Gene> introns=getIntrons(genes);
		genes.addAll(introns);
		
		Collection<Annotation> snps=new TreeSet<Annotation>();
		
		Map<String, Collection<Annotation>> probesByGene=new TreeMap<String, Collection<Annotation>>();
		
		for(Gene gene: genes) {
			String name=gene.getName().split("_")[0];
			if(!probesByGene.containsKey(name)) {probesByGene.put(name, new TreeSet<Annotation>());}
			
			Collection<SingleInterval> snpsByGene=getSNPs(gene, snpTree);
			Collection<Annotation> probesBySNP=getProbes(snpsByGene, probeLength, gene, snpTree.get(gene.getReferenceName()));
			
			Collection<Annotation> set=probesByGene.get(name);
			set.addAll(probesBySNP);
			
			snps.addAll(probesBySNP);
		}
		
		
		//write(probesByGene, save+".bed");
		//writeSequence(probesByGene, save+".fa", genomeSequenceMap, probeLength);
		writeProbes(probesByGene, save, genomeSequenceMap, probeLength);
	}
	
	
	public SNPProbeDesign(File geneFile, File snpFile, int probeLength, File genomeInputPath, int maxNum) throws NumberFormatException, IOException {
		this.maxNumberPerIntron=maxNum;
		Map<String, Sequence> genomeSequenceMap=FastaFileIOImpl.readFromFilesByName(genomeInputPath.listFiles());
		Map<String, IntervalTree<String>> snpTree=parseTree(snpFile);
		Collection<Gene> genes=BEDFileIO.loadRegionsFromFile(geneFile.getAbsolutePath());
		Collection<Gene> introns=getIntrons(genes);
		genes.addAll(introns);
		
		Collection<Annotation> snps=new TreeSet<Annotation>();
		
		Map<String, Collection<Annotation>> probesByGene=new TreeMap<String, Collection<Annotation>>();
		
		for(Gene gene: genes) {
			String name=gene.getName().split("_")[0];
			if(!probesByGene.containsKey(name)) {probesByGene.put(name, new TreeSet<Annotation>());}
			
			Collection<SingleInterval> snpsByGene=getSNPs(gene, snpTree);
			Collection<Annotation> probesBySNP=getProbes(snpsByGene, probeLength, gene, snpTree.get(gene.getReferenceName()));
			
			Collection<Annotation> set=probesByGene.get(name);
			set.addAll(probesBySNP);
			
			snps.addAll(probesBySNP);
		}
		
		
		//write(probesByGene, save+".bed");
		//writeSequence(probesByGene, save+".fa", genomeSequenceMap, probeLength);
		this.probes=writeProbes(probesByGene, genomeSequenceMap, probeLength);
	}
	
	
	private List<String> getProbes(){return this.probes;}
	
	private Collection<Gene> getIntrons(Collection<Gene> genes) {
		Collection<Gene> rtrn=new TreeSet<Gene>();
		
		for(Gene g: genes) {
			Collection<SingleInterval> introns=g.getIntrons();
			for(SingleInterval intron: introns) {
				Gene i=new Gene(intron);
				i.setName("intron:"+g.getName());
				rtrn.add(i);
			}
		}
		
		return rtrn;
	}


	private void writeProbes(Map<String, Collection<Annotation>> probesByGene, String save, Map<String,Sequence> genomeSeq, int probeLength) throws IOException {
		FileWriter writer=new FileWriter(save);
		
		int counter=0;
		
		for(String gene: probesByGene.keySet()) {
			Collection<Annotation> probes=probesByGene.get(gene);
			Map<Annotation, Sequence> probeSeq=filterProbes(probes, genomeSeq, probeLength);
			
			if(gene.contains("intron") && probeSeq.size()>maxNumberPerIntron){
				probeSeq=subsample(probeSeq, maxNumberPerIntron);
			}
			
			for(Annotation probe: probeSeq.keySet()) {
				Gene g=new Gene(probe);
				g.setOrientation(Strand.POSITIVE);
				Sequence seq=probeSeq.get(probe);
				String fullSeq=PP1_left+seq.getSequenceBases()+PP1_right;
				if(gene.contains("intron")) {
					fullSeq=PP2_left+seq.getSequenceBases()+PP2_right;
				}
				if(!rejectProbe(seq.getSequenceBases())) {
					writer.write(gene+"\t"+probe.getSingleInterval().toUCSCStrand()+"\t"+probe.getName()+"\t"+probe.getReferenceName()+"\t"+seq+"\t"+fullSeq+"\n");
					System.out.println(probe.toBED());
				}
			}
			counter++;
			if(counter%100==0) {System.err.println(counter+" "+probesByGene.size());}
		}
		writer.close();
	}
	
	private List<String> writeProbes(Map<String, Collection<Annotation>> probesByGene, Map<String,Sequence> genomeSeq, int probeLength) throws IOException {
		List<String> rtrn=new ArrayList<String>();
		
		int counter=0;
		
		for(String gene: probesByGene.keySet()) {
			Collection<Annotation> probes=probesByGene.get(gene);
			Map<Annotation, Sequence> probeSeq=filterProbes(probes, genomeSeq, probeLength);
			
			if(gene.contains("intron") && probeSeq.size()>maxNumberPerIntron){
				probeSeq=subsample(probeSeq, maxNumberPerIntron);
			}
			
			for(Annotation probe: probeSeq.keySet()) {
				Gene g=new Gene(probe);
				g.setOrientation(Strand.POSITIVE);
				Sequence seq=probeSeq.get(probe);
				String fullSeq=PP1_left+seq.getSequenceBases()+PP1_right;
				if(gene.contains("intron")) {
					fullSeq=PP2_left+seq.getSequenceBases()+PP2_right;
				}
				if(!rejectProbe(seq.getSequenceBases())) {
					rtrn.add(gene+"\t"+probe.getSingleInterval().toUCSCStrand()+"\t"+probe.getName()+"\t"+probe.getReferenceName()+"\t"+seq+"\t"+fullSeq);
				}
			}
			counter++;
			if(counter%100==0) {System.err.println(counter+" "+probesByGene.size());}
		}
		return rtrn;
	}
	
	private static boolean rejectProbe(String probe) {
		if(probe.contains("N")) {return true;}
		LowComplexityFilter lcf=new LowComplexityFilter();
		PolyBaseFilter pbf=new PolyBaseFilter("ACGT", 5, 5); //8,6
		PolyBaseFilter pbf2=new PolyBaseFilter("ACGT", 15, 12);
		PolyBaseFilter pbf3=new PolyBaseFilter("ACGT", 8, 6);
		//RepeatFilter rf=new RepeatFilter(0.07, true, true);
		return lcf.rejectSequence(probe) || pbf.rejectSequence(probe) || pbf2.rejectSequence(probe) || pbf3.rejectSequence(probe);
	}
	
	private Map<Annotation, Sequence> filterProbes(Collection<Annotation> probes, Map<String, Sequence> genomeSeq, int probeLength) {
		Map<Annotation, Sequence> rtrn=new TreeMap<Annotation, Sequence>();
		for(Annotation probe: probes) {
			Gene g=new Gene(probe);
			g.setOrientation(Strand.POSITIVE);
			Sequence seq=g.getSequence(genomeSeq.get(probe.getReferenceName()));
			if(seq.getLength()==probeLength && !seq.hasN()) {
				rtrn.put(probe, seq);
			}
		}
		return rtrn;
	}


	private Map<Annotation, Sequence> subsample(Map<Annotation, Sequence> probes, int num) {
		Map<Annotation, Sequence> rtrn=new TreeMap<Annotation, Sequence>();
		
		List<Annotation> list=new ArrayList<Annotation>();
		list.addAll(probes.keySet());
		
		for(int i=0; i<num; i++) {
			Annotation probe=Statistics.randomSample(list);
			rtrn.put(probe, probes.get(probe));
			list.remove(probe);
		}
		
		
		
		return rtrn;
	}


	private void writeSequence(Map<String, Collection<Annotation>> probesByGene, String save, Map<String,Sequence> genomeSeq, int probeLength) throws IOException {
		FileWriter writer=new FileWriter(save);
		for(String gene: probesByGene.keySet()) {
			for(Annotation probe: probesByGene.get(gene)) {
				Gene g=new Gene(probe);
				g.setOrientation(Strand.POSITIVE);
				Sequence seq=g.getSequence(genomeSeq.get(probe.getReferenceName()));
				if(seq.getLength()==probeLength && !seq.hasN()) {
					writer.write(gene+"\t"+probe.getName()+"\t"+probe.getReferenceName()+"\t"+seq+"\n");
				}
			}
		}
		writer.close();
	}


	private Collection<Annotation> getProbes(Collection<SingleInterval> snpsByGene, int probeLength, Gene gene, IntervalTree<String> snpTree) {
		Collection<Annotation> rtrn=new TreeSet<Annotation>();
		
		for(SingleInterval r: snpsByGene) {
			int add=probeLength/2;
			SingleInterval probe=new SingleInterval(r.getReferenceName(), r.getReferenceStartPosition()-add, r.getReferenceEndPosition()+add);
			SingleInterval upstream=new SingleInterval(r.getReferenceName(), r.getReferenceStartPosition()-probeLength, r.getReferenceEndPosition());
			SingleInterval downstream=new SingleInterval(r.getReferenceName(), r.getReferenceStartPosition(), r.getReferenceEndPosition()+probeLength);
			
			Annotation genomeProbe=convertToGenome(probe, gene);
			Annotation ugenomeProbe=convertToGenome(upstream, gene);
			Annotation dgenomeProbe=convertToGenome(downstream, gene);
			genomeProbe.setName(r.toUCSC()+"_centered");
			ugenomeProbe.setName(r.toUCSC()+"_upstream");
			dgenomeProbe.setName(r.toUCSC()+"_downstream");
			
			if(!overlapsSNP(ugenomeProbe, snpTree)) {rtrn.add(ugenomeProbe);}
			if(!overlapsSNP(dgenomeProbe, snpTree)) {rtrn.add(dgenomeProbe);}
			
			//rtrn.add(genomeProbe); //TODO get sequence and add variable nucleotide
			//rtrn.add(ugenomeProbe); //TODO get sequence and add variable nucleotide
			//rtrn.add(dgenomeProbe); //TODO get sequence and add variable nucleotide
		}
		
		return rtrn;
	}
	
	private boolean overlapsSNP(Annotation ugenomeProbe, IntervalTree<String> snpTree) {
		Iterator<SingleInterval> iter=ugenomeProbe.getBlocks();
		while(iter.hasNext()) {
			SingleInterval block=iter.next();
			if(snpTree.hasOverlappers(block.getReferenceStartPosition(), block.getReferenceEndPosition())) {return true;}
		}
		return false;
	}


	private Annotation convertToGenome(SingleInterval region, Annotation gene) {
		Annotation rtrn=gene.convertToReferenceSpace(region);
		rtrn.setOrientation(Strand.POSITIVE);
		//rtrn.setOrientation(gene.getOrientation());
		//SingleInterval rtrn= gene.convertToReferenceSpace(region).getSingleInterval();
		//rtrn.setOrientation(gene.getOrientation());
		
		return rtrn;
	}
	
	private Collection<SingleInterval> bin(Collection<SingleInterval> snpsByGene, int binSize) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		
		for(SingleInterval r: snpsByGene) {
			rtrn.addAll(r.allBins(binSize));
		}
		
		return rtrn;
	}

	private Collection<SingleInterval> getSNPs(Gene gene, Map<String, IntervalTree<String>> snpTree) {
		Collection<SingleInterval> rtrn=new TreeSet<SingleInterval>();
		if(snpTree.containsKey(gene.getReferenceName())) {
			IntervalTree<String> tree=snpTree.get(gene.getReferenceName());
			
			Iterator<Node<String>> iter=tree.overlappers(gene.getReferenceStartPosition(), gene.getReferenceEndPosition());
			while(iter.hasNext()) {
				Node<String> snp=iter.next();
				SingleInterval r=new SingleInterval(gene.getReferenceName(), snp.getStart(), snp.getEnd());
				r.setOrientation(gene.getOrientation());
				r.setName(snp.getValue());
				
				if(gene.overlapsExon(r)) {
					SingleInterval relative=getRelativeCoordinates(gene, r);
					rtrn.add(relative);
				}
			}
		}
		return rtrn;
	}

	public static SingleInterval getRelativeCoordinates(Annotation gene, SingleInterval feature){
		int featureStart=gene.getRelativePositionFrom5PrimeOfFeature(feature.getReferenceStartPosition());
		int featureEnd=gene.getRelativePositionFrom5PrimeOfFeature(feature.getReferenceEndPosition());
		
		SingleInterval rtrn=new SingleInterval(gene.getName(), featureStart, featureEnd);
		if(gene.getOrientation().equals(Strand.NEGATIVE)) {rtrn=new SingleInterval(gene.getName(), featureEnd, featureStart);}
		
		/*if(gene.getNumberOfBlocks()==1) {
			System.err.println(gene.getName()+" "+gene.toUCSC()+" "+feature.toUCSC());
			System.err.println(gene.getOrientation()+" "+ rtrn.toUCSC());
		}*/
		
		return rtrn;
	}
	
	private void write(Collection<SingleInterval> snpsByGene) {
		for(SingleInterval snp: snpsByGene) {
			System.out.println(snp.toShortBED()+"\t"+snp.getName());
		}
		
	}
	
	
	private void write(Map<String, Collection<Annotation>> snpsByGene, String save) throws IOException {
		FileWriter writer=new FileWriter(save);
		for(String gene: snpsByGene.keySet()) {
			for(Annotation snp: snpsByGene.get(gene)) {
				writer.write(snp.toBED()+"\n");
			}
		}
		writer.close();
	}

	private Map<String, IntervalTree<String>> parseTree(File snpFile) throws NumberFormatException, IOException {
		Map<String, IntervalTree<String>> rtrn=new TreeMap<String, IntervalTree<String>>();
		BufferedReader reader=new BufferedReader(new InputStreamReader(new FileInputStream(snpFile)));
		String nextLine;
		int counter=0;
		while ((nextLine = reader.readLine()) != null) {
			String[] tokens=nextLine.split("\t");
			SingleInterval region=new SingleInterval(tokens[0], new Integer(tokens[1]), new Integer(tokens[2]));
			String genotype=tokens[3];
			IntervalTree<String> tree;
			if(rtrn.containsKey(region.getReferenceName())){tree=rtrn.get(region.getReferenceName());}
			else{tree=new IntervalTree<String>();}
			tree.put(region.getReferenceStartPosition(), region.getReferenceEndPosition(), genotype);
			rtrn.put(region.getReferenceName(), tree);
			counter++;
			if(counter%1000000==0){System.err.println(counter);}
		}
		reader.close();
		return rtrn;
	}
	
	private static List<String> getUnique(List<String> probes1, List<String> probes2) {
		List<String> rtrn=new ArrayList<String>();
		Collection<String> seen=new TreeSet<String>();
		
		int counter=0;
		for(String p: probes1) {
			String[] tokens=p.split("\t");
			String seq=tokens[5];
			if(!seen.contains(seq)) {
				System.out.println(p);
				rtrn.add(p);
				seen.add(seq);
			}
			if(counter%1000==0) {System.err.println(counter);}
			counter++;
		}
		
		for(String p: probes2) {
			String[] tokens=p.split("\t");
			String seq=tokens[5];
			if(!seen.contains(seq)) {
				System.out.println(p);
				rtrn.add(p);
				seen.add(seq);
			}
			if(counter%1000==0) {System.err.println(counter);}
			counter++;
		}
		
		return rtrn;
	}
	
	public static void main(String[] args) throws NumberFormatException, IOException {
		if(args.length>3) {
			File genes=new File(args[0]);
			File snpFile1=new File(args[1]);
			File snpFile2=new File(args[2]);
			File genomeSeq=new File(args[3]);
			int probeLength=Integer.parseInt(args[4]);
			int maxNumIntronProbes=Integer.parseInt(args[5]);
			//int binSize=100;
			
			List<String> probes1=new SNPProbeDesign(genes, snpFile1, probeLength, genomeSeq, maxNumIntronProbes).getProbes();
			List<String> probes2=new SNPProbeDesign(genes, snpFile2, probeLength, genomeSeq, maxNumIntronProbes).getProbes();
			
			List<String> unique=getUnique(probes1, probes2);
			
		}
		else {System.err.println(usage);}
	}
	
	static String usage=" args[0]=genes (bed) \n args[1]=snp file1 \n args[2]=snp file2 \n args[3]=genomeSeq \n args[4]=probe length \n args[5]=maxNumberIntronProbes";
	
}
