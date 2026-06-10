package guttmanlab.core.probes;

// Standard Java imports for file I/O and data structures
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

// Custom imports for genomic annotations and data structures
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

    // Maximum number of probes allowed per intron (to avoid oversampling)
    int maxNumberPerIntron;

    // Primer pair sequences that will be added to the beginning and end of probes
    // PP1 is for exon probes, PP2 is for intron probes
    String PP1_left="TTTCGTCCGCGAGTGACCAG";
    String PP1_right="CAACGTCCATGTCGGGATGC";
    String PP2_left="TCAGGGCACGAGGACATTCG";
    String PP2_right="TCCGGCAAGATTGCTCTCCC";

    // List to store the final probe sequences
    List<String> probes;

    // Constructor 1: Generates probes and writes them to a file
    public SNPProbeDesign(File geneFile, File snpFile, int probeLength, String save, File genomeInputPath, int maxNum) throws NumberFormatException, IOException {
        this.maxNumberPerIntron = maxNum;

        // Read genome sequences from FASTA files into a map (chromosome -> sequence)
        Map<String, Sequence> genomeSequenceMap = FastaFileIOImpl.readFromFilesByName(genomeInputPath.listFiles());

        // Parse SNP file and create interval trees for efficient overlap queries
        Map<String, IntervalTree<String>> snpTree = parseTree(snpFile);

        // Load gene annotations from BED file
        Collection<Gene> genes = BEDFileIO.loadRegionsFromFile(geneFile.getAbsolutePath());

        // Extract intron regions from genes
        Collection<Gene> introns = getIntrons(genes);

        // Combine genes and introns for probe design
        genes.addAll(introns);

        // Collections to store SNPs and probes organized by gene
        Collection<Annotation> snps = new TreeSet<Annotation>();
        Map<String, Collection<Annotation>> probesByGene = new TreeMap<String, Collection<Annotation>>();

        // Process each gene/intron
        for(Gene gene: genes) {
            // Extract gene name (remove suffix after underscore)
            String name = gene.getName().split("_")[0];

            // Initialize probe collection for this gene if not exists
            if(!probesByGene.containsKey(name)) {
                probesByGene.put(name, new TreeSet<Annotation>());
            }

            // Find SNPs that overlap with this gene's exons
            Collection<SingleInterval> snpsByGene = getSNPs(gene, snpTree);

            // Generate probes around these SNPs
            Collection<Annotation> probesBySNP = getProbes(snpsByGene, probeLength, gene, snpTree.get(gene.getReferenceName()));

            // Add probes to gene-specific collection
            Collection<Annotation> set = probesByGene.get(name);
            set.addAll(probesBySNP);
            snps.addAll(probesBySNP);
        }

        // Write probes to output file with sequences and primer additions
        writeProbes(probesByGene, save, genomeSequenceMap, probeLength);
    }

    // Constructor 2: Generates probes and stores them in memory (doesn't write to file)
    public SNPProbeDesign(File geneFile, File snpFile, int probeLength, File genomeInputPath, int maxNum) throws NumberFormatException, IOException {
        this.maxNumberPerIntron = maxNum;

        // Same logic as constructor 1, but stores probes in memory
        Map<String, Sequence> genomeSequenceMap = FastaFileIOImpl.readFromFilesByName(genomeInputPath.listFiles());
        Map<String, IntervalTree<String>> snpTree = parseTree(snpFile);
        Collection<Gene> genes = BEDFileIO.loadRegionsFromFile(geneFile.getAbsolutePath());
        Collection<Gene> introns = getIntrons(genes);
        genes.addAll(introns);

        Collection<Annotation> snps = new TreeSet<Annotation>();
        Map<String, Collection<Annotation>> probesByGene = new TreeMap<String, Collection<Annotation>>();

        for(Gene gene: genes) {
            String name = gene.getName().split("_")[0];
            if(!probesByGene.containsKey(name)) {
                probesByGene.put(name, new TreeSet<Annotation>());
            }
            Collection<SingleInterval> snpsByGene = getSNPs(gene, snpTree);
            Collection<Annotation> probesBySNP = getProbes(snpsByGene, probeLength, gene, snpTree.get(gene.getReferenceName()));
            Collection<Annotation> set = probesByGene.get(name);
            set.addAll(probesBySNP);
            snps.addAll(probesBySNP);
        }

        // Store probes in instance variable instead of writing to file
        this.probes = writeProbes(probesByGene, genomeSequenceMap, probeLength);
    }

    // Getter method for probes
    private List<String> getProbes() {
        return this.probes;
    }

    // Extract intron regions from gene annotations
    private Collection<Gene> getIntrons(Collection<Gene> genes) {
        Collection<Gene> rtrn = new TreeSet<Gene>();

        for(Gene g: genes) {
            // Get intron intervals for each gene
            Collection<SingleInterval> introns = g.getIntrons();

            // Convert each intron to a Gene object with special naming
            for(SingleInterval intron: introns) {
                Gene i = new Gene(intron);
                i.setName("intron:" + g.getName());
                rtrn.add(i);
            }
        }
        return rtrn;
    }

    // Write probes to file with primer sequences added
    private void writeProbes(Map<String, Collection<Annotation>> probesByGene, String save, Map<String, Sequence> genomeSeq, int probeLength) throws IOException {
        FileWriter writer = new FileWriter(save);
        int counter = 0;

        for(String gene: probesByGene.keySet()) {
            Collection<Annotation> probes = probesByGene.get(gene);

            // Filter probes based on length and sequence quality
            Map<Annotation, Sequence> probeSeq = filterProbes(probes, genomeSeq, probeLength);

            // If too many intron probes, subsample to maxNumberPerIntron
            if(gene.contains("intron") && probeSeq.size() > maxNumberPerIntron) {
                probeSeq = subsample(probeSeq, maxNumberPerIntron);
            }

            // Write each probe with its sequence and primer-flanked version
            for(Annotation probe: probeSeq.keySet()) {
                Gene g = new Gene(probe);
                g.setOrientation(Strand.POSITIVE);
                Sequence seq = probeSeq.get(probe);

                // Add appropriate primers based on probe type
                String fullSeq = PP1_left + seq.getSequenceBases() + PP1_right;  // For exon probes
                if(gene.contains("intron")) {
                    fullSeq = PP2_left + seq.getSequenceBases() + PP2_right;     // For intron probes
                }

                // Only write probe if it passes quality filters
                if(!rejectProbe(seq.getSequenceBases())) {
                    writer.write(gene + "\t" + probe.getSingleInterval().toUCSCStrand() + "\t" + probe.getName() + "\t" + probe.getReferenceName() + "\t" + seq + "\t" + fullSeq + "\n");
                    System.out.println(probe.toBED());
                }

                counter++;
                if(counter % 100 == 0) {
                    System.err.println(counter + " " + probesByGene.size());
                }
            }
        }
        writer.close();
    }

    // Same as above but returns list instead of writing to file
    private List<String> writeProbes(Map<String, Collection<Annotation>> probesByGene, Map<String, Sequence> genomeSeq, int probeLength) throws IOException {
        List<String> rtrn = new ArrayList<String>();
        int counter = 0;

        for(String gene: probesByGene.keySet()) {
            Collection<Annotation> probes = probesByGene.get(gene);
            Map<Annotation, Sequence> probeSeq = filterProbes(probes, genomeSeq, probeLength);

            if(gene.contains("intron") && probeSeq.size() > maxNumberPerIntron) {
                probeSeq = subsample(probeSeq, maxNumberPerIntron);
            }

            for(Annotation probe: probeSeq.keySet()) {
                Gene g = new Gene(probe);
                g.setOrientation(Strand.POSITIVE);
                Sequence seq = probeSeq.get(probe);

                String fullSeq = PP1_left + seq.getSequenceBases() + PP1_right;
                if(gene.contains("intron")) {
                    fullSeq = PP2_left + seq.getSequenceBases() + PP2_right;
                }

                if(!rejectProbe(seq.getSequenceBases())) {
                    rtrn.add(gene + "\t" + probe.getSingleInterval().toUCSCStrand() + "\t" + probe.getName() + "\t" + probe.getReferenceName() + "\t" + seq + "\t" + fullSeq);
                }

                counter++;
                if(counter % 100 == 0) {
                    System.err.println(counter + " " + probesByGene.size());
                }
            }
        }
        return rtrn;
    }

    // Quality control: reject probes with poor sequence characteristics
    private static boolean rejectProbe(String probe) {
        // Reject if contains N bases
        if(probe.contains("N")) {
            return true;
        }

        // Initialize various sequence filters
        LowComplexityFilter lcf = new LowComplexityFilter();           // Filters low complexity sequences
        PolyBaseFilter pbf = new PolyBaseFilter("ACGT", 5, 5);        // Filters short polybase runs
        PolyBaseFilter pbf2 = new PolyBaseFilter("ACGT", 15, 12);     // Filters long polybase runs
        PolyBaseFilter pbf3 = new PolyBaseFilter("ACGT", 8, 6);       // Filters medium polybase runs

        // Return true if any filter rejects the sequence
        return lcf.rejectSequence(probe) || pbf.rejectSequence(probe) || pbf2.rejectSequence(probe) || pbf3.rejectSequence(probe);
    }

    // Filter probes based on length and sequence quality
    private Map<Annotation, Sequence> filterProbes(Collection<Annotation> probes, Map<String, Sequence> genomeSeq, int probeLength) {
        Map<Annotation, Sequence> rtrn = new TreeMap<Annotation, Sequence>();

        for(Annotation probe: probes) {
            Gene g = new Gene(probe);
            g.setOrientation(Strand.POSITIVE);

            // Extract sequence for this probe
            Sequence seq = g.getSequence(genomeSeq.get(probe.getReferenceName()));

            // Only keep probes with correct length and no N bases
            if(seq.getLength() == probeLength && !seq.hasN()) {
                rtrn.put(probe, seq);
            }
        }
        return rtrn;
    }

    // Randomly subsample probes if there are too many
    private Map<Annotation, Sequence> subsample(Map<Annotation, Sequence> probes, int num) {
        Map<Annotation, Sequence> rtrn = new TreeMap<Annotation, Sequence>();
        List<Annotation> list = new ArrayList<Annotation>();
        list.addAll(probes.keySet());

        // Randomly select 'num' probes
        for(int i = 0; i < num; i++) {
            int index = (int)(Math.random() * list.size());
            Annotation probe = list.get(index);
            rtrn.put(probe, probes.get(probe));
            list.remove(index);
        }
        return rtrn;
    }

    // Generate probe sequences around SNPs
    private Collection<Annotation> getProbes(Collection<SingleInterval> snpsByGene, int probeLength, Gene gene, IntervalTree<String> snpTree) {
        Collection<Annotation> rtrn = new TreeSet<Annotation>();

        for(SingleInterval r: snpsByGene) {
            int add = probeLength / 2;

            // Create three potential probe positions around each SNP:
            // 1. Centered on SNP
            SingleInterval probe = new SingleInterval(r.getReferenceName(), r.getReferenceStartPosition() - add, r.getReferenceEndPosition() + add);
            // 2. Upstream of SNP
            SingleInterval upstream = new SingleInterval(r.getReferenceName(), r.getReferenceStartPosition() - probeLength, r.getReferenceEndPosition());
            // 3. Downstream of SNP
            SingleInterval downstream = new SingleInterval(r.getReferenceName(), r.getReferenceStartPosition(), r.getReferenceEndPosition() + probeLength);

            // Convert to genome coordinates
            Annotation genomeProbe = convertToGenome(probe, gene);
            Annotation ugenomeProbe = convertToGenome(upstream, gene);
            Annotation dgenomeProbe = convertToGenome(downstream, gene);

            // Set names for identification
            genomeProbe.setName(r.toUCSC() + "_centered");
            ugenomeProbe.setName(r.toUCSC() + "_upstream");
            dgenomeProbe.setName(r.toUCSC() + "_downstream");

            // Only add probes that don't overlap with other SNPs
            if(!overlapsSNP(ugenomeProbe, snpTree)) {
                rtrn.add(ugenomeProbe);
            }
            if(!overlapsSNP(dgenomeProbe, snpTree)) {
                rtrn.add(dgenomeProbe);
            }
        }
        return rtrn;
    }

    // Check if a probe overlaps with any SNPs (to avoid multiple SNPs per probe)
    private boolean overlapsSNP(Annotation ugenomeProbe, IntervalTree<String> snpTree) {
        Iterator<SingleInterval> iter = ugenomeProbe.getBlocks();

        while(iter.hasNext()) {
            SingleInterval block = iter.next();
            // Check if any block of the probe overlaps with SNPs
            if(snpTree.hasOverlappers(block.getReferenceStartPosition(), block.getReferenceEndPosition())) {
                return true;
            }
        }
        return false;
    }

    // Convert relative coordinates to genome coordinates
    private Annotation convertToGenome(SingleInterval region, Annotation gene) {
        Annotation rtrn = gene.convertToReferenceSpace(region);
        rtrn.setOrientation(Strand.POSITIVE);
        return rtrn;
    }

    // Find SNPs that overlap with a gene's exons
    private Collection<SingleInterval> getSNPs(Gene gene, Map<String, IntervalTree<String>> snpTree) {
        Collection<SingleInterval> rtrn = new TreeSet<SingleInterval>();

        if(snpTree.containsKey(gene.getReferenceName())) {
            IntervalTree<String> tree = snpTree.get(gene.getReferenceName());

            // Find all SNPs that overlap with this gene
            Iterator<Node<String>> iter = tree.overlappers(gene.getReferenceStartPosition(), gene.getReferenceEndPosition());

            while(iter.hasNext()) {
                Node<String> snp = iter.next();
                SingleInterval r = new SingleInterval(gene.getReferenceName(), snp.getStart(), snp.getEnd());
                r.setOrientation(gene.getOrientation());
                r.setName(snp.getValue());

                // Only keep SNPs that overlap with exons (not introns)
                if(gene.overlapsExon(r)) {
                    SingleInterval relative = getRelativeCoordinates(gene, r);
                    rtrn.add(relative);
                }
            }
        }
        return rtrn;
    }

    // Convert genome coordinates to relative coordinates within a gene
    public static SingleInterval getRelativeCoordinates(Annotation gene, SingleInterval feature) {
        int featureStart = gene.getRelativePositionFrom5PrimeOfFeature(feature.getReferenceStartPosition());
        int featureEnd = gene.getRelativePositionFrom5PrimeOfFeature(feature.getReferenceEndPosition());

        SingleInterval rtrn = new SingleInterval(gene.getName(), featureStart, featureEnd);

        // Handle negative strand genes
        if(gene.getOrientation().equals(Strand.NEGATIVE)) {
            rtrn = new SingleInterval(gene.getName(), featureEnd, featureStart);
        }

        return rtrn;
    }

    // Parse SNP file into interval trees for efficient querying
    private Map<String, IntervalTree<String>> parseTree(File snpFile) throws NumberFormatException, IOException {
        Map<String, IntervalTree<String>> rtrn = new TreeMap<String, IntervalTree<String>>();
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(snpFile)));

        String nextLine;
        int counter = 0;

        while ((nextLine = reader.readLine()) != null) {
            String[] tokens = nextLine.split("\t");

            // Parse SNP coordinates and genotype
            SingleInterval region = new SingleInterval(tokens[0], new Integer(tokens[1]), new Integer(tokens[2]));
            String genotype = tokens[3];

            // Add to interval tree for this chromosome
            IntervalTree<String> tree;
            if(rtrn.containsKey(region.getReferenceName())) {
                tree = rtrn.get(region.getReferenceName());
            } else {
                tree = new IntervalTree<String>();
            }

            tree.put(region.getReferenceStartPosition(), region.getReferenceEndPosition(), genotype);
            rtrn.put(region.getReferenceName(), tree);

            counter++;
            if(counter % 1000000 == 0) {
                System.err.println(counter);
            }
        }

        reader.close();
        return rtrn;
    }

    // Get unique probes from two different probe sets (removes duplicates based on sequence)
    private static List<String> getUnique(List<String> probes1, List<String> probes2) {
        List<String> rtrn = new ArrayList<String>();
        Collection<String> seen = new TreeSet<String>();
        int counter = 0;

        // Process first probe set
        for(String p: probes1) {
            String[] tokens = p.split("\t");
            String seq = tokens[5];  // Full sequence with primers

            if(!seen.contains(seq)) {
                System.out.println(p);
                rtrn.add(p);
                seen.add(seq);
            }

            if(counter % 1000 == 0) {
                System.err.println(counter);
            }
            counter++;
        }

        // Process second probe set
        for(String p: probes2) {
            String[] tokens = p.split("\t");
            String seq = tokens[5];

            if(!seen.contains(seq)) {
                System.out.println(p);
                rtrn.add(p);
                seen.add(seq);
            }

            if(counter % 1000 == 0) {
                System.err.println(counter);
            }
            counter++;
        }

        return rtrn;
    }

    // Main method: processes two SNP files and merges unique probes
    public static void main(String[] args) throws NumberFormatException, IOException {
        if(args.length > 3) {
            File genes = new File(args[0]);              // Gene annotations (BED format)
            File snpFile1 = new File(args[1]);           // First SNP file
            File snpFile2 = new File(args[2]);           // Second SNP file
            File genomeSeq = new File(args[3]);          // Genome sequence files
            int probeLength = Integer.parseInt(args[4]);  // Desired probe length
            int maxNumIntronProbes = Integer.parseInt(args[5]); // Max intron probes per gene

            // Generate probes for both SNP files
            List<String> probes1 = new SNPProbeDesign(genes, snpFile1, probeLength, genomeSeq, maxNumIntronProbes).getProbes();
            List<String> probes2 = new SNPProbeDesign(genes, snpFile2, probeLength, genomeSeq, maxNumIntronProbes).getProbes();

            // Merge and get unique probes
            List<String> unique = getUnique(probes1, probes2);

        } else {
            System.err.println(usage);
        }
    }

    // Usage message
    static String usage = " args[0]=genes (bed) \n args[1]=snp file1 \n args[2]=snp file2 \n args[3]=genomeSeq \n args[4]=probe length \n args[5]=maxNumberIntronProbes";
}
