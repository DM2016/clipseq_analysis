#CLIP-seq data analysis walkthrough

Copyright &copy; 2015 Brian Sebastian Cole, PhD, Kristen Lynch, PhD, and the University of Pennsylvania, PhD

Permission is granted to copy, distribute, and/or modify this document under the terms of the GNU Free Documentation License, Version 1.3 or any later version published by the Free Software Foundation; with no Invariant Sections, with the Front-Cover Texts being TITLE and COPYRIGHT, and with no Back-Cover Texts. To view a copy of this license, visit:
https://www.gnu.org/licenses/fdl.html


##Table of contents
#####1. Overview and philosophy
#####2. Trimming and preprocessing raw CLIP-seq reads
#####3. Aligning CLIP-seq reads to a genome and/or transcriptome index
#####4. Removing duplicate alignments
#####5. Optional: pooling replicates
#####6. Calling peaks
#####7. Post-processing peaks
#####8. Motif enrichment analysis
#####9. Comparing sets of CLIP-seq peaks
#####10. RNA cartography

##Chapter 1. Overview and Philosophy

CLIP-seq, or <u>C</u>ross<u>l</u>inking and <u>i</u>mmuno<u>p</u>recipitation followed by <u>seq</u>uencing, is an experimental technique for identifying protein-RNA interaction sites <i>in vivo</i>.  The ability of this approach to allow investigators to assay protein-RNA interactions in living cells and tissues has provided valuable insights into the activity and potential mechanisms of dozens of RNA binding proteins with a diverse array of functions, including splicing factors, mRNA stability factors, miRNA silencing complexes, and many more.
Because RNA binding proteins control nearly all aspects of RNA metabolism, CLIP-seq is a vital component of the RNA biologist's toolchain. However, CLIP-seq requires extensive data analysis efforts to generate information from the raw sequence data it produces. The purpose of this walkthrough is to fully explain all of the steps in a complete CLIP-seq processing pipeline that begins with raw reads and ends with insights. Along the way, technical challenges arising from the many types of experimental designs commonly used in existing CLIP-seq datasets are discussed, including caveats and limitations.
The philosophy of this document is to explore the consequences of different decisions and parameters along the way.  There are several steps in the analysis of CLIP-seq data for which there is (as of yet) no objective standard or "best practice."  **There is more than one way to do it,** and the only way to truly understand the consequences of each action when faced with multiple choices is to explore them. This document will highlight a few key places in which there is no clear consensus on best practice, but will always recommend a specific approach.  Keep in mind that these recommendations may not work best for every CLIP-seq analysis. **Explore** and keep an open mind - that's what makes discovery possible. 


##Chapter 2. Trimming and preprocessing raw CLIP-seq reads

CLIP-seq reads are short oligonucleotide sequences generated from cDNAs prepared from CLIP RNA.  All pipelines that process CLIP-seq reads must begin with the sequences themselves, often in a format called FASTQ that contains the sequence information as well as a quality string that relates the confidence of the basecalls.  
The first step in most CLIP-seq analysis pipelines is to trim the raw reads from the 3' end to remove sequencing adaptors, which are oligonucleotides that were ligated to the 3' end of crosslinked RNA fragments during library preparation.  To do this, tools such as [cutadapt](https://code.google.com/p/cutadapt/) are commonly used.  Running cutadapt is easy, just pass it the sequence of the particular 3' adaptor used in the CLIP-seq library preparation.  
Let's say the adaptor sequence was 'CAGTGCAG' and the FASTQ file containing CLIP-seq reads is called raw_clipseq_reads.fastq on your computer: 

    cutadapt -a CAGTGCAG -o raw_clipseq_reads.fastq trimmed_clipseq_reads.fastq

This will generate a new output file called `trimmed_clipseq_reads.fastq` that also contains FASTQ formatted-data, but the 3' sequencing adaptors have been removed, along with any sequences that might have occurred downstream.
The purpose of removing sequencing adaptors is to recapture the cDNA sequence of the RNA that was originally crosslinked to the protein during irradiation with UV light.  Use of cutadapt or a similar tool to remove 3' adaptors is strongly recommended before aligning CLIP-seq reads to the genome and/or transcriptome under study.

>**Note:** Some analysts have also removed homopolymeric end runs and low-quality stretches at the 3' end of reads after removing 3' sequencing adaptors.  Cutadapt has an option for quality trimming, and aligners such as bowtie can use a 5' end-seeded alignment strategy that mitigates the deleterious effects of homopolymeric end runs. 

>**Note:** Some analysts "collapse" duplicated sequences after trimming and before alignment to the genome/transcriptome. The reason for this is to remove PCR duplicates.  Because sequencing errors can generate single nucleotide sequence differences (polymorphisms) between CLIP-seq reads that will still align to the same genomic coordinates, reomving duplicate alignment coordinates instead of duplicate sequences is recommended.  However, this is one of many examples of differences between pipelines in different analyses, and experimentation and exploration of the consequences of different approaches is the only way to be truly convinced of the validity of any given approach.


##Chapter 3. Aligning CLIP-seq reads to a genome and/or transcriptome index

Sequence data themselves are almost never interesting, instead, binding sites, or protein-RNA interaction sites, are the chief unit of interest.  Binding sites are not defined by the sequences within CLIP-seq reads, either raw or trimmed reads.  Instead, they are a higher-order unit of information that emerges from the overlap of aligned CLIP-seq reads. For this reason, alignment of CLIP-seq reads to a genome and/or transcriptome index is almost universally at the heart of CLIP-seq analysis pipelines.
RNA does not exist alone in cells: instead, it is almost always accompanied by proteins.  Proteins make contact with RNA throughout virtually all steps of RNA metabolism.  For this reason, a protein under study in a given CLIP-seq experiment might crosslink to RNA sequences from different steps of the RNA lifecycle.  Splicing regulatory proteins, for instance, may crosslink to pre-mRNA and mature mRNA in the same CLIP sample.  However, analysis that focuses on capturing these two types of interactions relies on separate alignment strategies.
The first choice that a CLIP-seq analyst must make is what index to align trimmed CLIP-seq reads to.  Aligning reads to an index generated from a processed transcriptome, e.g. refSeq or ENSEMBL annotations, can allow an analyst to examine protein-mRNA interactions directly.  However, potential interactions with intronic sequences in pre-mRNAs cannot be identified by this approach.  On the other hand, aligning reads directly to the genome will often disallow alignment of reads from spliced mRNA.  
For this reason, a balanced approach is recommended.  First, consider the protein under study.  Is it likely to interact mostly with processed mRNA, or pre-mRNA?  If it's the former, alignment to a transcriptome index might be the best approach as it directly assesses protein-mRNA interactions.  If it's the latter, alignment to the genome is almost certainly a better choice.

>**Note:** Alignment to the genome is recommended for CLIP-seq analysis. Protein-mRNA interactions across spliced junctions can still be captured by a second step alignment approach in which reads that do not linearly align to the genome index are realigned using an aligner that can introduce gaps into reads (a junction-capable aligner) such as bowtie2. Stay tuned.

Sequence alignment is at the heart of a huge diversity of bioinformatics analyses, so unsurprisingly there exists an embarassment of riches with respect to alignment software.  One of the most popular aligners is [bowtie]{http://bowtie-bio.sourceforge.net/index.shtml).  Bowtie is simple, easy-to-use, fast, efficient, and flexible.  There are many other choices out there, but bowtie offers the advantage of wide adoption and the accompanying support of a broad user base.  Additionally, bowtie and other tools availalbe by the authors are actively supported and under active development with a mature release cycle.  The rest of this tutorial will use bowtie.

Let's assume that the protein under study is hypothesized to interact with pre-mRNA, perhaps to influence pre-mRNA processing. In this case, alignment of trimmed CLIP-seq reads to the genome and not the transcriptome is recommended.  Using the trimmed reads example from above, the following bowtie invocation will generate alignments from trimmed CLIP-seq reads:

    bowtie -m 1 /path/to/genome/index/basename trimmed_clipseq_reads.fastq > clipseq_alignments.txt
    
The resulting output file contains alignment coordinates for the trimmed CLIP-seq reads.  In this case, the `-m 1` option specifies that only reads that align to one place in the genome index are reported.  This recommended policy prevents one trimmed read from generating multiple alignments (which distorts arithmetic such as alignment performance, albeit not intractably) and additionally does not report any alignment for which there was another possible alignment position.  This level of stringency is highly recommended.  Consult the excellent bowtie manual online or the embedded documentation by invoking the program with no arguments for more information.

>**Note:** Some reads will not align to the genome at all.  This could be because they are too short to generate unique alignment positions, because the entire read was removed by trimming, or because the read contains too many sequencing errors.  Additionally, it is possible that polymorphisms between the particular organism under study and the organism(s) from which the genome index was sequenced contribute to poor alignment of a given sequence (even if that sequence contained no sequencing errors).  It is highly recommended to compute the fraction of trimmed CLIP-seq reads that **could** be aligned (say, those with length of at least 8 nucleotides after trimming) to the total number of unique alignments reported by bowtie.

>**Note:** Some reads could not align linearly to the genome because they result from a spliced mRNA.  An approach wherein the reads that did not align linearly to the genome with a linear aligner such as bowtie are then realigned with a junction-capable aligner can properly map these reads.  However, strong caution must be taken in interpreting the results.  Keeping these junction-spanning alignments separate from linear alignments might be the best option.  The rest of this tutorial will assume linear alignments only, for example a CLIP-seq study of a pre-mRNA processing factor.

>**Note:** Many CLIP-seq experiments incorporate replicates.  Replicates should be trimmed and aligned separately.


##Chapter 4: Removing duplicate alignments

Many cycles of PCR are used in most CLIP-seq library preparations.  This results in the introduction of well-studied bias in PCR amplification with respect to fragment length and sequence composition.  For this reason, many CLIP-seq analysts recommend removing duplicated CLIP-seq reads to control for the bias inherent to PCR amplification.
As mentioned above, there are two reasonable approaches to removing PCR duplicates.  The first is to remove duplicated oligonucleotide sequences after trimming but **before** alignment.  The second is to remove duplicated alignment positions **after** alignment.

>**Note:** Removing duplicated alignment positions is recommended to control for sequencing errors.  Sequencing error rate is relatively high on modern sequencers, so high indeed, that two sequence reads resulting from identical oligonucleotides are expected to contain 1 polymorphism between them.  This means that removing duplicated sequences instead of duplicated alignments is more conservative.  An approach in which alignment positions are "collapsed" to one unique representative is recommended because it allows the alignment algorithm, for instance bowtie, to decide which sequence reads originated from which region of the genome or transcriptome.

There are three possible ways to collapse duplicated alignments.  The first is to remove any duplicated start **and** end coordinates, such that for a given start and end position on a given strand of a given chromosome in the genome, there can be at most one alignment.  A more conservative approach would be to remove any duplicated **5' coordinates.**

>**Note:** Removing duplicated alignments on the basis of 5' coordinate is recommended.  This controls for the increasing sequencing error rate at the 3' ends of reads.  Reads resulting from identical cDNAs display markedly higher entropy at the 3' ends as compared to the 5' ends.  This affects not only alignment, but also trimming and adaptor removal (e.g. with cutadapt).  For this reason, collapsing duplicated alignments on the basis of **5' coordinate** alone is a more conservative and recommended approach.

Duplicates should be removed separately for each replicate before pooling replicates, if pooling is the desired approach. Further discussion of replicates is provided below.

Once duplicate alignments have been removed (or "collapsed"), the alignments are ready to be searched for peaks.  The goal of peak calling is to identify significant interaction sites and discard noise.

>**Note:** Creating a genome browser session to visually inspect alignments before and after collapsing duplicates is highly recommended.  Bedgraph format is a useful visualization format for aligned CLIP-seq reads.  At this point in the pipeline, creating a genome browser session with bedgraph files of aligned reads before and after duplicate removal allows the analyst to visually inspect the effects of duplicate removal.  Addtionally, the visual signals in the bedgraph files can provide early insights into the performance of the CLIP-seq experiment.  For instance, if there is a "known" site of protein-RNA interaction for the protein under study, inspecting this site in the genome can confirm the ability of the CLIP-seq experiment to generate coverage at that site.

>**Note:** Visualizing CLIP-seq data in a genome browser is crucial.  CLIP-seq often generates diffuse signals with isolated regions of increased signal.  The level of this diffuseness in the signal depends on the biology of the protein as well as the sequencing depth.  There is really no substitute for visualization, and genome browsers are an essential resource for gaining understanding and intuition about a given CLIP-seq study.

The default bowtie output format is a tab-delimited flat text file that is highly amenable to collapsing duplicate alignments.  Alternatively, this tab-delimited output format can be easily converted to BED format, a standard 6-column tab-delimited format that contains genomic coordinates.  BED format can be directly visualized in genome browsers, and supports the powerful genomic arithmetic suite [bedtools](https://github.com/arq5x/bedtools).  Many aligners also support the SAM format, which can be converted to BED format.  The particular formatting choices are left up to the user, however, suites such as bedtools provide easy reformatting capabilities such that most formats can be reformatted into most other formats with a single command.


##Chapter 5: Optional: pooling replicates

If replicate CLIP-seq samples are used, combining alignments before calling peaks can increase the discovery power of the peak caller.  This approach is not always highly desired, for instance some analysts may wish to call peaks separately on each replicate (after collapsing duplicated alignments) and then apply a replication criterion.
One possible replication criterion would be to require that a peak be called in at least n - 1 replicates where n replicate libraries were created.  For instance, if three CLIP-seq libraries were prepared from the same sample group, calling peaks on each of three alignments and then discarding any peak that wasn't called in at least 2 of the 3 replicates is a valid approach to incorporate replicate information.  The chief challenge in this replicate approach is to consider the varying degrees of overlap that are possible.  For instance, two CLIP-seq peaks could overlap by only a small number of nucleotides, like two bell curves that only cross at the tails.  Is this instance considered a shared peak, or an unshared peak?  If it is considered shared, what is the shared peak? Is it the union of the two peaks, or the intersection?  These decisions have to be made if peaks are called separately on each replicate.

>**Note:** Calling peaks on pooled replicates is simpler as it considers the sample group instead of individual samples and therefore does require the analyst to decide upon a strategy to merge replicate peaks.  Replication information can then still be incorporated by requiring that CLIP-seq peaks have coverage from at least n - 1 replicates at the level of alignments.

Pooling replicates can be as easy as concatenating files.  In addition, some peak callers support multiple replicates and implicitly pool the replicates automatically.

>**Note:** Calling peaks on pooled replicates and then discarding peaks that had did not have CLIP-seq reads in all replicates (or, optionally, n - 1 replicates) is a hybrid approach that enforces a replicate requirement but avoids the tricky process of merging replicate peaks together.  Optionally, no replicate requirement could be used at all.  The particular choice of replicate requirement is left to the user, but an n - 1 replicate requirement strikes a balance between conservative and liberal approaches.


##Chapter 6: Calling peaks

Calling peaks is the process by which significant interaction sites are identified within CLIP-seq alignment signals.  This process is crucial in that it generates potential binding sites and discards noise.  The primary unit of interest in CLIP-seq experiments is often the peaks, not the aligned reads themselves.

A popular algorithm for defining CLIP-seq peaks was described in the literature and is implmented in the clipseq_analysis distribution.  Other popular implementations include [clipper](https://github.com/gpratt/clipper) and [pyicoteo](https://github.com/RegulatoryGenomicsUPF/pyicoteo).  All three implementations share underlying design features, and the particular choice of peak caller is left to the user.  The rest of this documentation will cover the [clipseq_anlaysis](https://github.com/bryketos/clipseq_analysis) distribution.

This algorithm identifies peaks by using a transcript-specific permutation algorithm.  This algorithm considers alignments within each transcript separately, setting thresholds for peaks that are specific to that transcript.  Transcript-specific analaysis is necessary due to the large dynamic range of RNA expression levels.  For this reason, a set of transcripts to be searched for peaks is required.  Two popular transcriptome annotations are refSeq and ENSEMBL.

>**Note:** ENSEMBL contains more transcripts and more isoforms of those transcripts than refSeq does.  While the particular choice of transcriptome annotation to search for peaks is left to the analyst, the refSeq transcriptome is a stringently curated and conservative annotation.  Additionally, refSeq annotations are available from the UCSC Genome Bioinformatics server in BED12 format, which is a convenient format for bioinformatic analysis.

Peaks should be called separately on each sample group where multiple sample groups exists.  For example, the same CLIP-seq experiment could have been performed in two different cell types, or the same cell type in two conditions.  To continue the above example, let's assume an experiment in which a single replicate of CLIP-seq reads was generated.  The reads were then preprocessed to remove sequencing adaptors and aligned to the genome.  Then, duplicate alignments were removed, generating a file called `collapsed_clipseq_alignments.bed`, which is a BED file containing the alignment coordinates of CLIP-seq reads after duplicates were removed.

Calling peaks requires only this BED file and a transcriptome annotation file.  Let's assume the transcriptome annotation was obtained from UCSC and is called `refseq.hg19.bed`, which exists in either BED or BED12 format.  Calling peaks is achieved via the following command:

    discover_peaks -a collapsed_clipseq_alignments.bed -t refseq.hg19.bed -o clipseq_peaks.bed
    
This command will search the entire transcriptome for peaks and generate a BED file containing the peaks called `clipseq_peaks.bed`.  This file contains the binding sites and discards the noise.  These peaks are then highly amenable to downstream analysis.

The `-a`, `-t`, and `-o` options are all required and respectively contain the CLIP-seq alignments to be searched for peaks, the transcriptome annotation to be searched for peaks, and the output file to write peaks to.  

Several additional options allow fine-grained control over the stringency of the algorithm and the verboseness of reporting.  Consult the documentation of the program (or the particular peak caller used) for more information.  However, the default settings are rational and need not be altered in most cases.

>**Note:** Visualization of this peaks file using a genome browser is highly recommended.  In addition to the BED file itself, reads within the BED file can be converted to bedgraph format using bedtools.  The resulting bedgraph files provide a useful visualization when loaded into a genome browser session also containing the aligned reads.  This allows an analyst to easily see which regions of the genome were called as peaks and which were discarded as noise.

If replicates CLIP-seq libraries were used, they can be provided as a comma-separated list:

    discover_peaks -a collapsed_clipseq_alignments_1.bed,collapsed_clipseq_alignments_2.bed -t refseq.hg19.bed -o clipseq_peaks.bed
    

##Chapter 7: Post-processing peaks

After peak calling, some CLIP-seq analysts wish to join neighboring peaks that occur within short distances of each other.  For instance, if two peaks occur next to each other in the genome (on the same strand) with only 9 nucleotides between, the analyst may which to consider this region as 1 peak instead of 2.  The particular choice of whether or not to join neighboring peaks is left up to the user.  

>**Note:** Joining neighboring CLIP-seq peaks is not required.  Additionally, for subsequent downstream analysis, the inclusion of the interpeak regions might affect motif enrichment results.  One might argue that regions of the genome that were not identified as peaks by the peak caller should not be included in the analysis - that is, that peaks should not be joined 


##Chapter 8: Motif enrichment analysis

Motif enrichment analysis from CLIP-seq peaks is a useful knowledge output from CLIP-seq experiments.  First, the motifs elicited by de novo motif enrichment of CLIP-seq peaks can be compared to prior knowledge obtained from biochemical studies such as SELEX.

An algorithm for motif enrichment analysis using CLIP-seq peaks is provided in the [clipseq_anlaysis](https://github.com/bryketos/clipseq_analysis) distribution.  It uses a permutation approach to assign Z-scores (enrichment scores) to kmers.  Kmers with high Z-scores are enriched within CLIP-seq peaks as compared to backgrounds obtained from permutation.

To identify enriched hexamers (kmers 6 nucleotides in length) in CLIP-seq peaks, the following command may be used:

    compute_z_scores clipseq_peaks.bed refseq.hg19.bed motifs.tsv
    
The command will use the same transcriptome annotation file used in peak calling to obtain background motif frequencies.  The resulting output file will contain Z scores and other useful information in a tab-delimited format.  This output file can be opened with spreadsheet software or a text editor.

>**Note:** Many users wish to generate a composite motif logo.  One approach to this is to take the top 20 most enriched kmers by Z score (highest Z score kmers) and align them against each other with a tool such as Clustal.  The resulting multiple sequence alignment can be used to generate a motif logo with tools such as [weblogo](http://weblogo.berkeley.edu/logo.cgi).  The multiple sequence alignment step is essential in this approach.  The resulting logo output conveys a consensus motif and is often an important visualization of motif enrichment results from CLIP-seq.  The choice of the number of kmers to include, in this example 20, is not set in stone and left to the discretion of the user.  Too few kmers may generate poor multiple sequence alignment results, but too many might introduce noise into the resulting logo.


##Chapter 9: Comparing sets of CLIP-seq peaks

CLIP-seq and its variants (PAR-CLIP and iCLIP) have found extensive adoption in scientific research, resulting in a plethora of publically accessible CLIP-seq datasets.  This affords the possibility of comparative analysis of CLIP-seq datasets.

Additionally, some experimental designs that incorporate multiple sample groups, for instance the same protein in two cell types or in two conditions of the same cell type, are amenable to comparative analysis.

>**Note:** Strong caution in interpretation of comparative CLIP-seq analysis is recommended.  Experience argues that CLIP-seq is a noisy experiment, with overlap between replicates of even monoclonal cell lines often highly variable within the same experiment.  Additionally, many scientists would argue that CLIP-seq is not a quantitative experiment in the sense that the number of reads within a peak is directly reflective of the intensity of the interaction at that locus.  Reasons for this argument against quantitative interpretation of CLIP-seq include not only high levels of variability between replicates, but also biases introduced during CLIP-seq library preparation (e.g. sequence specificity of RNase and ligase enzymes, PCR amplification bias) and biases introduced by bioinformatic analaysis (e.g. alignment bias due to unequal mappability of the search space.)  While comparative analysis of CLIP-seq experiments is possible using the following techniques, any results should be interpreted with the highest degree of caution.  Conversely, negative results must be interpreted via the maxim "absense of evidence is not evidence of absense."

Let's assume a block study design in which CLIP-seq libraries for a protein of interest were prepared in two cell states, such as treatment with an experimental drug versus mock treatment.  Let's assume these two sample groups were processed as outlined above, and the extent and significance of the overlap between CLIP-seq peaks between these two conditions is of interest.  If these two CLIP-seq peaks exist in files named `control_peaks.bed` and `treated_peaks.bed` for the CLIP-seq peaks resulting from untreated and treated sample groups, respectively, the following tool will compare their overlap:

    compare_clip_profiles -a control_peaks.bed -b treated_peaks.bed -t refgene.hg19.bed -o peaks_comparison.tsv
    
The result is an output file that contains the total number of peaks in the first file provided by the `-a` argument, the total count of these peaks that overlap any peaks in the second file provided by the `-b` argument, and the total number of overlaps that are expected based on permutation.  Importantly, this input is amenable to downstream statistical hypothesis testing.

The particular statistical analysis used is left to the discretion of the user.

##Chapter 10: RNA cartography

Integrative genomics is a field that has experienced rapid recent growth and seen wide adoption in the study of the function of RNA binding proteins that regulate pre-mRNA processing.  In particular, the study of alternative splicing regulated by an RNA binding protein can utilize an approach called RNA mapping wherein the fraction of regulated alternative splicing events that contain CLIP-seq peaks is compared to the fraction of unregulated events.  
Importantly, RNA maps provide insight into RNA binding protein occupancy within and around regulated splicing events.  Patterns of enrichment that are visualized in RNA maps can provide mechanistic hypotheses relating the position of binding to the outcome of splicing regulation.  For example, a splicing regulatory protein that displays marked enrichment within exons that it enhances could be evidence of a direct, exonic splicing enhancer function.

>**Note:** RNA maps are associative and do not imply causal relationships.  Results from integrative genomics analysis such as RNA maps should be treated as hypothesis-generating results, not hypothesis testing.  Additionally, experience argues that the overlap between binding and function is often low (10% or less).  These results must be interpreted with a grain of salt, and further functional validation with directed experiments is necessary.

RNA maps require prior establishment of sets of exons that respond to protein knockdown with altered inclusion.  While detailed coverage of functional studies of splicing regulators is beyond the scope of this document, techniques such as RNA-seq and splicing-sensitive microarrays can generate sets of exons that respond to knockdown of a splicing factor under study with significant alterations in inclusion levels.

Let's assume a BED file containing exons of interest exists and is named `interesting_exons.bed`.  These exons could be splicing targets of the RNA binding protein under study, for example.  To generate an RNA map for these exons using CLIP-seq peaks, the following command could be used:

    binding_mapper -f -p clipseq_peaks.bed -u upstream_map.tsv -d downstream_map.tsv interesting_exons.bed
    
This will generate two output files, one containing the fraction of exons bound at each nucleotide in a window that begins 250 nucleotides upstream of the exons' 3' splice sites and ending 50 nucleotides downstream of the 3' splice sites, therefore containing 250 nucleotides of proximal introns and 50 nucleotides of the exons.  Similarly a downstream ouptut file will be generated containing the 3' most 50 nucleotides of the exons and the 250 nucleotides of proximal introns downstream of the 5' splice sites.
These output files are tab-delimited text and are amenable to plotting with tools such as R, Python's matplotlib, gnuplot, and others.  Importantly, comparison of multiple series within a single plot allows the user to derive mechanistic hypotheses.  For example, enrichment of protein occupancy within exons that are enhanced by a splicing factor (and therfore respond to depletion of that splicing factor with significantly reduced inclusion levels) relative to a set of unresponsive exons (or to a set of repressed exons) could indicate a direct, exonic splicing enhancer role.

>**Note:** The fraction of a set of exons that are bound is a recommended output, enabled with the `-f` option to the `binding_mapper` program.  This mode treats the presence or absence of CLIP-seq peaks in a Boolean manner and therefore is not subject to the quantitative interpretation of CLIP-seq peaks.  This is in contrast to other possible approaches to RNA mapping, such as computing average number of CLIP-seq reads at each nucleotide position.  The goal of the peak caller is to classify, in a binary fashion, each nucleotide of the genome as either bound (contains a peak) or unbound (everything else).  With this binary classification of the entire genome, differences in peak height are discarded and RNA maps based on the fraction of intervals bound at a given position are easier to interpret.

>**Note:** RNA cartography is an evolving analysis.  Some researchers apply smoothing to RNA maps while others choose to display the raw data.  This field is changing rapidly and regarding RNA cartography as an associative visualization tool that is useful for exploratory data anlysis is encouraged.

While RNA maps have found wide adoption in the literature, several strong caveats exist.  First, other features such as splice site strengths and coassociated proteins are not displayed.  Second, multiple categories of exons could exist within the provided exon file that are regulated by distinct mechanisms.  These distinctions are lost in an RNA map and it is left to the analyst to further explore the underlying data.

>**Note:** Strong caution is urged in inferring causal relationships from RNA maps because features relevant to splicing outcomes are not included.  It is advised to use RNA maps as a discovery tool that forms a launchpad for further exploration and mechanistic/functinoal studies.

