
# CRISPR-A documentation
CRISPR-A is a computational pipeline that helps you from design to analysis of your genome editing experiments. CRISPR-A simulates gene editing outcomes of a target sequence, and assesses the quality of gene editing experiments using next generation sequencing (NGS) data. A web implementation is available at: <http://synbio.upf.edu/crispr-a/>.

#### Which can of data can you analysed with CRISPR-A?
Targeted Next Generation Sequencing data, preferentially from sequencing by synthesis like Illumina. 

#### Which can of experiments I can analyse?
- Double strand break (DSB) with CRISPR-based technologies
- Deletions produced with double nicking strategies
- Base editing (BE)
- Prime editing (PE)
- Homology directed repair (HDR)
- Predicted off-target sites 

In summary, any kind of edit small enough to be characterized through targeted sequencing 

#### Which is the recomended coverage?
- For clonal populations it is recommended to have at least 500 reads
- For cell bulks it is recommended to have at least 2000 reads

Even though, it will highly depend on the quality of the sample, complexity of the amplicon sequence, diversity of outcomes...

#### What does CRISPR-A analysis pipeline?
The analysis algorithm is composed of three mandatory steps: 

1) reads pre-processing for quality assessment, 

2) reads alignment against reference amplicon, and 

3) edit calling. 

There are other optional processes: UMI clustering, reference discovery, size bias correction, and noise subtraction based on an empirical model from negative control samples. 

This is a versatile tool that can be used in a broad list of genome editing experiments (Base Edinting, Prime Editing, CRISPR-based Doble Strand Break, Homology directed repair...) without the need of specifying much information. 

## CRISPR-A processing
#### Quality filtering
First, in case of paired-end sequencing reads, these are merged with PEAR to increase quality at the ends of the sequences. Then, input reads are filtered based on the quality score (phred33) in order to remove potentially false positive indels using fastq_quality_trimmer.
#### Adapter trimming
Next, adapters are searched and trimmed from the reads. The adapters are trimmed from the reads using FastQC and cutadapt. 
#### Read merging
If paired-end reads are provided, reads are merged using FLASh . This produces a single read for alignment to the amplicon sequence, and reduces sequencing errors that may be present at the end of sequencing reads.
#### Clustering by Uni-molecular identifiers (OPTIONAL)
Reads can be clustered by uni-molecular identifiers, if these barcodes are added in a first PCR amplification of just two cycles prior to standard PCR amplification and sample barcodes attachment. This process will increase precision by avoiding miss quantification of sequencing biases. When this process is run, cluster consensus sequences are used in the following steps. 
#### Self-referecens discovery (OPTIONAL)
If a genome is given and there is no amplicon sequence, the amplicon reference sequences is obtained by aligning the reads against the reference genome.  
#### Flip amplicon sequence into gRNA orientation
To make all results comparable, amplicon sequence is written in the same orientation as the protospacer prior to alignment. 
#### Alignment
The preprocessed reads are then aligned to the amplicon reference sequence with a global sequence alignment algorithm that has been optimized to take into account the biological knowladge of doble strand break (DBS) repair. 
#### Variant or edits calling
From the alignement, the sequences are classified and inspected infunction of found indels, delitions and insertions, and substitutions. In the case of deletions, in-frame and out-of-frame information is also retrived. 
#### Noise cleaning (OPTIONAL)
Noise subtraction based on an empirical model from negative control samples.
#### Visualization and analysis
Finally, a set of interactive plots are generated to visualize the position and type or outcomes within the amplicon sequence. Tables with pre-processing information and a IGV browser to inspect the alignment are also build. The de-noise processing of the data is also shown in a visual way. 

## Installation
#### Requirements
- nextflow version 21.10.6
- Singularity

To install nextflow, follow this instructions:
```
wget -qO- nextflow-21.10.6-all
chmod +x nextflow
```
Finally, you can move the nextflow file to a directory accessible by your $PATH variable to avoid writing the full path each time you need to run nextflow.
```
nextflow -version # check nextflow is correctly installe din the required version (21.10.6)
```
To install Singularity you can follow this instructions:
```
sudo apt-get update
sudo apt-get install -y build-essential libseccomp-dev pkg-config squashfs-tools cryptsetup curl wget git
export GOVERSION=1.17.3 OS=linux ARCH=amd64  # change this as you need
wget -O /tmp/go${GOVERSION}.${OS}-${ARCH}.tar.gz   https://dl.google.com/go/go${GOVERSION}.${OS}-${ARCH}.tar.gz
sudo tar -C /usr/local -xzf /tmp/go${GOVERSION}.${OS}-${ARCH}.tar.gz
echo 'export PATH=$PATH:/usr/local/go/bin' >> ~/.bashrc
source ~/.bashrc
curl -sSfL https://raw.githubusercontent.com/golangci/golangci-lint/master/install.sh | sh -s -- -b $(go env GOPATH)/bin v1.43.0
echo 'export PATH=$PATH:$(go env GOPATH)/bin' >> ~/.bashrc
source ~/.bashrc
git clone https://github.com/hpcng/singularity.git
cd singularity/
./mconfig
cd ./builddir
make
sudo make install
singularity --version # check singularity is installed
```

#### CRISPR-A nextflow pipeline
Clone the repository with CRISPR-A nextflow pipeline
```
git clone https://your-user@bitbucket.org/synbiolab/crispr-a_nextflow.git
cd crispr-a_nextflow
```
Check that all scripts of the pipeline are executable.
```
chmod -R a+x bin/*
```
#### Testing 
```
nextflow run crispr-a.nf -c nf_input-example.config
```
First time, it will take a while, since a Singularity image with all requirements will be downloaded. Next times you can determine the folder with the Singularity image with the variable $SINGULARITY_IMAGES. Paths of the config file are relative to ./crispr-a_nextflow. If you run it from anothe folder, you should change the paths. 

## CRISPR-A usage
This is a nextflow pipeline that can be run with the nextflow script and the following parameters:
- outDir: path to save all results

- input: path to fastq or fastq.gz files

- r2file: boolean; true when apir-end sequencing technologies have neen used and false for single-end sequencing samples

- inSize: -1

- umiClustering: true when UMI proocol is used and false for standard sequencing amplification

- umiType: "TV" if nanopore TV UMI is used. Otherwise sequnce has to be indicated. 

- minlen: UMI minimal length  (default: 55)

- maxlen: UMI maximal length (default: 57)

- id: minimal identity to cluster to UMI sequences together (default: 0.99)

- ubs:  UMI bin size; minimal number of sequences in a cluster to consider its consensus sequences in the following steps (default: 1)

- error: allowe error rate in UMI clustering (default: 0.2)

- simGE: boolean to indicate if simulated data has to be generated for the analysis. It should be false if data is provided for the analysis 

- simGEfilesPath: path with executable scripts used to simulate data (default: "./bin/SimGE")

- simGEid: name of the simulated samples (default: "sample_0")

- simGEsequence: target wild type sequence with which simulate multiple events of edition

- cutSite: relative position, from the beggining of simGEsequence, where cut sit should take place in the simulation

- referenceFasta: path where all reference single fasta files are found. The name of the file has to be composed by the name of the sample and terminate with fasta or fa.  When it is not 
used, emply files have to be provided

- indexHuman: path to human reference genome and index files

- indexMouse:  path to mouse reference genome and index files

- refOrganism: list with all samples and the reference organisme (possible values: "human", "mouse", "other")

- gRNAseq: list with sample name and used protospacer

- selfRef: list of booleans to indicate for each sample if amplicon reference has to be found in reference genome

- aligner: aligner to use (default: "minimap2")

- template_seq: path where all template reference (HDR, PE, BE..) single fasta files are found. The name of the file has to be composed by the name of the sample and terminate with fasta 
or fa. When it is not used, emply files have to be provided

- mock: boolean to indicate if negative control is provided for noise cleaning

- template: boolean to indicate if template sequence has to be used  = false

- spike: value "no" or "yes" to indicate if spike-in model has to be used to correct size sequencing bias  

- samplesNames_file: php list to connect ids with sample names (default: "samples.php") 

- protCutSite: list indicating the position of the cut relative to the protospacer sequence

## Troubleshooting
- Check nextflow version. If the version is not the correct one the pipeline will not work properly.
- Paths and inputs have to be given as required. Check everything is as it should be. 
- Is the singularity image with all software and packages downloaded correctly? 

## Reference
Sanvicente-García, M., García-Valiente, A., Jouide, S., Jaraba-Wallace, J., Bautista, E., Escobosa, M., Sánchez-Mejías, A., & Güell, M. (2023). CRISPR-Analytics (CRISPR-A): A platform for precise analytics and simulations for gene editing. PLOS Computational Biology, 19(5), e1011137.
