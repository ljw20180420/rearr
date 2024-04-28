# Download
```{bash}
git clone git@github.com:ljw20180420/sx_lcy.git
git clone https://github.com/ljw20180420/sx_lcy.git
```

# Install
```{bash}
cd sx_lcy
./install.sh
```

# Usage
```{bash}
rearr_run.sh path_to_fastq/fq[.gz] ref1 ref2 cut1 cut2 NGGCCNtype1 NGGCCNtype2
```

# Parameters
```{list}
ref1: reference string of locus 1
ref2: reference string of locus 2 (for single cut, ref2 = ref1)
cut1: cut point for ref1
cut2: cut point for ref2
NGGCCNtype1: ref1 is on the NGG strand or the CCN strand
NGGCCNtype2: ref2 is on the NGG strand or the CCN strand
```

# Output
```{list}
fastq.count: the count of duplicates in the fastq file
fastq.alg.cut1.cut2.ref1len: the alignments
fastq.table.cut1.cut2.ref1len: the summarization table for indel information
```

# TODO
```[tasklist]
- [ ] resemble indelphi
- [ ] implement tidymodels (need to install tidymodels in shiny rocker, which must not be installed through CRAN)
- [ ] merge shiny APP visualize alg file into single APP with multiple tabPanels in narbarPage
- [ ] explore base
- [ ] feature choose
- [ ] AI review
- [ ] PCR
- [X] Alignment browser with folded insertion lines
- [ ] CDN
- [ ] Single page webUI (reply results to the page instead of saving it on the server)
- [ ] Rewrite post_process in R
- [ ] Remove python code from draw_figures.Rmd
- [ ] Sequence heatmap, and more visualizations
- [ ] Make single target version a special case of screening case (merge normal mode with barcode mode)
- [ ] Modulerize native installation
- [ ] Substitute native installation into Dockerfiles
- [ ] Optimize DockerFiles by merge install commands
- [ ] Move sed/awk/perl codes embedded in bash scripts to individual files for easier debugging and reusing
- [ ] stage docker build process (say a celery image can be used for both flask image and worker image)
- [ ] Minimize docker image size, pre-build images and deployed to docker hub
- [ ] Add nginx
- [ ] Add waitress
- [ ] Add tls
- [ ] Install necessary dependencies to shiny server
- [ ] Organize the whole project to a more extensible form (easily add now function both natively, in docker, or online)
- [ ] Add 3D structure prediction shiny App
- [ ] Add large language model for DNA
- [ ] Use explicit base in shiny app microHomology
- [ ] Add shiny app to predict indel events
- [ ] Add shiny app for kpLogo
- [X] Avoid using system command in shiny app for better portability
- [ ] Deploy to JCloud
- [ ] Use probability language to inplement Gibbs sampling for predicting the frequencies of blunt end cleavage events
- [ ] Move all downstream analyses to shiny App
- [ ] Automatic document
- [ ] Build a docker with only binary program, awk, sed, perl, bash scripts
- [ ] Python program all run in shiny, so only need a python shiny docker
- [ ] R program all run in shiny, so only need rocker/shiny-verse image
```

# TODO (Long term)
```[tasklist]
- [ ] Use GNU autotools to install Rearrangement
- [ ] Hi-C apps
- [ ] Uniprot database
- [ ] RNA structure prediction
- [ ] Protein structure prediction (Alphafold2)
- [ ] DeepFri
- [ ] PDB structure prediction
- [ ] molecular dynamics simulation
- [ ] Call TADs
```