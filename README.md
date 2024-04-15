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

# Tools
Display the first readnum (default: 50) read's alignments in ANSI fomrat, the cut points and random insertions are aligned
```{bash}
rearr_view.sh fastq.alg.cut1.cut2.ref1len [readnum]
```
Generate the html report
```{bash}
rearr_render.r fastq.alg.cut1.cut2.ref1len
```

# TODO
```[tasklist]
- [ ] Alignment browser with folded insertion lines
- [ ] Shiny interactive plots
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
- [ ] Use probability language to inplement Gibbs sampling for predicting the frequencies of blunt end cleavage events
- [ ] Move all downstream analyses to shiny App
```

# TODO (Long term)
```[tasklist]
- [ ] Use GNU autotools to install Rearrangement
```