# Download
```{bash}
git clone git@github.com:ljw20180420/sx_lcy.git
git clone https://github.com/ljw20180420/sx_lcy.git
```

# Install
```{bash}
It is recommanded to use rearr in docker. If you want to install natively, cd into the project fold and execute ./install/sh core sx
```

# Usage
```{bash}
See rearrTest.sh
If you use docker, first login into a docker by ./loginWorker.sh dataPath. dataPath will be mounted to /app/data in docker. Then just use as native.
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
target.count: file of target and pair without duplicates
target.demultiplex: file after demultiplex
target.post: file ready to align
target.alg: alignments
```

# TODO
```[tasklist]
- [ ] add a link to genome in install.sh
- [ ] add a snakemake workflow
- [ ] add CDCI support by github action
- [ ] resemble indelphi
- [ ] implement tidymodels (need to install tidymodels in shiny rocker, which must not be installed through CRAN)
- [ ] merge shiny APP visualize alg file into single APP with multiple tabPanels in narbarPage
- [ ] explore base
- [ ] feature choose
- [ ] AI review
- [ ] PCR
- [ ] CDN
- [ ] Rewrite post_process in R
- [ ] Sequence heatmap, and more visualizations
- [ ] Make single target version a special case of screening case (merge normal mode with barcode mode)
- [ ] Modulerize native installation
- [ ] Substitute native installation into Dockerfiles
- [ ] Optimize DockerFiles by merge install commands
- [ ] Move sed/awk/perl codes embedded in bash scripts to individual files for easier debugging and reusing
- [ ] stage docker build process (say a celery image can be used for both flask image and worker image)
- [ ] Add tls
- [ ] Install necessary dependencies to shiny server
- [ ] Organize the whole project to a more extensible form (easily add now function both natively, in docker, or online)
- [ ] Add 3D structure prediction shiny App
- [ ] Add large language model for DNA
- [ ] Use explicit base in shiny app microHomology
- [ ] Add shiny app to predict indel events
- [ ] Add shiny app for kpLogo
- [ ] Deploy to JCloud
- [ ] Use probability language to inplement Gibbs sampling for predicting the frequencies of blunt end cleavage events
- [ ] Move all downstream analyses to shiny App
- [ ] Automatic document
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