# Install

## Conda

```console
$ conda install bioconda::rearr
```

## Container

Each bioconda package has a [biocontainer docker image](https://quay.io/repository/biocontainers/rearr) and a [galaxy singularity image](https://depot.galaxyproject.org/singularity). As far as I know, the easiest way to use these images is `apptainer`. Install `apptainer` by conda.
```console
conda install conda-forge::apptainer
```

To get an interactive shell environment,
```console
$ apptainer shell docker://quay.io/biocontainers/rearr:tag
Apptainer> rearrangement -h
```
You should replace `tag` by the latest working image tag found at [biocontainer docker image](https://quay.io/repository/biocontainers/rearr). For example,
```console
$ apptainer shell docker://quay.io/biocontainers/rearr:1.0.3--h9948957_0
```

To run commands non-interactively,
```
$ apptainer run docker://quay.io/biocontainers/rearr:tag rearrangement -h
```
This is handy when you want to build a pipeline depends on multiple containers.

Containers of `apptainer` share io and network with the host. The current working directory is where you invoke the `apptainer` command like `shell` or `run`.

## WebUI

The webUI is heavy, so not suitable to release as conda package. There is a prebuild [github docker package](https://ghcr.io/ljw20180420/rearr). To launch the webUI locally, either clone the repository
```console
$ git clone https://github.com/ljw20180420/rearr.git
```
or download the latest working [release](https://github.com/ljw20180420/rearr/releases). Then in the project folder, invoke
```console
$ ./compose.sh
```
The webUI is served at `http://ocalhost:80`. The first running of `./compose.sh` will pull necessary images including the prebuild [github docker package](https://ghcr.io/ljw20180420/rearr).

# Documentation

[Here](https://ljw20180420.github.io/rearr/).

# TODO
```[tasklist]
- [ ] Add homepage for qiangwulab.sjtu.edu.cn.
- [ ] improve BWT-SW and apply it to the sgRNA library demultiplex and genome-wide CRISPR
    - [ ] regex DFA and NFA
    - [ ] SIMD
    - [ ] prune resembling repeat alignment
    - [ ] statistic algorithm
    - [ ] decrease memory usage by record DNA in 2bit
    - [ ] load suffix array into memory
- [ ] list bad examples in benchmark
- [ ] Write core functions in c++
- [ ] Add assumption check.
- [ ] Rewrite pruned backtracking.
- [ ] Add a simulation for branch-and-bound backtracking.
- [ ] Resemble outputs of previous software
- [ ] Use small genome data to recover genome test.
- [ ] compare with divide-and-conquer
- [ ] modulerize demultiplex
- [ ] add manim to show the core algorithm of rearr
- [ ] add method to search genome-wide off-target for sequences not found in sgRNA libraries
- [ ] add benchmark to demultiplex
- [ ] github action to deploy to qiangwulab.sjtu.edu.cn
- [ ] add github wiki
- [ ] other CDCI support by github (git -> github cli -> github docs|skills|support|community)
- [ ] use functools LRU-cache to speed up python code
- [ ] add benchmark for SIQ: https://github.com/RobinVanSchendel/SIQ
- [ ] add benchmark for PEM-Q: https://github.com/liumz93/PEM-Q
- [ ] convert alg to sam
- [ ] Deploy to JCloud. Celery flower does not work properly on JCloud. Maybe permission problem.
- [ ] asgi is more advance than wsgi
- [ ] resemble indelphi to explicitly list sequences
- [ ] implement tidymodels (need to install tidymodels in shiny rocker, not install by apt)
- [ ] CDN
- [ ] Use explicit base in shiny app microHomology
- [ ] Use probability language to inplement Gibbs sampling for predicting the frequencies of blunt end cleavage events
```

# TODO (Long term)
```[tasklist]
- [ ] Use GNU autotools to install Rearrangement
```
