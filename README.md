# Download
https://github.com/ljw20180420/rearr/releases

# Install
## Native
It is recommanded to use `rearr` in docker. If you want to install natively, then in the project folder, execute
```bash
sudo ./install.md core sx
```

## Docker
1. Install docker engine [here](https://docs.docker.com/engine/install).
2. Install rootless docker [here](https://docs.docker.com/engine/security/rootless/#install).
3. Configure the proxy of docker daemon if necessary [here](https://docs.docker.com/engine/daemon/proxy).
4. Allow rootless docker to access loopback (e.g. local proxy listening at `localhost`) if necessary [here](https://forums.docker.com/t/no-longer-able-to-access-local-ips-in-rootless-docker-after-update/141890).
5. Login into a tempary container.
```bash
./loginWorker.md
```
6. Just use as natively.

## Setup web server
The web server is containerized by docker. In the project folder, execute
```bash
./compose.md
```

# Usage
See [`rearrTest.md`][`rearrTest.md`].

[`rearrTest.md`]: /sx_lcy/other/rearr-test/

# Documentation
[Here](https://ljw20180420.github.io/sx_lcy/).

# TODO
```[tasklist]
- [ ] Add assumption check.
- [ ] Rename the repo.
- [ ] Rewrite pruned backtracking.
- [ ] Add a simulation for branch-and-bound backtracking.
- [ ] Resemble outputs of previous software
- [ ] Use small genome data to recover genome test.
- [ ] do benchmark based on DWGSIM, Rabema, Seal
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
- [ ] DeepFri
```
