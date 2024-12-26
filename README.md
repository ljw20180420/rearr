# Download
https://github.com/ljw20180420/sx_lcy/releases

# Install
It is recommanded to use rearr in docker. If you want to install natively, cd into the project fold and execute `./install.md` core sx.

# Install rootless docker
First, install docker engine: https://docs.docker.com/engine/install
Then install rootless docker: https://docs.docker.com/engine/security/rootless/#install

# Usage
See `rearrTest.md`.
If you use docker, first login into docker.
```bash
./loginWorker.md
```
Then just use as native.

# Output
```{list}
target.count: file of target and pair without duplicates
target.demultiplex: file after demultiplex
target.post: file ready to align
target.alg: alignments
```

# Setup the docker based server
First, you need docker: https://docs.docker.com/engine/install.
Maybe you need to configure the proxy of docker daemon: https://docs.docker.com/engine/daemon/proxy.
To access host loopback in rootless docker: https://forums.docker.com/t/no-longer-able-to-access-local-ips-in-rootless-docker-after-update/141890.
Docker buildx does not respect the daemon proxy. One has to use system proxy, say
```bash
HTTPS_PROXY=socks5://127.0.0.1:1080 docker compose build
```
```bash
cd sx_lcy
./compose.md
```

# TODO
```[tasklist]
- [ ] modify flask
- [ ] modify shiny
- [ ] document all codes
- [ ] use github pages (classic) to host documents
- [ ] use Doxygen to generate documents from code (Automatic document)
- [ ] use github action to host github pages by run Doxygen
- [ ] add github wiki
- [ ] github action for build containers
- [ ] other CDCI support by github (git -> github cli -> github docs|skills|support|community)
- [ ] add benchmark for SIQ: https://github.com/RobinVanSchendel/SIQ
- [ ] convert alg to sam
- [ ] unittest
- [ ] kvm
- [ ] javascript -> html -> css -> npm -> typescript -> tamper monkey -> selenium
- [ ] vue -> uniapp -> tauri2
- [ ] Deploy to JCloud. Celery flower does not work properly on JCloud. Maybe permission problem.
- [ ] asgi is more advance than wsgi
- [ ] resemble indelphi
- [ ] implement tidymodels (need to install tidymodels in shiny rocker, which must not be installed through CRAN)
- [ ] CDN
- [ ] Add 3D structure prediction shiny App
- [ ] Add large language model for DNA
- [ ] Use explicit base in shiny app microHomology
- [ ] Use probability language to inplement Gibbs sampling for predicting the frequencies of blunt end cleavage events
```

# TODO (Long term)
```[tasklist]
- [ ] 4C normalization algorithm
- [ ] Use GNU autotools to install Rearrangement
- [ ] Hi-C apps
- [ ] DeepFri
- [ ] PDB structure prediction
- [ ] molecular dynamics simulation
- [ ] Call TADs
```