#!/bin/bash

shopt -s expand_aliases

alias ~~~=":<<'~~~bash'"

:<<'~~~bash'

# Usage
## Run test data
```bash
./rearrTest.md
```
## Common
`makeTarget` is the file you want to generate. The underly `make` engine use file extensions to determine which step to run, so the file extension matters. Depending on `makeTarget`, you may need to provide additional parameters and input files.

## Remove duplicates
If you want to remove duplicates for paired (or multiply paired) `fastq` files, then run
```bash
makeTarget=<dataDir>/<name>.noDup \
fastqFiles=<fq1>,<fq2>,... \
./rearrTest.md
```
For more details, see [`removeDuplicates.md`][`removeDuplicates.md`].

## Demultiplex
If you has a file `<dataDir>/<name>.noDup` output by [`removeDuplicates.md`][`removeDuplicates.md`], then run
```bash
makeTarget=<dataDir>/<name>.demultiplex \
spliterIndices=<path1>,<path2>,... \
minScores=<score1>,<scores2>,... \
./rearrTest.md
```
For more details, see [`demultiplex.md`][`demultiplex.md`].

## Align
If you has a file `<dataDir>/<name>.post` ready to align, then run
```bash
makeTarget=<dataDir>/<name>.alg \
refFile=<pathToRefFile> \
correctFile=<pathToCorrectFile> \
./rearrTest.md
```
This will generate the alignment file `<dataDir>/<name>.alg`. Note that the file ready to align must have extension `.post`. These following parameters with defaults control the chimeric alignment.
```bash
s0=-6
s1=4
s2=2
u=-3
v=-9
ru=0
rv=0
qu=0
qv=-5
```
Note that the defaults in this scripts override those in [`rearr`][`rearr`]. If you do not have `correctFile` and your `refFile=<path>.ref`, then you may just set `correctFile=<path>.correct`. `rearrTest.md` will generate `correctFile` with all corrections target `up`. For more details, see [`workFlow.mak`][`workFlow.mak`].

## Post-process by sx module
The output of [`demultiplex.md`][`demultiplex.md`] does not fit the input of [`rearr`][`rearr`]. The transformation between them is highly customized and changes from now and that. For Shi Xing's data, this is done by [`sxCutR2AdapterFilterCumulate.md`][`sxCutR2AdapterFilterCumulate.md`]. If you has a file `<dataDir>/<name>.demultiplex` output by [`demultiplex.md`][`demultiplex.md`], then just run
```bash
makeTarget=<dataDir>/<name>.post \
./rearrTest.md
```
Although the default `minToMapShear=30` works well for Shi Xing's data, you may modifty it as you like.
```bash
makeTarget=<dataDir>/<name>.post
minToMapShear=31 \
./rearrTest.md
```

## Extract spliter by sx module
If you has a csv file `<fullPathToCsvFile>` in the same format as Shi Xing, then you can extract spliters by [`sxExtractSpliter.md`][`sxExtractSpliter.md`].
```bash
makeTarget=<fullPathToCsvFile>.target.fa \
./rearrTest.md
```
Besides `<fullPathToCsvFile>.target.fa`, another file `<fullPathToCsvFile>.pair.fa` will be generated as well. This is because [`sxExtractSpliter.md`][`sxExtractSpliter.md`] always generate both `spliterIndices` simultaneously.

## Get reference by sx module
If you has a csv file `<fullPathToCsvFile>` in the same format as Shi Xing, then you can extract spliters by [`getSxCsvFileRef.md`][`getSxCsvFileRef.md`].
```bash
makeTarget=<fullPathToCsvFile>.ref \
genome=<pathToGenome> \
bowtie2index=<pathToGenomeIndex> \
./rearrTest.md
```
If you has run [`loginWorker.md`][`loginWorker.md`], and specify both `genome` and `bowtie2index`, then the default settings will work for you.
```bash
genome=genome/genome.fa
bowtie2index=genome/genome
```
So just run
```bash
makeTarget=<fullPathToCsvFile>.ref \
./rearrTest.md
```

## Full workflow
Assume that your design is the same as Shi Xing.
```
Then just run
```bash
makeTarget=<dataDir>/<name>.alg \
fastqFiles=<pathToR2>,<pathToR1> \
spliterIndices=<pathToCsvFile>.target.fa,<pathToCsvFile>.pair.fa \
genome=<pathToGenome> \
bowtie2index=<pathToGenomeIndex> \
refFile=<pathToCsvFile>.ref \
correctFile=<pathToCsvFile>.correct \
./rearrTest.md
```
`rearrTest.md` will run all steps above for you to generate `<dataDir>/<name>.alg`. You may also try to modify the following parameters with defaults for better results.
```bash
minScores=30,100
s0=-6
s1=4
s2=2
u=-3
v=-9
ru=0
rv=0
qu=0
qv=-5
minToMapShear=30
```

# Introduction
This scripts integrate all steps (remove duplicates, demultiplex, alignment and so on). However, why not just use the corresponding script to run each step? If you run test data for the second time, then `rearrTest.md` will not do anythin because the underlying `make` engine is smart enough to skip the updating of the outputs when no change is detected in the inputs. Thus, the reason to use `rearrTest.md` is that it may skip some duplicated computations for you. To actually run the test again, one need delete the previous results first.

Another reason to use `rearrTest.md` is that it hides two cumbersome steps from the user. For example, if `spliterIndices` used in [`demultiplex.md`][`demultiplex.md`] is not indexed by `bowtie2`, then `rearrTest.md` will do this silently. Also, if `correctFile` has the same path and name as `refFile` (see [`rearr`][`rearr`]), but with the file extension `.correct` instead of `.ref`, and `rearrTest.md` cannot find `correctFile` on the file system, then it will generate a default `correctFile` for you with all fields filled with `up`.

[`rearr`]: /sx_lcy/core/rearr/
[`removeDuplicates.md`]: /sx_lcy/core/remove-duplicates/
[`demultiplex.md`]: /sx_lcy/core/demultiplex/
[`workFlow.mak`]: /sx_lcy/other/rearr-test/work-flow/
[`sxCutR2AdapterFilterCumulate.md`]: /sx_lcy/sx/sx-cut-r2-adapter-filter-cumulate/
[`sxExtractSpliter.md`]: /sx_lcy/sx/sx-extract-spliter/
[`getSxCsvFileRef.md`]: /sx_lcy/sx/get-sx-csvfile-ref/
[`loginWorker.md`]: /sx_lcy/other/login-worker/

# Source
~~~bash
# The following parameters should be replaced.
makeTarget=${makeTarget:-test/rearr.alg}
fastqFiles=${fastqFiles:-test/A2-g1n-3.R2.fq.gz,test/A2-g1n-3.fq.gz}
spliterIndices=${spliterIndices:-sx/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.target.fa,sx/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.pair.fa}
minScores=${minScores:-30,100}

minToMapShear=${minToMapShear:-30}
refFile=${refFile:-sx/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.ref}
correctFile=${correctFile:-sx/csvfiles/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.correct}
ext1up=${ext1up:-50}
ext1down=${ext1down:-0}
ext2up=${ext2up:-10}
ext2down=${ext2down:-100}

# The following parameters are default in most cases.
genome=${genome:-genome/genome.fa}
bowtie2index=${bowtie2index:-genome/genome}
s0=${s0:--6}
s1=${s1:-4}
s2=${s2:-2}
u=${u:--3}
v=${v:--9}
ru=${ru:-0}
rv=${rv:-0}
qu=${qu:-0}
qv=${qv:--5}

make -f workFlow.mak $makeTarget fastqFiles=$fastqFiles spliterIndices=$spliterIndices minScores=$minScores genome=$genome bowtie2index=$bowtie2index refFile=$refFile correctFile=$correctFile s0=$s0 s1=$s1 s2=$s2 u=$u v=$v ru=$ru rv=$rv qu=$qu qv=$qv minToMapShear=$minToMapShear
~~~

~~~bash
alias ~~~=":" # This suppresses a warning and is not part of source.
~~~
