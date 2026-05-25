#!/bin/bash

shopt -s expand_aliases

alias ~~~=":<<'~~~bash'"

:<<'~~~bash'

# Usage

```shell
param1=value1 param2=value2 ... runWorkFlow.sh [options]
```

<pre class="mermaid">
---
title: runWorkflow.sh
---
flowchart TD
    REF[(
        <h1>reference_file</h1>
        <table>
            <tr>
                <th>start1</th>
                <th>ref1</th>
                <th>end1</th>
                <th>start2</th>
                <th>ref2</th>
                <th>end2</th>
            </tr>
            <tr>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
            </tr>
        </table>
    )]

    MK1[(marker1)] --> DM
    MK2[(marker2)] --> DM

    R1[(fastqR1)] --> RD[removeDuplicates.sh]
    R2[(faqstR2)] --> RD
    RD --> UNIQUE[(
        <h1>removeDuplicates_file</h1>
        <table>
            <tr>
                <th>R1</th>
                <th>R2</th>
                <th>#</th>
            </tr>
            <tr>
                <td>...</td>
                <td>...</td>
                <td>...</td>
            </tr>
        </table>
    )]

    UNIQUE --> DM[demultiplex.sh]
    SCORE[(
        <h1>minScores</h1>
        <table>
            <tr>
                <th>score1</th>
                <th>score2</th>
            </tr>
            <tr>
                <td>...</td>
                <td>...</td>
            </tr>
        </table>
    )]
    SCORE --> DM
    DM --> ONTARGET[(
        <h1>demultiplex_file</h1>
        <table>
            <tr>
                <th>R1</th>
                <th>R2</th>
                <th>#</th>
                <th>id</th>
                <th>rstart1</th>
                <th>rend1</th>
                <th>qstart1</th>
                <th>qend1</th>
                <th>rstart2</th>
                <th>rend2</th>
                <th>qstart2</th>
                <th>qend2</th>
            </tr>
            <tr>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
                <td>...</td>
            </tr>
        </table>
    )]

    ONTARGET --> QUERY[(
        <h1>input_file</h1>
        <table>
            <tr>
                <th>query</th>
                <th>#</th>
                <th>id</th>
            </tr>
            <tr>
                <td>...</td>
                <td>...</td>
                <td>...</td>
            </tr>
        </table>
    )] --> REARR[rearrangement]
    REF --> REARR
    REARR --> ALG[(
        <h1>rearrangement_file</h1>
        <table>
            <tr>
                <td>idx</td>
                <td>#</td>
                <td>score</td>
                <td>id</td>
            </tr>
            <tr>
                <td>ref1</td>
                <td>ref2</td>
            </tr>
            <tr>
                <td>query</td>
            </tr>
        </table>
    )]

    REF --> CMH[correct_micro_homology.awk]
    DIRECTION[(
        <h1>direction_file</h1>
        <table>
            <tr>
                <th>up/down</th>
            </tr>
            <tr>
                <td>...</td>
            </tr>
        </table>
    )] --> CMH
    ALG --> CMH
    CMH --> CORRECTED[(
        <h1>correct_micro_homology_file</h1>
        <table>
            <tr>
                <td>idx</td>
                <td>#</td>
                <td>score</td>
                <td>id</td>
                <td>udangle</td>
                <td>rstart1</td>
                <td>qstart1</td>
                <td>rend1</td>
                <td>qend1</td>
                <td>random</td>
                <td>rstart2</td>
                <td>qstart2</td>
                <td>rend2</td>
                <td>qend2</td>
                <td>ddangle</td>
                <td>cut1</td>
                <td>ref1+cut2</td>
            </tr>
            <tr>
                <td>ref1</td>
                <td>ref2</td>
            </tr>
            <tr>
                <td>query</td>
            </tr>
        </table>
    )]
</pre>

- `options` are passed to the underlying `make` calling. `makeTarget` is the file you want to generate. The underly `make` engine use file extensions to determine which step to run, so the file extension matters. Depending on `makeTarget`, you may need to provide additional parameters and input files.
- To remove duplicates for paired (or multiply paired) `fastq` files, run
```shell
makeTarget=removeDuplicates_file.noDup \
fastqFiles=fastqR1,fastqR2,... \
runWorkFlow.sh
```
For more details, see [`removeDuplicates.sh`][removeDuplicates.sh.md].
- To demultiplex `removeDuplicates_file.noDup`, run
```shell
makeTarget=demultiplex_file.demultiplex \
markerIndices=marker1,marker2,... \
minScores=score1,scores2,... \
runWorkFlow.sh
```
For more details, see [`demultiplex.sh`][demultiplex.sh.md]. If `marker` is not indexed by `bowtie2`, `runWorkFlow.sh` will index it silently
- To align `query.post` to reference and correct microhomology, run
```shell
makeTarget=correct_micro_homology_file.alg \
refFile=reference_file \
directionFile=direction_file \
runWorkFlow.sh
```
Chimeric alignment scores used by `rearrangement` can be set as follows.
```shell
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
For more details, see [core part of rearr][core.md]. If only `refFile` is provided, a default `directionFile=${refFile}.direct` will be created with all `up`. For more details, see [`workFlow.mak`][workFlow.mak.md].
- The output of [`demultiplex.sh`][demultiplex.sh.md] does not fit the input of [core part of rearr][core.md]. The transformation between them is highly dependent on the design of experiment and changes from now and that. By default, [`workFlow.mak`][workFlow.mak.md] use the first column of the output of [`demultiplex.sh`][demultiplex.sh.md] as query.
- To run the full workflow for the test data (the test data put fastqR2 before fastqR1),
```shell
makeTarget=correct_micro_homology_file.alg \
fastqFiles=fastqR2,fastqR1 \
markerIndices=pasmid_file.target.fa,pasmid_file.pair.fa \
refFile=plasmid_file.ref \
./runWorkFlow.sh
```
`runWorkFlow.sh` will run all steps above for you to generate `correct_micro_homology_file.alg`.
- `runWorkFlow.sh` use make engine, which skips the updating of the outputs if no change is detected in the inputs necesary to generate that output. This saves computations for you.

[core.md]: ./core.html
[removeDuplicates.sh.md]: ./auxilary/removeDuplicates.sh.html
[demultiplex.sh.md]: ./auxilary/demultiplex.sh.html
[workFlow.mak.md]: ./runWorkFlow/workFlow.mak.html

# Source
~~~bash
# The following parameters should be replaced.
makeTarget=${makeTarget:-test/test_work_flow/rearr.alg}
fastqFiles=${fastqFiles:-test/test_work_flow/A2-g1n-3.R2.fq.gz,test/test_work_flow/A2-g1n-3.fq.gz}
markerIndices=${markerIndices:-test/test_work_flow/target.fa,test/test_work_flow/pair.fa}
minScores=${minScores:-30,100}

refFile=${refFile:-test/test_work_flow/ref}
directionFile=${directionFile:-"${refFile}.direct"}

# The following parameters are default in most cases.
s0=${s0:--6}
s1=${s1:-4}
s2=${s2:-2}
u=${u:--3}
v=${v:--9}
ru=${ru:-0}
rv=${rv:-0}
qu=${qu:-0}
qv=${qv:--5}

default_prefix=${CONDA_PREFIX:-"${HOME}/.local"}
prefix=${prefix:-"${default_prefix}"}
if [ -f "workFlow.mak" ]
then
    make_file="workFlow.mak"
else
    make_file="${prefix}/share/rearr/workFlow.mak"
fi

make $@ -d -f "${make_file}" "${makeTarget}" \
    fastqFiles="${fastqFiles}" \
    markerIndices="${markerIndices}" \
    minScores="${minScores}" \
    refFile="${refFile}" \
    directionFile="${directionFile}" \
    s0="${s0}" \
    s1="${s1}" \
    s2="${s2}" \
    u="${u}" \
    v="${v}" \
    ru="${ru}" \
    rv="${rv}" \
    qu="${qu}" \
    qv="${qv}"
~~~

~~~bash
alias ~~~=":" # This suppresses a warning and is not part of source.
~~~
