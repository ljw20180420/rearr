#!/bin/bash

shopt -s expand_aliases

alias ~~~=":<<'~~~bash'"

:<<'~~~bash'

# Usage

```bash
param1=value1 param2=value2 ... runWorkFlow.sh [options]
```

<pre class="mermaid">
---
title: runWorkflow.sh
---
flowchart TD
    PF[(
        <h1>plasmid_file</h1>
        <table>
            <tr>
                <th>adapter#0040;20bp#0041; + sgRNA#0040;20bp#0041; + scaffold#0040;83/93bp#0041; + query#0040;44bp#0041; + 3bp + RCbarcode#0040;18bp#0041; + RCprimer#0040;21bp#0041;</th>
            </tr>
            <tr>
                <td>...</td>
            </tr>
        </table>
    )] --> GSPFR[getSxPlasmidFileRef.sh]
    GF[(<h1>genome_file</h1>)] --> GSPFR
    BI[(<h1>bowtie2index</h1>)] --> GSPFR
    EXT[(
        <h1>extentions</h1>
        <table>
            <tr>
                <td>ext1up</td>
                <td>ext1down</td>
                <td>ext2up</td>
                <td>ext2down</td>
            </tr>
        </table>
    )] --> GSPFR
    GSPFR --> REF[(
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

    PF --> SEM[sxExtractMarker.sh]
    SEM --> MK1[(
        <h1>stdout</h1>
        <table>
            <tr>
                <th>primer#0040;21bp#0041; + barcode#0040;18bp#0041;</th>
            </tr>
            <tr>
                <td>...</td>
            </tr>
        </table>
    )]
    SEM --> MK2[(
        <h1>fd3</h1>
        <table>
            <tr>
                <th>adapter#0040;20bp#0041; + sgRNA#0040;20bp#0041; + scaffold#0040;83/93bp#0041;</th>
            </tr>
            <tr>
                <td>...</td>
            </tr>
        </table>
    )]

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
    MK1 --> DM
    MK2 --> DM
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

    ONTARGET --> sxCRAFC[sxCutR2AdapterFilterCumulate.sh] --> QUERY[(
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
```bash
makeTarget=removeDuplicates_file.noDup \
fastqFiles=fastqR1,fastqR2,... \
runWorkFlow.sh
```
For more details, see [`removeDuplicates.sh`][removeDuplicates.sh.md].
- To demultiplex `removeDuplicates_file.noDup`, run
```bash
makeTarget=demultiplex_file.demultiplex \
markerIndices=marker1,marker2,... \
minScores=score1,scores2,... \
runWorkFlow.sh
```
For more details, see [`demultiplex.sh`][demultiplex.sh.md]. If `marker` is not indexed by `bowtie2`, `runWorkFlow.sh` will index it silently
- To align `query.post` to reference and correct microhomology, run
```bash
makeTarget=correct_micro_homology_file.alg \
refFile=reference_file \
directionFile=direction_file \
runWorkFlow.sh
```
Chimeric alignment scores used by `rearrangement` can be set as follows.
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
For more details, see [core part of rearr][core.md]. If only `refFile` is provided, a default `directionFile=${refFile}.direct` will be created with all `up`. For more details, see [`workFlow.mak`][workFlow.mak.md].
- The output of [`demultiplex.sh`][demultiplex.sh.md] does not fit the input of [core part of rearr][core.md]. The transformation between them is highly dependent on the design of experiment and changes from now and that. For out in-house data, this is done by [`sxCutR2AdapterFilterCumulate.sh`][sxCutR2AdapterFilterCumulate.sh.md] as follows.
```bash
makeTarget=query.post \
minToMapShear=30 \
./runWorkFlow.sh
```
- Our in-house data use plasmids in a `plasmid_file`. We extract demultiplex markers from those plasmids by [`sxExtractMarker.sh`][sxExtractMarker.sh.md].
```bash
makeTarget=plasmid_file.target.fa \
./runWorkFlow.sh
```
Besides `plasmid_file.target.fa` used as demutiplex marker for `R2`, another file `plasmid_file.pair.fa` will be generated as well used as demutiplex marker for `R1`.
- The `plasmid_file` also contain reference sequences (sgRNAs). These references are extract by [`getSxPlasmidFileRef.sh`][getSxPlasmidFileRef.sh.md].
```bash
makeTarget=plasmid_file.ref \
genome=genome_file \
bowtie2index=bowtie2index_prefix \
./runWorkFlow.sh
```
Our in-house data use hg19.
- To run the full workflow for our in-house data (our in-house data put fastqR2 before fastqR1),
```bash
makeTarget=correct_micro_homology_file.alg \
fastqFiles=fastqR2,fastqR1 \
markerIndices=pasmid_file.target.fa,pasmid_file.pair.fa \
genome=genome_file \
bowtie2index=bowtie2index_prefix \
refFile=plasmid_file.ref \
./runWorkFlow.sh
```
`runWorkFlow.sh` will run all steps above for you to generate `correct_micro_homology_file.alg`.
- `runWorkFlow.sh` use make engine, which skips the updating of the outputs if no change is detected in the inputs necesary to generate that output. This saves computations for you.

[core.md]: ./core.html
[removeDuplicates.sh.md]: ./auxilary/removeDuplicates.sh.html
[demultiplex.sh.md]: ./auxilary/demultiplex.sh.html
[workFlow.mak.md]: ./runWorkFlow/workFlow.mak.html
[sxCutR2AdapterFilterCumulate.sh.md]: ./sx/sxCutR2AdapterFilterCumulate.sh.html
[sxExtractMarker.sh.md]: ./sx/sxExtractMarker.sh.html
[getSxPlasmidFileRef.sh.md]: ./sx/getSxPlasmidFileRef.sh.html

# Source
~~~bash
# The following parameters should be replaced.
makeTarget=${makeTarget:-test/test_work_flow/rearr.alg}
fastqFiles=${fastqFiles:-test/test_work_flow/A2-g1n-3.R2.fq.gz,test/test_work_flow/A2-g1n-3.fq.gz}
markerIndices=${markerIndices:-test/test_work_flow/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.target.fa,test/test_work_flow/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.pair.fa}
minScores=${minScores:-30,100}

minToMapShear=${minToMapShear:-30}
refFile=${refFile:-test/test_work_flow/final_hgsgrna_libb_all_0811_NGG_scaffold_nor_G1.csv.ref}
directionFile=${directionFile:-"${refFile}.direct"}
ext1up=${ext1up:-50}
ext1down=${ext1down:-0}
ext2up=${ext2up:-10}
ext2down=${ext2down:-100}

# The following parameters are default in most cases.
genome=${genome:-"${GENOME}"}
bowtie2index=${bowtie2index:-"${BOWTIE2INDEX}"}
s0=${s0:--6}
s1=${s1:-4}
s2=${s2:-2}
u=${u:--3}
v=${v:--9}
ru=${ru:-0}
rv=${rv:-0}
qu=${qu:-0}
qv=${qv:--5}

if [ -f "workFlow.mak" ]
then
    make_file="workFlow.mak"
else
    make_file="$CONDA_PREFIX/share/rearr/workFlow.mak"
fi

make $@ -f "${make_file}" "${makeTarget}" \
    fastqFiles="${fastqFiles}" \
    markerIndices="${markerIndices}" \
    minScores="${minScores}" \
    genome="${genome}" \
    bowtie2index="${bowtie2index}" \
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
    qv="${qv}" \
    minToMapShear="${minToMapShear}"
~~~

~~~bash
alias ~~~=":" # This suppresses a warning and is not part of source.
~~~
