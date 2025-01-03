from .app import celeryApp
import subprocess
import os
import time
import shutil

@celeryApp.on_after_finalize.connect
def setup_periodic_tasks(sender, **kwargs):
    sender.add_periodic_task(86400, celeryClearTmp.s("tmp"), name='clear temporary file')

@celeryApp.task
def celeryClearTmp(tmpPath):
    now = time.time()
    for uuid in os.listdir(tmpPath):
        if uuid == '.gitignore':
            continue
        sessionDir = os.path.join(tmpPath, uuid)
        # save for a week
        if os.stat(sessionDir).st_mtime < now - 7 * 86400:
            shutil.rmtree(sessionDir)
    return "SUCCESS"

@celeryApp.task
def celeryRemoveDuplicates(inputFiles, rmDupFile):
    subprocess.run(f'''removeDuplicates.md {" ".join(inputFiles)} >{rmDupFile} ''', shell=True, executable="/bin/bash")
    return os.path.basename(rmDupFile)

@celeryApp.task
def celeryBuildSpliter(spliters):
    for spliter in spliters:
        subprocess.run(f'''bowtie2-build {spliter} {spliter}''', shell=True, executable="/bin/bash")
    return [f'{os.path.basename(spliter)}.{ext}' for spliter in spliters for ext in ['1.bt2', '2.bt2', '3.bt2', '4.bt2', 'rev.1.bt2', 'rev.2.bt2']]

@celeryApp.task
def celeryDemultiplex(rmDupFile, spliterIndices, minScores, demultiplexFile):
    subprocess.run(f'''spliterIndices={",".join(spliterIndices)} minScores={",".join([str(minScore) for minScore in minScores])} demultiplex.md {rmDupFile} >{demultiplexFile}''', shell=True, executable="/bin/bash")
    return os.path.basename(demultiplexFile)

@celeryApp.task
def celerySxPostProcess(demultiplexFile, minToMapShear, toMapFile):
    subprocess.run(f'''sxCutR2AdapterFilterCumulate.md {demultiplexFile} {minToMapShear} >{toMapFile}''', shell=True, executable="/bin/bash")
    return os.path.basename(toMapFile)

@celeryApp.task
def celeryRearrange(toMapFile, refFile, correctFile, s0, s1, s2, u, v, ru, rv, qu, qv, alignFile):
    subprocess.run(f'''rearrangement <{toMapFile} 3<{refFile} -s0 {s0} -s1 {s1} -s2 {s2} -u {u} -v {v} -ru {ru} -rv {rv} -qu {qu} -qv {qv} | gawk -f correct_micro_homology.awk -- {refFile} {correctFile} >{alignFile}''', shell=True, executable="/bin/bash")
    return os.path.basename(alignFile)

@celeryApp.task
def celeryDefaultCorrect(refFile, correctFile):
    with open(refFile, 'r') as rfd, open(correctFile, 'w') as cfd:
        for line in rfd:
            cfd.write('\t'.join(['up'] * (len(line.split('\t')) // 3 - 1)) + '\n')
    return os.path.basename(correctFile)

@celeryApp.task
def celeryGetReference(csvFile, genome, bowtie2index, ext1up, ext1down, ext2up, ext2down, refFile):
    subprocess.run(f'''getSxCsvFileRef.md {csvFile} {genome} {bowtie2index} {ext1up} {ext1down} {ext2up} {ext2down} >{refFile}''', shell=True, executable="/bin/bash")
    return [refFile]

@celeryApp.task
def celeryGetSpliters(csvFile, targetSpliter, pairSpliter):
    subprocess.run(f'''sxExtractSpliter.md {csvFile} >{targetSpliter} 3>{pairSpliter}''', shell=True, executable="/bin/bash")
    return [targetSpliter, pairSpliter]