from .app import celeryApp
import subprocess
import os
import time

@celeryApp.on_after_finalize.connect
def setup_periodic_tasks(sender, **kwargs):
    sender.add_periodic_task(86400, celeryClearTmp.s("tmp"), name='clear temporary file')

@celeryApp.task
def celeryClearTmp(tmpPath):
    now = time.time()
    for file in os.listdir(tmpPath):
        file = os.path.join(tmpPath, file)
        if os.stat(file).st_mtime < now - 7 * 86400:
            os.remove(file)
    return "success"

@celeryApp.task
def celeryUpload(uploadFiles):
    return uploadFiles

@celeryApp.task
def celeryRemoveDuplicates(inputFiles, rmDupFile):
    subprocess.run(f'''removeDuplicates.sh {" ".join(inputFiles)} >{rmDupFile} ''', shell=True, executable="/bin/bash")
    return [rmDupFile]

@celeryApp.task
def celeryBuildSpliter(spliter, spliterIndex):
    subprocess.run(f'''bowtie2-build {spliter} {spliterIndex}''', shell=True, executable="/bin/bash")
    return [f'{spliterIndex}.{ext}' for ext in ['1.bt2', '2.bt2', '3.bt2', '4.bt2', 'rev.1.bt2', 'rev.2.bt2']]

@celeryApp.task
def celeryDemultiplex(rmDupFile, targetSpliterIndex, pairSpliterIndex, minScoreTarget, minScorePair, demultiplexFile):
    subprocess.run(f'''demultiplex.sh {rmDupFile} {targetSpliterIndex} {pairSpliterIndex} {minScoreTarget} {minScorePair} >{demultiplexFile}''', shell=True, executable="/bin/bash")
    return [demultiplexFile]

@celeryApp.task
def celerySxPostProcess(demultiplexFile, minToMapShear, toAlignFile):
    subprocess.run(f'''sxCutR2AdapterFilterCumulate.sh {demultiplexFile} {minToMapShear} >{toAlignFile}''', shell=True, executable="/bin/bash")
    return [toAlignFile]

@celeryApp.task
def celeryRearrange(toAlignFile, refFile, s0, s1, s2, u, v, ru, rv, qu, qv, PAM1, PAM2, alignFile):
    subprocess.run(f'''rearrangement <{toAlignFile} 3<{refFile} -s0 {s0} -s1 {s1} -s2 {s2} -u {u} -v {v} -ru {ru} -rv {rv} -qu {qu} -qv {qv} | gawk -f correct_micro_homology.awk -- {refFile} {PAM1} {PAM2} >{alignFile}''', shell=True, executable="/bin/bash")
    return [alignFile]

@celeryApp.task
def celeryGetReference(csvFile, genome, bowtie2index, ext1up, ext1down, ext2up, ext2down, refFile):
    subprocess.run(f'''getSxCsvFileRef.sh {csvFile} {genome} {bowtie2index} {ext1up} {ext1down} {ext2up} {ext2down} >{refFile}''', shell=True, executable="/bin/bash")
    return [refFile]

@celeryApp.task
def celeryGetSpliters(csvFile, targetSpliter, pairSpliter):
    subprocess.run(f'''sxExtractSpliter.sh {csvFile} >{targetSpliter} 3>{pairSpliter}''', shell=True, executable="/bin/bash")
    return [targetSpliter, pairSpliter]