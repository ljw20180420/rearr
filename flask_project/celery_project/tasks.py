from .app import celeryApp
import os
import subprocess

@celeryApp.task(bind=True)
def celeryAlignReads(self, fileName=None, ref1=None, ref2=None, cut1=None, cut2=None, PAM1='NGG', PAM2='NGG', jobPath=None, exePath=None, downloadURL=None):
    savePath = os.path.join(jobPath, self.request.id)
    fastqFile = os.path.join(savePath, fileName) 
    cmd = f'''{os.path.join(exePath, "rearr_run.sh")} {fastqFile} {ref1} {ref2} {cut1} {cut2} {PAM1} {PAM2}'''
    with open(f"{fastqFile}.stdout", "w") as sdo, open(f"{fastqFile}.stderr", "w") as sde:
        _ = subprocess.run(cmd, stdout = sdo, stderr = sde, shell=True, executable="/bin/bash")
    _ = subprocess.run(f'''cd {jobPath}; zip -rm {savePath}.zip {self.request.id} -x {self.request.id}/{fileName}; rm -r {savePath}''', shell=True, executable="/bin/bash")
    return downloadURL