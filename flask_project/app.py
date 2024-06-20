#!/usr/bin/env python

import os, tempfile, io, zipfile, re
from flask import Flask, render_template, request, send_file, send_from_directory, redirect
from celery_project.tasks import celeryUpload, celeryRemoveDuplicates, celeryBuildSpliter, celeryDemultiplex, celerySxPostProcess, celeryRearrange, celeryGetReference, celeryGetSpliters
from celery.result import AsyncResult
from werkzeug.middleware.proxy_fix import ProxyFix

flaskApp = Flask(__name__, static_folder="vue_project/dist/assets", template_folder="vue_project/dist", static_url_path="/assets")
flaskApp.wsgi_app = ProxyFix(flaskApp.wsgi_app, x_for=1, x_proto=1, x_host=1, x_prefix=1)
flaskApp.config['MAX_CONTENT_LENGTH'] = 10 * 1024 * 1024 * 1024

@flaskApp.route('/favicon.ico')
def favicon():
    return send_from_directory(
        os.path.join(flaskApp.root_path, 'vue_project/dist'),
        'favicon.ico', mimetype='favicon.ico'
    )

@flaskApp.get("/")
def homePage():
    return render_template("index.html")

tmpFold = os.path.join(flaskApp.root_path, "tmp")

@flaskApp.put('/upload')
def upload():
    files = request.files.getlist('files[]')
    mat = re.search('\.(rev\.)?[1-4]\.bt2', files[0].filename)
    if not mat:
        uploadFiles = [os.path.relpath(tempfile.mkstemp(dir=tmpFold)[1], start=flaskApp.root_path) for file in files]
    else:
        pathBase = os.path.relpath(tempfile.mkstemp(dir=tmpFold)[1], start=flaskApp.root_path)
        os.remove(os.path.join(flaskApp.root_path, pathBase))
        uploadFiles = [pathBase + re.search('\.(rev\.)?[1-4]\.bt2', file.filename).group() for file in files]
    for file, uploadFile in zip(files, uploadFiles):
        file.save(os.path.join(flaskApp.root_path, uploadFile))
    result = celeryUpload.delay(uploadFiles)
    return {'taskId': result.id, 'value': [os.path.basename(uploadFile) for uploadFile in uploadFiles]}

@flaskApp.put('/runJob/removeDuplicates')
def removeDuplicates():
    inputFiles = [os.path.join("tmp", file) for file in [request.form["target file[]"], request.form["pair file[]"]]]
    rmDupFile = os.path.relpath(tempfile.mkstemp(dir=tmpFold)[1], start=flaskApp.root_path)
    result = celeryRemoveDuplicates.delay(inputFiles, rmDupFile)
    return [{'taskId': result.id, 'name': 'file without duplicates', 'value': [os.path.basename(rmDupFile)]}]

@flaskApp.put('/runJob/buildSpliter/<string:spliterFile>')
def buildSpliter(spliterFile):
    spliter = os.path.join("tmp", request.form[f"{spliterFile}[]"])
    spliterIndex = os.path.relpath(tempfile.mkstemp(dir=tmpFold)[1], start=flaskApp.root_path)
    os.remove(os.path.join(flaskApp.root_path, spliterIndex))
    result = celeryBuildSpliter.delay(spliter, spliterIndex)
    return [{'taskId': result.id, 'name': f'{spliterFile} index', 'value': [f'{os.path.basename(spliterIndex)}.{ext}' for ext in ['1.bt2', '2.bt2', '3.bt2', '4.bt2', 'rev.1.bt2', 'rev.2.bt2']]}]

@flaskApp.put('/runJob/demultiplex')
def demultiplex():
    rmDupFile = os.path.join("tmp", request.form["file without duplicates[]"])
    try:
        targetSpliterIndex = os.path.join("tmp", request.form["target spliter index[]"])
        mat = re.search('\.(rev\.)?[1-4]\.bt2', targetSpliterIndex)
        targetSpliterIndex = targetSpliterIndex[:mat.span()[0]]
        minScoreTarget = request.form["minimal alignment score of target spliter"]
    except:
        targetSpliterIndex = ""
        minScoreTarget = ""
    try:
        pairSpliterIndex = os.path.join("tmp", request.form["pair spliter index[]"])
        mat = re.search('\.(rev\.)?[1-4]\.bt2', pairSpliterIndex)
        pairSpliterIndex = pairSpliterIndex[:mat.span()[0]]
        minScorePair = request.form["minimal alignment score of pair spliter"]
    except:
        pairSpliterIndex = ""
        minScorePair = ""
    if not targetSpliterIndex and not pairSpliterIndex:
        raise Exception("At least one of targetSpliter and pairSpliter should be specified")
    demultiplexFile = os.path.relpath(tempfile.mkstemp(dir=tmpFold)[1], start=flaskApp.root_path)
    result = celeryDemultiplex.delay(rmDupFile, targetSpliterIndex, pairSpliterIndex, minScoreTarget, minScorePair, demultiplexFile)
    return [{'taskId': result.id, 'name': 'demultiplex file', 'value': [os.path.basename(demultiplexFile)]}]

@flaskApp.put('/runJob/sxPostProcess')
def sxPostProcess():
    demultiplexFile = os.path.join("tmp", request.form["demultiplex file[]"])
    minToMapShear = request.form["minimal base number after remove 5' spliter and 3' adapter from target"]
    toAlignFile = os.path.relpath(tempfile.mkstemp(dir=tmpFold)[1], start=flaskApp.root_path)
    result = celerySxPostProcess.delay(demultiplexFile, minToMapShear, toAlignFile)
    return [{'taskId': result.id, 'name': 'file of reads to align', 'value': [os.path.basename(toAlignFile)]}]

@flaskApp.put('/runJob/rearrange')
def rearrange():
    toAlignFile = os.path.join("tmp", request.form["file of reads to align[]"])
    refFile = os.path.join("tmp", request.form["file of reference[]"])
    s0 = request.form["mismatching score"]
    s1 = request.form["matching score for non-extension reference part"]
    s2 = request.form["matching score for extension reference part"]
    u = request.form["gap-extending penalty"]
    v = request.form["gap-opening penalty"]
    ru = request.form["gap-extending penalty for unaligned reference end"]
    rv = request.form["gap-opening penalty for unaligned reference end"]
    qu = request.form["gap-extending penalty for unaligned query part"]
    qv = request.form["gap-opening penalty for unaligned query part"]
    PAM1 = request.form["PAM1"]
    PAM2 = request.form["PAM2"]
    alignFile = os.path.relpath(tempfile.mkstemp(dir=tmpFold)[1], start=flaskApp.root_path)
    result = celeryRearrange.delay(toAlignFile, refFile, s0, s1, s2, u, v, ru, rv, qu, qv, PAM1, PAM2, alignFile)
    return [{'taskId': result.id, 'name': 'alignments', 'value': [os.path.basename(alignFile)]}]

@flaskApp.put('/runJob/indexGenome')
def indexGenome():
    genome = os.path.join("tmp", request.form["genome[]"])
    bowtie2index = os.path.relpath(tempfile.mkstemp(dir=tmpFold)[1], start=flaskApp.root_path)
    os.remove(os.path.join(flaskApp.root_path, bowtie2index))
    result = celeryBuildSpliter.delay(genome, bowtie2index)
    return [{'taskId': result.id, 'name': 'genome index', 'value': [f'{os.path.basename(bowtie2index)}.{ext}' for ext in ['1.bt2', '2.bt2', '3.bt2', '4.bt2', 'rev.1.bt2', 'rev.2.bt2']]}]

@flaskApp.put('/runJob/getReference')
def getReference():
    csvFile = os.path.join("tmp", request.form["csvfile[]"])
    genome = os.path.join("tmp", request.form["genome[]"])
    bowtie2index = os.path.join("tmp", request.form["genome index[]"])
    mat = re.search('\.(rev\.)?[1-4]\.bt2', bowtie2index)
    bowtie2index = bowtie2index[:mat.span()[0]]
    ext1up = request.form["cleavage 1 extend upstream"]
    ext1down = request.form["cleavage 1 extend downstream"]
    ext2up = request.form["cleavage 2 extend upstream"]
    ext2down = request.form["cleavage 2 extend downstream"]
    refFile = os.path.relpath(tempfile.mkstemp(dir=tmpFold)[1], start=flaskApp.root_path)
    result = celeryGetReference.delay(csvFile, genome, bowtie2index, ext1up, ext1down, ext2up, ext2down, refFile)
    return [{'taskId': result.id, 'name': 'file of reference', 'value': [os.path.basename(refFile)]}]

@flaskApp.put('/runJob/getSpliters')
def getSpliters():
    csvFile = os.path.join("tmp", request.form["csvfile[]"])
    targetSpliter = os.path.relpath(tempfile.mkstemp(dir=tmpFold)[1], start=flaskApp.root_path)
    os.remove(os.path.join(flaskApp.root_path, targetSpliter))
    pairSpliter = targetSpliter + ".pair"
    targetSpliter = targetSpliter + ".target"
    result = celeryGetSpliters.delay(csvFile, targetSpliter, pairSpliter)
    return [
        {'taskId': result.id, 'name': 'target spliter', 'value': [os.path.basename(targetSpliter)]},
        {'taskId': result.id, 'name': 'pair spliter', 'value': [os.path.basename(pairSpliter)]}
    ]

@flaskApp.get("/download/<string:taskId>")
def download(taskId):
    result = AsyncResult(taskId)
    if result.status == 'SUCCESS':
        stream = io.BytesIO()
        with zipfile.ZipFile(stream, 'w') as zf:
            for file in result.get():
                zf.write(os.path.join(flaskApp.root_path, file), os.path.basename(file))
        stream.seek(0)
        download_name = os.path.basename(result.get()[0]).split(".")[0]
        return send_file(
            stream,
            as_attachment=True,
            download_name=f'{download_name}.zip'
        )
    return 'Accepted', 202

@flaskApp.get("/inspect/<string:taskId>")
def inspect(taskId):
    result = AsyncResult(taskId)
    return result.status

if __name__ == "__main__":
    flaskApp.run()
