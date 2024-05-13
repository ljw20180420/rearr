#!/usr/bin/env python

import os, tempfile
from flask import Flask, render_template, request, send_file, send_from_directory, redirect
from celery_project.tasks import celeryRemoveDuplicates, celeryBuildSpliter, celeryDemultiplex, celeryRearrange
from celery.result import AsyncResult
from werkzeug.middleware.proxy_fix import ProxyFix

flaskApp = Flask(__name__, static_folder="vue_project/dist/assets", template_folder="vue_project/dist", static_url_path="/assets")
flaskApp.wsgi_app = ProxyFix(flaskApp.wsgi_app, x_for=1, x_proto=1, x_host=1, x_prefix=1)

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

@flaskApp.put('/uploadFile')
def uploadFile():
    tmpFile = tempfile.mkstemp(dir=tmpFold)[1]
    request.files['file[]'].save(tmpFile)
    return os.path.basename(tmpFile)

@flaskApp.put('/runJob/removeDuplicates')
def removeDuplicates():
    inputFiles = [os.path.join("tmp", file) for file in [request.form["target file"], request.form["pair file"]]]
    rmDupFile = os.path.relpath(tempfile.mkstemp(dir=tmpFold)[1], start=flaskApp.root_path)
    result = celeryRemoveDuplicates.delay(inputFiles, rmDupFile)
    return {'taskId': result.id, 'fileName': os.path.basename(rmDupFile)}

@flaskApp.put('/runJob/buildSpliter/<string:spliterFile>')
def buildSpliter(spliterFile):
    spliter = os.path.join("tmp", request.form[spliterFile])
    spliterIndex = os.path.relpath(tempfile.mkstemp(dir=tmpFold)[1], start=flaskApp.root_path)
    bowtie2buildLog = os.path.relpath(tempfile.mkstemp(dir=tmpFold)[1], start=flaskApp.root_path)
    result = celeryBuildSpliter.delay(spliter, spliterIndex, bowtie2buildLog)
    return {'taskId': result.id, 'fileName': os.path.basename(spliterIndex)}

@flaskApp.put('/runJob/demultiplex')
def demultiplex():
    rmDupFile = os.path.join("tmp", request.form["file without duplicates"])
    targetSpliterIndex = os.path.join("tmp", request.form["target spliter index"])
    pairSpliterIndex = os.path.join("tmp", request.form["pair spliter index"])
    minScoreTarget = request.form["minimal alignment score of target spliter"]
    minScorePair = request.form["minimal alignment score of pair spliter"]
    demultiplexFile = os.path.relpath(tempfile.mkstemp(dir=tmpFold)[1], start=flaskApp.root_path)
    result = celeryDemultiplex.delay(rmDupFile, targetSpliterIndex, pairSpliterIndex, minScoreTarget, minScorePair, demultiplexFile)
    return {'taskId': result.id, 'fileName': os.path.basename(demultiplexFile)}

@flaskApp.put('/runJob/rearrange')
def rearrange():
    toAlignFile = os.path.join("tmp", request.form["file of reads to align"])
    refFile = os.path.join("tmp", request.form["file of reference"])
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
    return {'taskId': result.id, 'fileName': os.path.basename(alignFile)}

@flaskApp.get("/download/<string:taskId>")
def downloadResult(taskId):
    result = AsyncResult(taskId)
    if result.status == 'SUCCESS':
        return send_file(result.get(), as_attachment=True)
    return 'Accepted', 202

@flaskApp.get("/flower")
def flower():
    return redirect('http://www.rearr.xyz:5555', code=301)

@flaskApp.get("/shiny")
def shiny():
    return redirect('http://www.rearr.xyz:3838', code=301)

if __name__ == "__main__":
    flaskApp.run()