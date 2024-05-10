#!/usr/bin/env python

import os, tempfile
from flask import Flask, render_template, request, send_file, send_from_directory, redirect
from celery_project.tasks import celeryRemoveDuplicates, celeryBowtie2build, celeryDemultiplex, celeryRearrange
from celery.result import AsyncResult

# change Jinja2 delimiter to solve the comflicts with vue3
class VueFlask(Flask):
    jinja_options = Flask.jinja_options.copy()
    jinja_options.update(dict(
        block_start_string='{b{s',
        block_end_string='b}e}',
        variable_start_string='{v{s',
        variable_end_string='v}e}',
        comment_start_string='{c{s',
        comment_end_string='c}e}',
    ))

flaskApp = VueFlask(__name__, static_folder="vue_project/dist/assets", template_folder="vue_project/dist", static_url_path="/assets")

@flaskApp.route('/favicon.ico')
def favicon():
    return send_from_directory(
        os.path.join(flaskApp.root_path, 'vue_project/dist'),
        'favicon.ico', mimetype='favicon.ico'
    )

@flaskApp.get("/")
def homePage():
    return render_template("index.html")

fileType2Name = {}
name2Value = {}
tmpFold = os.path.join(flaskApp.root_path, "tmp")
os.makedirs(tmpFold, exist_ok=True)

@flaskApp.put('/uploadFile/<string:fileType>')
def uploadFile(fileType):
    global fileType2Name
    tmpFile = tempfile.mkstemp(dir=tmpFold)[1]
    fileType2Name[fileType] = os.path.relpath(tmpFile, start=flaskApp.root_path)
    request.files['file[]'].save(tmpFile)
    return 'OK', 200

@flaskApp.put('/setValues')
def setValues():
    global name2Value
    for key, val in request.form.items():
        name2Value[key] = val
    return 'OK', 200

@flaskApp.get('/runJob/removeDup')
def removeDup():
    global fileType2Name
    fileType2Name["rmDupFile"] = os.path.relpath(tempfile.mkstemp(dir=tmpFold)[1], start=flaskApp.root_path)
    result = celeryRemoveDuplicates.delay([fileType2Name["targetFile"], fileType2Name["pairFile"]], fileType2Name["rmDupFile"])
    return result.id

@flaskApp.get('/runJob/bowtie2build')
def bowtie2build():
    global fileType2Name
    fileType2Name["targetSpliterIndex"] = os.path.relpath(tempfile.mkstemp(dir=tmpFold)[1], start=flaskApp.root_path)
    fileType2Name["pairSpliterIndex"] = os.path.relpath(tempfile.mkstemp(dir=tmpFold)[1], start=flaskApp.root_path)
    fileType2Name["bowtie2buildLog"] = os.path.relpath(tempfile.mkstemp(dir=tmpFold)[1], start=flaskApp.root_path)
    result = celeryBowtie2build.delay(fileType2Name["targetSpliter"], fileType2Name["pairSpliter"], fileType2Name["targetSpliterIndex"], fileType2Name["pairSpliterIndex"], fileType2Name["bowtie2buildLog"])
    return result.id

@flaskApp.get('/runJob/demultiplex')
def demultiplex():
    global fileType2Name
    fileType2Name["demultiplexFile"] = os.path.relpath(tempfile.mkstemp(dir=tmpFold)[1], start=flaskApp.root_path)
    result = celeryDemultiplex.delay(fileType2Name["rmDupFile"], fileType2Name["targetSpliterIndex"], fileType2Name["pairSpliterIndex"], name2Value["minScoreTarget"], name2Value["minScorePair"], fileType2Name["demultiplexFile"])
    return result.id

@flaskApp.get('/runJob/rearrange')
def rearrange():
    global fileType2Name
    fileType2Name["alignFile"] = os.path.relpath(tempfile.mkstemp(dir=tmpFold)[1], start=flaskApp.root_path)
    result = celeryRearrange.delay(fileType2Name["toAlignFile"], fileType2Name["refFile"], name2Value["s0"], name2Value["s1"], name2Value["s2"], name2Value["u"], name2Value["v"], name2Value["ru"], name2Value["rv"], name2Value["qu"], name2Value["qv"], name2Value["PAM1"], name2Value["PAM2"], fileType2Name["alignFile"])
    return result.id

@flaskApp.get("/download/<string:taskId>")
def downloadResult(taskId):
    result = AsyncResult(taskId)
    if result.status == 'SUCCESS':
        return send_file(result.get(), as_attachment=True)
    return 'Accepted', 202

@flaskApp.get("/flower")
def flower():
    return redirect('http://localhost:5555', code=301)

@flaskApp.get("/shiny")
def shiny():
    return redirect('http://localhost:3838', code=301)

if __name__ == "__main__":
    # flaskApp.run()
    flaskApp.run(debug=True)