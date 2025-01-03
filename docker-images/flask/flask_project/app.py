#!/usr/bin/env python

import os, tempfile, io, zipfile, re, uuid, shutil, pathlib
from flask import Flask, render_template, send_file, send_from_directory, request, session
from celery_project.tasks import celeryRemoveDuplicates, celeryBuildSpliter, celeryDemultiplex, celerySxPostProcess, celeryRearrange, celeryDefaultCorrect, celeryGetReference, celeryGetSpliters
from celery.result import AsyncResult
# from werkzeug.middleware.proxy_fix import ProxyFix

flaskApp = Flask(__name__, static_folder="vue_project/dist/assets", template_folder="vue_project/dist", static_url_path="/assets")
# flaskApp.wsgi_app = ProxyFix(flaskApp.wsgi_app, x_for=1, x_proto=1, x_host=1, x_prefix=1)
flaskApp.config['MAX_CONTENT_LENGTH'] = 10 * 1024 * 1024 * 1024
flaskApp.secret_key = b'913d1c26d46f82f119662371bc71efa4de6902a8ee2a378a6f342b8eb39b2d52'

@flaskApp.route('/favicon.ico')
def favicon():
    return send_from_directory(
        os.path.join(flaskApp.root_path, 'vue_project/dist'),
        'favicon.ico',
        mimetype='favicon.ico'
    ), 200

@flaskApp.get("/")
def homePage():
    if 'dir' not in session:
        session['dir'] = os.path.join(flaskApp.root_path, 'tmp', str(uuid.uuid4()))
    os.makedirs(session['dir'], exist_ok=True)
    return render_template("index.html")

@flaskApp.route('/stop')
def stop():
    if os.path.exists(session['dir']):
        shutil.rmtree(session['dir'])
    return "Stop", 200

@flaskApp.put('/upload')
def upload():
    file = request.files['file[]']
    file.save(os.path.join(session['dir'], file.filename))
    return file.filename, 200

# This api is accessed by href. Firefox will add trailing slash for cached pages. To avoid using strict_slashes=False, add a slash to this api.
@flaskApp.get("/download/<string:filename>/")
def download(filename):
    try:
        return send_file(
            os.path.join(session['dir'], filename),
            as_attachment=True
        ), 200
    except FileNotFoundError:
        return "Accepted", 202

# This api is accessed by href. Firefox will add trailing slash for cached pages. To avoid using strict_slashes=False, add a slash to this api.
@flaskApp.get("/inspect/<string:taskId>/")
def inspect(taskId):
    result = AsyncResult(taskId)
    return result.status, 200

@flaskApp.put('/runJob/removeDuplicates')
def removeDuplicates():
    fastqFiles = [
        os.path.join(session['dir'], file['value'])
        for file in request.get_json()['.fastq files']
    ]
    rmDupFile = os.path.join(session['dir'], 'rearr.noDup')
    result = celeryRemoveDuplicates.delay(
        fastqFiles,
        rmDupFile
    )
    return {
      '.noDup file': {
        'taskId': result.id,
        'value': 'rearr.noDup'
      }
    }

@flaskApp.put('/runJob/buildSpliter')
def buildSpliter():
    spliters = [
        os.path.join(session['dir'], file['value'])
        for file in request.get_json()['.fasta files']
    ]
    result = celeryBuildSpliter.delay(spliters)
    return {
        'auxiliaries': [
            {
                'spliterIndex': {
                    ext: {
                        'taskId': result.id,
                        'value': f'{os.path.basename(spliter)}.{ext}'
                    }
                    for ext in ['1.bt2', '2.bt2', '3.bt2', '4.bt2', 'rev.1.bt2', 'rev.2.bt2']
                }
            }
            for spliter in spliters
        ]
    }

@flaskApp.put('/runJob/demultiplex')
def demultiplex():
    json = request.get_json()
    rmDupFile = os.path.join(session['dir'], json['.noDup file']['value'])
    spliterIndices = [
        os.path.join(session['dir'], auxiliary['spliterIndex']['1.bt2']['value'][:-6])
        for auxiliary in json['auxiliaries']
    ]
    minScores = [
        str(auxiliary['minScore']['value'])
        for auxiliary in json['auxiliaries']
    ]
    demultiplexFile = os.path.splitext(rmDupFile)[0] + '.demultiplex'
    result = celeryDemultiplex.delay(
        rmDupFile,
        spliterIndices,
        minScores,
        demultiplexFile
    )
    return {
        '.demultiplex file': {
            'taskId': result.id,
            'value': 'rearr.demultiplex'
        }
    }

@flaskApp.put('/runJob/sxPostProcess')
def sxPostProcess():
    json = request.get_json()
    demultiplexFile = os.path.join(session['dir'], json['.demultiplex file']['value'])
    toMapFile = os.path.splitext(demultiplexFile)[0] + '.post'
    result = celerySxPostProcess.delay(
        demultiplexFile,
        json['minimal base number']['value'],
        toMapFile
    )
    return {
        '.post file': {
            'taskId': result.id,
            'value': os.path.basename(toMapFile)
        }
    }

@flaskApp.put('/runJob/rearrange')
def rearrange():
    json = request.get_json()
    toMapFile = os.path.join(session['dir'], json['.post file']['value'])
    refFile = os.path.join(session['dir'], json['.ref file']['value'])
    correctFile = os.path.join(session['dir'], json['.correct file']['value'])
    alignFile = os.path.splitext(toMapFile)[0] + '.alg'
    result = celeryRearrange.delay(
        toMapFile,
        refFile,
        correctFile,
        json['align scores']['s0']['value'],
        json['align scores']['s1']['value'],
        json['align scores']['s2']['value'],
        json['align scores']['u']['value'],
        json['align scores']['v']['value'],
        json['align scores']['ru']['value'],
        json['align scores']['rv']['value'],
        json['align scores']['qu']['value'],
        json['align scores']['qv']['value'],
        alignFile
    )
    return {
        '.alg file': {
            'taskId': result.id,
            'value': os.path.basename(alignFile)
        }
    }

@flaskApp.put('/runJob/defaultCorrect')
def defaultCorrect():
    refFile = os.path.join(session['dir'], request.get_json()['.ref file']['value'])
    correctFile = os.path.splitext(refFile)[0] + '.correct'
    result = celeryDefaultCorrect.delay(refFile, correctFile)
    return {
        '.correct file': {
            'taskId': result.id,
            'value': os.path.basename(correctFile)
        }
    }
    

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

if __name__ == "__main__":
    flaskApp.run()
