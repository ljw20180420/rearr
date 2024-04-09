#!/usr/bin/env python

import os, sys
from flask import render_template, request, send_file, send_from_directory, redirect, url_for, jsonify
from werkzeug.utils import secure_filename
from app import VueFlask
from celery import uuid
from celery_project.app import celery_init_app
from celery_project.tasks import celeryAlignReads
from celery.result import AsyncResult

flaskApp = VueFlask(__name__, static_folder="vue_project/dist/assets", template_folder="vue_project/dist", static_url_path="/assets")

flaskApp.config.from_mapping(
    CELERY=dict(
        broker_url="amqp://localhost",
        result_backend="rpc://localhost",
        include=['celery_project.tasks'],
        broker_connection_retry_on_startup=True,
    ),
)

celeryApp = celery_init_app(flaskApp)

@flaskApp.route('/favicon.ico')
def favicon():
    return send_from_directory(os.path.join(flaskApp.root_path, 'vue_project/dist'),
                               'favicon.ico', mimetype='favicon.ico')

@flaskApp.get("/")
def homePage():
    return render_template("index.html")

@flaskApp.put('/align')
def alignReads():
    # flask cannot serialize FileStorage, so FileStorage cannot be passed to celery. The compromise is to apply a synchronic saving.
    task_id = uuid()
    file = request.files['files[]']
    fileName = secure_filename(file.filename)
    savePath = os.path.join(flaskApp.root_path, "jobs", task_id)
    os.makedirs(savePath, exist_ok=True)
    file.save(os.path.join(savePath, fileName))

    result = celeryAlignReads.apply_async(kwargs={
            'fileName': fileName,
            'ref1': request.form['ref1'],
            'ref2': request.form['ref2'],
            'cut1': request.form['cut1'],
            'cut2': request.form['cut2'],
            'PAM1': request.form['PAM1'],
            'PAM2': request.form['PAM2'],
            'jobPath': os.path.join(flaskApp.root_path, "jobs"),
            'exePath': os.path.join(flaskApp.root_path, ".."),
            'downloadURL': url_for('downloadResult', task_id=task_id, _external=True, schema=None)
        }, task_id=task_id)
    return f'''{task_id} {result.status}'''

@flaskApp.get("/flower")
def flower():
    return redirect('localhost:5555', code=301)

@flaskApp.get("/inspect/<string:task_id>")
def inspect(task_id):
    return AsyncResult(task_id).status

@flaskApp.get("/download/<string:task_id>")
def downloadResult(task_id):
    zipFile = os.path.join(flaskApp.root_path, "jobs", f'''{task_id}.zip''')
    return send_file(zipFile, as_attachment=True)

if __name__ == "__main__":
    # flaskApp.run()
    flaskApp.run(debug=True)