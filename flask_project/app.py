#!/usr/bin/env python

import os
from flask import Flask, render_template, request, send_file, send_from_directory, redirect, url_for
from werkzeug.utils import secure_filename
# from celery_project.app import celeryApp
from celery_project.tasks import celeryAlignReads
from celery.result import AsyncResult
from celery import uuid

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
    os.makedirs(os.path.join(flaskApp.root_path, "jobs", task_id), exist_ok=True)
    file.save(os.path.join(flaskApp.root_path, "jobs", task_id, fileName))

    result = celeryAlignReads.apply_async(kwargs={
            'fileName': fileName,
            'ref1': request.form['ref1'],
            'ref2': request.form['ref2'],
            'cut1': request.form['cut1'],
            'cut2': request.form['cut2'],
            'PAM1': request.form['PAM1'],
            'PAM2': request.form['PAM2'],
            'jobPath': "jobs",
            'exePath': "..",
            'downloadURL': url_for('downloadResult', task_id=task_id, _external=True),
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