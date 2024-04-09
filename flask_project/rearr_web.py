#!/usr/bin/env python

import os, uuid, subprocess, sys
from flask import render_template, request, send_file, send_from_directory
from werkzeug.utils import secure_filename
from flaskApp import CustomFlask
from celery_project.celery import celery_init_app

flaskApp = CustomFlask(__name__, static_folder="../vue_project/dist/assets", template_folder="../vue_project/dist", static_url_path="/assets")

flaskApp.config.from_mapping(
    CELERY=dict(
        broker_url="amqp://",
        result_backend="rpc://",
    ),
)

celeryApp = celery_init_app(flaskApp)

@flaskApp.route('/favicon.ico')
def favicon():
    return send_from_directory(os.path.join(flaskApp.root_path, '../vue_project/dist'),
                               'favicon.ico', mimetype='favicon.ico')

@flaskApp.get("/")
def home_page():
    return render_template("index.html")

@flaskApp.put('/upload')
def uploadFile():
    file = request.files['files[]']
    uUiD = uuid.uuid4().hex
    savepath = os.path.join(flaskApp.root_path, "jobs", uUiD)
    os.makedirs(savepath, exist_ok=True)
    fastqfile = os.path.join(savepath, secure_filename(file.filename)) 
    file.save(fastqfile)
    cmd = f'''{os.path.join(flaskApp.root_path, "../rearr_run.sh")} {fastqfile} {request.form['ref1']} {request.form['ref2']} {request.form['cut1']} {request.form['cut2']} {request.form['PAM1']} {request.form['PAM2']}'''
    with open(f"{fastqfile}.stdout", "w") as sdo, open(f"{fastqfile}.stderr", "w") as sde:
        _ = subprocess.run(cmd, stdout = sdo, stderr = sde, shell=True, executable="/bin/bash")
    _ = subprocess.run(f'''cd {os.path.join(flaskApp.root_path, "jobs")}; zip -rm {savepath}.zip {uUiD} -x {uUiD}/{file.filename}; rm -r {savepath}''', shell=True, executable="/bin/bash")
    return send_file(f'''{savepath}.zip''', as_attachment=True)

if __name__ == "__main__":
    # flaskApp.run()
    flaskApp.run(debug=True)