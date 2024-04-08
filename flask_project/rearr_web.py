#!/usr/bin/env python

import os, uuid, subprocess, sys
from flask import Flask, render_template, request, send_file, jsonify
from werkzeug.utils import secure_filename

# change Jinja2 delimiter to solve the comflicts with vue3
class CustomFlask(Flask):
    jinja_options = Flask.jinja_options.copy()
    jinja_options.update(dict(
        block_start_string='{b{s',
        block_end_string='b}e}',
        variable_start_string='{v{s',
        variable_end_string='v}e}',
        comment_start_string='{c{s',
        comment_end_string='c}e}',
    ))

app = CustomFlask(__name__, static_folder="../vue_project/dist/assets", template_folder="../vue-project/dist", static_url_path="/assets")

@app.route("/")
def home_page():
    return render_template("index.html")

@app.route('/upload', methods = ["PUT"])
def uploadFile():
    file = request.files['files[]']
    savepath = os.path.join(app.root_path, "jobs", uuid.uuid4().hex)
    os.makedirs(savepath, exist_ok=True)
    fastqfile = os.path.join(savepath, secure_filename(file.filename)) 
    file.save(fastqfile)
    cmd = f'''{os.path.join(app.root_path, "../rearr_run.sh")} {fastqfile} {request.form['ref1']} {request.form['ref2']} {request.form['cut1']} {request.form['cut2']} {request.form['PAM1']} {request.form['PAM2']}'''
    with open(f"{fastqfile}.stdout", "w") as sdo, open(f"{fastqfile}.stderr", "w") as sde:
        _ = subprocess.run(cmd, stdout = sdo, stderr = sde, shell=True, executable="/bin/bash")
    _ = subprocess.run(f'''zip -rj {savepath}.zip {savepath}''', shell=True, executable="/bin/bash")
    return send_file(f'''{savepath}.zip''', as_attachment=True)

if __name__ == "__main__":
    print(app.root_path)
    # app.run()
    app.run(debug=True)