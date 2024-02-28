import os, shutil, subprocess
from flask import Flask, render_template, request, send_file, jsonify
from werkzeug.utils import secure_filename
# from flask_autoindex import AutoIndex

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

app = CustomFlask(__name__)
# AutoIndex(app)   

@app.route("/")
def home_page():
    return render_template("rearrangement.html")

@app.route('/upload', methods = ["PUT"])
def uploadFile():
    try:
        savepath = os.path.join(app.root_path, "users", request.form['username'], request.form['jobname'])
        os.makedirs(savepath, exist_ok=True)
        for file in request.files.getlist("files[]"):
            file.save(os.path.join(savepath, secure_filename(file.filename)))
        return "success"
    except Exception as err:
        return str(err)
    
def displayDir(links, rootpath, subpath):
    for fsterm in os.listdir(os.path.join(rootpath, "users", subpath)):
        if os.path.isdir(os.path.join(rootpath, "users", subpath, fsterm)):
            links[fsterm + "/"] = os.path.join("/users", subpath, fsterm)
        else:
            links[fsterm] = os.path.join("/users", subpath, fsterm)
    return jsonify(links)

@app.route('/users')
def initBrowser():
    return displayDir({}, app.root_path, "")

@app.route('/users/<path:subpath>', methods=["GET", "DELETE"])
def webBrowser(subpath):
    subpath = subpath.rstrip("/")
    fullpath = os.path.join(app.root_path, "users", subpath)
    if request.method == "GET":
        if os.path.isfile(fullpath):
            return send_file(fullpath, as_attachment=True)
        else:
            return displayDir({"../": os.path.dirname(os.path.join("/users", subpath))}, app.root_path, subpath)
    else:
        try:
            if os.path.isdir(fullpath):
                try:
                    shutil.rmtree(fullpath)
                except:
                    os.unlink(fullpath)
            else:
                os.remove(fullpath)
            return "success"
        except Exception as err:
            return str(err)
        
@app.route('/symlink', methods = ["PUT"])
def createSymlink():
    try:
        savepath = os.path.join(app.root_path, "users", request.form['username'], request.form['jobname'])
        os.makedirs(savepath, exist_ok=True)
        symlink = request.form['symlink']
        os.symlink(symlink, os.path.join(savepath, os.path.basename(symlink)))
        return "success"
    except Exception as err:
        return str(err)
    
@app.route('/run/<path:subpath>', methods = ["PUT"])
def runRearr(subpath):
    try:
        project_path = os.path.dirname(app.root_path)
        fastqfile = os.path.join(app.root_path, subpath)
        ref1, ref2, cut1, cut2, NGGCCNtype1, NGGCCNtype2 = request.form['ref1'], request.form['ref2'], request.form['cut1'], request.form['cut1'], request.form['NGGCCNtype1'], request.form['NGGCCNtype2']
        cmd = f'''{os.path.join(project_path, "rearr_run.sh")} {fastqfile} {ref1} {ref2} {cut1} {cut2} {NGGCCNtype1} {NGGCCNtype2}'''
        with open(f"{fastqfile}.stdout", "w") as sdo, open(f"{fastqfile}.stderr", "w") as sde:
            results = subprocess.run(cmd, stdout = sdo, stderr = sde, shell=True, executable="/bin/bash")
        return str(results)
    except Exception as err:
        return str(err)
    
if __name__ == "__main__":
    app.run()