import os, shutil, subprocess, sys
import concurrent.futures
from flask import Flask, render_template, request, send_file, jsonify
from werkzeug.utils import secure_filename
# from flask_autoindex import AutoIndex
# from flaskwebgui import FlaskUI

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

# cpulimit = max(os.cpu_count() * 3 // 4, 2)
cpulimit = 1
executor = concurrent.futures.ThreadPoolExecutor(max_workers=cpulimit)
sep = "--xy9b92jiKcrG--"

app = CustomFlask(__name__)
# AutoIndex(app)
# firefox = subprocess.run("which firefox", capture_output=True, shell=True).stdout.decode().rstrip()
# ui = FlaskUI(app=app, server="flask", browser_path=firefox, width=500, height=500)


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
    runnings = os.listdir(os.path.join(app.root_path, "running_query"))
    for running in runnings: # debug
        sys.stderr.write(running + "\n") # debug
    for fsterm in os.listdir(os.path.join(rootpath, "users", subpath)):
        if os.path.isdir(os.path.join(rootpath, "users", subpath, fsterm)):
            links[fsterm + "/"] = os.path.join("/users", subpath, fsterm)
        else:
            if fsterm.endswith(".fq.gz") or fsterm.endswith(".fastq.gz") or fsterm.endswith(".fq") or fsterm.endswith(".fastq"):
                runname = os.path.join("users", subpath, fsterm).replace("/", ".")
                state = "new"
                for qstate in ["running", "wait", "finish", "fail"]:
                    sys.stderr.write(f"{runname}.{qstate}\n") # debug
                    if f"{runname}.{qstate}" in runnings:
                        state = qstate
                        break
                links[fsterm] = os.path.join("/users", subpath, fsterm) + sep + state
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

def start_run(runname, fastqfile, cmd):
    os.remove(os.path.join(app.root_path, "running_query", f"{runname}.wait"))
    with open(os.path.join(app.root_path, "running_query", f"{runname}.running"), "w") as fd:
        pass
    try:
        with open(f"{fastqfile}.stdout", "w") as sdo, open(f"{fastqfile}.stderr", "w") as sde:
            _ = subprocess.run(cmd, stdout = sdo, stderr = sde, shell=True, executable="/bin/bash")
        with open(os.path.join(app.root_path, "running_query", f"{runname}.finish"), "w") as fd:
            pass
    except:
        with open(os.path.join(app.root_path, "running_query", f"{runname}.fail"), "w") as fd:
            pass
    os.remove(os.path.join(app.root_path, "running_query", f"{runname}.running"))

def run_controller(subpath, fastqfile, cmd):
    try:
        runname = subpath.replace("/",".")
        runnings = os.listdir(os.path.join(app.root_path, "running_query"))
        if f"{runname}.running" in runnings:
            raise Exception("job is already running")
        elif f"{runname}.wait" in runnings:
            raise Exception("job is already in query")
        else:
            for qstate in ["fail", "finish"]:
                if f"{runname}.{qstate}" in runnings:
                    os.remove(os.path.join(app.root_path, "running_query", f"{runname}.{qstate}"))
        with open(os.path.join(app.root_path, "running_query", f"{runname}.wait"), "w") as fd:
            pass
        executor.submit(start_run, runname, fastqfile, cmd)
        return "your job is added to queue"
    except Exception as err:
        with open(os.path.join(app.root_path, "running_query", f"{runname}.fail"), "w") as fd:
            pass
        raise err

@app.route('/run/<path:subpath>', methods = ["PUT"])
def runRearr(subpath):
    try:
        project_path = os.path.dirname(app.root_path)
        fastqfile = os.path.join(app.root_path, subpath)
        ref1, ref2, cut1, cut2, NGGCCNtype1, NGGCCNtype2 = request.form['ref1'], request.form['ref2'], request.form['cut1'], request.form['cut1'], request.form['NGGCCNtype1'], request.form['NGGCCNtype2']
        cmd = f'''{os.path.join(project_path, "rearr_run.sh")} {fastqfile} {ref1} {ref2} {cut1} {cut2} {NGGCCNtype1} {NGGCCNtype2}'''
        results = run_controller(subpath, fastqfile, cmd)
        return str(results)
    except Exception as err:
        return str(err)



if __name__ == "__main__":
    app.run()
    # app.run(debug=True)
    # if len(sys.argv) > 1 and sys.argv[1] == "UI":
    #     ui.run()
    # else:
    #     app.run()