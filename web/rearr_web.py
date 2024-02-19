import os, sys
from flask import Flask, flash, render_template, request, redirect, session, url_for
from werkzeug.utils import secure_filename
from markupsafe import escape

app = Flask(__name__, root_path=os.path.dirname(__file__))
app.secret_key = 'super secret key'
app.config['SESSION_TYPE'] = 'filesystem'

# session.init_app(app)

# app.debug = True

@app.route("/")
def home_page():
    return render_template("rearrangement.html")

@app.route('/upload', methods=['POST'])
def upload_file():
    if request.method == "POST":
        try:
            user = request.form['user']
            job = request.form['job']
            os.makedirs(os.path.join(app.root_path, "users", user, job), exist_ok=True)
            for file in request.files.getlist("file"):
                file.save(os.path.join(app.root_path, "users", user, job, secure_filename(file.filename)))            
            return f"success"
        except Exception as err:
            return f"str(err)"

# @app.route("/<user>")
# def user(user):
#     return f"Hello, {escape(user)}!"

# @app.route("/<user>/<job>")
# def job(user, job):
#     return f"Hello, {escape(job)} of {escape(user)}!"