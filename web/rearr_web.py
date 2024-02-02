from flask import Flask, render_template, request

app = Flask(__name__)

@app.route("/")
def home_page():
    return render_template("rearrangement.html")

@app.route("/<user>")
def user():
    return "<p>Hello, World!</p>"

