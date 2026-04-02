#!/usr/bin/env python

import argparse
import os
import pathlib
import shutil
import uuid

import waitress
from celery.result import AsyncResult
from flask import (
    Flask,
    render_template,
    request,
    send_file,
    send_from_directory,
    session,
)

from .tasks import (
    celeryApp,
    celeryBuildMarker,
    celeryDefaultDirection,
    celeryDemultiplex,
    celeryPostProcess,
    celeryRearrange,
    celeryRemoveDuplicates,
)

parser = argparse.ArgumentParser()
parser.add_argument("-b", "--broker", type=str, required=True, help="TEXT")
parser.add_argument("--result-backend", type=str, required=True, help="TEXT")
args = parser.parse_args()

# set_default make default_app = celeryApp. This is necessary because in flask function, _tls.current_app = None, which causes celery.current_app to be default_app. For more details, see https://stackoverflow.com/questions/26527214/why-celery-current-app-refers-the-default-instance-inside-flask-view-functions.
celeryApp.set_default()

celeryApp.conf.update(
    broker_url=args.broker,
    result_backend=args.result_backend,
)

flaskApp = Flask(
    import_name="flask_app",
    static_folder="vue_project/dist/assets",
    template_folder="vue_project/dist",
)


flaskApp.config.update(
    MAX_CONTENT_LENGTH=10 * 1024 * 1024 * 1024,
    SECRET_KEY=b"913d1c26d46f82f119662371bc71efa4de6902a8ee2a378a6f342b8eb39b2d52",
)


@flaskApp.route("/favicon.ico")
def favicon():
    return (
        send_from_directory(
            pathlib.Path(flaskApp.import_name) / "vue_project/dist",
            "favicon.ico",
            mimetype="favicon.ico",
        ),
        200,
    )


@flaskApp.get("/")
def homePage():
    if "uuid" not in session:
        session["uuid"] = uuid.uuid4().hex
    session_path = os.path.join(flaskApp.import_name, "tmp", session["uuid"])
    os.makedirs(session_path, exist_ok=True)
    return render_template("index.html")


@flaskApp.route("/stop")
def stop():
    session_path = os.path.join(flaskApp.import_name, "tmp", session["uuid"])
    if os.path.exists(session_path):
        shutil.rmtree(session_path)
    return "STOP", 200


@flaskApp.put("/upload")
def upload():
    file = request.files["file[]"]
    file_path = os.path.join(
        flaskApp.import_name, "tmp", session["uuid"], file.filename
    )
    file.save(file_path)
    return file.filename, 200


# This api is accessed by href. Firefox will add trailing slash for cached pages. To avoid using strict_slashes=False, add a slash to this api.
@flaskApp.get("/download/<string:filename>/")
def download(filename):
    file_path = os.path.join("tmp", session["uuid"], filename)
    try:
        return (
            send_file(file_path, as_attachment=True),
            200,
        )
    except FileNotFoundError:
        return "ACCEPTED", 202


# This api is accessed by href. Firefox will add trailing slash for cached pages. To avoid using strict_slashes=False, add a slash to this api.
@flaskApp.get("/inspect/<string:taskId>/")
def inspect(taskId):
    result = AsyncResult(taskId)
    return result.status, 200


@flaskApp.put("/runJob/removeDuplicates")
def removeDuplicates():
    session_path = os.path.join(flaskApp.import_name, "tmp", session["uuid"])
    fastqFiles = [
        os.path.join(session_path, file["value"])
        for file in request.get_json()[".fastq files"]
    ]
    rmDupFile = os.path.join(session_path, "rearr.noDup")
    result = celeryRemoveDuplicates.delay(fastqFiles, rmDupFile)
    return {".noDup file": {"taskId": result.id, "value": "rearr.noDup"}}


@flaskApp.put("/runJob/buildMarker")
def buildMarker():
    session_path = os.path.join(flaskApp.import_name, "tmp", session["uuid"])
    markers = [
        os.path.join(session_path, file["value"])
        for file in request.get_json()[".fasta files"]
    ]
    result = celeryBuildMarker.delay(markers)
    return {
        "auxiliaries": [
            {
                "markerIndex": {
                    ext: {
                        "taskId": result.id,
                        "value": f"{os.path.basename(marker)}.{ext}",
                    }
                    for ext in [
                        "1.bt2",
                        "2.bt2",
                        "3.bt2",
                        "4.bt2",
                        "rev.1.bt2",
                        "rev.2.bt2",
                    ]
                }
            }
            for marker in markers
        ]
    }


@flaskApp.put("/runJob/demultiplex")
def demultiplex():
    json = request.get_json()
    session_path = os.path.join(flaskApp.import_name, "tmp", session["uuid"])
    rmDupFile = os.path.join(session_path, json[".noDup file"]["value"])
    markerIndices = [
        os.path.join(session_path, auxiliary["markerIndex"]["1.bt2"]["value"][:-6])
        for auxiliary in json["auxiliaries"]
    ]
    minScores = [
        str(auxiliary["minScore"]["value"]) for auxiliary in json["auxiliaries"]
    ]
    demultiplexFile = os.path.splitext(rmDupFile)[0] + ".demultiplex"
    result = celeryDemultiplex.delay(
        rmDupFile, markerIndices, minScores, demultiplexFile
    )
    return {
        ".demultiplex file": {
            "taskId": result.id,
            "value": os.path.basename(demultiplexFile),
        }
    }


@flaskApp.put("/runJob/postProcess")
def postProcess():
    json = request.get_json()
    demultiplexFile = os.path.join(
        flaskApp.import_name, "tmp", session["uuid"], json[".demultiplex file"]["value"]
    )
    toMapFile = os.path.splitext(demultiplexFile)[0] + ".post"
    result = celeryPostProcess.delay(demultiplexFile, toMapFile)
    return {".post file": {"taskId": result.id, "value": os.path.basename(toMapFile)}}


@flaskApp.put("/runJob/rearrange")
def rearrange():
    json = request.get_json()
    session_path = os.path.join(flaskApp.import_name, "tmp", session["uuid"])
    toMapFile = os.path.join(session_path, json[".post file"]["value"])
    refFile = os.path.join(session_path, json[".ref file"]["value"])
    directionFile = os.path.join(session_path, json[".direct file"]["value"])
    alignFile = os.path.splitext(toMapFile)[0] + ".alg"
    result = celeryRearrange.delay(
        toMapFile,
        refFile,
        directionFile,
        json["align scores"]["s0"]["value"],
        json["align scores"]["s1"]["value"],
        json["align scores"]["s2"]["value"],
        json["align scores"]["u"]["value"],
        json["align scores"]["v"]["value"],
        json["align scores"]["ru"]["value"],
        json["align scores"]["rv"]["value"],
        json["align scores"]["qu"]["value"],
        json["align scores"]["qv"]["value"],
        alignFile,
    )
    return {".alg file": {"taskId": result.id, "value": os.path.basename(alignFile)}}


@flaskApp.put("/runJob/defaultDirection")
def defaultDirection():
    refFile = os.path.join(
        flaskApp.import_name,
        "tmp",
        session["uuid"],
        request.get_json()[".ref file"]["value"],
    )
    directionFile = os.path.splitext(refFile)[0] + ".direct"
    result = celeryDefaultDirection.delay(refFile, directionFile)
    return {
        ".direct file": {"taskId": result.id, "value": os.path.basename(directionFile)}
    }


@flaskApp.put("/runJob/indexGenome")
def indexGenome():
    genomeFile = os.path.join(
        flaskApp.import_name,
        "tmp",
        session["uuid"],
        request.get_json()["genome file"]["value"],
    )
    result = celeryBuildMarker.delay([genomeFile])
    return {
        "genome index": {
            ext: {"taskId": result.id, "value": f"{os.path.basename(genomeFile)}.{ext}"}
            for ext in ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]
        }
    }


if __name__ == "__main__":
    waitress.serve(
        flaskApp,
        url_prefix="/rearr",
        max_request_body_size=flaskApp.config["MAX_CONTENT_LENGTH"],
    )
