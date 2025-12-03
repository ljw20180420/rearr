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
    celeryRearrange,
    celeryRemoveDuplicates,
    celerySxGetMarkers,
    celerySxGetReference,
    celerySxPostProcess,
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
    if "tempdir" not in session:
        session["tempdir"] = (
            pathlib.Path(flaskApp.import_name) / "tmp" / uuid.uuid4().hex
        ).as_posix()
    os.makedirs(session["tempdir"], exist_ok=True)
    return render_template("index.html")


@flaskApp.route("/stop")
def stop():
    if os.path.exists(session["tempdir"]):
        shutil.rmtree(session["tempdir"])
    return "STOP", 200


@flaskApp.put("/upload")
def upload():
    file = request.files["file[]"]
    file.save(pathlib.Path(session["tempdir"]) / file.filename)
    return file.filename, 200


# This api is accessed by href. Firefox will add trailing slash for cached pages. To avoid using strict_slashes=False, add a slash to this api.
@flaskApp.get("/download/<string:filename>/")
def download(filename):
    try:
        return (
            send_file(pathlib.Path(session["tempdir"]) / filename, as_attachment=True),
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
    fastqFiles = [
        (pathlib.Path(session["tempdir"]) / file["value"]).as_posix()
        for file in request.get_json()[".fastq files"]
    ]
    rmDupFile = (pathlib.Path(session["tempdir"]) / "rearr.noDup").as_posix()
    result = celeryRemoveDuplicates.delay(fastqFiles, rmDupFile)
    return {".noDup file": {"taskId": result.id, "value": "rearr.noDup"}}


@flaskApp.put("/runJob/buildMarker")
def buildMarker():
    markers = [
        (pathlib.Path(session["tempdir"]) / file["value"]).as_posix()
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
    rmDupFile = (
        pathlib.Path(session["tempdir"]) / json[".noDup file"]["value"]
    ).as_posix()
    markerIndices = [
        (
            pathlib.Path(session["tempdir"])
            / auxiliary["markerIndex"]["1.bt2"]["value"][:-6]
        ).as_posix()
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


@flaskApp.put("/runJob/sxPostProcess")
def sxPostProcess():
    json = request.get_json()
    demultiplexFile = (
        pathlib.Path(session["tempdir"]) / json[".demultiplex file"]["value"]
    ).as_posix()
    toMapFile = os.path.splitext(demultiplexFile)[0] + ".post"
    result = celerySxPostProcess.delay(
        demultiplexFile, json["minimal base number"]["value"], toMapFile
    )
    return {".post file": {"taskId": result.id, "value": os.path.basename(toMapFile)}}


@flaskApp.put("/runJob/rearrange")
def rearrange():
    json = request.get_json()
    toMapFile = (
        pathlib.Path(session["tempdir"]) / json[".post file"]["value"]
    ).as_posix()
    refFile = (pathlib.Path(session["tempdir"]) / json[".ref file"]["value"]).as_posix()
    directionFile = (
        pathlib.Path(session["tempdir"]) / json[".direct file"]["value"]
    ).as_posix()
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
    refFile = (
        pathlib.Path(session["tempdir"]) / request.get_json()[".ref file"]["value"]
    ).as_posix()
    directionFile = os.path.splitext(refFile)[0] + ".direct"
    result = celeryDefaultDirection.delay(refFile, directionFile)
    return {
        ".direct file": {"taskId": result.id, "value": os.path.basename(directionFile)}
    }


@flaskApp.put("/runJob/indexGenome")
def indexGenome():
    genomeFile = (
        pathlib.Path(session["tempdir"]) / request.get_json()["genome file"]["value"]
    ).as_posix()
    result = celeryBuildMarker.delay([genomeFile])
    return {
        "genome index": {
            ext: {"taskId": result.id, "value": f"{os.path.basename(genomeFile)}.{ext}"}
            for ext in ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]
        }
    }


@flaskApp.put("/runJob/sxGetReference")
def sxGetReference():
    json = request.get_json()
    plasmid_file = (
        pathlib.Path(session["tempdir"]) / json[".csv file"]["value"]
    ).ax_posix()
    bowtie2index = (
        pathlib.Path(session["tempdir"]) / json["genome index"]["1.bt2"]["value"][:-6]
    ).as_posix()
    genome = (
        pathlib.Path(session["tempdir"]) / json["genome file"]["value"]
    ).as_posix()
    ext1up = json["extensions"]["cut1 upstream"]["value"]
    ext1down = json["extensions"]["cut1 downstream"]["value"]
    ext2up = json["extensions"]["cut2 upstream"]["value"]
    ext2down = json["extensions"]["cut2 downstream"]["value"]
    refFile = f"{plasmid_file}.ref"
    result = celerySxGetReference.delay(
        plasmid_file, genome, bowtie2index, ext1up, ext1down, ext2up, ext2down, refFile
    )
    return {".ref file": {"taskId": result.id, "value": os.path.basename(refFile)}}


@flaskApp.put("/runJob/sxGetMarkers")
def getSxMarkers():
    json = request.get_json()
    plasmid_file = (
        pathlib.Path(session["tempdir"]) / json[".csv file"]["value"]
    ).as_posix()
    targetMarker = f"{plasmid_file}.target.fa"
    pairMarker = f"{plasmid_file}.pair.fa"
    result = celerySxGetMarkers.delay(plasmid_file, targetMarker, pairMarker)
    return {
        ".fasta files": [
            {"taskId": result.id, "value": os.path.basename(targetMarker)},
            {"taskId": result.id, "value": os.path.basename(pairMarker)},
        ]
    }


if __name__ == "__main__":
    waitress.serve(
        flaskApp,
        url_prefix="/workflow",
        max_request_body_size=flaskApp.config["MAX_CONTENT_LENGTH"],
    )
