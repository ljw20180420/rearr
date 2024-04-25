from shiny import App, req, render, ui

import subprocess
import os
import tempfile
from hic2cool import hic2cool_convert

app_ui = ui.page_fixed(
    ui.input_file(id = "upload", label = "upload file"),
    ui.input_select(id = "format", label = "which format to trans", choices = [".mcool", ".hic"], selected = ".mcool"),
    ui.download_button(id = "download", label = "download file"),
    title="Convert .hic to .mcool",
)

def server(input, output, session):
    @render.download
    def download():
        req(input.upload())
        inExt = os.path.splitext(input.upload()[0]["name"])[1]
        outExt = input.format()
        req(inExt != outExt)
        tempdir = tempfile.mkdtemp()
        helperPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "../helpers")
        inputFile = input.upload()[0]['datapath']
        outputFile = os.path.join(tempdir, f'''{os.path.splitext(input.upload()[0]["name"])[0]}{outExt}''')
        if inExt == ".allValidPairs" and outExt == ".hic":
            subprocess.check_output(f'''{os.path.join(helperPath, "hicpro2juicebox.sh")} -i {inputFile} -o {tempdir} -g {os.path.join(helperPath, "hg19.chrom.sizes")} -j {os.path.join(helperPath, "juicer_tools_1.22.01.jar")}''', shell=True)
            os.rename(os.path.join(tempdir, f"{os.path.basename(inputFile)}.hic"), outputFile)
        elif inExt == ".allValidPairs" and outExt == ".mcool":
            raise Exception("not support yet")
        elif inExt == ".hic" and outExt == ".mcool":
            hic2cool_convert(inputFile, outputFile, resolution = 0)
        elif inExt == ".mcool" and outExt == ".hic":
            raise Exception("not support yet")
        return outputFile

app = App(app_ui, server)