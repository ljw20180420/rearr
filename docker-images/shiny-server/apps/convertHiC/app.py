from shiny import App, req, render, ui

import subprocess
import os
import tempfile
from hic2cool import hic2cool_convert
import cooler
import cooltools

def load_rename_add_normVec_cov_tot(file, minDis = 0, force = False):
    chr_dict = {'X' : 'chrX', 'Y' : 'chrY', 'MT' : 'chrM', 'M' : 'chrM'}
    for i in range(1, 23):
        chr_dict[f'{i}'] = f'chr{i}'

    for path in cooler.fileops.list_coolers(file):
        clr = cooler.Cooler(f'{file}::{path}')
        cooler.rename_chroms(clr, chr_dict)
        if force and ('weight' in clr.bins().columns):
            with clr.open('r+') as f:
                del f['bins']['weight']
        if 'weight' not in clr.bins().columns:
            ctmp = cooler.Cooler(f'{file}::{cooler.fileops.list_coolers(file)[-1]}')
            btmp = ctmp.bins()[:]
            ptmp = ctmp.pixels()[:]
            totalreads = sum(ptmp.loc[(ptmp.bin2_id-ptmp.bin1_id>=minDis) & (btmp.chrom[ptmp.bin1_id].reset_index(drop=True)==btmp.chrom[ptmp.bin2_id].reset_index(drop=True)), 'count'])
            with clr.open('r+') as f:
                f["bins"].create_dataset("weight", data=1./f["bins/SCALE"][:]/(totalreads**0.5)*1e4, compression="gzip", compression_opts=6)
        if 'cov_cis_raw' not in clr.bins().columns:
            cooltools.coverage(clr, ignore_diags=2, store=True)

app_ui = ui.page_navbar(
    ui.nav_panel(
        ui.tooltip(
            ".allValidPairs to .hic",
            "convert .allValidPairs to .hic"
        ),
        ui.tooltip(
            ui.input_file(id = "uploadAllValidPairs", label = ".allValidPairs"),
            "upload .allValidPairs file"
        ),
        ui.tooltip(
            ui.download_button(id = "downloadHic", label = ".hic"),
            "download .hic file"
        )
    ),
    ui.nav_panel(
        ui.tooltip(
            ".hic to .mcool",
            "convert .hic to .mcool"
        ),
        ui.tooltip(
            ui.input_file(id = "uploadHic", label = ".hic"),
            "upload .hic file"
        ),
        ui.tooltip(
            ui.download_button(id = "downloadMcool", label = ".mcool"),
            "download .mcool file"
        )
    )
)

def server(input, output, session):
    @render.download
    def downloadHic():
        req(input.uploadAllValidPairs())
        tempdir = tempfile.mkdtemp()
        helperPath = os.path.join(os.path.dirname(os.path.realpath(__file__)), "helpers")
        inputFile = input.uploadAllValidPairs()[0]['datapath']
        outputFile = os.path.join(tempdir, f'''{os.path.splitext(input.uploadAllValidPairs()[0]["name"])[0]}.hic''')
        subprocess.check_output(f'''{os.path.join(helperPath, "hicpro2juicebox.sh")} -i {inputFile} -o {tempdir} -g {os.path.join(helperPath, "hg19.chrom.sizes")} -j {os.path.join(helperPath, "juicer_tools_1.22.01.jar")}''', shell=True)
        os.rename(os.path.join(tempdir, f"{os.path.basename(inputFile)}.hic"), outputFile)
        return outputFile
    
    @render.download
    def downloadMcool():
        req(input.uploadHic())
        tempdir = tempfile.mkdtemp()
        inputFile = input.uploadHic()[0]['datapath']
        outputFile = os.path.join(tempdir, f'''{os.path.splitext(input.uploadHic()[0]["name"])[0]}.mcool''')
        hic2cool_convert(inputFile, outputFile, resolution = 0)
        load_rename_add_normVec_cov_tot(outputFile)
        return outputFile

app = App(app_ui, server)