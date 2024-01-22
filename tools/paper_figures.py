import matplotlib.pyplot, numpy, matplotlib.animation, matplotlib.patches, matplotlib.transforms, matplotlib.patches, matplotlib.path, matplotlib.offsetbox, matplotlib.image, PIL.Image

#####################################
# stagger cut
#####################################

fig, axs = matplotlib.pyplot.subplot_mosaic([["(a)", "(an)"], ["(b)", "(bn)"], ["(c)", "(cn)"]], width_ratios=[0.9,0.1], gridspec_kw={"hspace": 0, "wspace": 0}, figsize=[8,1.3])
DNA = numpy.random.choice(['A', 'C', 'G', 'T'], 60)
DNA[34:36] = "G"
DNA = "".join(DNA)
rcDNA = DNA.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c").upper()
upscissor = numpy.asarray(PIL.Image.open('pictures/Scissors2.png').rotate(90))
downscissor = numpy.asarray(PIL.Image.open('pictures/Scissors2.png').rotate(-90))
for label, ax in axs.items():
    ax.set_xlim([-1,1])
    ax.set_ylim([-1,1])
    ax.set_axis_off()
    if label == "(a)":
        unit = 1.05/33
        ax.text(0, 0, f''' {DNA[:33]+"   "+DNA[36:]} \n {rcDNA} ''', va='center', ha='center', fontfamily='courier new')
        ax.text(0, 0, f''' {" "*33+DNA[33:36]+" "*24} \n''', va='center', ha='center', fontweight="bold", fontfamily='courier new')
        ax.add_artist(matplotlib.offsetbox.AnnotationBbox(matplotlib.offsetbox.OffsetImage(upscissor, zoom=0.015), (0, -0.5), frameon=False))
        ax.add_artist(matplotlib.offsetbox.AnnotationBbox(matplotlib.offsetbox.OffsetImage(downscissor, zoom=0.015), (-unit, 0.5), frameon=False))
        ax.text(-1, 1, label, va="top", ha="left", fontfamily="arial")
    elif label == "(an)":
        ax.text(-1, 0, "wild type", va='center', ha='left', fontfamily='arial')
    elif label == "(b)":
        ax.text(0, 0, f''' {DNA[:29]+" "+DNA[29:33]+"   "+DNA[36:]}\n {rcDNA[:30]+" "+rcDNA[30:]}''', va='center', ha='center', fontfamily='courier new')
        ax.text(0, 0, f''' {" "*34+DNA[33:36]+" "*24}\n''', va='center', ha='center', fontweight="bold", fontfamily='courier new')
        ax.text(-1, 1, label, va="top", ha="left", fontfamily="arial")
    elif label == "(bn)":
        ax.text(-1, 0, "DSB", va='center', ha='left', fontfamily='arial')
    elif label == "(c)":
        ax.text(0, 0, f''' {DNA[:29]+" "+DNA[29:33]+"   "+DNA[36:]}\n {rcDNA[:30]+" "+rcDNA[30:]}''', va='center', ha='center', fontfamily='courier new')
        ax.text(0, 0, f''' {" "*34+DNA[33:36]+" "*24}\n''', va='center', ha='center', fontweight="bold", fontfamily='courier new')
        ax.text(0, 0, f''' {" "*29+DNA[29]+" "*31}\n {" "*30+rcDNA[29]+" "*30}''', va='center', ha='center', color="red", fontfamily='courier new')
        ax.text(-1, 1, label, va="top", ha="left", fontfamily="arial")
    elif label == "(cn)":
        ax.text(-1, 0, "1bp templated insertion", va='center', ha='left', fontfamily='arial')
    fig.savefig("paper_figures/stagger_cut.png", bbox_inches="tight", pad_inches=0, dpi=300)

############################################################
# predictable insertion, deletion, random insertion, templated insertion, templated deletion, small indel disaster, microhomology
############################################################

fig, axs = matplotlib.pyplot.subplot_mosaic([["(a)", "(an)"], ["(b)", "(bn)"], ["(c)", "(cn)"], ["(d)", "(dn)"]], width_ratios=[0.9,0.1], gridspec_kw={"hspace": 0, "wspace": 0}, figsize=[8,8])
for label, ax in axs.items():
    ax.set_xlim([-1,1])
    ax.set_ylim([-1,1])
    ax.set_axis_off()
    if label == "(a)":
        DNA = numpy.random.choice(['A', 'C', 'G', 'T'], 60)
        DNA[34:36] = "G"
        DNA = "".join(DNA)
        riDNA = "".join(numpy.random.choice(['A', 'C', 'G', 'T'], 3))
        while riDNA[0] == DNA[28] or riDNA[-1] == DNA[31]:
            riDNA = "".join(numpy.random.choice(['A', 'C', 'G', 'T'], 3))
        ax.text(0, 0, f''' {DNA[:33]+"   "+DNA[36:]} \n {" "*28+DNA[28:32]+" "*28} \n {" "*10+DNA[10:28]+riDNA+"-"+DNA[32:50]+" "*10} \n {" "*10+DNA[10:28]+" "*32} \n {" "*28+riDNA+" "*29} \n {" "*32+DNA[32:50]+" "*10} \n {" "*28+DNA[28:30]+" "*30} \n {" "*30+DNA[30:32]+" "*28} ''', va='center', ha='center', fontfamily='courier new', linespacing=1.1)
        ax.text(0, 0, f''' {" "*33+DNA[33:36]+" "*24} \n\n\n\n\n\n\n''', va='center', ha='center', fontfamily='courier new', linespacing=1.1, fontweight="bold")
        ax.text(-1, 1, label, va="top", ha="left", fontfamily="arial")
        ax.plot([0,0],[-1,1],c="grey",ls="--")
    elif label == "(an)":
        ax.text(-1, 0, f''' ref \n deletion \n read \n block1 \n random insertion \n block2 \n upstream deletion \n downstream deletion ''', fontfamily="arial", va='center', ha='left', linespacing=1)
    elif label == "(b)":
        DNA = numpy.random.choice(['A', 'C', 'G', 'T'], 60)
        DNA[34:36] = "G"
        DNA = "".join(DNA)
        riDNA = "".join(numpy.random.choice(['A', 'C', 'G', 'T'], 3))
        while riDNA[0] == DNA[29] or riDNA[-1] == DNA[25]:
            riDNA = "".join(numpy.random.choice(['A', 'C', 'G', 'T'], 3))
        ax.text(0, 0, f''' {DNA[:33]+"   "+DNA[36:]} \n {" "*26+DNA[26:29]+" "*31} \n {" "*10+DNA[10:29]+" "*31} \n {" "*29+riDNA+" "*28} \n {" "*26+DNA[26:50]+" "*10} \n {" "*29+DNA[29:30]+" "*30} \n {" "*26+DNA[26:30]+" "*30} ''', va='center', ha='center', fontfamily='courier new', linespacing=1.1)
        ax.text(0, 0, f''' {" "*33+DNA[33:36]+" "*24} \n\n\n\n\n\n''', va='center', ha='center', fontfamily='courier new', linespacing=1.1, fontweight="bold")
        ax.text(-1, 1, label, va="top", ha="left", fontfamily="arial")
        ax.plot([0,0],[-1,1],c="grey",ls="--")
    elif label == "(bn)":
        ax.text(-1, 0, f''' ref \n templated insertion \n block1 \n random insertion \n block2 \n upstream deletion \n downstream templated insertion ''', fontfamily="arial", va='center', ha='left', linespacing=1)
    elif label == "(c)":
        DNA = numpy.random.choice(['A', 'C', 'G', 'T'], 60)
        DNA[34:36] = 'G'
        DNA = "".join(DNA)
        ax.text(0, 0, f''' {DNA[:33]+"   "+DNA[36:]} \n {" "*10+DNA[10:23]+"-"+DNA[24:50]+" "*10} \n {" "*10+DNA[10:23]+"-"+DNA[24:30]+" "*30} \n {" "*30+DNA[30:50]+" "*10} \n {" "*10+DNA[10:23]+" "*37} \n {" "*24+DNA[24:50]+" "*10} ''', va='center', ha='center', fontfamily='courier new', linespacing=1.1)
        ax.text(0, 0, f''' {" "*33+DNA[33:36]+" "*24} \n\n\n\n\n''', va='center', ha='center', fontfamily='courier new', linespacing=1.1, fontweight="bold")
        ax.text(-1, 1, label, va="top", ha="left", fontfamily="arial")
        ax.plot([0,0],[-1,1],c="grey",ls="--")
    elif label == "(cn)":
        ax.text(-1, 0, f''' ref \n read \n block1 \n block2 \n block1 (false) \n block2 (false) ''', fontfamily="arial", va='center', ha='left', linespacing=1)
    elif label == "(d)":
        DNA = numpy.random.choice(['A', 'C', 'G', 'T'], 60)
        while DNA[29] == "G" or DNA[25] == DNA[31]:
            DNA = numpy.random.choice(['A', 'C', 'G', 'T'], 60)
        DNA[34:36] = "G"
        DNA = "".join(DNA)
        micro = "".join(numpy.random.choice(['A', 'C', 'G', 'T'], 2)) + "G"
        DNA = DNA[:26] + micro + DNA[29:32] + micro + DNA[35:]
        ax.text(0, 0, f''' {DNA[:26]+"   "+DNA[29:32]+"    "+DNA[36:]} \n {" "*10+DNA[10:26]+" "*6+DNA[32:50]+" "*10} \n {" "*10+DNA[10:27]+" "*6+DNA[33:50]+" "*10} \n {" "*10+DNA[10:28]+" "*6+DNA[34:50]+" "*10} \n {" "*10+DNA[10:29]+" "*6+DNA[35:50]+" "*10} ''', va='center', ha='center', fontfamily='courier new', linespacing=1.1)
        ax.text(0, 0, f''' {" "*26+DNA[26:29]+"   "+DNA[32:33]+" "*27} \n\n\n\n''', va='center', ha='center', fontfamily='courier new', linespacing=1.1, color="red")
        ax.text(0, 0, f''' {" "*33+DNA[33:35]+" "*25} \n\n\n\n''', va='center', ha='center', fontfamily='courier new', linespacing=1.1, color="red", fontweight="bold")
        ax.text(0, 0, f''' {" "*35+DNA[35:36]+" "*24} \n\n\n\n''', va='center', ha='center', fontfamily='courier new', linespacing=1.1, fontweight="bold")
        ax.text(-1, 1, label, va="top", ha="left", fontfamily="arial")
        ax.plot([0,0],[-1,1],c="grey",ls="--")
    elif label == "(dn)":
        ax.text(-1, 0, f''' ref \n align1 \n align2 \n align3 \n align4 ''', fontfamily="arial", va='center', ha='left', linespacing=1)
    fig.savefig("paper_figures/indel_concept_disaster_microhomology.png", bbox_inches="tight", pad_inches=0, dpi=300)

##########################################################
# indel to replacement and templated insertion and double cut defects
##########################################################

fig, axs = matplotlib.pyplot.subplot_mosaic([["(a)", "(an)"], ["(b)", "(bn)"], ["(c)", "(cn)"]], width_ratios=[0.9,0.1], gridspec_kw={"hspace": 0, "wspace": 0}, figsize=[8,5])
DNA = numpy.random.choice(['A', 'C', 'G', 'T'], 60)
DNA[34:36] = "G"
DNA = "".join(DNA)
rep = "".join(numpy.random.choice(['A', 'C', 'G', 'T'], 3))
while rep[0] == DNA[24] or rep[1] == DNA[25] or rep[2] == DNA[26] or rep[0] == DNA[30]:
    rep = "".join(numpy.random.choice(['A', 'C', 'G', 'T'], 3))
DNA2 = "".join(numpy.random.choice(['A', 'C', 'G', 'T'], 60))
while DNA2[27] == DNA[27] or DNA2[28] == DNA[28] or DNA2[29] == DNA[29]:
    DNA2 = "".join(numpy.random.choice(['A', 'C', 'G', 'T'], 60))

for label, ax in axs.items():
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_axis_off()
    if label == "(a)":
        ax.set_title("indel", y=0.8, fontfamily="arial")
        ax.text(0, 0, f''' {DNA[:33]+"   "+DNA[36:]} \n {DNA[:24]+" "*36} \n {" "*24+rep+" "*33} \n {" "*27+DNA[27:]} \n {DNA[:24]+rep+DNA[27:]} ''', va='center', ha='center', fontfamily='courier new', linespacing=1.1)
        ax.text(0, 0, f''' {" "*33+DNA[33:36]+" "*24} \n\n\n\n''', va='center', ha='center', fontfamily='courier new', linespacing=1.1, fontweight="bold")
        ax.text(-1, 1, label, va="top", ha="left", fontfamily="arial")
        ax.plot([0, 0], [-1, 1], c="grey", ls="--")
    elif label == "(an)":
        ax.text(-1, 0, f''' reference \n block1 \n random insertion \n block2 \n replacement (false) ''', fontfamily="arial", va='center', ha='left', linespacing=1)
    elif label == "(b)":
        ax.set_title("templated insertion", y=0.8, fontfamily="arial")
        ax.text(0, 0, f''' {DNA[:33]+"   "+DNA[36:]} \n {DNA[:30]+" "*30} \n {" "*24+rep+" "*33} \n {" "*27+DNA[27:]} \n {" "*24+rep+DNA[27:30]+" "*30} \n {" "*30+DNA[30:]} ''', va='center', ha='center', fontfamily='courier new', linespacing=1.1)
        ax.text(0, 0, f''' {" "*33+DNA[33:36]+" "*24} \n\n\n\n\n''', va='center', ha='center', fontfamily='courier new', linespacing=1.1, fontweight="bold")
        ax.text(-1, 1, label, va="top", ha="left", fontfamily="arial")
        ax.plot([0, 0], [-1, 1], c="grey", ls="--")
    elif label == "(bn)":
        ax.text(-1, 0, f''' reference \n block1 \n random insertion \n block2 \n random insertion (false) \n block2 (false) ''', fontfamily="arial", va='center', ha='left', linespacing=1)
    elif label == "(c)":
        ax.set_title("double cut", y=0.8, fontfamily="arial")
        ax.text(0, 0, f''' {DNA[:33]+"   "+DNA[36:]} \n {DNA2} \n {DNA[:24]+" "*36} \n {" "*27+DNA[27:30]+" "*30} \n {" "*30+DNA2[30:]} \n {DNA[:30]+DNA2[30:]} \n {" "*27+DNA[27:30]+DNA2[30:]} ''', va='center', ha='center', fontfamily='courier new', linespacing=1.1)
        ax.text(0, 0, f''' {" "*33+DNA[33:36]+" "*24} \n\n\n\n\n\n''', va='center', ha='center', fontfamily='courier new', linespacing=1.1, fontweight="bold")
        ax.text(-1, 1, label, va="top", ha="left", fontfamily="arial")
        ax.plot([0, 0], [-1, 1], c="grey", ls="--")
    elif label == "(cn)":
        ax.text(-1, 0, f''' locus1 \n locus2 \n block1 \n random insertion \n block2 \n locus1+2 \n block2 (false) ''', fontfamily="arial", va='center', ha='left', linespacing=1)
fig.savefig("paper_figures/defects.png", bbox_inches="tight", pad_inches=0, dpi=300)

########################################
# double cut
########################################

fig, axs = matplotlib.pyplot.subplot_mosaic([["(a)", "(an)"], ["(b)", "(bn)"], ["(c)", "(cn)"], ["(d)", "(dn)"]], width_ratios=[0.9,0.1], height_ratios=[0.2, 0.2, 0.2, 0.4], gridspec_kw={"hspace": 0, "wspace": 0}, figsize=(8, 3))

DNA = numpy.random.choice(['A', 'C', 'G', 'T'], 60)
DNA[24:26] = "G"
DNA[44:46] = "G"
DNA = "".join(DNA)
cDNA = DNA.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c").upper()

unit=1.04
for label, ax in axs.items():
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_axis_off()
    if label == "(a)":
        ax.text(0, 0, f''' {DNA[:23]+"   "+DNA[26:43]+"   "+DNA[46:]} \n {cDNA} ''', va='center', ha='center', fontfamily='courier new', linespacing=1.1)
        ax.text(0, 0, f''' {" "*23+DNA[23:26]+" "*17+DNA[43:46]+" "*14} \n''', va='center', ha='center', fontfamily='courier new', fontweight="bold", linespacing=1.1)
        ax.text(-1, 1, label, va="top", ha="left", fontfamily="arial")
        ax.plot([-0.3, -0.3], [-1,1], [0.3, 0.3], [-1,1], c="grey", ls="--")
    elif label == "(an)":
        ax.text(-1, 0, f''' wild type ''', fontfamily="arial", va='center', ha='left', linespacing=1)
    elif label == "(b)":
        ax.text(0, 0, f''' {DNA[:20]+"-"*20+DNA[40:43]+"   "+DNA[46:]} \n {cDNA[:20]+"-"*20+cDNA[40:]} ''', va='center', ha='center', fontfamily='courier new', linespacing=1.1)
        ax.text(0, 0, f''' {" "*43+DNA[43:46]+" "*14} \n''', va='center', ha='center', fontfamily='courier new', fontweight="bold", linespacing=1.1)
        ax.text(-1, 1, label, va="top", ha="left", fontfamily="arial")
        ax.plot([-0.3, -0.3], [-1,1], [0.3, 0.3], [-1,1], c="grey", ls="--")
    elif label == "(bn)":
        ax.text(-1, 0, f''' deletion ''', fontfamily="arial", va='center', ha='left', linespacing=1)
    elif label == "(c)":
        ax.text(0, 0, f''' {DNA[:20]+cDNA[39:19:-1]+DNA[40:43]+"   "+DNA[46:]} \n {cDNA[:20]+DNA[39:25:-1]+"   "+DNA[22:19:-1]+cDNA[40:]} ''', va='center', ha='center', fontfamily='courier new', linespacing=1.1)
        ax.text(0, 0, f''' {" "*43+DNA[43:46]+" "*14} \n {" "*34+DNA[25:22:-1]+" "*23} ''', va='center', ha='center', fontfamily='courier new', fontweight="bold", linespacing=1.1)
        ax.text(-1, 1, label, va="top", ha="left", fontfamily="arial")
        ax.plot([-0.3, -0.3], [-1,1], [0.3, 0.3], [-1,1], c="grey", ls="--")
    elif label == "(cn)":
        ax.text(-1, 0, f''' inversion ''', fontfamily="arial", va='center', ha='left', linespacing=1)
    elif label == "(d)":
        ax.text(0, 0.5, f''' {DNA[:23]+"   "+DNA[26:40]+" "*20} \n {cDNA[:40]+" "*20} ''', va='center', ha='center', fontfamily='courier new', linespacing=1.1)
        ax.text(0, 0.5, f''' {" "*23+DNA[23:26]+" "*34} \n''', va='center', ha='center', fontfamily='courier new', fontweight="bold", linespacing=1.1)
        ax.text(0, -0.5, f''' {" "*20+DNA[20:23]+"   "+DNA[26:43]+"   "+DNA[46:]} \n {" "*20+cDNA[20:]} ''', va='center', ha='center', fontfamily='courier new', linespacing=1.1)
        ax.text(0, -0.5, f''' {" "*23+DNA[23:26]+" "*17+DNA[43:46]+" "*14} \n''', va='center', ha='center', fontfamily='courier new', fontweight="bold", linespacing=1.1)
        ax.text(-1, 1, label, va="top", ha="left", fontfamily="arial")
        ax.plot([-0.3, -0.3], [-1,1], [0.3, 0.3], [-1,1], c="grey", ls="--")
        Hshift, Vshift = 1, 0.3
        path = matplotlib.path.Path([(-0.3, -0.5), (-0.3-Hshift, -0.5 + Vshift), (0.3+Hshift, 0.5-Vshift), (0.3, 0.5)], [1, 4, 4, 4])
        ax.add_patch(matplotlib.patches.PathPatch(path, color="grey", fill=False, ls="--"))
    elif label == "(dn)":
        ax.text(-1, 0, f''' duplication ''', fontfamily="arial", va='center', ha='left', linespacing=1)
fig.savefig("paper_figures/double_cut_rearrangement.png", bbox_inches="tight", pad_inches=0, dpi=300)

