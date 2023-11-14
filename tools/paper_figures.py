#!/usr/bin/env python

import matplotlib.pyplot, numpy, matplotlib.animation, matplotlib.patches

#####################################
# stagger cut
#####################################

DNA = numpy.random.choice(['A', 'C', 'G', 'T'], 60)
DNA[34:36] = "G"
DNA = "".join(DNA)
rcDNA = DNA.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c").upper()[::-1]
uDNA, dDNA, urcDNA, drcDNA = DNA[:29], DNA[29:], rcDNA[:30], rcDNA[30:]
PAM = dDNA[4:7]
dDNA = dDNA[:4] + "   " + dDNA[7:]

unit = 0.9
fig, ax = matplotlib.pyplot.subplots(figsize=[8,2])
uDNAtext = ax.text(-unit, 0, uDNA, verticalalignment='bottom', horizontalalignment='right', fontfamily='courier new')
dDNAtext = ax.text(-unit, 0, dDNA, verticalalignment='bottom', horizontalalignment='left', fontfamily='courier new')
PAMtext = ax.text(3 * unit, 0, PAM, verticalalignment='bottom', horizontalalignment='left', fontfamily='courier new', color = "blue")
urcDNAtext = ax.text(0, 0, urcDNA, verticalalignment='top', horizontalalignment='right', fontfamily='courier new')
drcDNAtext = ax.text(0, 0, drcDNA, verticalalignment='top', horizontalalignment='left', fontfamily='courier new')
ecut = ax.plot([0, 0], [-1, 1], c = 'grey', ls = '--')
cuttext = ax.text(-unit, 0, "|", verticalalignment='bottom', horizontalalignment='center', alpha = 0)
rccuttext = ax.text(0, 0, "|", verticalalignment='top', horizontalalignment='center', alpha = 0)
abasetext = ax.text(-unit, 0, dDNA[0], verticalalignment='bottom', horizontalalignment='left', fontfamily='courier new', alpha = 0, color = "red")
arcbasetext = ax.text(0, 0, urcDNA[-1], verticalalignment='top', horizontalalignment='left', fontfamily='courier new', alpha = 0, color = "red")
ax.set_xlim([-33,33])
ax.set_ylim([-1,1])
ax.set_axis_off()

reso = 10
def update(frame):
    if frame <= reso:
        cuttext.set_alpha(frame/reso)
        rccuttext.set_alpha(frame/reso)
    if frame == reso + reso // 2:
        dDNAtext.set_x(0)
        PAMtext.set_x(4 * unit)
        drcDNAtext.set_x(unit)
        ecut[0].set_xdata([unit, unit])
    if frame > 2*reso and frame <= 3*reso:
        abasetext.set_alpha((frame-2*reso)/reso)
        arcbasetext.set_alpha((frame-2*reso)/reso)
    if frame > 3*reso and frame <= 4*reso:
        cuttext.set_alpha((4*reso-frame)/reso)
        rccuttext.set_alpha((4*reso-frame)/reso)
    return dDNAtext, PAMtext, drcDNAtext, ecut, cuttext, rccuttext, abasetext, arcbasetext

ani = matplotlib.animation.FuncAnimation(fig=fig, func=update, frames=6*reso, interval=120)
ani.save(filename="paper_figures/stagger_cut.gif", writer="pillow")

##################################
# templated insertion
##################################

DNA = numpy.random.choice(['A', 'C', 'G', 'T'], 60)
DNA[34:36] = "G"
DNA = "".join(DNA)
insuDNA = DNA[10:30]
insdDNA = DNA[27:50]
insrDNA = "".join(numpy.random.choice(['A', 'C', 'G', 'T'], 3))
while insrDNA[0] == DNA[30] or insrDNA[-1] == DNA[26]:
    insrDNA = "".join(numpy.random.choice(['A', 'C', 'G', 'T'], 3))
uhDNA, uaDNA, tDNA, daDNA, dhDNA = DNA[:10], DNA[10:27], DNA[27:30], DNA[30:50], DNA[50:]

unit = 0.9
fig, ax = matplotlib.pyplot.subplots(figsize=[8,2])
ax.text(-20 * unit, 0.5, uhDNA, verticalalignment='bottom', horizontalalignment='right', fontfamily='courier new')
ax.text(-3 * unit, 0.5, uaDNA, verticalalignment='bottom', horizontalalignment='right', fontfamily='courier new', color = "red")
ax.text(0, 0.5, tDNA, verticalalignment='bottom', horizontalalignment='right', fontfamily='courier new', color = "orange")
ax.text(0, 0.5, daDNA, verticalalignment='bottom', horizontalalignment='left', fontfamily='courier new', color = "green")
ax.text(20 * unit, 0.5, dhDNA, verticalalignment='bottom', horizontalalignment='left', fontfamily='courier new')
ax.add_patch(matplotlib.patches.Rectangle([3*unit, 0.5], 3*unit, 0.17, angle=0.0, fill=False, color="black"))
ax.text(-3 * unit, -0.5, insuDNA, verticalalignment='top', horizontalalignment='right', fontfamily='courier new', color = "red")
ax.text(0, -0.5, insrDNA, verticalalignment='top', horizontalalignment='right', fontfamily='courier new', color = "purple")
ax.text(0, -0.5, insdDNA, verticalalignment='top', horizontalalignment='left', fontfamily='courier new', color="green")
ax.plot([-23 * unit, -20 * unit], [-0.5, 0.5], [-3 * unit, 0], [-0.5, 0.5], [0, -3 * unit], [-0.5, 0.5], [23 * unit, 20 * unit], [-0.5, 0.5], c = "black", ls = "--")
ax.text(-1.5 * unit, 1, "templated insertion", verticalalignment='bottom', horizontalalignment='center', fontsize=10, bbox={"boxstyle": "square", "fill": False, "pad": 0}, fontfamily='arial')
ax.arrow(-1.5 * unit, 1, 0, -0.25, head_width = 0.5, head_length = 0.05, color="black")
ax.text(-1.5 * unit, -1, "random insertion", verticalalignment='top', horizontalalignment='center', fontsize=10, bbox={"boxstyle": "square", "fill": False, "pad": 0}, fontfamily='arial')
ax.arrow(-1.5 * unit, -1, 0, 0.3, head_width = 0.5, head_length = 0.05, color="black")
ax.plot([0,0], [0.3, 0.85], c="grey", ls="--")
ax.set_xlim([-33,33])
ax.set_ylim([-1,1])
ax.set_axis_off()
fig.savefig("paper_figures/templated_insertion.pdf")

##################################
# deletion
##################################

DNA = numpy.random.choice(['A', 'C', 'G', 'T'], 60)
DNA[34:36] = "G"
DNA = "".join(DNA)
insuDNA = DNA[10:29]
insdDNA = DNA[32:50]
insrDNA = "".join(numpy.random.choice(['A', 'C', 'G', 'T'], 3))
while insrDNA[0] == DNA[29] or insrDNA[-1] == DNA[31]:
    insrDNA = "".join(numpy.random.choice(['A', 'C', 'G', 'T'], 3))
uhDNA, uaDNA, tDNA, daDNA, dhDNA = DNA[:10], DNA[10:29], DNA[29:32], DNA[32:50], DNA[50:]

unit = 0.9
fig, ax = matplotlib.pyplot.subplots(figsize=[8,2])
ax.text(-20 * unit, 0.5, uhDNA, verticalalignment='bottom', horizontalalignment='right', fontfamily='courier new')
ax.text(-1 * unit, 0.5, uaDNA, verticalalignment='bottom', horizontalalignment='right', fontfamily='courier new', color = "red")
ax.text(-1 * unit, 0.5, tDNA, verticalalignment='bottom', horizontalalignment='left', fontfamily='courier new', color = "orange")
ax.text(2 * unit, 0.5, daDNA, verticalalignment='bottom', horizontalalignment='left', fontfamily='courier new', color = "green")
ax.text(20 * unit, 0.5, dhDNA, verticalalignment='bottom', horizontalalignment='left', fontfamily='courier new')
ax.add_patch(matplotlib.patches.Rectangle([3*unit, 0.5], 3*unit, 0.17, angle=0.0, fill=False, color="black"))
ax.text(-1 * unit, -0.5, insuDNA, verticalalignment='top', horizontalalignment='right', fontfamily='courier new', color = "red")
ax.text(-1 * unit, -0.5, insrDNA, verticalalignment='top', horizontalalignment='left', fontfamily='courier new', color = "purple")
ax.text(2 * unit, -0.5, insdDNA, verticalalignment='top', horizontalalignment='left', fontfamily='courier new', color="green")
ax.plot([-20 * unit, -20 * unit], [-0.5, 0.5], [-1 * unit, -1 * unit], [-0.5, 0.5], [2 * unit, 2 * unit], [-0.5, 0.5], [20 * unit, 20 * unit], [-0.5, 0.5], c = "black", ls = "--")
ax.text(0.5 * unit, 1, "deletion", verticalalignment='bottom', horizontalalignment='center', fontsize=10, bbox={"boxstyle": "square", "fill": False, "pad": 0}, fontfamily='arial')
ax.arrow(0.5 * unit, 1, 0, -0.25, head_width = 0.5, head_length = 0.05, color="black")
ax.text(0.5 * unit, -1, "random insertion", verticalalignment='top', horizontalalignment='center', fontsize=10, bbox={"boxstyle": "square", "fill": False, "pad": 0}, fontfamily='arial')
ax.arrow(0.5 * unit, -1, 0, 0.3, head_width = 0.5, head_length = 0.05, color="black")
ax.plot([0,0], [0.3, 0.85], c="grey", ls="--")
ax.set_xlim([-33,33])
ax.set_ylim([-1,1])
ax.set_axis_off()
fig.savefig("paper_figures/deletion.pdf")

##################################
# microhomology
##################################

DNA = numpy.random.choice(['A', 'C', 'G', 'T'], 60)
while DNA[29] == "G" or DNA[25] == DNA[31]:
    DNA = numpy.random.choice(['A', 'C', 'G', 'T'], 60)
DNA[34:36] = "G"
DNA = "".join(DNA)
micro = "".join(numpy.random.choice(['A', 'C', 'G', 'T'], 2)) + "G"
umDNA, imDNA, dmDNA = DNA[:26], DNA[29:32], DNA[35:]

unit = 0.9
fig, ax = matplotlib.pyplot.subplots(figsize=[8,2])
ax.text(-4 * unit, 0.5, umDNA, verticalalignment='bottom', horizontalalignment='right', fontfamily='courier new')
ax.text(-4 * unit, 0.5, micro, verticalalignment='bottom', horizontalalignment='left', fontfamily='courier new', color='red')
ax.text(-1 * unit, 0.5, imDNA, verticalalignment='bottom', horizontalalignment='left', fontfamily='courier new')
ax.text(2 * unit, 0.5, micro, verticalalignment='bottom', horizontalalignment='left', fontfamily='courier new', color='red')
ax.text(5 * unit, 0.5, dmDNA, verticalalignment='bottom', horizontalalignment='left', fontfamily='courier new')
ax.add_patch(matplotlib.patches.Rectangle([3*unit, 0.5], 3*unit, 0.17, angle=0.0, fill=False, color="black"))
ax.text(-1 * unit, -0.5, umDNA, verticalalignment='top', horizontalalignment='right', fontfamily='courier new')
ax.text(-1 * unit, -0.5, micro, verticalalignment='top', horizontalalignment='left', fontfamily='courier new', color='red')
ax.text(2 * unit, -0.5, dmDNA, verticalalignment='top', horizontalalignment='left', fontfamily='courier new')
ax.plot([-27 * unit, -30 * unit], [-0.5, 0.5], [27 * unit, 30 * unit], [-0.5, 0.5], ls="--", c="black")
intralines = ax.plot([-unit, -4 * unit], [-0.5, 0.5], [-unit, 2 * unit], [-0.5, 0.5], ls="--", c="black")
ax.set_xlim([-33,33])
ax.set_ylim([-1,1])
ax.set_axis_off()

reso = 1
def update_micro(frame):
    if frame == reso:
        intralines[0].set_xdata([0, -3 * unit])
        intralines[1].set_xdata([0, 3 * unit])
    if frame == 2 * reso:
        intralines[0].set_xdata([unit, -2 * unit])
        intralines[1].set_xdata([unit, 4 * unit])
    if frame == 3 * reso:
        intralines[0].set_xdata([2 * unit, -unit])
        intralines[1].set_xdata([2 * unit, 5 * unit])
    return intralines

ani = matplotlib.animation.FuncAnimation(fig=fig, func=update_micro, frames=4*reso, interval=1200)
ani.save(filename="paper_figures/microhomology.gif", writer="pillow")