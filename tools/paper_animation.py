import matplotlib.pyplot, numpy, matplotlib.animation, matplotlib.patches, matplotlib.transforms, matplotlib.patches, matplotlib.path

#####################################
# stagger cut (animation)
#####################################

DNA = numpy.random.choice(['A', 'C', 'G', 'T'], 60)
DNA[34:36] = "G"
DNA = "".join(DNA)
rcDNA = DNA.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c").upper()[::-1]
uDNA, dDNA, urcDNA, drcDNA = DNA[:29], DNA[29:], rcDNA[:30], rcDNA[30:]

unit = 1.04
fig, ax = matplotlib.pyplot.subplots(figsize=[5.2,0.5])
ax.set_position([0,0,1,1])
uDNAtext = ax.text(-unit, 0, uDNA, va='bottom', ha='right', fontfamily='courier new')
dDNAtext = ax.text(-unit, 0, dDNA, va='bottom', ha='left', fontfamily='courier new')
PAMbox = ax.add_patch(matplotlib.patches.Rectangle([3*unit, 0.1], 3*unit, 0.4, angle=0.0, fill=False, color="black"))
urcDNAtext = ax.text(0, 0, urcDNA, va='top', ha='right', fontfamily='courier new')
drcDNAtext = ax.text(0, 0, drcDNA, va='top', ha='left', fontfamily='courier new')
ecut = ax.plot([0, 0], [-1, 1], c = 'grey', ls = '--')
cuttext = ax.text(-unit, 0, "|", va='bottom', ha='center', alpha = 0)
rccuttext = ax.text(0, 0, "|", va='top', ha='center', alpha = 0)
abasetext = ax.text(-unit, 0, dDNA[0], va='bottom', ha='left', fontfamily='courier new', alpha = 0, color = "red")
arcbasetext = ax.text(0, 0, urcDNA[-1], va='top', ha='left', fontfamily='courier new', alpha = 0, color = "red")
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
        PAMbox.set_x(4 * unit)
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
# microhomology (animation)
##################################

DNA = numpy.random.choice(['A', 'C', 'G', 'T'], 60)
while DNA[29] == "G" or DNA[25] == DNA[31]:
    DNA = numpy.random.choice(['A', 'C', 'G', 'T'], 60)
DNA[34:36] = "G"
DNA = "".join(DNA)
micro = "".join(numpy.random.choice(['A', 'C', 'G', 'T'], 2)) + "G"
umDNA, imDNA, dmDNA = DNA[:26], DNA[29:32], DNA[35:]

unit = 1.04
fig, ax = matplotlib.pyplot.subplots(figsize=[5.2,0.5])
ax.set_position([0,0,1,1])
ax.text(-4 * unit, 0.5, umDNA, va='bottom', ha='right', fontfamily='courier new')
ax.text(-4 * unit, 0.5, micro, va='bottom', ha='left', fontfamily='courier new', color='red')
ax.text(-1 * unit, 0.5, imDNA, va='bottom', ha='left', fontfamily='courier new')
ax.text(2 * unit, 0.5, micro, va='bottom', ha='left', fontfamily='courier new', color='red')
ax.text(5 * unit, 0.5, dmDNA, va='bottom', ha='left', fontfamily='courier new')
ax.add_patch(matplotlib.patches.Rectangle([3*unit, 0.55], 3*unit, 0.45, angle=0.0, fill=False, color="black"))
ax.plot([0, 0], [-1, 1], c = 'grey', ls = '--')
ax.text(-1 * unit, -0.5, umDNA, va='top', ha='right', fontfamily='courier new')
ax.text(-1 * unit, -0.5, micro, va='top', ha='left', fontfamily='courier new', color='red')
ax.text(2 * unit, -0.5, dmDNA, va='top', ha='left', fontfamily='courier new')
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