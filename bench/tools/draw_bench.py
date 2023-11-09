#!/usr/bin/env python

import sys, pandas, matplotlib.pyplot, seaborn, numpy, os
readnum = int(sys.argv[1])
difffiles = sys.argv[2:]

# readnum = 50
# difffiles = ["rearr.diff", "CRISPResso.diff", "CrisprVariants.diff", "amplican.diff", "ampliconDIVider.diff"]

pandict = {}
pandict["index"] = [i+1 for i in range(readnum)] * len(difffiles)
pandict["program"] = numpy.repeat([os.path.basename(difffile) for difffile in difffiles], readnum).tolist()
pandict["diff"] = []
for i in range(len(difffiles)):
    diffs = [None] * readnum
    with open(difffiles[i], "r") as fd:
        for line in fd:
            query, ld, rd, md = line.split("\t")
            diffs[int(query[3:]) - 1] = abs(float(ld)) + abs(float(rd)) + abs(float(md))
    pandict["diff"].extend(diffs)

diffdata = pandas.DataFrame(pandict)

f, ax = matplotlib.pyplot.subplots()
seaborn.scatterplot(data = diffdata, x = "index", y = "diff", hue = "program", style = "program", alpha = 0.5, ax = ax)
f.tight_layout()
f.savefig(f'bench/diff_scatter.pdf')

f, ax = matplotlib.pyplot.subplots()
seaborn.violinplot(data = diffdata, x = "program", y = "diff", ax = ax)
ax.tick_params(axis='x', rotation=10)
f.tight_layout()
f.savefig(f'bench/diff_violinplot.pdf')