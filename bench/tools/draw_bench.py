#!/usr/bin/env python

import sys, pandas, matplotlib.pyplot, seaborn, numpy, os
readnum = int(sys.argv[1])
difffiles = sys.argv[2:]
workdir = os.path.dirname(difffiles[0])

pandict = {}
pandict["index"] = numpy.repeat([i+1 for i in range(readnum)] * len(difffiles), 4).tolist()
pandict["program"] = numpy.repeat([os.path.basename(os.path.splitext(difffile)[0]) for difffile in difffiles], 4 * readnum).tolist()
pandict["end"] = ["us2", "ds1", "uw2", "dw1"] * readnum * len(difffiles)
pandict["diff"] = []
for i in range(len(difffiles)):
    diffs = [None] * 4 * readnum
    with open(difffiles[i], "r") as fd:
        for line in fd:
            query, us2, ds1, uw2, dw1 = line.split("\t")
            pos = 4 * (int(query[3:]) - 1)
            diffs[pos], diffs[pos + 1], diffs[pos + 2], diffs[pos + 3] = abs(float(us2)), abs(float(ds1)), abs(float(uw2)), abs(float(dw1))
    pandict["diff"].extend(diffs)

diffdata = pandas.DataFrame(pandict)

f, ax = matplotlib.pyplot.subplots()
seaborn.scatterplot(data = diffdata, x = "index", y = "diff", hue = "program", style = "program", alpha = 0.5, ax = ax)
f.tight_layout()
f.savefig(os.path.join(workdir, 'diff_scatter.png'))

f, ax = matplotlib.pyplot.subplots()
seaborn.violinplot(data = diffdata, x = "program", y = "diff", ax = ax)
ax.tick_params(axis='x', rotation=10)
f.tight_layout()
f.savefig(os.path.join(workdir, 'diff_violinplot.png'))