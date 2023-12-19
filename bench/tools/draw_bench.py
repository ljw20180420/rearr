#!/usr/bin/env python

import sys, pandas, matplotlib.pyplot, seaborn, numpy, os
import plotnine
from plotnine import ggplot, aes, facet_grid, facet_wrap, geom_point, geom_smooth, geom_violin, ggsave

# readnum = int(sys.argv[1])
# difffiles = sys.argv[2:]


readnum = 1000
difffiles = ["bench/single/rearr.diff", "bench/single/CRISPResso.diff", "bench/single/CrisprVariants.diff", "bench/single/amplican.diff", "bench/single/ampliconDIVider.diff", "bench/single/CRISPR-GRANT.diff", "bench/single/ZhangFeng.diff", "bench/single/SelfTarget.diff"]

workdir = os.path.dirname(difffiles[0])
pandict = {}
pandict["index"] = numpy.repeat([i+1 for i in range(readnum)] * len(difffiles), 4).tolist()
programs = [os.path.basename(os.path.splitext(difffile)[0]) for difffile in difffiles]
pandict["program"] = numpy.repeat(programs, 4 * readnum).tolist()
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

ggfig = ggplot(diffdata, aes("index", "diff")) + geom_point(alpha=0.2, color="blue") + geom_smooth() + facet_grid("end ~ program", margins=True, shrink=False)
ggsave(ggfig, os.path.join(workdir, 'diff_scatter.png'), width = 16, height = 12)

ggfig = ggplot(diffdata, aes("program", "diff")) + geom_violin() + facet_grid(". ~ end", margins=True, shrink=False)
ggsave(ggfig, os.path.join(workdir, 'diff_violinplot.png'), width = 24, height = 6)