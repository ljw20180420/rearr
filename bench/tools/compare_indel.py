import sys, numpy

def NDsegment_distance(indel1, indel2, range1, range2):
    range = [range1[0] - range2[1], range1[1] - range2[0]]
    target = numpy.median(indel2 - indel1)
    actual = range[1] if target > range[1] else range[0] if target < range[0] else target
    return indel1 + actual - indel2

def get_range(indel, ref1, ref2):
    range = [0, 0]
    while indel[0] + range[0] - 1 >= 0 and indel[1] + range[0] - 1 >= 0 and ref1[indel[0] + range[0] - 1] == ref2[indel[1] + range[0] - 1]:
        range[0] = range[0] - 1
    while indel[0] + range[1] < len(ref1) and indel[1] + range[1] < len(ref2) and ref1[indel[0] + range[1]] == ref2[indel[1] + range[1]]:
        range[1] = range[1] + 1
    return range

def indel_dis(indel1, indel2, ref1, ref2):
    range1, range2 = get_range(indel1, ref1, ref2), get_range(indel2, ref1, ref2)
    return NDsegment_distance(indel1, indel2, range1, range2)

_, ref1, ref2, randomseq, mode, reflen, probability, readnum, program = sys.argv

with open(randomseq, "r") as rd:
    usertime, systime, realtime, memory = ["*"] * 4
    for line in sys.stdin:
        _, _, _, _, _, usertime, systime, realtime, memory, name, c2 = line.rstrip().split("\t", 10)
        for randline in rd:
            query, refl1, refr1, queryl1, queryr1 = randline.rstrip().split("\t")
            if query == name:
                break
            sys.stdout.write(f"{mode}\t{reflen}\t{probability}\t{readnum}\t{program}\t{usertime}\t{systime}\t{realtime}\t{memory}\t{query}\t*\t*\t*\t*\n")
        mindis = numpy.inf
        indel1 = numpy.array([int(refl1), int(refr1), int(queryl1), int(queryr1)])
        c2 = c2.split("\t")
        for i in range(0, len(c2), 4):
            indel2 = numpy.array([int(c2[i]), int(c2[i + 1]), int(c2[i + 2]), int(c2[i + 3])])
            disv = indel_dis(indel1, indel2, ref1, ref2)
            dis = numpy.sum(numpy.abs(disv))
            if dis < mindis:
                mindis = dis
                minindel2 = indel2
        sys.stdout.write(f"{mode}\t{reflen}\t{probability}\t{readnum}\t{program}\t{usertime}\t{systime}\t{realtime}\t{memory}\t{query}\t{minindel2[0]}\t{minindel2[1]}\t{minindel2[2]}\t{minindel2[3]}\n")
    for randline in rd:
        query, refl1, refr1, queryl1, queryr1 = randline.rstrip().split("\t")
        sys.stdout.write(f"{mode}\t{reflen}\t{probability}\t{readnum}\t{program}\t{usertime}\t{systime}\t{realtime}\t{memory}\t{query}\t*\t*\t*\t*\n")