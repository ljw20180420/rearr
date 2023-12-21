#!/usr/bin/env python
import numpy, sys

def generate_random_DNA(length):
    # usage: generate_random_DNA of length
    return "".join(numpy.random.choice(["A","C","G","T"], length))

def SNP_DNA(DNA, probability):
    DNA = list(DNA)
    for i in range(len(DNA)):
        if numpy.random.rand() < probability:
            DNA[i] = numpy.random.choice(["A","C","G","T"])
    return "".join(DNA)

def indel_DNA(DNA, probability):
    inslens = numpy.random.negative_binomial(1, 1 - probability, len(DNA) + 1)
    insDNAs = generate_random_DNA(numpy.sum(inslens))
    DNAnew, start = [], 0
    for i in range(len(DNA)):
        end = start + inslens[i]
        DNAnew.append(insDNAs[start:end])
        if numpy.random.rand() > probability:
            DNAnew.append(DNA[i])
        start = end
    DNAnew.append(insDNAs[start:])
    return "".join(DNAnew)

if __name__ == "__main__":
    ref1, ref2, probability, num = sys.argv[1], sys.argv[2], float(sys.argv[3]), int(sys.argv[4])
    for i in range(num):
        lpos = len(ref1) // 2 + numpy.random.randint(-len(ref1) // 20, len(ref1) // 20)
        rpos = len(ref2) // 2 + numpy.random.randint(-len(ref2) // 20, len(ref2) // 20)
        readleft, readright = ref1[len(ref1) // 4: lpos], ref2[rpos : len(ref2) // 4 * 3]
        # rand_ins1 = generate_random_DNA(numpy.random.randint(len(ref1) // 20))
        rand_ins1 = ""
        readleft_mut = SNP_DNA(indel_DNA(readleft, probability), probability)
        rand_ins2 = generate_random_DNA(numpy.random.randint((len(ref1) + len(ref2)) // 40))
        readright_mut = SNP_DNA(indel_DNA(readright, probability), probability)
        # rand_ins3 = generate_random_DNA(numpy.random.randint(len(ref2) // 20))
        rand_ins3 = ""
        sys.stdout.write(f"{rand_ins1}{readleft_mut}{rand_ins2}{readright_mut}{rand_ins3}\tseq{i + 1}\t{lpos}\t{rpos}\t{len(rand_ins1)+len(readleft_mut)}\t{len(rand_ins1)+len(readleft_mut)+len(rand_ins2)}\n")


