#!/usr/bin/env python
import numpy, sys, os, io

def generate_random_DNA(length):
    # usage: generate_random_DNA of length
    return "".join(numpy.random.choice(["A","C","G","T"], length))

def SNP_DNA(DNA, probability):
    DNAnew = ""
    for i in range(len(DNA)):
        DNAnew += numpy.random.choice(["A","C","G","T"]) if numpy.random.rand() < probability else DNA[i]
    return DNAnew

def indel_DNA(DNA, probability):
    DNAnew = ""
    for i in range(len(DNA)):
        DNAnew += generate_random_DNA(numpy.random.negative_binomial(1, 1 - probability))
        if numpy.random.rand() > probability:
            DNAnew += DNA[i]
    DNAnew += generate_random_DNA(numpy.random.negative_binomial(1, 1 - probability))
    return DNAnew

if __name__ == "__main__":
    ref, probability, num = sys.argv[1], float(sys.argv[2]), int(sys.argv[3])
    for i in range(num):
        lpos = len(ref) // 2 + numpy.random.randint(-len(ref) // 20, len(ref) // 20)
        rpos = len(ref) // 2 + numpy.random.randint(-len(ref) // 20, len(ref) // 20)
        readleft, readright = ref[len(ref) // 4: lpos], ref[rpos : len(ref) // 4 *3]
        rand_ins1 = generate_random_DNA(numpy.random.randint(len(ref) // 20))
        readleft_mut = SNP_DNA(indel_DNA(readleft, probability), probability)
        rand_ins2 = generate_random_DNA(numpy.random.randint(len(ref) // 20))
        readright_mut = SNP_DNA(indel_DNA(readright, probability), probability)
        rand_ins3 = generate_random_DNA(numpy.random.randint(len(ref) // 20))
        sys.stdout.write(f"{rand_ins1}{readleft_mut}{rand_ins2}{readright_mut}{rand_ins3}\tseq{i + 1}\t{lpos}\t{rpos}\t{len(rand_ins1)+len(readleft_mut)}\t{len(rand_ins1)+len(readleft_mut)+len(rand_ins2)}\n")