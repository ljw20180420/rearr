#ifndef PARSER_H
#define PARSER_H

#include <bits/stdc++.h>

void print_help(int argc, char **argv)
{
    for(int i=1; i<argc; ++i)
        if(!strcmp(argv[i],"--help") || !strcmp(argv[i],"-h"))
        {
            std::cout << "###Basic Usage\n"
            << "rearrangement -file input_file -ref_file reference_file -ALIGN_MAX 1 -THR_NUM 6\n"

            << "\n### Parameters\n"
            << "\n# Help Related\n"
            << "-help: Display help.\n"

            << "\n# Input Files\n"
            << "-ref_file: Reference file.\n"

            << "\n# Aligning Parameters\n"
            << "-s0: Mismatching score. (default: -3)\n"
            << "-s1: Matching score for non-extension reference part. (default: +1)\n"
            << "-s2: Matching score for extension reference part. (default: +1)\n"
            << "-u: Gap-extending penalty. (default: -2)\n"
            << "-v: Gap-opening penalty. (default: -5)\n"
            << "-ru: Gap-extending penalty for unaligned reference end. (default: 0)\n"
            << "-rv: Gap-opening penalty for unaligned reference end. (default: 0)\n"
            << "-qu: Gap-extending penalty for unaligned query part. (default: 0)\n"
            << "-qv: Gap-opening penalty for unaligned query part. (default: 0)\n"
            << "-alg_type: Method to scoring unaligned parts between aligned segments in query sequence and reference (-alg_type local|imbed|contain|local_imbed). (default: local)\n"

            << "\n# Output Options\n"
            << "-ALIGN_MAX: The maximally reserved number of best alignments for each read (Each read may have several best alignments to the reference). (default: 5)\n";
            exit(0);
        }
    return;
}

struct Command_content
{
    // Input Files
    std::string ref_file="";
    // Aligning Parameters
    int s0=-3; // mismatch score
    int s1=1; // match score for non-extension part of reference
    int s2=1; // match score for extension part of reference
    int u=-2; // gap extension
    int v=-5; // gap open
    int ru=0; // reference end gap extension
    int rv=0; // reference end gap open
    int qu=0; // query unaligned gap extension
    int qv=0; // query unaligned gap open
    std::string alg_type="local";
    // Output Options
    int ALIGN_MAX=5;
};

Command_content command(int argc, char **argv)
{
    Command_content cc;

    for(int i=1; i<argc-1; ++i)
    {
        // Input Files
        if(!strcmp(argv[i],"-ref_file"))
            cc.ref_file=argv[i+1];
        // Aligning Parameters
        if(!strcmp(argv[i],"-s0"))
            cc.s0=atoi(argv[i+1]);
        if(!strcmp(argv[i],"-s1"))
            cc.s1=atoi(argv[i+1]);
        if(!strcmp(argv[i],"-s2"))
            cc.s2=atoi(argv[i+1]);
        if(!strcmp(argv[i],"-u"))
            cc.u=atoi(argv[i+1]);
        if(!strcmp(argv[i],"-v"))
            cc.v=atoi(argv[i+1]);
        if(!strcmp(argv[i],"-ru"))
            cc.ru=atoi(argv[i+1]);
        if(!strcmp(argv[i],"-rv"))
            cc.rv=atoi(argv[i+1]);
        if(!strcmp(argv[i],"-qu"))
            cc.qu=atoi(argv[i+1]);
        if(!strcmp(argv[i],"-qv"))
            cc.qv=atoi(argv[i+1]);
        if(!strcmp(argv[i],"-alg_type"))
            cc.alg_type=argv[i+1];
        // Output Options
        if(!strcmp(argv[i],"-ALIGN_MAX"))
            cc.ALIGN_MAX=atoi(argv[i+1]);
    }
    
    if (cc.ref_file.empty())
    {
        std::cerr << "-ref_file is not specified\n";
        exit(-1);
    }
    
    return cc;
}

#endif