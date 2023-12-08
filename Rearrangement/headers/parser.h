#ifndef PARSER_H
#define PARSER_H

#include <bits/stdc++.h>

typedef int32_t SCORETYPE;

void print_help(int argc, char **argv)
{
    bool print_help = false;
    if (argc <= 1)
        print_help = true;
    else
        for(int i=1; i<argc; ++i)
            if(!strcmp(argv[i],"--help") || !strcmp(argv[i],"-help") || !strcmp(argv[i],"-h"))
            {
                print_help = true;
                break;
            }
    if (print_help)
    {
        std::cout << "###Basic Usage\n"
        << "rearrangement <input_file 3<reference_file\n"

        << "\n### Parameters\n"
        << "-h, -help, --help: Display help.\n"
        << "# Aligning Parameters\n"
        << "-s0: Mismatching score. (default: -3)\n"
        << "-s1: Matching score for non-extension reference part. (default: +1)\n"
        << "-s2: Matching score for extension reference part. (default: +1)\n"
        << "-u: Gap-extending penalty. (default: -2)\n"
        << "-v: Gap-opening penalty. (default: -5)\n"
        << "-ru: Gap-extending penalty for unaligned reference end. (default: 0)\n"
        << "-rv: Gap-opening penalty for unaligned reference end. (default: 0)\n"
        << "-qu: Gap-extending penalty for unaligned query part. (default: 0)\n"
        << "-qv: Gap-opening penalty for unaligned query part. (default: 0)\n";
        exit(EXIT_SUCCESS);
    }
    return;
}

struct Command_content
{
    // Aligning Parameters
    SCORETYPE s0=-3; // mismatch score
    SCORETYPE s1=1; // match score for non-extension part of reference
    SCORETYPE s2=1; // match score for extension part of reference
    SCORETYPE u=-2; // gap extension
    SCORETYPE v=-5; // gap open
    SCORETYPE ru=0; // reference end gap extension
    SCORETYPE rv=0; // reference end gap open
    SCORETYPE qu=0; // query unaligned gap extension
    SCORETYPE qv=0; // query unaligned gap open
};

Command_content command(int argc, char **argv)
{
    Command_content cc;

    for(size_t i=1; i<argc-1; ++i)
    {
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
    }
    
    return cc;
}

#endif