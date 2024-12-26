/**
 * @file parser.h
 * @author Jingwei Li (ljw2017@sjtu.edu.cn)
 * @brief Handle command-line arguments.
 * @date 2024-12-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef PARSER_H
#define PARSER_H

#include <bits/stdc++.h>

/** @brief Alignment scores are of type int32_t. */
typedef int32_t SCORETYPE;

/**
 * @brief Print help if --help, -help, or -h is passed through command-line.
 * 
 * @param[in] argc argc of main function.
 * @param[in] argv argv of main function.
 */
void print_help(
    int argc,
    char **argv
) {
    bool print_help = false;
    for(int i=1; i<argc; ++i) {
        if(!strcmp(argv[i],"--help") || !strcmp(argv[i],"-help") || !strcmp(argv[i],"-h")) {
            print_help = true;
            break;
        }
    }
    if (print_help) {
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
        << "-ru: Gap-extending penalty for unaligned reference ends. (default: 0)\n"
        << "-rv: Gap-opening penalty for unaligned reference ends. (default: 0)\n"
        << "-qu: Gap-extending penalty for unaligned query parts. (default: 0)\n"
        << "-qv: Gap-opening penalty for unaligned query parts. (default: 0)\n";
        exit(EXIT_SUCCESS);
    }
    return;
}

/**
 * @brief Structure holding command-line arguments.
 * 
 */
struct Command_content {
    /// @brief Mismatching score.
    SCORETYPE s0=-3;
    /// @brief Matching score for non-extension reference part.
    SCORETYPE s1=1;
    /// @brief Matching score for extension reference part.
    SCORETYPE s2=1;
    /// @brief Gap-extending penalty.
    SCORETYPE u=-2;
    /// @brief Gap-opening penalty.
    SCORETYPE v=-5;
    /// @brief Gap-extending penalty for unaligned reference ends.
    SCORETYPE ru=0;
    /// @brief Gap-opening penalty for unaligned reference ends.
    SCORETYPE rv=0;
    /// @brief Gap-extending penalty for unaligned query parts.
    SCORETYPE qu=0;
    /// @brief Gap-opening penalty for unaligned query parts.
    SCORETYPE qv=0;
};

/**
 * @brief Parse and save command-line arguments in returned Command_content struct.
 * 
 * @param argc argc of main function
 * @param argv argv of main function
 * @return Command_content 
 */
Command_content command(
    int argc,
    char **argv
) {
    Command_content cc;

    for(size_t i=1; i<argc-1; ++i) {
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