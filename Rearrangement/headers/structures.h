#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <string>
#include <vector>

struct Back
{
    int max_val;
    std::vector<Back*> biters;
};

#endif