#ifndef LOADER_H
#define LOADER_H

#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <algorithm>

void load_ref(std::string& ref_file, std::string& x, std::vector<int>& S, std::vector<int>& left_exp, std::vector<int>& right_exp)
{
    S.push_back(0);
    std::ifstream fin(ref_file);
    int le_tmp, re_tmp;
    std::string str_tmp;
    while(fin >> re_tmp >> str_tmp >> le_tmp)
    {
        std::transform(str_tmp.begin(), str_tmp.end(), str_tmp.begin(), toupper);
        str_tmp.front()=tolower(str_tmp.front());
        str_tmp.back()=tolower(str_tmp.back());
        S.push_back(S.back()+str_tmp.size());
        x.append(str_tmp);
        left_exp.push_back(le_tmp);
        right_exp.push_back(re_tmp);
    }
    fin.close();
}

#endif