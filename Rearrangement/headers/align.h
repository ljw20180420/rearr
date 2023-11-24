#ifndef ALIGN_H
#define ALIGN_H

#include "TD_array.h"
#include "structures.h"
#include <tuple>

void ori_cmp(Back& back, int new_val, Back* address)
{
    if(back.max_val<=new_val)
    {
        if(back.max_val<new_val)
        {
            back.max_val=new_val;
            back.biters.clear();
        }
        back.biters.push_back(address);
    }
}

std::tuple<std::vector<int>,std::vector<std::string>,std::vector<std::string>> column_wise(std::vector<int>::iterator ve_b, std::vector<int>::iterator ue_b, int row, int v, int u, std::vector<int>::iterator tvf_b, std::vector<int>::iterator tuf_b, TD_array<int>& gamma, std::string::iterator x_b, const std::string &o, std::vector<int>::iterator S_b, std::vector<int>::iterator S_e, bool head, bool tail, bool ceil, bool floor, int ALIGN_MAX, std::map<char,int>& nt2int, TD_array<Back>& EFG, std::vector<int>& left_exp, std::vector<int>& right_exp)
{
    EFG(0,0,2).max_val=0;
    for(int j=0; j<S_e-S_b; ++j)
        for(int s=*(S_b+j)+1; s<=*(S_b+j+1); ++s)
        {
            if(ceil)
                EFG(s,0,2).max_val=EFG(s-1,0,2).max_val+((s==*(S_b+j)+1)?(*(tvf_b+j)):(*(tuf_b+j)));
            else
                EFG(s,0,2).max_val = EFG(s-1,0,2).max_val + ((s==*(S_b+j)+1) ? (v) : (u));
            EFG(s,0,0).max_val=std::numeric_limits<int>::min()/2;
        }
    for(int w = 1; w <= o.size(); ++w)
    {
        EFG(0,w,0).max_val=EFG(0,w,2).max_val=EFG(0,w-1,2).max_val+((w==1 && head)?(*ve_b):(*ue_b));
        for(int j=0; j<S_e-S_b; ++j)
        {
            if(floor)
                EFG(*(S_b+j+1),w,1).max_val=std::numeric_limits<int>::min()/2;
            for(int s=*(S_b+j)+1; s<=*(S_b+j+1); ++s)
            {
                // std::cerr << *(S_b+j) << '\t' << s << '\t' << right_exp[j] << '\t' << left_exp[j] << '\n'; // debug
                int gamma_slice = ((s - *(S_b+j)) > right_exp[j] && (s - *(S_b+j)) <= left_exp[j]) ? 0 : 1;
                if(floor)
                    ori_cmp(EFG(*(S_b+j+1),w,1), EFG(s-1,w,2).max_val+*(tvf_b+j)+(*(S_b+j+1)-s)*(*(tuf_b+j)), &EFG(s-1,w,2));
                if(s!=*(S_b+j+1) || !floor)
                {
                    EFG(s,w,1).max_val=std::numeric_limits<int>::min()/2;
                    ori_cmp(EFG(s,w,1), EFG(s-1,w,2).max_val + v, &EFG(s-1,w,2));
                    if(s!=*(S_b+j)+1)
                        ori_cmp(EFG(s,w,1), EFG(s-1,w,1).max_val + u, &EFG(s-1,w,1));
                    if(ceil)
                        ori_cmp(EFG(s,w,1), EFG(*(S_b+j),w,2).max_val+*(tvf_b+j)+(s-*(S_b+j)-1)*(*(tuf_b+j)), &EFG(*(S_b+j),w,2));
                }
                EFG(s,w,0).max_val=std::numeric_limits<int>::min()/2;
                if(!tail && s==*S_e)
                    ori_cmp(EFG(s,w,0), EFG(s,w-1,2).max_val+*(ue_b+s), &EFG(s,w-1,2));
                else
                {
                    ori_cmp(EFG(s,w,0), EFG(s,w-1,2).max_val+*(ve_b+s), &EFG(s,w-1,2));
                    if(EFG(s,w-1,2).max_val!=EFG(s,w-1,0).max_val || *(ve_b+s)!=*(ue_b+s))
                        ori_cmp(EFG(s,w,0), EFG(s,w-1,0).max_val+*(ue_b+s), &EFG(s,w-1,0));
                }
                EFG(s,w,2).max_val=std::numeric_limits<int>::min()/2;
                ori_cmp(EFG(s,w,2), EFG(s-1,w-1,2).max_val+gamma(nt2int[*(x_b+s-1)],nt2int[o[w-1]],gamma_slice), &EFG(s-1,w-1,2));
                ori_cmp(EFG(s,w,2), EFG(s,w,1).max_val, &EFG(s,w,1));
                ori_cmp(EFG(s,w,2), EFG(s,w,0).max_val, &EFG(s,w,0));
            }
        }
    }
    std::vector<std::vector<std::pair<int,int>>> rivets;
    std::vector<Back*> biters;
    rivets.push_back(std::vector<std::pair<int,int>>());
    biters.push_back(&EFG(*S_e, o.size(), 2));
    int i=0;
    while(i<biters.size())
    {
        rivets[i].push_back(std::make_pair(EFG.grow(biters[i]),EFG.gcol(biters[i])));
        if(rivets[i].back().first==0 || rivets[i].back().second==0)
        {
            rivets[i].push_back(std::make_pair(0,0));
            ++i;
        }
        else
        {
            std::vector<Back*>::iterator iter=biters[i]->biters.begin();
            std::vector<Back*>::iterator end_iter=biters[i]->biters.end();
            biters[i]=*iter;
            ++iter;
            while(biters.size()<ALIGN_MAX && iter!=end_iter)
            {
                rivets.push_back(rivets[i]);
                biters.push_back(*iter);
                ++iter;
            }
        }
    }
    std::vector<int> max_scores;
    std::vector<std::string> references, queries;
    for(int i=0; i<rivets.size(); ++i)
    {
        references.emplace_back();
        queries.emplace_back();
        max_scores.push_back(EFG(*S_e, o.size(), 2).max_val);
        for(std::vector<std::pair<int,int>>::reverse_iterator iter=rivets[i].rbegin(); iter!=rivets[i].rend()-1; ++iter)
        {
            if(iter->first==(iter+1)->first)
            {
                references.back().append(std::string((iter+1)->second-iter->second,'-'));
                queries.back().append(o.begin() + iter->second, o.begin() + (iter+1)->second);
            }
            else if(iter->second==(iter+1)->second)
            {
                references.back().append(x_b+iter->first,x_b+(iter+1)->first);
                queries.back().append(std::string((iter+1)->first-iter->first,'-'));
            }
            else
            {
                references.back().append(x_b+iter->first,x_b+(iter+1)->first);
                queries.back().append(o.begin() + iter->second, o.begin() + (iter+1)->second);
            }
        }
    }
    return std::make_tuple(max_scores, references, queries);
}

void wapper_column_wise(std::string x, std::vector<int> S, std::map<char,int> &nt2int, TD_array<int> &gamma, int v, int u, std::vector<int> &ve, std::vector<int> &ue, std::vector<int> &tvf, std::vector<int> &tuf, TD_array<Back> &EFG, int ALIGN_MAX, const std::string &o, size_t index, size_t num, std::vector<int> left_exp, std::vector<int> right_exp)
{
    std::vector<int> max_scores;
    std::vector<std::string> references, queries;
    std::tie(max_scores, references, queries) = column_wise(ve.begin(), ue.begin(), S.size()-1, v, u, tvf.begin(), tuf.begin(), gamma, x.begin(), o, S.begin(), S.end()-1, true, true, true, true, ALIGN_MAX, nt2int, EFG, left_exp, right_exp);
    for(int j=0; j<max_scores.size(); ++j)
    {
        std::cout << index << '\t' << num << '\t' << max_scores[j] << '\n' << references[j] << '\n' << queries[j] << '\n';
    }
}

#endif
