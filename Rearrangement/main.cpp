#include "headers/loader.h"
#include "headers/parser.h"
#include "headers/align.h"
#include <limits>

int main(int argc, char **argv)
{
    print_help(argc, argv);
    
    Command_content cc=command(argc, argv);
    
    std::string x;
    std::vector<int> S, left_exp, right_exp;
    load_ref(cc.ref_file, x, S, left_exp, right_exp);
    if(x.empty())
    {
        std::cerr << "empty reference loaded\n";
        return EXIT_FAILURE;
    }

    int max_len = 0;
    std::string o;
    size_t num;

    std::map<char,int> nt2int{{'N',0},{'A',1},{'C',2},{'G',3},{'T',4},{'n',0},{'a',1},{'c',2},{'g',3},{'t',4}};
    TD_array<int> gamma(5, 5, 2, cc.s0);
    gamma(1, 1, 0) = gamma(2, 2, 0) = gamma(3, 3, 0) = gamma(4, 4, 0) = cc.s1;
    gamma(1, 1, 1) = gamma(2, 2, 1) = gamma(3, 3, 1) = gamma(4, 4, 1) = cc.s2;
    std::vector<int> ve(S.back()+1, cc.v), ue(S.back()+1, cc.u), tvf(S.size()-1, cc.v), tuf(S.size()-1, cc.u);
    for(int j=0; j<S.size(); ++j)
    {
        if(cc.alg_type=="local" || cc.alg_type=="contain")
        {
            ve[S[j]] = cc.qv;
            ue[S[j]] = cc.qu;
        }
        else if(cc.alg_type == "local_imbed" && j != 0 && j != S.size() - 1)
        {
            ve[S[j]] = cc.qv;
            ue[S[j]] = cc.qu;
        }
        if(j!=S.size()-1 && (cc.alg_type=="local" || cc.alg_type=="local_imbed" || cc.alg_type=="imbed"))
        {
            tvf[j] = cc.rv;
            tuf[j] = cc.ru;
        }
    }
    TD_array<Back> EFG;
    for(size_t index=1; std::cin >> o >> num; ++index)
    {
        bool new_max = false;
        if (max_len)
        {
            while (max_len < o.size())
            {
                max_len *= 2;
                new_max = true;
            }
        }
        else
        {
            max_len = o.size();
            new_max = true;
        }
        if (new_max)
        {
            EFG.resize(S.back() + 1, max_len + 1, 3);
        }
    
        wapper_column_wise(x, S, nt2int, gamma, cc.v, cc.u, ve, ue, tvf, tuf, EFG, cc.ALIGN_MAX, o, index, num, left_exp, right_exp);
    }

    return 0;
}

































