#include "headers/parser.h"
#include "headers/align.h"

int main(int argc, char **argv)
{
    print_help(argc, argv);
    
    Command_content cc=command(argc, argv);
    
    std::deque<std::deque<std::string>> refss;
    std::deque<std::deque<size_t>> upper_boundariess, down_boundariess;
    size_t max_ref_sz = 0;
    FILE *ref_fd = fdopen(3, "r");
    char buffer_ref[128];
    while (fgets(buffer_ref, 128, ref_fd))
    {
        std::stringstream line;
        while (strcspn(buffer_ref, "\n") == strlen(buffer_ref))
        {
                    line << buffer_ref;
                    fgets(buffer_ref, 128, ref_fd);
        }
        buffer_ref[strcspn(buffer_ref, "\n")] = '\0';
        line << buffer_ref;
        refss.emplace_back();
        upper_boundariess.emplace_back();
        down_boundariess.emplace_back();
        size_t ub, db;
        std::string ref;
        while (line >> ub >> ref >> db)
        {
            upper_boundariess.back().push_back(ub);
            refss.back().push_back(ref);
            down_boundariess.back().push_back(db);
            std::transform(refss.back().back().begin(), refss.back().back().end(), refss.back().back().begin(), toupper);
            max_ref_sz = refss.back().back().size() > max_ref_sz ? refss.back().back().size() : max_ref_sz;
        }
    }
    if(!max_ref_sz)
    {
        std::cerr << "empty reference loaded\n";
        return EXIT_FAILURE;
    }
    fclose(ref_fd);
    size_t max_seg_sz = max_ref_sz / simd_sz + 1;
    std::vector<SCORETYPE> gr = {0, cc.rv};
    for (size_t s = 2; s < max_seg_sz * simd_sz; ++s)
        gr.push_back(gr.back() + cc.ru);
    vsimd *Es = new vsimd[max_seg_sz], *Gs = new vsimd[max_seg_sz], *Gps = new vsimd[max_seg_sz];
    std::vector<std::vector<vsimd *>> grpss, grmss;
    std::vector<std::vector<vsimd **>> gammass;
    for (size_t r = 0; r < refss.size(); ++r)
    {
        grpss.emplace_back();
        grmss.emplace_back();
        gammass.emplace_back();
        for (size_t i = 0; i < refss[r].size(); ++i)
        {
            size_t seg_sz = refss[r][i].size() / simd_sz + 1;
            grpss.back().push_back(new vsimd[seg_sz]);
            grmss.back().push_back(new vsimd[seg_sz]);
            for (size_t j = 0; j < seg_sz; ++j)
                for (size_t k = 0; k < simd_sz; ++k)
                {
                    size_t s = seg_sz * k + j;
                    grpss.back()[i][j][k] = s > 0 ? (cc.rv + (s - 1) * cc.ru) : 0;
                    grmss.back()[i][j][k] = s == refss[r][i].size() ? 0 : (cc.rv + (refss[r][i].size() - s - 1) * cc.ru);
                }

            gammass.back().push_back(new vsimd *[seg_sz]);
            for (size_t j = 0; j < seg_sz; ++j)
            {
                gammass.back()[i][j] = new vsimd[5];
                for (size_t k = 0; k < simd_sz; ++k)
                {
                    size_t s = k * seg_sz + j; 
                    const char *bases = "ACGTN";
                    for (size_t l = 0; l < 5; ++l)
                    {
                        if (s >= refss[r][i].size())
                            gammass.back()[i][j][l][k] = -inf;
                        else
                        {
                            if (refss[r][i][s] != 'N' && refss[r][i][s] == bases[l])
                            {
                                if (s >= upper_boundariess[r][i] &&  s < down_boundariess[r][i])
                                    gammass.back()[i][j][l][k] = cc.s1;
                                else
                                    gammass.back()[i][j][l][k] = cc.s2;
                            }
                            else
                                gammass.back()[i][j][l][k] = cc.s0;
                        }
                    }
                }
            }
        }
    }

    size_t Omax = 0;
    std::string O;
    size_t count, ref_id;
    std::vector<SCORETYPE> As, Bs, Cs;
    std::vector<size_t> Ds;
    for(size_t index=1; std::cin >> O >> count >> ref_id; ++index)
    {
        std::transform(O.begin(), O.end(), O.begin(), toupper);
        bool new_max = false;
        if (index == 1)
        {
            Omax = O.size();
            new_max = true;
        }
        while (Omax < O.size())
        {
            Omax *= 2;
            new_max = true;
        }
        if (new_max)
        {
            As.resize((Omax + 1) * (refss[ref_id].size() + 1));
            As[0] = 0;
            As[1] = cc.qv;
            for (size_t i = 2; i <= Omax; ++i)
                As[i] = As[i - 1] + cc.qu;
            Bs.resize((Omax + 1) * (refss[ref_id].size()));
            Cs.resize((Omax + 1) * (refss[ref_id].size()));
            Ds.resize((Omax + 1) * (refss[ref_id].size()));
        }
    
        for (size_t i = 0; i < refss[ref_id].size(); ++i)
            cross_align(As.data() + i * (Omax + 1), As.data() + (i + 1) * (Omax + 1), Bs.data() + i * (Omax + 1), Cs.data() + i * (Omax + 1), Ds.data() + i * (Omax + 1), Es, Gs, Gps, refss[ref_id][i].data(), refss[ref_id][i].size(), O.data(), O.size(), cc.u, cc.v, grpss[ref_id][i], grmss[ref_id][i], gammass[ref_id][i], cc.qu, cc.qv, cc.s0, cc.s1, cc.s2, upper_boundariess[ref_id][i], down_boundariess[ref_id][i]);

        std::deque<std::string> reflines, querylines, random_parts;
        size_t w = O.size();
        for (size_t i = 0; i < refss[ref_id].size(); ++i)
        {
            size_t j = refss[ref_id].size() - i;
            std::pair<size_t, size_t> swp = NodeTrack(w, As.data() + j * (Omax + 1), Bs.data() + (j - 1) * (Omax + 1), Cs.data() + (j - 1) * (Omax + 1), Ds.data() + (j - 1) * (Omax + 1), cc.qv);
            random_parts.push_front(O.substr(swp.second, w - swp.second));
            reflines.emplace_front();
            querylines.emplace_front();
            w = EdgeTrack(reflines.front(), querylines.front(), swp.first, swp.second, As.data() + (j - 1) * (Omax + 1), As.data() + j * (Omax + 1), refss[ref_id][j - 1].data(), refss[ref_id][j - 1].size(), O.data(), O.size(), cc.u, cc.v, gr.data(), cc.qv, cc.s0, cc.s1, cc.s2, upper_boundariess[ref_id][j - 1], down_boundariess[ref_id][j - 1]);
        }
        random_parts.push_front(O.substr(0, w));

        std::cout << index << '\t' << count << '\t' << As[refss[ref_id].size() * (Omax + 1) + O.size()] << '\t' << ref_id << '\n';
        for (size_t i = 0; i < reflines.size(); ++i)
            std::cout << std::string(random_parts[i].size(), '-') << reflines[i];
        std::cout << std::string(random_parts.back().size(), '-') << '\n';
        for (size_t i = 0; i < querylines.size(); ++i)
            std::cout << random_parts[i] << querylines[i];
        std::cout << random_parts.back() << '\n';
    }

    return 0;
}

































