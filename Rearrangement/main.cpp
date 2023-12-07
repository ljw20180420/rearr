#include "headers/parser.h"
#include "headers/align.h"

int main(int argc, char **argv)
{
    print_help(argc, argv);
    
    Command_content cc=command(argc, argv);
    
    std::deque<std::string> refs;
    std::deque<uint32_t> upper_boundaries, down_boundaries;
    uint32_t max_ref_sz = 0;
    FILE *ref_fd = fdopen(3, "r");
    char buffer_ref[128];
    for (size_t i = 0; fgets(buffer_ref, 128, ref_fd); ++i)
    {
        switch (i % 3)
        {
            case 0:
                upper_boundaries.push_back(std::atoi(buffer_ref));
                break;
            case 1:
                refs.emplace_back();
                while (strcspn(buffer_ref, "\n") == strlen(buffer_ref))
                {
                    refs.back().append(buffer_ref);
                    fgets(buffer_ref, 128, ref_fd);
                }
                buffer_ref[strcspn(buffer_ref, "\n")] = '\0';
                refs.back().append(buffer_ref);
                std::transform(refs.back().begin(), refs.back().end(), refs.back().begin(), toupper);
                max_ref_sz = refs.back().size() > max_ref_sz ? refs.back().size() : max_ref_sz;
                break;
            case 2:
                down_boundaries.push_back(std::atoi(buffer_ref));
                break;
        }
    }
    if(!max_ref_sz)
    {
        std::cerr << "empty reference loaded\n";
        return EXIT_FAILURE;
    }
    fclose(ref_fd);

    std::vector<int32_t> Es(max_ref_sz + 1), Gs(max_ref_sz + 1), Gps(max_ref_sz + 1), gr(max_ref_sz + 1);
    gr[0] = 0;
    gr[1] = cc.rv;
    for (size_t i = 2; i <= max_ref_sz; ++i)
        gr[i] = gr[i - 1] + cc.ru;

    uint32_t Omax = 0;
    std::string O;
    size_t count;
    std::vector<int32_t> As, Bs, Cs;
    std::vector<uint32_t> Ds;
    for(size_t index=1; std::cin >> O >> count; ++index)
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
            As.resize((Omax + 1) * (refs.size() + 1));
            As[0] = 0;
            As[1] = cc.qv;
            for (size_t i = 2; i <= Omax; ++i)
                As[i] = As[i - 1] + cc.qu;
            Bs.resize((Omax + 1) * (refs.size()));
            Cs.resize((Omax + 1) * (refs.size()));
            Ds.resize((Omax + 1) * (refs.size()));
        }
    
        for (size_t i = 0; i < refs.size(); ++i)
        {
            cross_align(As.data() + i * (Omax + 1), As.data() + (i + 1) * (Omax + 1), Bs.data() + i * (Omax + 1), Cs.data() + i * (Omax + 1), Ds.data() + i * (Omax + 1), Es.data(), Gs.data(), Gps.data(), refs[i].data(), refs[i].size(), O.data(), O.size(), cc.u, cc.v, gr.data(), cc.qu, cc.qv, cc.s0, cc.s1, cc.s2, upper_boundaries[i], down_boundaries[i]);
        }

        std::deque<std::string> reflines, querylines, random_parts;
        uint32_t w = O.size();
        for (size_t i = 0; i < refs.size(); ++i)
        {
            size_t j = refs.size() - i;
            std::pair<uint32_t, uint32_t> swp = NodeTrack(w, As.data() + j * (Omax + 1), Bs.data() + (j - 1) * (Omax + 1), Cs.data() + (j - 1) * (Omax + 1), Ds.data() + (j - 1) * (Omax + 1), cc.qv);
            random_parts.push_front(O.substr(swp.second, w - swp.second));
            reflines.emplace_front();
            querylines.emplace_front();
            w = EdgeTrack(reflines.front(), querylines.front(), swp.first, swp.second, As.data() + (j - 1) * (Omax + 1), As.data() + j * (Omax + 1), refs[j - 1].data(), refs[j - 1].size(), O.data(), O.size(), cc.u, cc.v, gr.data(), cc.qv, cc.s0, cc.s1, cc.s2, upper_boundaries[j - 1], down_boundaries[j - 1]);
        }
        random_parts.push_front(O.substr(0, w));

        std::cout << index << '\t' << count << '\t' << As[refs.size() * (Omax + 1) + O.size()] << '\n';
        for (size_t i = 0; i < reflines.size(); ++i)
            std::cout << std::string(random_parts[i].size(), '-') << reflines[i];
        std::cout << std::string(random_parts.back().size(), '-') << '\n';
        for (size_t i = 0; i < querylines.size(); ++i)
            std::cout << random_parts[i] << querylines[i];
        std::cout << random_parts.back() << '\n';
    }

    return 0;
}

































