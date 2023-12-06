#ifndef ALIGN_H
#define ALIGN_H

#include <bits/stdc++.h>

const static int32_t inf = std::numeric_limits<int32_t>::max() / 2;

void cross_align(int32_t *tAs, int32_t *hAs, int32_t *hBs, int32_t *hCs, uint32_t *hDs, int32_t *Es, int32_t *Gs, int32_t *Gps, const char *ref, const uint32_t ref_sz, const char *O, const uint32_t O_sz, const int32_t u, const int32_t v, const int32_t *gr, const int32_t qu, const int32_t qv, const int32_t s0, const int32_t s1, const int32_t s2, const uint32_t upper_boundary, const uint32_t down_boundary)
{
    for (uint32_t s = 0; s <= ref_sz; ++s)
    {
        Es[s] = -inf;
        Gs[s] = tAs[0] + gr[s];
    }
    hCs[0] = Gs[0] + gr[ref_sz];
    hDs[0] = 0;
    hBs[0] = -inf;

    for (uint32_t w = 1; w <= O_sz; ++w)
    {
        std::swap(Gs, Gps);
        for (uint32_t s = 0; s <= ref_sz; ++s)
        {
            Es[s] = std::max(Es[s] + u, Gps[s] + v);
            Gs[s] = std::max(tAs[w] + gr[s], Es[s]);
        }
        hCs[w] = Gs[0] + gr[ref_sz];
        hDs[w] = 0;
        int32_t F = -inf;
        for (uint32_t s = 0; s < ref_sz; ++s)
        {
            F = std::max(F + u, Gs[s] + v);
            Gs[s + 1] = std::max(Gs[s + 1], F);
            int32_t gamma;
            if (ref[s] != O[w-1])
                gamma = s0;
            else if (s < upper_boundary || s >= down_boundary)
                gamma = s2;
            else
                gamma = s1;
            Gs[s + 1] = std::max(Gs[s + 1], Gps[s] + gamma);
            if (Gs[s + 1] + gr[ref_sz - s - 1] > hCs[w])
            {
                hCs[w] = Gs[s + 1] + gr[ref_sz - s - 1];
                hDs[w] = s + 1;
            }
        }
        hBs[w] = std::max(hBs[w - 1] + qu, hCs[w - 1] + qv);
        hAs[w] = std::max(hBs[w], hCs[w]);
    }
}

uint32_t EdgeTrack(std::string &refline, std::string &queryline, const uint32_t s, const uint32_t w, const int32_t *tAs, const int32_t *hAs, const char *ref, const uint32_t ref_sz, const char *O, const uint32_t O_sz, const int32_t u, const int32_t v, const int32_t *gr, const int32_t qv, const int32_t s0, const int32_t s1, const int32_t s2, const uint32_t upper_boundary, const uint32_t down_boundary)
{
    std::vector<int32_t> Es, Fs, Gs;
    std::vector<uint32_t> shifts, wends;
    shifts.push_back(0);
    wends.push_back(w);
    Gs.push_back(hAs[w] - gr[ref_sz - s]);
    Fs.push_back(inf);
    Es.push_back(inf);
    for (uint32_t i = 1; ; ++i)
    {
        shifts.push_back(Gs.size());
        wends.push_back(wends[i - 1]);
        for (uint32_t j = 0; j <= wends[i];)
        {
            int32_t E, F, G, gamma;
            if (j == 0)
                E = inf;
            else
                E = std::min(Es[shifts[i] + j - 1], Gs[shifts[i] + j - 1]) - u;
            uint32_t pj = wends[i - 1] - wends[i] + j;
            if (pj >= shifts[i] - shifts[i - 1])
                F = inf;
            else
                F = std::min(Fs[shifts[i - 1] + pj], Gs[shifts[i - 1] + pj]) - u;
            G = std::min(E, F) - v + u;
            uint32_t ss = s - i, ww = wends[i] - j;
            if (pj > 0 && pj - 1 < shifts[i] - shifts[i - 1])
            {
                if (ref[ss] != O[ww])
                    gamma = s0;
                else if (ss < upper_boundary || ss >= down_boundary)
                    gamma = s2;
                else
                    gamma = s1;
                G = std::min(G, Gs[shifts[i - 1] + pj - 1] - gamma);
            }
            bool legal_G = (G + gr[ref_sz - ss] <= hAs[ww]), legal_E = (E - v + u + gr[ref_sz - ss] <= hAs[ww] + std::max(int32_t(0), u - qv));
            if (j > 0 || legal_G || legal_E)
            {
                Gs.push_back(G);
                Fs.push_back(F);
                Es.push_back(E);
                ++j;
            }
            else
                --wends[i];
            if (!legal_G && !legal_E && wends[i - 1] - ww >= shifts[i] - shifts[i - 1])
                break;

            if (G - gr[ss] == tAs[ww])
            {
                uint32_t www = ww;
                refline.append(ref, ss);
                queryline.append(ss, '-');
                char type = 'G';
                while (ss != s || ww != w)
                {
                    uint32_t ii = s - ss, jj = wends[ii] - ww, pjj;
                    if (ii > 0)
                        pjj = wends[ii - 1] - ww;
                    switch (type)
                    {
                        case 'G':
                            if (ii > 0 && pjj > 0 && pjj - 1 < shifts[ii] - shifts[ii - 1])
                            {
                                int32_t gg;
                                if (ref[ss] != O[ww])
                                    gg = s0;
                                else if (ss < upper_boundary || ss >= down_boundary)
                                    gg = s2;
                                else
                                    gg = s1;
                                if (Gs[shifts[ii] + jj] == Gs[shifts[ii - 1] + pjj - 1] - gg)
                                {
                                    refline += ref[ss++];
                                    queryline += O[ww++];
                                    break;
                                }
                            }
                            if (Gs[shifts[ii] + jj] == Es[shifts[ii] + jj] - v + u)
                                type = 'E';
                            else
                                type = 'F';
                            break;
                        case 'F':
                            if (Gs[shifts[ii - 1] + pjj] <= Fs[shifts[ii - 1] + pjj])
                                type = 'G';
                            refline += ref[ss++];
                            queryline += '-';
                            break;
                        case 'E':
                            if (Gs[shifts[ii] + jj - 1] <= Es[shifts[ii] + jj - 1])
                                type = 'G';
                            queryline += O[ww++];
                            refline += '-';
                            break;
                    }
                }
                refline.append(ref + s, ref_sz - s);
                queryline.append(ref_sz - s, '-');
                size_t pos = refline.find_first_not_of('-');
                refline[pos] = std::tolower(refline[pos]);
                pos = refline.find_last_not_of('-');
                refline[pos] = std::tolower(refline[pos]);
                return www;
            }
        }
    }
}

std::pair<uint32_t, uint32_t> NodeTrack(const uint32_t w, const int32_t *hAs, const int32_t *hBs, const int32_t *hCs, const uint32_t *hDs, const int32_t qv)
{
    uint32_t ww = w;
    char type = 'A';
    while (true)
    {
        switch (type)
        {
            case 'A':
                if (hAs[ww] == hCs[ww])
                    return std::make_pair(hDs[ww], ww);
                type = 'B';
                break;
            case 'B':
                if (hBs[ww] == hAs[ww - 1] + qv)
                    type = 'A';
                --ww;
                break;
        }
    }
}

#endif
