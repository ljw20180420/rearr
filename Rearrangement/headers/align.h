#ifndef ALIGN_H
#define ALIGN_H

#include <bits/stdc++.h>

const static int32_t inf = std::numeric_limits<int32_t>::max() / 2;
typedef int32_t v8si __attribute__ ((vector_size (32)));

void cross_align(int32_t *tAs, int32_t *hAs, int32_t *hBs, int32_t *hCs, uint32_t *hDs, v8si *Es, v8si *Gs, v8si *Gps, const char *ref, const uint32_t ref_sz, const char *O, const uint32_t O_sz, const int32_t u, const int32_t v, const v8si *grp, const v8si *grm, v8si **gamma, const int32_t qu, const int32_t qv, const int32_t s0, const int32_t s1, const int32_t s2, const uint32_t upper_boundary, const uint32_t down_boundary)
{
    v8si inf8 = {-inf, -inf, -inf, -inf, -inf, -inf, -inf, -inf};
    uint32_t seg_sz = ref_sz / 8 + 1;
    for (uint32_t i = 0; i < seg_sz; ++i)
    {
        Es[i] = inf8;
        Gs[i] = tAs[0] + grp[i];
    }
    hCs[0] = Gs[0][0] + grp[seg_sz - 1][ref_sz - (seg_sz - 1) * 8];
    hDs[0] = 0;
    hBs[0] = -inf;

    v8si umask = {7, 0, 1, 2, 3, 4, 5, 6};
    std::map<char, uint8_t> base2int = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}, {'N', 4}};
    for (uint32_t w = 1; w <= O_sz; ++w)
    {
        uint8_t baseint = base2int[O[w - 1]];
        std::swap(Gs, Gps);
        for (uint32_t i = 0; i < seg_sz; ++i)
        {
            Es[i] = Es[i] + u > Gps[i] + v ? Es[i] + u : Gps[i] + v; 
            Gs[i] = tAs[w] + grp[i] > Es[i] ? tAs[w] + grp[i] : Es[i];
        }
        hCs[w] = Gs[0][0] + grp[seg_sz - 1][ref_sz - (seg_sz - 1) * 8];
        hDs[w] = 0;

        v8si G = __builtin_shuffle(Gps[seg_sz - 1] + gamma[seg_sz - 1][baseint], umask);
        Gs[0] = G > Gs[0] ? G : Gs[0];
        v8si F = Gs[0] + v;
        for (uint32_t i = 1; i < seg_sz; ++i)
        {
            Gs[i] = F > Gs[i] ? F : Gs[i];
            G = Gps[i - 1] + gamma[i - 1][baseint];
            Gs[i] = G > Gs[i] ? G : Gs[i];
            F = Gs[i] + v > F + u ? Gs[i] + v : F + u;
        }

        F = __builtin_shuffle(F, umask);
        F[0] = -inf; 
        for (uint32_t i = 0; ;)
        {
            bool FisLarge = false;
            for (uint32_t k = 0; k < 8; ++k)
            {
                if (F[k] + u > Gs[i][k] + v)
                {
                    FisLarge = true;
                    break;
                }
            }
            if (!FisLarge)
                break;
            Gs[i] = F > Gs[i] ? F : Gs[i]; 
            F = F + u;
            if (++i == seg_sz)
            {
                F = __builtin_shuffle(F, umask);
                F[0] = -inf;
                i = 0;
            }
        }
        v8si hC = Gs[0], hD, hDnow = {0, 1, 2, 3, 4, 5, 6, 7};
        hDnow *= int32_t(seg_sz);
        hD = hDnow;
        for (uint32_t i = 1; i < seg_sz; ++i)
        {
            hDnow += 1;
            v8si cmp = Gs[i] > hC;
            hC = cmp ? Gs[i] : hC;
            hD = cmp ? hDnow : hD;
        }
        hCs[w] = hC[0];
        hDs[w] = hD[0];
        for (uint32_t k = 1; k < 8; ++k)
            if (hC[k] > hCs[w])
            {
                hCs[w] = hC[k];
                hDs[w] = hD[k];
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
                if (ref[ss] == 'N' || O[ww] == 'N' || ref[ss] != O[ww])
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
                                if (ref[ss] == 'N' || O[ww] == 'N' || ref[ss] != O[ww])
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

std::pair<uint32_t, uint32_t> NodeTrack(uint32_t w, const int32_t *hAs, const int32_t *hBs, const int32_t *hCs, const uint32_t *hDs, const int32_t qv)
{
    char type = 'A';
    while (true)
    {
        switch (type)
        {
            case 'A':
                if (hAs[w] == hCs[w])
                    return std::make_pair(hDs[w], w);
                type = 'B';
                break;
            case 'B':
                if (hBs[w] == hAs[w - 1] + qv)
                    type = 'A';
                --w;
                break;
        }
    }
}

#endif
