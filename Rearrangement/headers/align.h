#ifndef ALIGN_H
#define ALIGN_H

#include "parser.h"

const static SCORETYPE inf = std::numeric_limits<SCORETYPE>::max() / 2;
typedef SCORETYPE vsimd __attribute__ ((vector_size (16)));
const static int simd_sz = sizeof(vsimd) / sizeof(SCORETYPE);

void cross_align(SCORETYPE *tAs, SCORETYPE *hAs, SCORETYPE *hBs, SCORETYPE *hCs, size_t *hDs, vsimd *Es, vsimd *Gs, vsimd *Gps, const char *ref, const size_t ref_sz, const char *O, const size_t O_sz, const SCORETYPE u, const SCORETYPE v, const vsimd *grp, const vsimd *grm, vsimd **gamma, const SCORETYPE qu, const SCORETYPE qv, const SCORETYPE s0, const SCORETYPE s1, const SCORETYPE s2, const size_t upper_boundary, const size_t down_boundary)
{
    vsimd infsimd;
    for (size_t k = 0; k < simd_sz; ++k)
        infsimd[k] = -inf;
    size_t seg_sz = ref_sz / simd_sz + 1;
    for (size_t i = 0; i < seg_sz; ++i)
    {
        Es[i] = infsimd;
        Gs[i] = tAs[0] + grp[i];
    }
    hAs[0] = hCs[0] = Gs[0][0] + grp[seg_sz - 1][ref_sz - (seg_sz - 1) * simd_sz];
    hDs[0] = 0;
    hBs[0] = -inf;

    vsimd umask;
    umask[0] = simd_sz - 1;
    for (size_t k = 1; k < simd_sz; ++k)
        umask[k] = k - 1;
    std::map<char, uint8_t> base2int = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}, {'N', 4}};
    for (size_t w = 1; w <= O_sz; ++w)
    {
        uint8_t baseint = base2int[O[w - 1]];
        std::swap(Gs, Gps);
        for (size_t i = 0; i < seg_sz; ++i)
        {
            Es[i] = (Es[i] + u > Gps[i] + v) ? (Es[i] + u) : (Gps[i] + v); 
            Gs[i] = (tAs[w] + grp[i] > Es[i]) ? (tAs[w] + grp[i]) : Es[i];
        }
        hCs[w] = Gs[0][0] + grp[seg_sz - 1][ref_sz - (seg_sz - 1) * simd_sz];
        hDs[w] = 0;

        vsimd G = __builtin_shuffle(Gps[seg_sz - 1] + gamma[seg_sz - 1][baseint], umask);
        Gs[0] = G > Gs[0] ? G : Gs[0];
        vsimd F = Gs[0] + v;
        for (size_t i = 1; i < seg_sz; ++i)
        {
            Gs[i] = F > Gs[i] ? F : Gs[i];
            G = Gps[i - 1] + gamma[i - 1][baseint];
            Gs[i] = G > Gs[i] ? G : Gs[i];
            F = (Gs[i] + v > F + u) ? (Gs[i] + v) : (F + u);
        }

        F = __builtin_shuffle(F, umask);
        F[0] = -inf; 
        for (size_t i = 0; ;)
        {
            bool FisLarge = false;
            for (size_t k = 0; k < simd_sz; ++k)
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
        vsimd hC = Gs[0], hD, hDnow;
        for (size_t k = 0; k < simd_sz; ++k)
            hDnow[k] = k * seg_sz;
        hD = hDnow;
        for (size_t i = 1; i < seg_sz; ++i)
        {
            hDnow += 1;
            vsimd cmp = Gs[i] > hC;
            hC = cmp ? Gs[i] : hC;
            hD = cmp ? hDnow : hD;
        }

        hCs[w] = hC[0];
        hDs[w] = hD[0];
        for (size_t k = 1; k < simd_sz; ++k)
            if (hC[k] > hCs[w])
            {
                hCs[w] = hC[k];
                hDs[w] = hD[k];
            }
        hBs[w] = std::max(hBs[w - 1] + qu, hCs[w - 1] + qv);
        hAs[w] = std::max(hBs[w], hCs[w]);
    }
}

size_t EdgeTrack(std::string &refline, std::string &queryline, const size_t s, const size_t w, const SCORETYPE *tAs, const SCORETYPE *hAs, const char *ref, const size_t ref_sz, const char *O, const size_t O_sz, const SCORETYPE u, const SCORETYPE v, const SCORETYPE *gr, const SCORETYPE qv, const SCORETYPE s0, const SCORETYPE s1, const SCORETYPE s2, const size_t upper_boundary, const size_t down_boundary)
{
    std::vector<SCORETYPE> Es, Fs, Gs;
    std::vector<size_t> shifts, wends;
    shifts.push_back(0);
    wends.push_back(w);
    Gs.push_back(hAs[w] - gr[ref_sz - s]);
    Fs.push_back(inf);
    Es.push_back(inf);
    for (size_t i = 1; ; ++i)
    {
        shifts.push_back(Gs.size());
        wends.push_back(wends[i - 1]);
        for (size_t j = 0; j <= wends[i];)
        {
            SCORETYPE E, F, G, gamma;
            if (j == 0)
                E = inf;
            else
                E = std::min(Es[shifts[i] + j - 1], Gs[shifts[i] + j - 1]) - u;
            size_t pj = wends[i - 1] - wends[i] + j;
            if (pj >= shifts[i] - shifts[i - 1])
                F = inf;
            else
                F = std::min(Fs[shifts[i - 1] + pj], Gs[shifts[i - 1] + pj]) - u;
            G = std::min(E, F) - v + u;
            size_t ss = s - i, ww = wends[i] - j;
            if (pj > 0 && pj - 1 < shifts[i] - shifts[i - 1])
            {
                if (ref[ss] == 'N' || O[ww] == 'N' || ref[ss] != O[ww])
                    gamma = s0;
                else if (ss < upper_boundary || ss >= down_boundary)
                    gamma = s2;
                else
                    gamma = s1;
                G = (G < Gs[shifts[i - 1] + pj - 1] - gamma) ? G : (Gs[shifts[i - 1] + pj - 1] - gamma);
            }
            bool legal_G = (G + gr[ref_sz - ss] <= hAs[ww]), legal_E = (E - v + u + gr[ref_sz - ss] <= hAs[ww] + (qv > u ? 0 : (u - qv)));
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
                size_t www = ww;
                refline.append(ref, ss);
                queryline.append(ss, '-');
                char type = 'G';
                while (ss != s || ww != w)
                {
                    size_t ii = s - ss, jj = wends[ii] - ww, pjj;
                    if (ii > 0)
                        pjj = wends[ii - 1] - ww;
                    switch (type)
                    {
                        case 'G':
                            if (ii > 0 && pjj > 0 && pjj - 1 < shifts[ii] - shifts[ii - 1])
                            {
                                SCORETYPE gg;
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

std::pair<size_t, size_t> NodeTrack(size_t w, const SCORETYPE *hAs, const SCORETYPE *hBs, const SCORETYPE *hCs, const size_t *hDs, const SCORETYPE qv)
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
