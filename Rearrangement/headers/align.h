/**
 * @file align.h
 * @author Jingwei Li (ljw2017@sjtu.edu.cn)
 * @brief Align sequences chimerically.
 * 
 * @details The aligning consists of two parts. The forward dynamic programming gets alignment scores. The backtracking process gets the optimal alignment corresponding the highest score.
 * @date 2024-12-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef ALIGN_H
#define ALIGN_H

#include "parser.h"

/// @brief Define infinity as half of SCORETYPE upper limit. This prevents overflow.
const static SCORETYPE inf = std::numeric_limits<SCORETYPE>::max() / 2;
/// @brief simd vector of 16 bytes.
typedef SCORETYPE vsimd __attribute__ ((vector_size (16)));
const static int simd_sz = sizeof(vsimd) / sizeof(SCORETYPE);

/**
 * @brief Forward dynamic programming recording scores.
 * 
 * @details Traditional DP method record a full score table. cross_align generalize Farrar's SIMD implementation of SW algorithm to chimeric case. cross_align only records one row of scores. This saves memory. The backtracking process in EdgeTrack and NodeTrack is achieved by shearing backtracking branches.
 * 
 * @param[in] tAs Scores up to the alignment of the previous reference.
 * @param[out] hAs Scores up to the alignment of the current reference.
 * @param[out] hBs Auxiliary term to handle affine gap score of unaligned query part.
 * @param[out] hCs Auxiliary term to handle affine gap score of unaligned query part.
 * @param[out] hDs Record the PD source of hCs.
 * @param Es Scores of partial alignments end with a deletion in query.
 * @param Gs Scores of partial alignments end with a (mis)match.
 * @param Gps Auxiliary term recording previous Gs. Handle column-wise DP.
 * @param ref_sz Current reference size.
 * @param O Query sequence.
 * @param O_sz Query size
 * @param u
 * @param v 
 * @param grp Pre-calculated gap penalty for reference start.
 * @param grm Pre-calculated gap penalty for reference end.
 * @param gamma (Mis)matching score matrix.
 * @param qu 
 * @param qv 
 * @param s0 
 * @param s1 
 * @param s2 
 * 
 * @see For s0, s1, s2, u, v, ru, rv, qu, qv, refer to Command_content.
 * @see EdgeTrack
 * @see NodeTrack
 */
void cross_align(
    SCORETYPE *tAs,
    SCORETYPE *hAs,
    SCORETYPE *hBs,
    SCORETYPE *hCs,
    size_t *hDs,
    vsimd *Es,
    vsimd *Gs,
    vsimd *Gps,
    const size_t ref_sz,
    const char *O,
    const size_t O_sz,
    const SCORETYPE u,
    const SCORETYPE v,
    const vsimd *grp,
    const vsimd *grm,
    const vsimd **gamma,
    const SCORETYPE qu,
    const SCORETYPE qv,
    const SCORETYPE s0,
    const SCORETYPE s1,
    const SCORETYPE s2
) {
    /// @brief Assign -inf (SCORETYPE) to vsimd is not allowed. Initialize infsimd with -inf. Then initialize Es with infsimd.
    vsimd infsimd;
    for (size_t k = 0; k < simd_sz; ++k) {
        infsimd[k] = -inf;
    }

    /// @brief Initialize values for w = 0.
    size_t seg_sz = ref_sz / simd_sz + 1;
    for (size_t i = 0; i < seg_sz; ++i) {
        Es[i] = infsimd;
        Gs[i] = tAs[0] + grp[i];
    }
    hAs[0] = hCs[0] = Gs[0][0] + grp[seg_sz - 1][ref_sz - (seg_sz - 1) * simd_sz];
    hDs[0] = 0;
    hBs[0] = -inf;

    /// @brief Initialize umask to shift F and G in Farrar's SIMD implementation.
    vsimd umask;
    umask[0] = simd_sz - 1;
    for (size_t k = 1; k < simd_sz; ++k) {
        umask[k] = k - 1;
    }

    std::map<char, uint8_t> base2int = {{'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}, {'N', 4}};
    for (size_t w = 1; w <= O_sz; ++w) {
        uint8_t baseint = base2int[O[w - 1]];
        /// @brief Move previous Gs to Gps.
        std::swap(Gs, Gps);
        /// @brief Initialize values for current w.
        for (size_t i = 0; i < seg_sz; ++i) {
            Es[i] = (Es[i] + u > Gps[i] + v) ? (Es[i] + u) : (Gps[i] + v); 
            Gs[i] = (tAs[w] + grp[i] > Es[i]) ? (tAs[w] + grp[i]) : Es[i];
        }
        hCs[w] = Gs[0][0] + grp[seg_sz - 1][ref_sz - (seg_sz - 1) * simd_sz];
        hDs[w] = 0;

        /// @brief Update the last G. Apply Farrar's shift on G.
        vsimd G = __builtin_shuffle(Gps[seg_sz - 1] + gamma[seg_sz - 1][baseint], umask);
        Gs[0] = G > Gs[0] ? G : Gs[0];
        vsimd F = Gs[0] + v;
        for (size_t i = 1; i < seg_sz; ++i) {
            Gs[i] = F > Gs[i] ? F : Gs[i];
            G = Gps[i - 1] + gamma[i - 1][baseint];
            Gs[i] = G > Gs[i] ? G : Gs[i];
            F = (Gs[i] + v > F + u) ? (Gs[i] + v) : (F + u);
        }

        /// @brief Update F to its true value like Farrar.
        F = __builtin_shuffle(F, umask);
        F[0] = -inf; 
        for (size_t i = 0; ;) {
            bool FisLarge = false;
            for (size_t k = 0; k < simd_sz; ++k) {
                if (F[k] + u > Gs[i][k] + v) {
                    FisLarge = true;
                    break;
                }
            }
            if (!FisLarge) break;
            Gs[i] = F > Gs[i] ? F : Gs[i]; 
            F = F + u;
            if (++i == seg_sz) {
                F = __builtin_shuffle(F, umask);
                F[0] = -inf;
                i = 0;
            }
        }

        /// @brief hC records the highest score for each simd group. hD records the source of hC.
        vsimd hC = Gs[0], hD, hDnow;
        for (size_t k = 0; k < simd_sz; ++k) {
            hDnow[k] = k * seg_sz;
        }
        hD = hDnow;
        for (size_t i = 1; i < seg_sz; ++i) {
            hDnow += 1;
            vsimd cmp = Gs[i] > hC;
            hC = cmp ? Gs[i] : hC;
            hD = cmp ? hDnow : hD;
        }

        /// @brief Summary all simd group to get the final hCs and hDs for current w.
        hCs[w] = hC[0];
        hDs[w] = hD[0];
        for (size_t k = 1; k < simd_sz; ++k) {
            if (hC[k] > hCs[w]) {
                hCs[w] = hC[k];
                hDs[w] = hD[k];
            }
        }

        /// @brief Update hBs. Then hAs.
        hBs[w] = std::max(hBs[w - 1] + qu, hCs[w - 1] + qv);
        hAs[w] = std::max(hBs[w], hCs[w]);
    }
}

/**
 * @brief Track G(s, w) to its source of A[w'].
 * 
 * @details Instead of the full DP table, only tAs and hAs are recorded. At each step, EdgeTrack tries all possible sources. The score of each source is compared with tAs and hAs. For a real source, the score will be the highest score. To shear branch, EdgeTrack take use of two principles based on the definition of highest score. Firstly, the source score cannot provide a higher score to hAs. Secondly, it cannot require a higher score from tAs.
 * @param[out] refline Reference line of the alignment. It may contain gaps.
 * @param[out] queryline Query line of the alignment. It may contain gaps.
 * @param[in] s Backtrack start position in reference. 
 * @param[in] w Backtrack start position in query.
 * @param[in] tAs 
 * @param[in] hAs 
 * @param[in] ref Tracked reference.
 * @param ref_sz Reference size.
 * @param O Query sequence.
 * @param u 
 * @param v 
 * @param gr Pre-calculated gap penalty for reference start.
 * @param qv 
 * @param s0 
 * @param s1 
 * @param s2 
 * @param upper_boundary Match upstream to upper_boundary has bonus s2 instead of s1.
 * @param down_boundary Match downstream to down_boundary has bonus s2 instead of s1.
 * @return size_t 
 * 
 * @see For s0, s1, s2, u, v, qv, refer to Command_content.
 * @see For tAs, hAs, refer to cross_align.
 */
size_t EdgeTrack(
    std::string &refline,
    std::string &queryline,
    const size_t s,
    const size_t w,
    const SCORETYPE *tAs,
    const SCORETYPE *hAs,
    const char *ref,
    const size_t ref_sz,
    const char *O,
    const SCORETYPE u,
    const SCORETYPE v,
    const SCORETYPE *gr,
    const SCORETYPE qv,
    const SCORETYPE s0,
    const SCORETYPE s1,
    const SCORETYPE s2,
    const size_t upper_boundary,
    const size_t down_boundary
) {
    std::vector<SCORETYPE> Es, Fs, Gs;
    std::vector<size_t> shifts, wends;
    shifts.push_back(0);
    wends.push_back(w);
    Gs.push_back(hAs[w] - gr[ref_sz - s]);
    /// @brief Handle the special case that current reference is skipped completely.
    if (Gs[0] - gr[s] == tAs[w]) {
        refline.append(ref, ref_sz);
        refline[0] = std::tolower(refline[0]);
        refline[ref_sz - 1] = std::tolower(refline[ref_sz - 1]);
        queryline.append(ref_sz, '-');
        return w;
    }
    /// @brief Handle the row s.
    Fs.push_back(inf);
    Es.push_back(inf);
    for (size_t j = 1; j <= w; ++j) {
        /// @brief Get scores of possible sources at (s, w - j).
        SCORETYPE E = std::min(Es[j - 1], Gs[j - 1]) - u;
        SCORETYPE F = inf;
        SCORETYPE G = std::min(E, F) - v + u;
        /// @brief Check whether source scores satisfy hAs restriction.
        bool legal_G = (G + gr[ref_sz - s] <= hAs[w - j]);
        bool legal_E = (E - v + u + gr[ref_sz - s] <= hAs[w - j] + (qv > u ? 0 : (u - qv)));
        /// @brief If not, then hAs restriction is impossible for further upstream position in query, so break.
        if (!legal_G && !legal_E) {
            break;
        }
        /// @brief If so, then record EFG for (s, w - j).
        Gs.push_back(G);
        Fs.push_back(F);
        Es.push_back(E);
        /// @brief If G has a source in tAs, then the track is over.
        if (G - gr[s] == tAs[w - j]) {
            refline.append(ref, s);
            queryline.append(s, '-');
            refline.append(j, '-');
            queryline.append(O + w - j, j);
            refline.append(ref + s, ref_sz - s);
            queryline.append(ref_sz - s, '-');
            return w - j;
        }
    }
    /// @brief Track row s - i.
    for (size_t i = 1; ; ++i) {
        shifts.push_back(Gs.size());
        wends.push_back(wends[i - 1]);
        for (size_t j = 0; j <= wends[i];) {
            /// @brief Get scores of possible sources at (s - i, wends[i] - j).
            SCORETYPE E, F, G, gamma;
            if (j == 0) {
                E = inf;
            } else {
                E = std::min(Es[shifts[i] + j - 1], Gs[shifts[i] + j - 1]) - u;   
            }
            size_t pj = wends[i - 1] - wends[i] + j;
            if (pj >= shifts[i] - shifts[i - 1]) {
                F = inf;
            } else {
                F = std::min(Fs[shifts[i - 1] + pj], Gs[shifts[i - 1] + pj]) - u;
            }
            G = std::min(E, F) - v + u;
            size_t ss = s - i, ww = wends[i] - j;
            if (pj > 0 && pj - 1 < shifts[i] - shifts[i - 1]) {
                if (ref[ss] == 'N' || O[ww] == 'N' || ref[ss] != O[ww]) {
                    gamma = s0;
                } else if (ss < upper_boundary || ss >= down_boundary) {
                    gamma = s2;
                } else {
                    gamma = s1;
                }
                G = (G < Gs[shifts[i - 1] + pj - 1] - gamma) ? G : (Gs[shifts[i - 1] + pj - 1] - gamma);
            }
            /// @brief Check whether source scores satisfy hAs restriction. Determine the range.
            bool legal_G = (G + gr[ref_sz - ss] <= hAs[ww]);
            bool legal_E = (E - v + u + gr[ref_sz - ss] <= hAs[ww] + (qv > u ? 0 : (u - qv)));
            if (j > 0 || legal_G || legal_E) {
                Gs.push_back(G);
                Fs.push_back(F);
                Es.push_back(E);
                ++j;
            } else {
                --wends[i];
            }
            if (!legal_G && !legal_E && wends[i - 1] - ww >= shifts[i] - shifts[i - 1]) break;

            /// @brief If G has a source in tAs, then the track is over.
            if (G - gr[ss] == tAs[ww]) {
                size_t www = ww;
                /// @brief Track the actual alignment from (ss, ww) to (s, w).
                refline.append(ref, ss);
                queryline.append(ss, '-');
                char type = 'G';
                while (ss != s || ww != w) {
                    size_t ii = s - ss, jj = wends[ii] - ww, pjj;
                    if (ii > 0) pjj = wends[ii - 1] - ww;
                    switch (type) {
                        case 'G':
                            if (ii > 0 && pjj > 0 && pjj - 1 < shifts[ii] - shifts[ii - 1]) {
                                SCORETYPE gg;
                                if (ref[ss] == 'N' || O[ww] == 'N' || ref[ss] != O[ww]) {
                                    gg = s0;
                                } else if (ss < upper_boundary || ss >= down_boundary) {
                                    gg = s2;
                                } else {
                                    gg = s1;
                                }
                                if (Gs[shifts[ii] + jj] == Gs[shifts[ii - 1] + pjj - 1] - gg) {
                                    refline += ref[ss++];
                                    queryline += O[ww++];
                                    break;
                                }
                            }
                            if (Gs[shifts[ii] + jj] == Es[shifts[ii] + jj] - v + u) {
                                type = 'E';
                            } else {
                                type = 'F';
                            }
                            break;
                        case 'F':
                            if (Gs[shifts[ii - 1] + pjj] <= Fs[shifts[ii - 1] + pjj]) type = 'G';
                            refline += ref[ss++];
                            queryline += '-';
                            break;
                        case 'E':
                            if (Gs[shifts[ii] + jj - 1] <= Es[shifts[ii] + jj - 1]) type = 'G';
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

/**
 * @brief Backtrack from hAs until finding its source at hCs.
 * 
 * @details hAs[w] comes from either hBs[w] or hCs[w]. hBs[w] comes from either hBs[w-1] or hCs[w-1]. Because hBs[0] = -inf, hAs[w] must come from hCs[w'] for some w' <= w. NodeTrack finds this w', and returns (hDs[w'], w'). hDs[w'] saves the source of hCs[w'], i.e. hCs[w'] = Gs[hDs[w']].
 * @param[in] w The current iteration column.
 * @param[in] hAs 
 * @param[in] hBs 
 * @param[in] hCs 
 * @param[in] hDs 
 * @param[in] qv 
 * @return std::pair<size_t, size_t> 
 * @see For hAs, hBs, hCs, hDs, refer to cross_align. For qv, refer to Command_content.
 */
std::pair<size_t, size_t> NodeTrack(
    size_t w,
    const SCORETYPE *hAs,
    const SCORETYPE *hBs,
    const SCORETYPE *hCs,
    const size_t *hDs,
    const SCORETYPE qv
) {
    char type = 'A';
    while (true) {
        switch (type) {
            case 'A':
                if (hAs[w] == hCs[w]) return std::make_pair(hDs[w], w);
                type = 'B';
                break;
            case 'B':
                if (hBs[w] == hAs[w - 1] + qv) type = 'A';
                --w;
                break;
        }
    }
}

#endif
