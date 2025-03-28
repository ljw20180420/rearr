#!/usr/bin/env python

from manim import *
import numpy as np

rng = np.random.default_rng(63036)


def generate_random_DNA(length, rng):
    # usage: generate_random_DNA of length
    return "".join(rng.choice(["A", "C", "G", "T"], length))


def SNP_DNA(DNA, probability, rng):
    DNA = list(DNA)
    for i in range(len(DNA)):
        if rng.random() < probability:
            DNA[i] = rng.choice(["A", "C", "G", "T"])
    return "".join(DNA)


def indel_DNA(DNA, probability, rng):
    inslens = rng.negative_binomial(1, 1 - probability, len(DNA) + 1)
    insDNAs = generate_random_DNA(np.sum(inslens), rng)
    DNAnew, start = [], 0
    for i in range(len(DNA)):
        end = start + inslens[i]
        DNAnew.append(insDNAs[start:end])
        if rng.random() > probability:
            DNAnew.append(DNA[i])
        start = end
    DNAnew.append(insDNAs[start:])
    return "".join(DNAnew)


ref = generate_random_DNA(18, rng)
ref = ref[:13] + "GG" + ref[15:]
region1 = indel_DNA(SNP_DNA(ref[2:-2], 0.05, rng), 0.05, rng)
region2 = indel_DNA(SNP_DNA(ref[2:-2], 0.05, rng), 0.05, rng)
query = (
    generate_random_DNA(2, rng)
    + region1
    + generate_random_DNA(32 - len(region1) - len(region2), rng)
    + region2
    + generate_random_DNA(2, rng)
)


def initialize_number_plane(ref, query):
    number_plane = NumberPlane(
        x_range=[0, len(query) + 0.1],
        y_range=[0, len(ref) + 0.1],
        x_length=12,
        y_length=6,
        y_axis_config={"rotation": -np.pi / 2},
        tips=False,
    )
    x_ax, y_ax = number_plane.get_axes()
    x_ax.add_labels({i + 1 / 2: base for i, base in enumerate(query)}, direction=UP)
    y_ax.add_labels({i + 1 / 2: base for i, base in enumerate(ref)}, direction=LEFT)
    return number_plane


def initialize_node_scores(qu, qv, C_scores, number_plane):
    B_scores = [None for j in range(len(query) + 1)]
    A_scores = [None for j in range(len(query) + 1)]
    for j in range(0, len(query) + 1):
        coor = number_plane.p2c(C_scores[j].get_center())

        B_scores[j] = (
            Integer(
                -999999
                if j == 0
                else max(B_scores[j - 1].number + qu, C_scores[j - 1].number + qv)
            )
            .move_to(number_plane.c2p(coor[0] - 1 / 2, coor[1] + 1 / 2, 0))
            .scale(1 / 4)
        )
        A_scores[j] = (
            Integer(max(B_scores[j].number, C_scores[j].number))
            .move_to(number_plane.c2p(coor[0], coor[1] + 1 / 2, 0))
            .scale(1 / 4)
        )
    return A_scores, B_scores


def initialize_edge_scores(s0, s1, u, v, ru, rv, ref, query, A_scores, number_plane):
    G_scores = [[None for j in range(len(query) + 1)] for i in range(len(ref) + 1)]
    M_scores = [[None for j in range(len(query) + 1)] for i in range(len(ref) + 1)]
    R_scores = [[None for j in range(len(query) + 1)] for i in range(len(ref) + 1)]
    E_scores = [[None for j in range(len(query) + 1)] for i in range(len(ref) + 1)]
    F_scores = [[None for j in range(len(query) + 1)] for i in range(len(ref) + 1)]
    C_scores = [None for j in range(len(query) + 1)]

    # fill R_scores (restart score)
    for j in range(len(query) + 1):
        R_scores[0][j] = (
            Integer(A_scores[j].number)
            .move_to(number_plane.c2p(j + 1 / 4, 1 / 4, 0))
            .scale(1 / 4)
        )
        for i in range(1, len(ref) + 1):
            gap_score = rv if i == 1 else ru
            R_scores[i][j] = (
                Integer(R_scores[i - 1][j].number + gap_score)
                .move_to(number_plane.c2p(j + 1 / 4, i + 1 / 4, 0))
                .scale(1 / 4)
            )

    # fill boundary M_scores ([mis]match score)
    for j in range(len(query) + 1):
        M_scores[0][j] = (
            Integer(-999999)
            .move_to(number_plane.c2p(j - 1 / 4, -1 / 4, 0))
            .scale(1 / 4)
        )
    for i in range(1, len(ref) + 1):
        M_scores[i][0] = (
            Integer(-999999)
            .move_to(number_plane.c2p(-1 / 4, i - 1 / 4, 0))
            .scale(1 / 4)
        )

    # fill boundary E_scores (horizontal gap score) at y-axis
    for i in range(len(ref) + 1):
        E_scores[i][0] = (
            Integer(-999999)
            .move_to(number_plane.c2p(-1 / 4, i + 1 / 4, 0))
            .scale(1 / 4)
        )

    # fill boundary F_scores (vertical gap score) at x-axis
    for j in range(len(query) + 1):
        F_scores[0][j] = (
            Integer(-999999)
            .move_to(number_plane.c2p(j + 1 / 4, -1 / 4, 0))
            .scale(1 / 4)
        )

    # fill corner score
    G_scores[0][0] = (
        Integer(
            max(
                R_scores[0][0].number,
                M_scores[0][0].number,
                F_scores[0][0].number,
                E_scores[0][0].number,
            )
        )
        .move_to(number_plane.c2p(1 / 4, 1 / 4, 0))
        .scale(1 / 4)
    )

    # fill boundary E_scores and G_scores at x-axis
    for j in range(1, len(query) + 1):
        E_scores[0][j] = (
            Integer(max(E_scores[0][j - 1].number + u, G_scores[0][j - 1].number + v))
            .move_to(number_plane.c2p(j - 1 / 4, 1 / 4, 0))
            .scale(1 / 4)
        )
        G_scores[0][j] = (
            Integer(
                max(
                    R_scores[0][j].number,
                    M_scores[0][j].number,
                    F_scores[0][j].number,
                    E_scores[0][j].number,
                )
            )
            .move_to(number_plane.c2p(j + 1 / 4, 1 / 4, 0))
            .scale(1 / 4)
        )

    # fill boundary F_scores and G_scores at y-axis
    for i in range(1, len(ref) + 1):
        F_scores[i][0] = (
            Integer(max(F_scores[i - 1][0].number + u, G_scores[i - 1][0].number + v))
            .move_to(number_plane.c2p(1 / 4, i - 1 / 4, 0))
            .scale(1 / 4)
        )
        G_scores[i][0] = (
            Integer(
                max(
                    R_scores[i][0].number,
                    M_scores[i][0].number,
                    F_scores[i][0].number,
                    E_scores[i][0].number,
                )
            )
            .move_to(number_plane.c2p(1 / 4, i + 1 / 4, 0))
            .scale(1 / 4)
        )

    # fill main body E_scores, F_scores, M_scores, G_scores
    for i in range(1, len(ref) + 1):
        for j in range(1, len(query) + 1):
            E_scores[i][j] = (
                Integer(
                    max(E_scores[i][j - 1].number + u, G_scores[i][j - 1].number + v)
                )
                .move_to(number_plane.c2p(j - 1 / 4, i + 1 / 4, 0))
                .scale(1 / 4)
            )
            F_scores[i][j] = (
                Integer(
                    max(F_scores[i - 1][j].number + u, G_scores[i - 1][j].number + v)
                )
                .move_to(number_plane.c2p(j + 1 / 4, i - 1 / 4, 0))
                .scale(1 / 4)
            )
            M_scores[i][j] = (
                Integer(
                    G_scores[i - 1][j - 1].number
                    + (s0 if ref[i - 1] != query[j - 1] else s1)
                )
                .move_to(number_plane.c2p(j - 1 / 4, i - 1 / 4, 0))
                .scale(1 / 4)
            )
            G_scores[i][j] = (
                Integer(
                    max(
                        R_scores[i][j].number,
                        M_scores[i][j].number,
                        F_scores[i][j].number,
                        E_scores[i][j].number,
                    )
                )
                .move_to(number_plane.c2p(j + 1 / 4, i + 1 / 4, 0))
                .scale(1 / 4)
            )

    # fill C_scores
    C_scores, C_tracks = list(), list()
    for j in range(len(query) + 1):
        max_G_score = -np.inf
        for i in range(len(ref) + 1):
            if i == len(ref):
                gap_score = 0
            else:
                gap_score = (len(ref) - i - 1) * ru + rv
            if G_scores[i][j].number + gap_score > max_G_score:
                max_G_score = G_scores[i][j].number
                C_track = i
        C_scores.append(
            Integer(max_G_score)
            .move_to(number_plane.c2p(j + 1 / 4, len(ref) + 3 / 4, 0))
            .scale(1 / 4)
        )
        C_tracks.append(C_track)

    return G_scores, R_scores, M_scores, E_scores, F_scores, C_scores, C_tracks


def initialize_arrows(ref, query, number_plane):
    E_arrows = [[None for j in range(len(query) + 1)] for i in range(len(ref) + 1)]
    for i in range(len(ref) + 1):
        for j in range(1, len(query) + 1):
            E_arrows[i][j] = Arrow(
                start=number_plane.c2p(j - 1, i, 0),
                end=number_plane.c2p(j, i, 0),
                color=RED,
            )
    F_arrows = [[None for j in range(len(query) + 1)] for i in range(len(ref) + 1)]
    for i in range(1, len(ref) + 1):
        for j in range(len(query) + 1):
            F_arrows[i][j] = Arrow(
                start=number_plane.c2p(j, i - 1, 0),
                end=number_plane.c2p(j, i, 0),
                color=RED,
            )
    G_arrows = [[None for j in range(len(query) + 1)] for i in range(len(ref) + 1)]
    for i in range(1, len(ref) + 1):
        for j in range(1, len(query) + 1):
            G_arrows[i][j] = Arrow(
                start=number_plane.c2p(j - 1, i - 1, 0),
                end=number_plane.c2p(j, i, 0),
                color=RED,
            )

    T_arrows = [[None for j in range(len(query) + 1)] for i in range(len(ref) + 1)]
    for i in range(len(ref) + 1):
        for j in range(len(query) + 1):
            T_arrows[i][j] = Arrow(
                start=number_plane.c2p(j, -1, 0),
                end=number_plane.c2p(j, i, 0),
                color=RED,
            )

    H_arrows = [[None for j in range(len(query) + 1)] for i in range(len(ref) + 1)]
    for i in range(len(ref) + 1):
        for j in range(len(query) + 1):
            H_arrows[i][j] = Arrow(
                start=number_plane.c2p(j, i, 0),
                end=number_plane.c2p(j, len(ref) + 1, 0),
                color=RED,
            )

    N_arrows = [None for j in range(len(query) + 1)]
    for j in range(1, len(query) + 1):
        N_arrows[j] = Arrow(
            start=number_plane.c2p(j - 1, len(ref) + 1, 0),
            end=number_plane.c2p(j, len(ref) + 1, 0),
            color=RED,
        )

    return E_arrows, F_arrows, G_arrows, T_arrows, H_arrows, N_arrows


def track_node(qv, A_scores, B_scores, C_scores, j, modes):
    if modes[-1] == "C":
        modes.append("G")
        return j, modes
    elif modes[-1] == "A":
        if A_scores[j].number == C_scores[j].number:
            modes.append("C")
        else:
            modes.append("B")
        return track_node(qv, A_scores, B_scores, C_scores, j, modes)
    else:
        if B_scores[j].number == C_scores[j - 1].number + qv:
            modes.append("C")
        else:
            modes.append("B")
        return track_node(qv, A_scores, B_scores, C_scores, j - 1, modes)


def track_edge(v, E_scores, F_scores, G_scores, R_scores, i, j, modes):
    if modes[-1] == "G":
        if G_scores[i][j].number == R_scores[i][j].number:
            modes.append("A")
            return j, modes
        elif G_scores[i][j].number == E_scores[i][j].number:
            modes.append("E")
            return track_edge(v, E_scores, F_scores, G_scores, R_scores, i, j, modes)
        elif G_scores[i][j].number == F_scores[i][j].number:
            modes.append("F")
            return track_edge(v, E_scores, F_scores, G_scores, R_scores, i, j, modes)
        else:
            modes.append("G")
            return track_edge(
                v,
                E_scores,
                F_scores,
                G_scores,
                R_scores,
                i - 1,
                j - 1,
                modes,
            )
    elif modes[-1] == "E":
        if E_scores[i][j].number == G_scores[i][j - 1].number + v:
            modes.append("G")
        else:
            modes.append("E")
        return track_edge(v, E_scores, F_scores, G_scores, R_scores, i, j - 1, modes)
    elif modes[-1] == "F":
        if F_scores[i][j].number == G_scores[i - 1][j].number + v:
            modes.append("G")
        else:
            modes.append("F")
        return track_edge(v, E_scores, F_scores, G_scores, R_scores, i - 1, j, modes)


def track_arrow(i, j, C_tracks, mode1, mode2):
    if mode1 == "A":
        return i, j
    if mode1 == "B":
        return i, j - 1
    if mode1 == "C":
        return C_tracks[j], j
    if mode1 == "G":
        if mode2 == "A":
            return -1, j
        if mode2 == "G":
            return i - 1, j - 1
        return i, j
    if mode1 == "E":
        return i, j - 1
    return i - 1, j


def track_arrows(i, j, C_tracks, modes):
    path = [(i, j)]
    for k in range(len(modes) - 1):
        path.append(
            track_arrow(path[-1][0], path[-1][1], C_tracks, modes[k], modes[k + 1])
        )
    return path


class ChimericAlignment(Scene):
    def construct(self):
        self.next_section("initialize sequence and DP table")
        number_plane = initialize_number_plane(ref, query)
        self.play(Create(number_plane))
        self.wait()

        self.next_section("initialize boundary scores")
        C_scores0 = [None for i in range(len(query) + 1)]
        for j in range(len(query) + 1):
            C_scores0[j] = (
                Integer(0 if j == 0 else -999999)
                .move_to(number_plane.c2p(j + 1 / 4, -5 / 4, 0))
                .scale(1 / 4)
            )

        A_scores0, B_scores0 = initialize_node_scores(0, 0, C_scores0, number_plane)
        G_scores1, R_scores1, M_scores1, E_scores1, F_scores1, C_scores1, C_tracks1 = (
            initialize_edge_scores(
                -3, 1, -5, -5, 0, 0, ref, query, A_scores0, number_plane
            )
        )
        A_scores1, B_scores1 = initialize_node_scores(0, 0, C_scores1, number_plane)

        E_arrows, F_arrows, G_arrows, T_arrows, H_arrows, N_arrows = initialize_arrows(
            ref, query, number_plane
        )

        self.play(FadeIn(G_scores1[0][0]))
        for j in range(1, 4):
            self.play(
                GrowArrow(E_arrows[0][j]),
                FadeIn(R_scores1[0][j]),
                FadeIn(E_scores1[0][j]),
            )
            self.play(
                FadeOut(E_arrows[0][j]),
                FadeOut(R_scores1[0][j]),
                FadeOut(E_scores1[0][j]),
                FadeIn(G_scores1[0][j]),
            )
        self.play([FadeIn(G_scores1[0][j]) for j in range(4, len(query) + 1)])
        self.play([FadeIn(G_scores1[i][0]) for i in range(1, len(ref) + 1)])

        self.next_section("fill scores of DP table")
        for j in range(1, 4):
            self.play(
                FadeIn(R_scores1[1][j]),
                GrowArrow(E_arrows[1][j]),
                GrowArrow(F_arrows[1][j]),
                GrowArrow(G_arrows[1][j]),
                FadeIn(E_scores1[1][j]),
                FadeIn(F_scores1[1][j]),
                FadeIn(M_scores1[1][j]),
            )
            self.play(
                FadeOut(R_scores1[1][j]),
                FadeOut(E_arrows[1][j]),
                FadeOut(F_arrows[1][j]),
                FadeOut(G_arrows[1][j]),
                FadeOut(E_scores1[1][j]),
                FadeOut(F_scores1[1][j]),
                FadeOut(M_scores1[1][j]),
                FadeIn(G_scores1[1][j]),
            )
        self.play([FadeIn(G_scores1[1][j]) for j in range(4, len(query) + 1)])
        self.play(
            [
                FadeIn(G_scores1[i][j])
                for i in range(2, len(ref) + 1)
                for j in range(len(query) + 1)
            ]
        )

        self.next_section("track best alignment")
        j_start, _ = track_node(0, A_scores1, B_scores1, C_scores1, len(query), ["A"])
        j_end, modes = track_edge(
            -5,
            E_scores1,
            F_scores1,
            G_scores1,
            R_scores1,
            C_tracks1[j_start],
            j_start,
            ["G"],
        )
        path1 = track_arrows(C_tracks1[j_start], j_start, C_tracks1, modes[:-1])
        path_arrow1 = []
        assert path1[-1][1] == j_end, "j_end is not consistent"
        for k in range(1, len(path1)):
            if path1[k][1] < path1[k - 1][1] or path1[k][0] < path1[k - 1][0]:
                path_arrow1.append(
                    Arrow(
                        start=number_plane.c2p(path1[k - 1][1], path1[k - 1][0], 0),
                        end=number_plane.c2p(path1[k][1], path1[k][0], 0),
                        color=GREEN,
                    )
                )
                self.play(
                    GrowArrow(path_arrow1[-1]),
                )

        self.next_section("calculate node scores")
        for j in range(3):
            self.play(GrowArrow(H_arrows[0][j]), FadeIn(C_scores1[j]))
            self.play(FadeOut(H_arrows[0][j]))
        for j in range(3, len(query) + 1):
            self.play(FadeIn(C_scores1[j]))
        for j in range(1, 4):
            self.play(GrowArrow(N_arrows[j]), FadeIn(B_scores1[j]))
            self.play(
                FadeOut(N_arrows[j]),
                FadeIn(A_scores1[j]),
            )
        self.play([FadeIn(A_scores1[j]) for j in range(4, len(query) + 1)])

        self.next_section("alignment second")
        self.play(
            [
                A_scores1[j].animate.move_to(number_plane.c2p(j + 1 / 4, -3 / 4))
                for j in range(len(query) + 1)
            ]
            + [
                FadeOut(G_scores1[i][j])
                for i in range(len(ref) + 1)
                for j in range(len(query) + 1)
            ]
            + [FadeOut(B_scores1[j]) for j in range(1, len(query) + 1)]
            + [FadeOut(C_scores1[j]) for j in range(len(query) + 1)]
        )

        G_scores2, R_scores2, M_scores2, E_scores2, F_scores2, C_scores2, C_tracks2 = (
            initialize_edge_scores(
                -3, 1, -5, -5, 0, 0, ref, query, A_scores1, number_plane
            )
        )
        A_scores2, B_scores2 = initialize_node_scores(0, 0, C_scores2, number_plane)

        self.play(GrowArrow(T_arrows[0][0]), FadeIn(R_scores2[0][0]))
        self.play(
            FadeOut(T_arrows[0][0]), FadeOut(R_scores2[0][0]), FadeIn(G_scores2[0][0])
        )
        for j in range(1, 4):
            self.play(
                GrowArrow(E_arrows[0][j]),
                GrowArrow(T_arrows[0][j]),
                FadeIn(R_scores2[0][j]),
                FadeIn(E_scores2[0][j]),
            )
            self.play(
                FadeOut(E_arrows[0][j]),
                FadeOut(T_arrows[0][j]),
                FadeOut(R_scores2[0][j]),
                FadeOut(E_scores2[0][j]),
                FadeIn(G_scores2[0][j]),
            )
        self.play([FadeIn(G_scores2[0][j]) for j in range(4, len(query) + 1)])
        self.play([FadeIn(G_scores2[i][0]) for i in range(1, len(ref) + 1)])

        for j in range(1, 4):
            self.play(
                GrowArrow(T_arrows[1][j]),
                GrowArrow(E_arrows[1][j]),
                GrowArrow(F_arrows[1][j]),
                GrowArrow(G_arrows[1][j]),
                FadeIn(R_scores2[1][j]),
                FadeIn(E_scores2[1][j]),
                FadeIn(F_scores2[1][j]),
                FadeIn(M_scores2[1][j]),
            )
            self.play(
                FadeOut(T_arrows[1][j]),
                FadeOut(E_arrows[1][j]),
                FadeOut(F_arrows[1][j]),
                FadeOut(G_arrows[1][j]),
                FadeOut(R_scores2[1][j]),
                FadeOut(E_scores2[1][j]),
                FadeOut(F_scores2[1][j]),
                FadeOut(M_scores2[1][j]),
                FadeIn(G_scores2[1][j]),
            )
        self.play([FadeIn(G_scores2[1][j]) for j in range(4, len(query) + 1)])
        self.play(
            [
                FadeIn(G_scores2[i][j])
                for i in range(2, len(ref) + 1)
                for j in range(len(query) + 1)
            ]
        )

        self.next_section("track second")
        j_start, _ = track_node(0, A_scores2, B_scores2, C_scores2, len(query), ["A"])
        j_end, modes = track_edge(
            -5,
            E_scores2,
            F_scores2,
            G_scores2,
            R_scores2,
            C_tracks2[j_start],
            j_start,
            ["G"],
        )
        path2 = track_arrows(C_tracks2[j_start], j_start, C_tracks2, modes[:-1])
        path_arrow2 = []
        assert path2[-1][1] == j_end, "j_end is not consistent"
        for k in range(1, len(path2)):
            if path2[k][1] < path2[k - 1][1] or path2[k][0] < path2[k - 1][0]:
                path_arrow2.append(
                    Arrow(
                        start=number_plane.c2p(path2[k - 1][1], path2[k - 1][0], 0),
                        end=number_plane.c2p(path2[k][1], path2[k][0], 0),
                        color=GREEN,
                    )
                )
                self.play(
                    GrowArrow(path_arrow2[-1]),
                )

        self.next_section("alignment more")

        self.next_section("track more")

        self.next_section("extend to affine gap score")

        self.next_section("score unaligned query")

        self.next_section("score unaligned reference")

        self.next_section("empty DP table while still track (fail)")

        self.next_section("fill reverse DP table")

        self.next_section(
            "recall forward DP table, show that forward score is lower than reverse score, but always coincide on the best alignment path"
        )

        self.next_section("show that reverse DP table can be truncated")

        self.next_section("show one-liner forward DP")

        self.next_section("show classic SIMD")

        self.next_section("generalize SIMD to chimeric alignment")
