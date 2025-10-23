#!/usr/bin/env python

import manim
from manim_dna import ManimDNA
from Bio import Seq


class Stagger(manim.Scene):
    def construct(self):
        manimDNA = ManimDNA()
        watson = manimDNA.generate_random_DNA(20)
        watson = watson[:14] + "GG" + watson[16:]
        crick = Seq.Seq(watson).complement().__str__()

        watson, crick = manimDNA.display_DNA(
            [watson, crick], [-10, -10], [1, 0], self, color=[manim.WHITE] * 2
        )
        manimDNA.change_color(watson[13:16], manim.ManimColor("#FF0000"), self)
        manimDNA.move_DNA(watson[9:] + crick[10:], manim.RIGHT, 3, self)
        _, crick_i = manimDNA.display_DNA(
            [watson[9].text, crick[9].text],
            [-1, 2],
            [1, 0],
            self,
            [manim.ManimColor("#00FF00")] * 2,
        )
        manimDNA.move_DNA(watson[9:] + crick[10:] + crick_i, manim.LEFT, 2, self)


class MMEJ(manim.Scene):
    def construct(self):
        self.next_section("Resection")

        manimDNA = ManimDNA()
        watson = "TACGTGGTGCTCGGGGGGCT"
        crick = Seq.Seq(watson).complement().__str__()

        watson, crick = manimDNA.display_DNA(
            [watson, crick], [-10, -10], [1, 0], self, color=[manim.WHITE] * 2
        )

        manimDNA.change_color(watson[13:16], manim.ManimColor("#FF0000"), self)
        manimDNA.move_DNA(watson[10:] + crick[10:], manim.RIGHT, 3, self)

        for i in range(7):
            manimDNA.delete_DNA([watson[10 + i], crick[9 - i]], self, run_time=0.2)

        self.next_section("MMEJ")

        manimDNA.move_DNA(watson[-3:] + crick[10:], manim.DOWN, 1, self)
        manimDNA.move_DNA(watson[-3:] + crick[10:], manim.LEFT, 10, self)

        sticks = manimDNA.display_DNA("|| |", -5, 0, self)
        manimDNA.delete_DNA([watson[9], crick[10]], self, run_time=0.2)
        manimDNA.delete_DNA([watson[8], crick[11], sticks[-1]], self, run_time=0.2)
        manimDNA.delete_DNA([watson[7]], self, run_time=0.2)

        _, waston_sticks, crick_fills, crick_sticks = manimDNA.display_DNA(
            ["GGG", "|||", "CA", "||"],
            [-3, -3, -7, -7],
            [1, 0, -1, 0],
            self,
            [
                manim.ManimColor("#00FF00"),
                manim.WHITE,
                manim.ManimColor("#00FF00"),
                manim.WHITE,
            ],
        )

        manimDNA.delete_DNA(sticks[:2] + waston_sticks + crick_sticks, self)
        manimDNA.move_DNA(watson[-3:] + crick[12:] + crick_fills, manim.UP, 1, self)


class UNILATERAL(manim.Scene):
    def construct(self):
        self.next_section("Resection")

        manimDNA = ManimDNA()
        watson = "TACGTGGTGCTCGGGGGGCT"
        crick = Seq.Seq(watson).complement().__str__()

        watson, crick = manimDNA.display_DNA(
            [watson, crick], [-10, -10], [1, 0], self, color=[manim.WHITE] * 2
        )

        manimDNA.change_color(watson[13:16], manim.ManimColor("#FF0000"), self)
        manimDNA.move_DNA(watson[10:] + crick[10:], manim.RIGHT, 3, self)

        for i in range(7):
            manimDNA.delete_DNA([watson[10 + i], crick[9 - i]], self, run_time=0.2)

        self.next_section("Unilateral")

        manimDNA.move_DNA(watson[-3:] + crick[10:], manim.DOWN, 1, self)
        manimDNA.move_DNA(watson[-3:] + crick[10:], manim.LEFT, 10, self)

        sticks = manimDNA.display_DNA("|| |", -5, 0, self)

        manimDNA.move_DNA(watson[-3:] + crick[-3:], manim.UP, 1, self)
        manimDNA.delete_DNA(crick[10:17] + sticks, self)

        manimDNA.display_DNA("CACCACG", -7, 0, self, manim.ManimColor("#00FF00"))


class UNILATERAL2(manim.Scene):
    def construct(self):
        self.next_section("Resection")

        manimDNA = ManimDNA()
        watson = "TACGTGGTGCTCGGGGGGCT"
        crick = Seq.Seq(watson).complement().__str__()

        watson, crick = manimDNA.display_DNA(
            [watson, crick], [-10, -10], [1, 0], self, color=[manim.WHITE] * 2
        )

        manimDNA.change_color(watson[13:16], manim.ManimColor("#FF0000"), self)
        manimDNA.move_DNA(watson[10:] + crick[10:], manim.RIGHT, 3, self)

        for i in range(7):
            manimDNA.delete_DNA([watson[10 + i], crick[9 - i]], self, run_time=0.2)

        self.next_section("Unilateral")

        manimDNA.move_DNA(watson[-3:] + crick[10:], manim.DOWN, 1, self)
        manimDNA.move_DNA(watson[-3:] + crick[10:], manim.LEFT, 10, self)

        sticks = manimDNA.display_DNA("|| |", -5, 0, self)

        manimDNA.move_DNA(watson[:3] + crick[:3], manim.DOWN, 1, self)
        manimDNA.delete_DNA(watson[3:10] + sticks, self)

        manimDNA.display_DNA("TCGGGGG", -7, 0, self, manim.ManimColor("#00FF00"))


class CHIMERIC(manim.Scene):
    def construct(self):
        manimDNA = ManimDNA()

        watson = manimDNA.display_DNA(
            "TACGTGGTGCTCGGGGGGCT", -10, 1, self, color=manim.WHITE
        )

        repair_type = (
            manim.Text("blunt end")
            .move_to(manimDNA.number_plane.c2p(-15, 8, 0))
            .scale(0.5)
        )
        self.play(manim.FadeIn(repair_type))

        query1, query2 = manimDNA.display_DNA(
            ["TACGTGGTGC", "TCGGGGGGCT"],
            [-10, 0],
            [0, 0],
            self,
            color=[manim.WHITE] * 2,
        )

        self.play(manim.FadeOut(repair_type))
        repair_type = (
            manim.Text("unilateral")
            .move_to(manimDNA.number_plane.c2p(-15, 8, 0))
            .scale(0.5)
        )
        self.play(manim.FadeIn(repair_type))

        manimDNA.delete_DNA(query1[3:], self)
        self.play([manim.FadeIn(base) for base in query1[3:]])
        manimDNA.delete_DNA(query2[:-3], self)
        self.play([manim.FadeIn(base) for base in query2[:-3]])

        self.play(manim.FadeOut(repair_type))
        repair_type = (
            manim.Text("MMEJ").move_to(manimDNA.number_plane.c2p(-15, 8, 0)).scale(0.5)
        )
        self.play(manim.FadeIn(repair_type))

        manimDNA.delete_DNA(query1[-3:] + query2[:4], self)
        self.play([manim.FadeIn(base) for base in query1[-3:] + query2[:4]])

        self.play(manim.FadeOut(repair_type))
        repair_type = (
            manim.Text("stagger")
            .move_to(manimDNA.number_plane.c2p(-15, 8, 0))
            .scale(0.5)
        )
        self.play(manim.FadeIn(repair_type))

        manimDNA.move_DNA(query2, manim.DOWN, 1, self)
        query2 = manimDNA.display_DNA("C", -1, -1, self) + query2

        manimDNA.move_DNA(query2, manim.DOWN, 1, self)

        self.play(
            manim.FadeIn(
                manim.Text("random insertion")
                .move_to(manimDNA.number_plane.c2p(-1, -1, 0))
                .scale(0.5)
            )
        )


class DOUBLECLEAVAGE(manim.Scene):
    def construct(self):
        manimDNA = ManimDNA()
        watson = manimDNA.generate_random_DNA(30)
        watson = watson[:4] + "CC" + watson[6:24] + "GG" + watson[26:]
        crick = Seq.Seq(watson).complement().__str__()

        watson, crick = manimDNA.display_DNA(
            [watson, crick], [-15, -15], [1, 0], self, color=[manim.WHITE] * 2
        )
        manimDNA.change_color(
            watson[23:26] + crick[4:7], manim.ManimColor("#FF0000"), self
        )

        manimDNA.move_DNA(watson[:10] + crick[:10], manim.LEFT, 1, self)
        manimDNA.move_DNA(watson[20:] + crick[20:], manim.RIGHT, 1, self)
        manimDNA.delete_DNA(watson[10:20] + crick, self)

        manimDNA.move_DNA(watson[:10], manim.RIGHT, 6, self)
        manimDNA.move_DNA(watson[20:], manim.LEFT, 6, self)

        query1 = "".join([mobj.text for mobj in watson[:8]]) + " " + watson[9].text
        query2 = "".join([mobj.text for mobj in watson[20:]])
        query1, query2 = manimDNA.display_DNA(
            [query1, query2], [-10, 0], [0, 0], self, color=[manim.WHITE] * 2
        )

        manimDNA.move_DNA(watson[:10] + query1, manim.LEFT, 5, self)
        manimDNA.move_DNA(watson[20:] + query2, manim.RIGHT, 5, self)
        self.play([manim.FadeIn(mobj) for mobj in watson[10:20]])

        manimDNA.move_DNA(query1[-1], manim.DOWN, 3, self)
        self.play(
            manim.FadeIn(
                manim.Text("random insertion")
                .move_to(manimDNA.number_plane.c2p(0, -3, 0))
                .scale(0.5)
            )
        )


class ScoreCrossCleavage(manim.Scene):
    def construct(self):
        manimDNA = ManimDNA()
        watson = manimDNA.generate_random_DNA(20)
        watson = watson[:14] + "GG" + watson[16:]

        watson = manimDNA.display_DNA(watson, -10, 1, self, color=manim.WHITE)
        manimDNA.change_color(watson[13:16], manim.ManimColor("#FF0000"), self)

        query11 = "".join([mobj.text for mobj in watson[:4]])
        query12 = "".join([mobj.text for mobj in watson[5:10]])
        query2 = "".join([mobj.text for mobj in watson[10:]])

        query11, query12, query2 = manimDNA.display_DNA(
            [query11, query12, query2], [-10, -5, 0], [0, 0, 0], self, [manim.WHITE] * 3
        )
        manimDNA.move_DNA(query12 + query2, manim.DOWN, 1, self)
        manimDNA.move_DNA(query12, manim.UP, 1, self)


class MicroHomology(manim.Scene):
    def construct(self):
        manimDNA = ManimDNA()
        watson = manimDNA.generate_random_DNA(20)
        watson = watson[:10] + watson[7:10] + watson[13] + "GG" + watson[16:]

        watson = manimDNA.display_DNA(watson, -10, 1, self, color=manim.WHITE)
        manimDNA.change_color(watson[13:16], manim.ManimColor("#FF0000"), self)

        query1 = "".join([mobj.text for mobj in watson[:7]])
        queryM = "".join([mobj.text for mobj in watson[7:10]])
        query2 = "".join([mobj.text for mobj in watson[13:]])

        query1, queryM, query2 = manimDNA.display_DNA(
            [query1, queryM, query2],
            [-10, -3, 3],
            [0, 0, 0],
            self,
            color=[manim.WHITE, manim.ManimColor("#00FF00"), manim.WHITE],
        )
        for _ in range(2):
            for i in range(3):
                manimDNA.move_DNA(queryM[2 - i], manim.RIGHT, 3, self)
            for i in range(3):
                manimDNA.move_DNA(queryM[i], manim.LEFT, 3, self)


class MainScene(manim.Scene):
    def construct(self):
        stagger = Stagger()
        stagger.construct()
        mmej = MMEJ()
        mmej.construct()
