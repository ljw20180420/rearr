#!/usr/bin/env python

import manim
from manim_dna import ManimDNA, TIME_UNIT
from Bio import Seq


class Stagger(manim.Scene):
    def construct(self, scene: manim.Scene = None):
        if scene is None:
            scene = self

        scene.next_section("Stagger")

        manimDNA = ManimDNA()
        watson = manimDNA.generate_random_DNA(20)
        watson = watson[:14] + "GG" + watson[16:]
        crick = Seq.Seq(watson).complement().__str__()

        watson, crick = manimDNA.display_DNA(
            [watson, crick], [-10, -10], [1, 0], scene, color=[manim.WHITE] * 2
        )
        manimDNA.change_color(watson[13:16], manim.ManimColor("#FF0000"), scene)
        manimDNA.move_DNA(watson[9:] + crick[10:], manim.RIGHT, 3, scene)
        _, crick_i = manimDNA.display_DNA(
            [watson[9].text, crick[9].text],
            [-1, 2],
            [1, 0],
            scene,
            [manim.ManimColor("#00FF00")] * 2,
        )
        manimDNA.move_DNA(watson[9:] + crick[10:] + crick_i, manim.LEFT, 2, scene)


class MMEJ(manim.Scene):
    def construct(self, scene: manim.Scene = None):
        if scene is None:
            scene = self

        scene.next_section("MMEJ")

        manimDNA = ManimDNA()
        watson = "TACGTGGTGCTCGGGGGGCT"
        crick = Seq.Seq(watson).complement().__str__()

        watson, crick = manimDNA.display_DNA(
            [watson, crick], [-10, -10], [1, 0], scene, color=[manim.WHITE] * 2
        )

        manimDNA.change_color(watson[13:16], manim.ManimColor("#FF0000"), scene)
        manimDNA.move_DNA(watson[10:] + crick[10:], manim.RIGHT, 3, scene)

        for i in range(7):
            manimDNA.delete_DNA([watson[10 + i], crick[9 - i]], scene, run_time=0.2)

        manimDNA.move_DNA(watson[-3:] + crick[10:], manim.DOWN, 1, scene)
        manimDNA.move_DNA(watson[-3:] + crick[10:], manim.LEFT, 10, scene)

        sticks = manimDNA.display_DNA("|| |", -5, 0, scene)
        manimDNA.delete_DNA([watson[9], crick[10]], scene, run_time=0.2)
        manimDNA.delete_DNA([watson[8], crick[11], sticks[-1]], scene, run_time=0.2)
        manimDNA.delete_DNA([watson[7]], scene, run_time=0.2)

        _, waston_sticks, crick_fills, crick_sticks = manimDNA.display_DNA(
            ["GGG", "|||", "CA", "||"],
            [-3, -3, -7, -7],
            [1, 0, -1, 0],
            scene,
            [
                manim.ManimColor("#00FF00"),
                manim.WHITE,
                manim.ManimColor("#00FF00"),
                manim.WHITE,
            ],
        )

        manimDNA.delete_DNA(sticks[:2] + waston_sticks + crick_sticks, scene)
        manimDNA.move_DNA(watson[-3:] + crick[12:] + crick_fills, manim.UP, 1, scene)


class Unilateral(manim.Scene):
    def construct(self, scene: manim.Scene = None):
        if scene is None:
            scene = self

        scene.next_section("Unilateral")

        manimDNA = ManimDNA()
        watson = "TACGTGGTGCTCGGGGGGCT"
        crick = Seq.Seq(watson).complement().__str__()

        watson, crick = manimDNA.display_DNA(
            [watson, crick], [-10, -10], [1, 0], scene, color=[manim.WHITE] * 2
        )

        manimDNA.change_color(watson[13:16], manim.ManimColor("#FF0000"), scene)
        manimDNA.move_DNA(watson[10:] + crick[10:], manim.RIGHT, 3, scene)

        for i in range(7):
            manimDNA.delete_DNA([watson[10 + i], crick[9 - i]], scene, run_time=0.2)

        manimDNA.move_DNA(watson[-3:] + crick[10:], manim.DOWN, 1, scene)
        manimDNA.move_DNA(watson[-3:] + crick[10:], manim.LEFT, 10, scene)

        sticks = manimDNA.display_DNA("|| |", -5, 0, scene)

        manimDNA.move_DNA(watson[-3:] + crick[-3:], manim.UP, 1, scene)
        manimDNA.delete_DNA(crick[10:17] + sticks, scene)

        manimDNA.display_DNA("CACCACG", -7, 0, scene, manim.ManimColor("#00FF00"))


class Unilateral2(manim.Scene):
    def construct(self, scene: manim.Scene = None):
        if scene is None:
            scene = self

        scene.next_section("Unilateral2")

        manimDNA = ManimDNA()
        watson = "TACGTGGTGCTCGGGGGGCT"
        crick = Seq.Seq(watson).complement().__str__()

        watson, crick = manimDNA.display_DNA(
            [watson, crick], [-10, -10], [1, 0], scene, color=[manim.WHITE] * 2
        )

        manimDNA.change_color(watson[13:16], manim.ManimColor("#FF0000"), scene)
        manimDNA.move_DNA(watson[10:] + crick[10:], manim.RIGHT, 3, scene)

        for i in range(7):
            manimDNA.delete_DNA([watson[10 + i], crick[9 - i]], scene, run_time=0.2)

        manimDNA.move_DNA(watson[-3:] + crick[10:], manim.DOWN, 1, scene)
        manimDNA.move_DNA(watson[-3:] + crick[10:], manim.LEFT, 10, scene)

        sticks = manimDNA.display_DNA("|| |", -5, 0, scene)

        manimDNA.move_DNA(watson[:3] + crick[:3], manim.DOWN, 1, scene)
        manimDNA.delete_DNA(watson[3:10] + sticks, scene)

        manimDNA.display_DNA("TCGGGGG", -7, 0, scene, manim.ManimColor("#00FF00"))


class Chimeric(manim.Scene):
    def construct(self, scene: manim.Scene = None):
        if scene is None:
            scene = self

        scene.next_section("Chimeric")

        manimDNA = ManimDNA()

        watson = manimDNA.display_DNA(
            "TACGTGGTGCTCGGGGGGCT", -10, 1, scene, color=manim.WHITE
        )

        repair_type = (
            manim.Text("blunt end")
            .move_to(manimDNA.number_plane.c2p(-15, 8, 0))
            .scale(0.5)
        )
        scene.play(manim.FadeIn(repair_type), run_time=TIME_UNIT)

        query1, query2 = manimDNA.display_DNA(
            ["TACGTGGTGC", "TCGGGGGGCT"],
            [-10, 0],
            [0, 0],
            scene,
            color=[manim.WHITE] * 2,
        )

        scene.play(manim.FadeOut(repair_type), run_time=TIME_UNIT)
        repair_type = (
            manim.Text("unilateral")
            .move_to(manimDNA.number_plane.c2p(-15, 8, 0))
            .scale(0.5)
        )
        scene.play(manim.FadeIn(repair_type), run_time=TIME_UNIT)

        manimDNA.delete_DNA(query1[3:], scene)
        scene.play([manim.FadeIn(base) for base in query1[3:]], run_time=TIME_UNIT)
        manimDNA.delete_DNA(query2[:-3], scene)
        scene.play([manim.FadeIn(base) for base in query2[:-3]], run_time=TIME_UNIT)

        scene.play(manim.FadeOut(repair_type), run_time=TIME_UNIT)
        repair_type = (
            manim.Text("MMEJ").move_to(manimDNA.number_plane.c2p(-15, 8, 0)).scale(0.5)
        )
        scene.play(manim.FadeIn(repair_type), run_time=TIME_UNIT)

        manimDNA.delete_DNA(query1[-3:] + query2[:4], scene)
        scene.play(
            [manim.FadeIn(base) for base in query1[-3:] + query2[:4]],
            run_time=TIME_UNIT,
        )

        scene.play(manim.FadeOut(repair_type), run_time=TIME_UNIT)
        repair_type = (
            manim.Text("stagger")
            .move_to(manimDNA.number_plane.c2p(-15, 8, 0))
            .scale(0.5)
        )
        scene.play(manim.FadeIn(repair_type), run_time=TIME_UNIT)

        manimDNA.move_DNA(query2, manim.DOWN, 1, scene)
        query2 = manimDNA.display_DNA("C", -1, -1, scene) + query2

        manimDNA.move_DNA(query2, manim.DOWN, 1, scene)

        scene.play(
            manim.FadeIn(
                manim.Text("random insertion")
                .move_to(manimDNA.number_plane.c2p(-1, -1, 0))
                .scale(0.5)
            ),
            run_time=TIME_UNIT,
        )


class DoubleCleavage(manim.Scene):
    def construct(self, scene: manim.Scene = None):
        if scene is None:
            scene = self

        scene.next_section("DoubleCleavage")

        manimDNA = ManimDNA()
        watson = manimDNA.generate_random_DNA(30)
        watson = watson[:4] + "CC" + watson[6:24] + "GG" + watson[26:]
        crick = Seq.Seq(watson).complement().__str__()

        watson, crick = manimDNA.display_DNA(
            [watson, crick], [-15, -15], [1, 0], scene, color=[manim.WHITE] * 2
        )
        manimDNA.change_color(
            watson[23:26] + crick[4:7], manim.ManimColor("#FF0000"), scene
        )

        manimDNA.move_DNA(watson[:10] + crick[:10], manim.LEFT, 1, scene)
        manimDNA.move_DNA(watson[20:] + crick[20:], manim.RIGHT, 1, scene)
        manimDNA.delete_DNA(watson[10:20] + crick, scene)

        manimDNA.move_DNA(watson[:10], manim.RIGHT, 6, scene)
        manimDNA.move_DNA(watson[20:], manim.LEFT, 6, scene)

        query1 = "".join([mobj.text for mobj in watson[:8]]) + " " + watson[9].text
        query2 = "".join([mobj.text for mobj in watson[20:]])
        query1, query2 = manimDNA.display_DNA(
            [query1, query2], [-10, 0], [0, 0], scene, color=[manim.WHITE] * 2
        )

        manimDNA.move_DNA(watson[:10] + query1, manim.LEFT, 5, scene)
        manimDNA.move_DNA(watson[20:] + query2, manim.RIGHT, 5, scene)
        scene.play([manim.FadeIn(mobj) for mobj in watson[10:20]], run_time=TIME_UNIT)

        manimDNA.move_DNA(query1[-1], manim.DOWN, 3, scene)
        scene.play(
            manim.FadeIn(
                manim.Text("random insertion")
                .move_to(manimDNA.number_plane.c2p(0, -3, 0))
                .scale(0.5)
            ),
            run_time=TIME_UNIT,
        )


class ScoreCrossCleavage(manim.Scene):
    def construct(self, scene: manim.Scene = None):
        if scene is None:
            scene = self

        scene.next_section("ScoreCrossCleavage")

        manimDNA = ManimDNA()
        watson = manimDNA.generate_random_DNA(20)
        watson = watson[:14] + "GG" + watson[16:]

        watson = manimDNA.display_DNA(watson, -10, 1, scene, color=manim.WHITE)
        manimDNA.change_color(watson[13:16], manim.ManimColor("#FF0000"), scene)

        query11 = "".join([mobj.text for mobj in watson[:4]])
        query12 = "".join([mobj.text for mobj in watson[5:10]])
        query2 = "".join([mobj.text for mobj in watson[10:]])

        query11, query12, query2 = manimDNA.display_DNA(
            [query11, query12, query2],
            [-10, -5, 0],
            [0, 0, 0],
            scene,
            [manim.WHITE] * 3,
        )
        manimDNA.move_DNA(query12 + query2, manim.DOWN, 1, scene)
        manimDNA.move_DNA(query12, manim.UP, 1, scene)


class MicroHomology(manim.Scene):
    def construct(self, scene: manim.Scene = None):
        if scene is None:
            scene = self

        scene.next_section("MicroHomology")

        manimDNA = ManimDNA()
        watson = manimDNA.generate_random_DNA(20)
        watson = watson[:10] + watson[7:10] + watson[13] + "GG" + watson[16:]

        watson = manimDNA.display_DNA(watson, -10, 1, scene, color=manim.WHITE)
        manimDNA.change_color(watson[13:16], manim.ManimColor("#FF0000"), scene)

        query1 = "".join([mobj.text for mobj in watson[:7]])
        queryM = "".join([mobj.text for mobj in watson[7:10]])
        query2 = "".join([mobj.text for mobj in watson[13:]])

        query1, queryM, query2 = manimDNA.display_DNA(
            [query1, queryM, query2],
            [-10, -3, 3],
            [0, 0, 0],
            scene,
            color=[manim.WHITE, manim.ManimColor("#00FF00"), manim.WHITE],
        )
        for _ in range(2):
            for i in range(3):
                manimDNA.move_DNA(queryM[2 - i], manim.RIGHT, 3, scene)
            for i in range(3):
                manimDNA.move_DNA(queryM[i], manim.LEFT, 3, scene)


def fade_out(scene: manim.Scene):
    animations = []
    for mobject in scene.mobjects:
        animations.append(manim.FadeOut(mobject))
    scene.play(*animations, run_time=TIME_UNIT)


class MainScene(manim.Scene):
    def construct(self):
        manimDNA = ManimDNA()
        for cls_name in [
            "Stagger",
            "MMEJ",
            "Unilateral",
            "Unilateral2",
            "Chimeric",
            "DoubleCleavage",
            "ScoreCrossCleavage",
            "MicroHomology",
        ]:
            self.play(
                manim.FadeIn(
                    manim.Text(cls_name)
                    .move_to(manimDNA.number_plane.c2p(15, 8, 0))
                    .scale(0.5)
                ),
                run_time=TIME_UNIT,
            )
            obj = eval(cls_name)()
            obj.construct(self)
            fade_out(self)
