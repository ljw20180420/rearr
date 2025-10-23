#!/usr/bin/env python

from manim import *
import random
import math
import numpy as np

# config.pixel_width = 1920
# config.pixel_height = 1080
config.frame_rate = 30
config.save_sections = True

random.seed(1989)
ref_len = 40
rec_epsilon = 1e-6

def random_refs():
    ref = [random.choice('CGTA') for _ in range(ref_len)]
    ref[ref_len // 2 + 4: ref_len // 2 + 6] = ["G", "G"]
    RCmap = {"A": "T", "C": "G", "G": "C", "T": "A"}
    refRC = [RCmap[nuc] for nuc in ref]
    return ref, refRC

ref, refRC = random_refs()

def get_mh_matrix():
    mh_matrix = np.array(ref[len(ref)//2:]).reshape(-1, 1) == np.array(ref[:len(ref)//2]).reshape(1, -1)
    mh_matrix = mh_matrix.astype(int)
    for i in range(1, mh_matrix.shape[0]):
        for j in range(1, mh_matrix.shape[1]):
            if mh_matrix[i, j] > 0:
                mh_matrix[i, j] += mh_matrix[i - 1, j - 1]
                if i == mh_matrix.shape[0] - 1 or j == mh_matrix.shape[1] - 1:
                    for k in range(mh_matrix[i, j]):
                        mh_matrix[i - k, j - k] = mh_matrix[i, j]
            else:
                for k in range(mh_matrix[i - 1, j - 1]):
                    mh_matrix[i - 1 - k, j - 1 - k] = mh_matrix[i - 1, j - 1]
    return mh_matrix

mh_matrix = get_mh_matrix()

from manim import *

class NHEJ(Scene):
    class CrossTransform(Animation):
        def __init__(self, rectangle: Rectangle, new_end3: float, new_end5: float, one_or_two: int, **kwargs) -> None:
            super().__init__(rectangle,  **kwargs)
            self.one_or_two = one_or_two
            if self.one_or_two == 1:
                self.old_end3 = rectangle.get_right()[0]
                self.old_end5 = rectangle.get_left()[0]
            else:
                self.old_end3 = rectangle.get_top()[1]
                self.old_end5 = rectangle.get_bottom()[1]
            if rectangle.get_color() == RED:
                self.old_end3, self.old_end5 = self.old_end5, self.old_end3
            self.new_end3 = new_end3
            self.new_end5 = new_end5

        def interpolate_mobject(self, alpha: float) -> None:
            end3 = self.old_end3 + alpha * (self.new_end3 - self.old_end3)
            end5 = self.old_end5 + alpha * (self.new_end5 - self.old_end5)
            self.mobject.move_to([(end3 + end5) / 2 if self.one_or_two == 1 else -3.1, (end3 + end5) / 2 if self.one_or_two == 2 else 0, 0])
            self.mobject.set_color(YELLOW if abs(end3 - end5) < 2 * rec_epsilon else GREEN if end3 > end5 else RED)
            if self.one_or_two == 1:
                self.mobject.rescale_to_fit(length=max(abs(end3 - end5), rec_epsilon), dim=0, stretch=True)
            else:
                self.mobject.rescale_to_fit(length=max(abs(end3 - end5), rec_epsilon), dim=1, stretch=True)

    class NucletideColor(Animation):
        def __init__(self, nucMob: Text, old_end: float, new_end: float, one_or_two: int, **kwargs):
            super().__init__(nucMob,  **kwargs)
            self.one_or_two = one_or_two
            self.old_end = old_end
            self.new_end = new_end

        def interpolate_mobject(self, alpha: float) -> None:
            end = self.old_end + alpha * (self.new_end - self.old_end)
            if self.one_or_two == 1:
                if self.mobject.get_right()[0] > end:
                    self.mobject.set_color(BLACK)
                else:
                    self.mobject.set_color(WHITE)
            else:
                if self.mobject.get_top()[1] > end:
                    self.mobject.set_color(BLACK)
                else:
                    self.mobject.set_color(WHITE)


    def __init__(self):
        super().__init__()
        self.csize = 0.09 # nucletide size
        self.cdis_factor = 0.3 # nucletide distance factor
        self.buff_sz = 0.7 # nucletide buffer to boundary
        self.mh_thres = 99999 # initialize a large threshold to prevent displaying any micro-homology
        self.mh_markers = {} # empty micro-homology markers
        self.mh_marker_sz = 0.05 # micro-homology marker size
        self.top = 4
        self.bottom = -4
        self.left = -7
        self.right = 7

    def get_ends(self):
        end13 = self.rec1.get_right()[0]
        end15 = self.rec1.get_left()[0]
        if self.rec1.get_color() == RED:
            end13, end15 = end15, end13
        end23 = self.rec2.get_top()[1]
        end25 = self.rec2.get_bottom()[1]
        if self.rec2.get_color() == RED:
            end23, end25 = end25, end23
        
        return end13, end15, end23, end25

        

    def initial_ends(self):
        self.end15 = self.pos1thres + ref_len // 2 * self.pos1unit
        self.end13 = self.pos1thres + ref_len // 2 * self.pos1unit
        self.end25 = self.pos2thres + ref_len // 2 * self.pos2unit
        self.end23 = self.pos2thres + ref_len // 2 * self.pos2unit

        self.rec1 = Rectangle(
            color=YELLOW,
            height=8, width=rec_epsilon, fill_opacity=0.2
        ).move_to([(self.end13 + self.end15) / 2, 0, 0])

        self.rec2 = Rectangle(
            color=YELLOW,
            height=rec_epsilon, width=8, fill_opacity=0.2
        ).move_to([-3.1, (self.end23 + self.end25) / 2, 0])

        def update_Unilateral1(mob):
            end13, end15, end23, end25 = self.get_ends()
            mob.move_to([end15, end23, 0])
            mob.set_color(BLUE if end13 > end15 and end23 > end25 else GREY)

        def update_Unilateral2(mob):
            end13, end15, end23, end25 = self.get_ends()
            mob.move_to([end13, end25, 0])
            mob.set_color(BLUE if end13 > end15 and end23 > end25 else GREY)

        def update_NHEJ_name(mob, target, direction):
            mob.next_to(target, direction * 0.2)
            mob.set_color(target.get_color())

        self.NHEJ = Circle(0.05, color=BLUE, fill_opacity=1).move_to([max(self.end13, self.end15), max(self.end23, self.end25), 0])
        self.NHEJ.add_updater(lambda mob: mob.move_to([self.rec1.get_right()[0], self.rec2.get_top()[1], 0]))
        self.NHEJ_name = Text("NHEJ", font="Arial", font_size=10, color=BLUE)
        self.NHEJ_name.add_updater(lambda mob: update_NHEJ_name(mob, self.NHEJ, RIGHT))

        self.Unilateral1 = Circle(0.05, color=GREY, fill_opacity=1).move_to([self.end15, self.end23, 0])
        self.Unilateral1.add_updater(update_Unilateral1)
        self.Unilateral1_name = Text("Unilateral1", font="Arial", font_size=10, color=GREY)
        self.Unilateral1_name.add_updater(lambda mob: update_NHEJ_name(mob, self.Unilateral1, DOWN))

        self.Unilateral2 = Circle(0.05, color=GREY, fill_opacity=1).move_to([self.end13, self.end25, 0])
        self.Unilateral2.add_updater(update_Unilateral2)
        self.Unilateral2_name = Text("Unilateral2", font="Arial", font_size=10, color=GREY)
        self.Unilateral2_name.add_updater(lambda mob: update_NHEJ_name(mob, self.Unilateral2, UP))

        self.play(
            FadeIn(self.rec1),
            FadeIn(self.rec2),
            *[FadeToColor(self.nucMobs[i], color=BLACK) for i in range(ref_len // 2, ref_len)],
            *[FadeToColor(self.nucRcMobs[i], color=BLACK) for i in range(ref_len // 2, ref_len)],
            *[FadeToColor(self.nucMobs2[i],  color=BLACK) for i in range(ref_len // 2)],
            *[FadeToColor(self.nucRcMobs2[i],  color=BLACK) for i in range(ref_len // 2)],
            FadeIn(self.NHEJ),
            FadeIn(self.Unilateral1),
            FadeIn(self.Unilateral2),
            FadeIn(self.NHEJ_name),
            FadeIn(self.Unilateral1_name),
            FadeIn(self.Unilateral2_name),
        )

    def get_idx1(self, end1):
        return round((end1 - self.pos1thres) / self.pos1unit)
    
    def get_idx2(self, end2):
        return round((end2 - self.pos2thres) / self.pos2unit)

    def update_ends(self, new_end13, new_end15, new_end23, new_end25, run_time=1, wait_time=1):
        old_end13, old_end15, old_end23, old_end25 = self.get_ends()

        self.play(
            self.CrossTransform(self.rec1, new_end13, new_end15, 1),
            self.CrossTransform(self.rec2, new_end23, new_end25, 2),
            *[self.NucletideColor(nucMob, old_end13, new_end13, 1) for nucMob in self.nucMobs],
            *[self.NucletideColor(nucMob, old_end15, new_end15, 1) for nucMob in self.nucRcMobs],
            *[self.NucletideColor(nucMob, old_end23, new_end23, 2) for nucMob in self.nucRcMobs2],
            *[self.NucletideColor(nucMob, old_end25, new_end25, 2) for nucMob in self.nucMobs2],
            run_time=run_time,
            rate_func=linear
        )
        self.wait(wait_time)

    def draw_horizon_ref(self):
        self.nucMobs = [Text(nuc, font="Monospace", height=self.csize, width=self.csize) for nuc in ref]
        self.nucRcMobs = [Text(nuc, font="Monospace", height=self.csize, width=self.csize) for nuc in refRC]
        self.nucMobs[0].to_edge(UP, buff=0.2)
        self.nucMobs[0].to_edge(LEFT, buff=self.buff_sz)
        self.nucRcMobs[0].next_to(self.nucMobs[0], DOWN * self.cdis_factor)
        self.add(self.nucMobs[0])
        self.add(self.nucRcMobs[0])
        for i in range(len(self.nucMobs) - 1):
            self.nucMobs[i + 1].next_to(self.nucMobs[i], RIGHT * self.cdis_factor)
            self.nucRcMobs[i + 1].next_to(self.nucRcMobs[i], RIGHT * self.cdis_factor)

        self.play(
            *[FadeIn(nucMob) for nucMob in self.nucMobs],
            *[FadeIn(nucRcMob) for nucRcMob in self.nucRcMobs]
        )

        self.pos1unit = self.nucMobs[1].get_center()[0] - self.nucMobs[0].get_center()[0]
        self.pos1thres = self.nucMobs[0].get_center()[0] - self.pos1unit / 2

    def draw_vertial_ref(self):
        self.nucMobs2 = [nucMob.copy() for nucMob in self.nucMobs]
        self.nucRcMobs2 = [nucMob.copy() for nucMob in self.nucRcMobs]
        for i in range(len(self.nucRcMobs2)):
            self.nucRcMobs2[i].generate_target()
            self.nucRcMobs2[i].target.rotate(-math.pi/2)
            if i == 0:
                self.nucRcMobs2[i].target.to_edge(LEFT, buff=0.2).to_edge(UP, buff=self.buff_sz)
            else:
                self.nucRcMobs2[i].target.next_to(self.nucRcMobs2[i - 1].target, DOWN * self.cdis_factor)

        for i in range(len(self.nucMobs2)):
            self.nucMobs2[i].generate_target()
            self.nucMobs2[i].target.rotate(-math.pi/2)
            if i == 0:
                self.nucMobs2[i].target.next_to(self.nucRcMobs2[i].target, RIGHT * self.cdis_factor)
            else:
                self.nucMobs2[i].target.next_to(self.nucMobs2[i - 1].target, DOWN * self.cdis_factor)

        self.play(
            *[MoveToTarget(nucMob) for nucMob in self.nucMobs2],
            *[MoveToTarget(nucMob) for nucMob in self.nucRcMobs2]
        )

        self.pos2unit = self.nucMobs2[1].get_center()[1] - self.nucMobs2[0].get_center()[1]
        self.pos2thres = self.nucMobs2[0].get_center()[1] - self.pos2unit / 2

    def add_5p_3p(self):
        self.play(
            FadeIn(Text("5'", font="Monospace", height=self.csize, width=self.csize).next_to(self.nucMobs[0], LEFT * self.cdis_factor)),
            FadeIn(Text("3'", font="Monospace", height=self.csize, width=self.csize).next_to(self.nucRcMobs[0], LEFT * self.cdis_factor)),
            FadeIn(Text("5'", font="Monospace", height=self.csize, width=self.csize).rotate(-math.pi/2).next_to(self.nucMobs2[0], UP * self.cdis_factor)),
            FadeIn(Text("3'", font="Monospace", height=self.csize, width=self.csize).rotate(-math.pi/2).next_to(self.nucRcMobs2[0], UP * self.cdis_factor)),
            FadeIn(Text("3'", font="Monospace", height=self.csize, width=self.csize).next_to(self.nucMobs[-1], RIGHT * self.cdis_factor)),
            FadeIn(Text("5'", font="Monospace", height=self.csize, width=self.csize).next_to(self.nucRcMobs[-1], RIGHT * self.cdis_factor)),
            FadeIn(Text("3'", font="Monospace", height=self.csize, width=self.csize).rotate(-math.pi/2).next_to(self.nucMobs2[-1], DOWN * self.cdis_factor)),
            FadeIn(Text("5'", font="Monospace", height=self.csize, width=self.csize).rotate(-math.pi/2).next_to(self.nucRcMobs2[-1], DOWN * self.cdis_factor)),
        )

    def add_cut(self):
        cut1 = (self.nucMobs[ref_len // 2 - 1].get_center()[0] + self.nucMobs[ref_len // 2].get_center()[0]) / 2
        top1 = self.nucMobs[0].get_top()[1]
        bottom1 = self.nucRcMobs[0].get_bottom()[1]
        cut2 = (self.nucMobs2[ref_len // 2 - 1].get_center()[1] + self.nucMobs2[ref_len // 2].get_center()[1]) / 2
        left2 = self.nucRcMobs2[0].get_left()[0]
        right2 = self.nucMobs2[0].get_right()[0]
        cut1line = Line(start = [cut1, top1, 0], end = [cut1, bottom1, 0], color = ORANGE)
        cut2line = Line(start = [left2, cut2, 0], end = [right2, cut2, 0], color = ORANGE)
        PAM1 = Rectangle(color = ORANGE, height = top1 - bottom1, width = self.pos1unit * 3).move_to([cut1 + 4.5 * self.pos1unit, (top1 + bottom1) / 2, 0])
        PAM2 = Rectangle(color = ORANGE, width = right2 - left2, height = -self.pos2unit * 3).move_to([(left2 + right2) / 2, cut2 + 4.5 * self.pos2unit, 0])
        self.play(
            FadeIn(cut1line),
            FadeIn(cut2line),
            FadeIn(PAM1),
            FadeIn(PAM2)
        )

    def update_micro_homology(self, new_mh_thres):
        if self.mh_thres > new_mh_thres:
            fade_ins = []
            for i in range(1, mh_matrix.shape[0]):
                for j in range(1, mh_matrix.shape[1]):
                    if mh_matrix[i, j] >= self.mh_thres or mh_matrix[i, j] < new_mh_thres:
                        continue
                    if (i, j) not in self.mh_markers:
                        self.mh_markers[(i, j)] = Circle(radius=self.mh_marker_sz, fill_opacity=1).align_to(self.nucMobs2[mh_matrix.shape[1] + i], UP).align_to(self.nucMobs[j], RIGHT)
                    fade_ins.append(FadeIn(self.mh_markers[(i, j)]))
            self.play(*fade_ins)
        else:
            fade_outs = []
            for i in range(1, mh_matrix.shape[0]):
                for j in range(1, mh_matrix.shape[1]):
                    if mh_matrix[i, j] < self.mh_thres or mh_matrix[i, j] >= new_mh_thres:
                        continue
                    fade_outs.append(FadeOut(self.mh_markers[(i, j)]))
            self.play(*fade_outs)
        
        self.mh_thres = new_mh_thres

    def get_pos1(self, pos):
        return self.pos1thres + pos * self.pos1unit
    
    def get_pos2(self, pos):
        return self.pos2thres + pos * self.pos2unit

    def explain_NHEJ(self):
        self.next_section("NHEJ")
        self.play(FadeOut(self.title))
        self.remove(self.title)
        self.title = Text("NHEJ", font="Arial", font_size=30).move_to([4, 0, 0])
        self.play(FadeIn(self.title))

        for stagger in range(1, 4):
            self.next_section(f"NHEJ stagger {stagger}")
            self.update_ends(
                new_end13 = self.get_pos1(ref_len // 2 - stagger),
                new_end15 = self.get_pos1(ref_len // 2),
                new_end23 = self.get_pos2(ref_len // 2),
                new_end25 = self.get_pos2(ref_len // 2 - stagger),
            )

    def explain_Unilateral(self):
        self.next_section("Unilateral")
        self.play(FadeOut(self.title))
        self.remove(self.title)
        self.title = Text("Unilateral", font="Arial", font_size=30).move_to([4, 0, 0])
        self.play(FadeIn(self.title))

        self.update_ends(
            new_end13 = self.get_pos1(ref_len // 2),
            new_end15 = self.get_pos1(ref_len // 2),
            new_end23 = self.get_pos2(ref_len // 2),
            new_end25 = self.get_pos2(ref_len // 2),
        )

        self.next_section("Unilateral_start")
        del_len = 18
        self.update_ends(
            new_end13 = self.get_pos1(ref_len // 2),
            new_end15 = self.get_pos1(ref_len // 2 - del_len),
            new_end23 = self.get_pos2(ref_len // 2),
            new_end25 = self.get_pos2(ref_len // 2 + del_len),
            run_time = 5
        )

    def explain_MMEJ(self):
        self.next_section("MMEJ")
        self.play(FadeOut(self.title))
        self.remove(self.title)
        self.title = Text("MMEJ", font="Arial", font_size=30).move_to([4, 0, 0])
        self.play(FadeIn(self.title))

    def construct(self):
        self.title = Text("CRISPR", font="Arial", font_size=30).move_to([4, 0, 0])
        self.play(FadeIn(self.title))
        self.draw_horizon_ref()
        self.draw_vertial_ref()
        self.add_5p_3p()
        self.add_cut()        
        self.initial_ends()

        self.explain_NHEJ()

        self.explain_Unilateral()

        self.next_section("micro-homology")
        self.update_micro_homology(3)

        self.explain_MMEJ()
