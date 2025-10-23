import manim
import numpy as np

TIME_UNIT = 2


class ManimDNA:
    def __init__(self):
        self.rng = np.random.default_rng(63036)
        self.number_plane = manim.NumberPlane(
            x_range=[-18, 18],
            y_range=[-9, 9],
            x_length=12,
            y_length=6,
        )
        manim.Create(self.number_plane)

    def generate_random_DNA(self, length):
        # usage: generate_random_DNA of length
        return "".join(self.rng.choice(["A", "C", "G", "T"], length))

    def SNP_DNA(self, DNA, probability):
        DNA = list(DNA)
        for i in range(len(DNA)):
            if self.rng.random() < probability:
                DNA[i] = self.rng.choice(["A", "C", "G", "T"])
        return "".join(DNA)

    def indel_DNA(self, DNA, probability):
        inslens = self.rng.negative_binomial(1, 1 - probability, len(DNA) + 1)
        insDNAs = self.generate_random_DNA(np.sum(inslens))
        DNAnew, start = [], 0
        for i in range(len(DNA)):
            end = start + inslens[i]
            DNAnew.append(insDNAs[start:end])
            if self.rng.random() > probability:
                DNAnew.append(DNA[i])
            start = end
        DNAnew.append(insDNAs[start:])
        return "".join(DNAnew)

    def display_DNA(
        self,
        DNA: str | list[str],
        x: float | list[float],
        y: float | list[float],
        scene: manim.Scene,
        color: manim.ManimColor | list[manim.ManimColor] = manim.WHITE,
    ) -> list[list[manim.Mobject]]:
        if isinstance(DNA, str):
            DNA = [DNA]
            x = [x]
            y = [y]
            color = [color]

        objs_list = []
        for dd, xx, yy, cc in zip(DNA, x, y, color):
            objs = []
            for i, base in enumerate(dd):
                objs.append(
                    manim.Text(base, color=cc)
                    .move_to(self.number_plane.c2p(xx + i, yy, 0))
                    .scale(3 / 5)
                )
            objs_list.append(objs)

        scene.play(
            [manim.FadeIn(obj) for objs in objs_list for obj in objs],
            run_time=TIME_UNIT,
        )
        if len(objs_list) == 1:
            return objs_list[0]
        return objs_list

    def move_DNA(
        self,
        DNA: list[manim.Text],
        direction: np.ndarray,
        distance: float,
        scene: manim.Scene,
    ) -> None:
        moved_DNA = []
        for base in DNA:
            x, y = self.number_plane.p2c(base.get_center())
            x += direction[0] * distance
            y += direction[1] * distance
            moved_DNA.append(base.animate.move_to(self.number_plane.c2p(x, y, 0)))
        scene.play(moved_DNA, run_time=TIME_UNIT)

    def delete_DNA(
        self,
        DNA: list[manim.Text],
        scene: manim.Scene,
        run_time: float = 1.0,
    ) -> None:
        scene.play([manim.FadeOut(base) for base in DNA], run_time=run_time * TIME_UNIT)

    def change_color(
        self,
        DNA: list[manim.Text],
        color: manim.ManimColor,
        scene: manim.Scene,
    ) -> None:
        scene.play(
            [manim.FadeToColor(base, color=color) for base in DNA], run_time=TIME_UNIT
        )
