---
title: About
---

# Introduction

This is the github page for a chimeric alignment engine [`rearr`][rearr]. Although [`rearr`][rearr] can be used for any chimeric alignment job in theory, [`rearr`][rearr] is mainly used to analyze CRISPR editing data.

The CRISPR editing output is highly chimeric due to the following reasons.

- Interplay of stagger cleavage, MMEJ and unilateral.
- Random insertion with or without template.
- Unlisted/Unknown machanisms.

<iframe src="./assets/videos/MainScene.mp4" width="854" height="480"></iframe>
<!-- <video controls width="854" height="480">
<source src="/assets/videos/MainScene.mp4" type="video/mp4">
</video> -->

Most current tools rely on non-chimeric alignment engines. [`rearr`][rearr] implements a chimeric alignment engine by applying the Smith-Waterman algorithm twice.

<iframe src="./assets/videos/Forward.mp4" width="854" height="480"></iframe>
<!-- <video controls width="854" height="480">
<source src="/assets/videos/Forward.mp4" type="video/mp4">
</video> -->

[`rearr`][rearr] is very efficient in both time and space.

- [`rearr`][rearr] generalizes [Farrar's SIMD](https://doi.org/10.1093/bioinformatics/btl582) to the chimeric alignment.
- [`rearr`][rearr] implements a novel branch-and-bound backtrack process to get the optimal alignment in linear space. This method is more efficient in time compared with the famous linear-space implementation based on [divide and conquer](https://doi.org/10.1093/bioinformatics/4.1.11).

# Accessibility

See the project README.md [`rearr`][rearr].

[rearr]: https://github.com/ljw20180420/rearr
