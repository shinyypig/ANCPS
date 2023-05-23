# A New Translation Matching Method Based on Autocorrelated Normalized Cross-Power Spectrum

![license](https://img.shields.io/github/license/shinyypig/ancps)

This repository includes two test scripts, the code of the methods evaluated in the manuscript, and a sub-image of the HongKong dataset. It should be noted that all files were written in MATLAB R2020a.

Fulltext of the paper can be found [here](./ANCPS.pdf).

## Test scripts

- ANCPS_vs_NCPS.m: This script compares the phase of ANCPS and of NCPS under the noise-free and noisy situations. It will plot the figures of Fig.1 and Fig.2 in the manuscript.

- anti_noise_exp.m: This script runs the anti-noise experiment on the third sub-image in Fig.5 in the manuscript, and all parameters are set to the same as that in the III.A.2 section of the manuscript. This script will plot the three metrics: the mean value of e, the max value of e, and the standard deviation of e.

## Usage of ANCPS_itr.m

The file ANCPS_itr.m locates in the methods folder, which includes the function of the proposed method in the manuscript.

> offset = ANCPS_itr(im1, im2, num)

Input

- im1: the reference image
- im2: the image to be matched
- num: the number of iterations

Output

- offset: a 2\*num matrix, which includes the estimation of each iteration

## BibTex

```bib
@article{zhu2021a,
	author = {Zhu, Liangliang and Geng, Xiurui},
	journal = {IEEE Transactions on Geoscience and Remote Sensing},
	number = {8},
	year = {2021},
	pages = {6956--6968},
	publisher = {IEEE},
	title = {A {New} {Translation} {Matching} {Method} {Based} on {Autocorrelated} {Normalized} {Cross}-{Power} {Spectrum}},
	volume = {59},
}
```
