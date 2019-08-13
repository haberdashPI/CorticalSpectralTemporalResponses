# CorticalSpectralTemporalResponses

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Build Status](https://travis-ci.org/haberdashPI/CorticalSpectralTemporalResponses.svg?branch=master)](https://travis-ci.org/haberdashPI/CorticalSpectralTemporalResponses)
[![codecov](https://codecov.io/gh/haberdashPI/CorticalSpectralTemporalResponses/branch/master/graph/badge.svg)](https://codecov.io/gh/haberdashPI/CorticalSpectralTemporalResponses)

A somewhat minimal implementation of the auditory spectrogram (`Audiospect`)
and cortical model (`cortical`, `scalefitler`, `ratefilter`) as described in the following paper:

[Chi, T., Ru, P., & Shamma, S. A. (2005). Multiresolution spectrotemporal
analysis of complex sounds. The Journal of the Acoustical Society of America,
118(2), 887–906.](http://doi.org/10.1121/1.1945807)

Refer to the documentation of each function (ala ? at the REPL) for more information.

You can find the original MATLAB implementation of these models
[here](https://isr.umd.edu/Labs/NSL/Software.htm).

## Status

These functions are quite stable; I'm using them extensively in another project.
I have not put much work into documenting the interface yet, but it is relatively
straightforward to use. Take a look at `test/runtests.jl` for examples.

## TODO

- Document function interfaces (double-check the examples)
- And a few example uses, and usage with `PlotAxes` for visualization.
