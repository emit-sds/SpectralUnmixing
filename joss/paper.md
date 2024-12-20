---
title: 'SpectralUnmixing: A general Julia package for unmixing spectroscopy data'
tags:
  - Julia
  - imaging spectroscopy
  - EMIT
  - remote sensing
  - Earth science
authors:
  - name: Philip G. Brodrick
    orcid: 0000-0001-9497-7661
    affiliation: "1"
  - name: Francisco Ochoa
    orcid: 0000-0003-3363-2496
    affiliation: "2"
  - name: Gregory Okin
    orcid: 0000-0002-0484-3537
    affiliation: "2"
  - name: Vatsal A. Jhalani
    orcid: 0000-0003-0866-0858
    affiliation: "1"
  - name: Winston Olson-Duvall
    orcid: 0000-0002-4210-0283
    affiliation: "1"
  - name: Sarah R. Lundeen
    affiliation: "1"
  - name: David R. Thompson
    orcid: 0000-0003-1100-7550
    affiliation: "1"
  - name: Robert O. Green
    orcid: 0000-0001-9447-3076
    affiliation: "1"
affiliations:
 - name: Jet Propulsion Laboratory, California Institute of Technology, USA
   index: 1
 - name: University of California, Los Angeles, USA
   index: 2
date: 15 July 2022
bibliography: paper.bib
---

# Summary

On the Earth's surface, mixtures are the norm rather than the exception.
Taken to the limit, virtually every surface can be split into multiple constituents.
Spectral unmixing is the remote sensing retrieval process that attempts to quantitatively determine the relative fractions of various components that make up a surface based on optical data.
Imaging spectroscopy, in particular, has demonstrated the capacity for robust fractional retrievals across a wide range of domains, including mineralogical maps [@yan2010minerals; @combe2008analysis], urban land cover [@myint2009modelling], and vegetation [@okin2001practical].

Spectral unmixing is typically performed under the assumption of linear mixtures.
Some set of candidate 'endmembers' - constituents with known absolute (often pure) quantities and spectral signatures - are provided, and each pixel within an image is linearly unmixed with these endmembers to retrieve the relative contribution of each constituent.
In essence, this reduces to a simple linear algebra inversion where the known reference library can be inverted and multiplied by the observed reflectances to produce mixture fractions.
However, details arise regarding the nature of the selection of endmembers, with strategies ranging from dimensionality reduction of endmember 'classes' [@roberts1998mapping], bootstrapping [@asner2000biogeophysical], combinatorial selection [@roberts1998mapping; @franke2009hierarchical], and spectral brightness normalization [@asner2000biogeophysical].
The exact matrix inversion strategy to use is also an open and problem-specific decision, with candidates ranging from direct algebraic inversion to a constrained and regularized optimization [@hastie2015statistical].

# Statement of need

`SpectralUnmixing` is a one-stop-shop Julia package for all types of spectral unmixing strategies, focused on imaging spectroscopy data.
The code was designed for NASA's Earth Surface Mineral Dust Source Investigation (EMIT) [@emit2020] mission when, during the algorithm design phase, it became evident that there did not yet exist an optimized codebase that provided parallelized and flexible spectral unmixing strategies.
In particular, no central framework exists in which different unmixing strategies could be efficiently tested against one another.
`SpectralUnmixing` addresses this issue by drawing on a vast amount of existing literature and bringing the various proposed strategies together into a single codebase.
The framework also leverages a rich set of Julia packages for linear algebra, optimization, remote sensing data IO and more, resulting in a highly flexible and scalable unmixing package.
`SpectralUnmixing` has already supported the development of new unmixing strategies by accelerating the process of combining different components of the overall unmixing problem in novel ways [@OchoaQuantifying]
The package's scalability and breadth will allow it to continue providing this kind of coupled flexibility and operational capacity into the future.

![Example usage of `SpectralUnmixing`: fractional cover from an EMIT spectral image.  A) Red-green-blue (RGB) image of sample EMIT reflectance data observed near Sacramento, CA, USA with zoom-ins around labeled points (left).  B) Fractional cover output of `SpectralUnmixing` driver script `unmix.jl` on EMIT reflectance image with zoom-ins around labeled points (right). The RGB values in the fractional cover correspond to fractions of endmember library classes: non-photosynthetic vegetation (npv), photosynthetic vegetation (pv), and soil, respectively. C) EMIT reflectance spectra of sample points labeled in A) and B), chosen to each have high fractions of each of the 3 classes. Each spectra is colored by the RGB value corresponding to their class fraction. Here, `unmix.jl` was run using the included example endmember library and with the arguments `--normalization brightness --mode sma-best --n_mc 20 --num_endmembers 30`\label{fig:unmixing}](fig.pdf)

While `SpectralUnmixing` was originally created to be used operationally for the EMIT mission, it was also designed to be quickly adaptable for different researchers' needs.
The package currently supports the industry standard ENVI file format for unmixing raster maps, with support for more formats planned for the future.
The code is fully documented, including an example notebook, and features a script front-end that allows for arguments to be easily passed in and different options and datasets to be coupled together for rapid testing.
The package ultimately aims to benefit students, educators, professional researchers, and operational missions alike.

# Acknowledgements

EMIT is supported by the National Aeronautics and Space Administration Earth Venture Instrument program, USA, under the Earth Science Division of the Science Mission Directorate.
A portion of this research was carried out at the Jet Propulsion Laboratory, California Institute of Technology, under a contract with the National Aeronautics and Space Administration, USA (80NM0018D0004).
We acknowledge the support and assistance of NASA’s International Space Station Program, USA. Any use of trade, firm, or product names is for descriptive purposes only and does not imply endorsement by the U.S. Government. Copyright 2024 California Institute of Technology. All rights reserved. US Government Support Acknowledged.

# References