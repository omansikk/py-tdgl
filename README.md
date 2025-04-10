# pyTDGL
Original version: [py-TDGL](https://github.com/loganbvh/py-tdgl/)

### Citing `pyTDGL`

`pyTDGL` is described in the following paper:

>*pyTDGL: Time-dependent Ginzburg-Landau in Python*, Computer Physics Communications **291**, 108799 (2023), DOI: [10.1016/j.cpc.2023.108799](https://doi.org/10.1016/j.cpc.2023.108799).

If you use `pyTDGL` in your research, please cite the paper linked above.

    % BibTeX citation
    @article{
        Bishop-Van_Horn2023-wr,
        title    = "{pyTDGL}: Time-dependent {Ginzburg-Landau} in Python",
        author   = "Bishop-Van Horn, Logan",
        journal  = "Comput. Phys. Commun.",
        volume   =  291,
        pages    = "108799",
        month    =  may,
        year     =  2023,
        url      = "http://dx.doi.org/10.1016/j.cpc.2023.108799",
        issn     = "0010-4655",
        doi      = "10.1016/j.cpc.2023.108799"
    }

### Acknowledgments

Parts of this package have been adapted from [`SuperDetectorPy`](https://github.com/afsa/super-detector-py), a GitHub repo authored by [Mattias Jönsson](https://github.com/afsa). Both `SuperDetectorPy` and `py-tdgl` are released under the open-source MIT License. If you use either package in an academic publication or similar, please consider citing the following in addition to the `pyTDGL` paper:

- Mattias Jönsson, Theory for superconducting few-photon detectors (Doctoral dissertation), KTH Royal Institute of Technology (2022) ([Link](http://urn.kb.se/resolve?urn=urn:nbn:se:kth:diva-312132))
- Mattias Jönsson, Robert Vedin, Samuel Gyger, James A. Sutton, Stephan Steinhauer, Val Zwiller, Mats Wallin, Jack Lidmar, Current crowding in nanoscale superconductors within the Ginzburg-Landau model, Phys. Rev. Applied 17, 064046 (2022) ([Link](https://journals.aps.org/prapplied/abstract/10.1103/PhysRevApplied.17.064046))

The user interface is adapted from [`SuperScreen`](https://github.com/loganbvh/superscreen).

# This fork
This fork of the py-TDGL package is slightly modified in order to allow thermoelectric quasiparticle currents $\mathbf{J} = \sigma S \nabla T$, and was used in the article [Thermoelectric AC Josephson effect](https://arXiv.org/). The figures 2-4 of the article can be reproduced with the scripts in the `scripts` folder.
