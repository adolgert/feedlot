Feedlot
=======

A continuous time discrete event simulation of a cattle feedlot, written in C++.

- [Installation](src/install.md)

- [Usage](src/usage.md)

This prototype is a demonstration project for rapid construction of complex epidemiological simulations.
It is a command-line program which produces timelines for spread of FMDV within
feedlots of up to 100,000 cattle, separated into separate pens.

 - Individual disease progress follows susceptible, latent, infected, and recovered states
   as described by Mardones, et al.
 - Spread of FMDV among individuals within a pen is modeled with an exact equivalent
   to frequency-dependent SEIR, with a user-specified hazard rate for infection.
 - Spread from a pen to neighboring pens is also frequency-dependent where there
   is fenceline contact.
 - As an example of more complicated types of disease spread, a pen rider
   can carry infection to neighboring pens.

If you have any questions, please do not hesitate to email me.
