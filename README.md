# IDAES-GTEP

<!-- BEGIN Status Badges -->
[![GitHub CI](https://github.com/IDAES/idaes-gtep/actions/workflows/test_pr_and_main.yml/badge.svg?branch=main&event=push)](https://github.com/IDAES/idaes-gtep/actions/workflows/test_pr_and_main.yml)
[![Documentation Status](https://readthedocs.org/projects/idaes-gtep/badge/?version=latest)](http://idaes-gtep.readthedocs.org/en/latest/)
<!-- END Status Badges -->

The IDAES Generation and Transmission Expansion Planning (GTEP) package provides a [Pyomo](https://github.com/Pyomo/pyomo)-based implementation of a modular, flexible, Generalized Disjunctive Programming (GDP) formulation for power infrastructure planning problems.  This formulation is designed with the following goals in mind:

- Abstract GTEP modeling away from any particular case study or fixed modeling assumptions (e.g., technologies, temporal resolution, spatial resolution, policy implications, etc.)
- Admit flexible decision sets and heterogeneous parameterization
- Allow high-level modeling options to be understood easily, chosen modularly, and changed rapidly

You can find the latest documentation [here](https://idaes-gtep.readthedocs.io).