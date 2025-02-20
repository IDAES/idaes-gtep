# Generation and Transmission Expansion Planning with Reliability
# IDAES project
# author: Kyle Skolfield and Seolhee Cho
# date: 01/04/2024
# Model available at http://www.optimization-online.org/DB_FILE/2017/08/6162.pdf

import gtep_result_save


def reliability_data():

    d = {}

    if not gtep_result_save.results_sets:
        raise ValueError("solve_expansion_model must be run first")

    d.update(gtep_result_save.results_sets)

    return d
