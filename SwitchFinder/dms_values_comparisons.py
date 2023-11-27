import numpy as np
from scipy.stats import mannwhitneyu

import sys


import SwitchFinder.classes as classes
import SwitchFinder.dotbracket_comparisons as dotbracket_comparisons
import SwitchFinder.utils as utils


def get_shape_values_for_constant_ACs(inp_fragment, shape_profile):
    inp_fragment.get_differential_ACs()
    constant_paired_values = shape_profile[inp_fragment.major_conf.constant_paired_ACs]
    constant_paired_values = constant_paired_values[constant_paired_values >= 0]
    constant_unpaired_values = shape_profile[inp_fragment.major_conf.constant_unpaired_ACs]
    constant_unpaired_values = constant_unpaired_values[constant_unpaired_values >= 0]
    return constant_paired_values, constant_unpaired_values


def get_shape_values_for_variable_ACs(inp_fragment, shape_profile):
    inp_fragment.get_differential_ACs()
    unpaired_major_values = shape_profile[inp_fragment.major_conf.differentially_unpaired_ACs]
    unpaired_major_values = unpaired_major_values[unpaired_major_values >= 0]
    unpaired_second_values = shape_profile[inp_fragment.second_conf.differentially_unpaired_ACs]
    unpaired_second_values = unpaired_second_values[unpaired_second_values >= 0]
    return unpaired_major_values, unpaired_second_values


def compare_shape_values_paired_unpaired_ACs(constant_paired_values, constant_unpaired_values,
                                             alternative = 'greater'):
    stat, pv = mannwhitneyu(constant_paired_values,
                            constant_unpaired_values,
                            alternative = alternative)
    return pv

