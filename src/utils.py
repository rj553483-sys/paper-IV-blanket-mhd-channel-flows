"""
Utility functions for Shercliff-Hartmann MHD boundary layer calculations
and pressure drop scaling in liquid metal divertor channels.
"""

import numpy as np

def calc_hartmann_number(B0, a, sigma, mu):
    """
    Calculate Hartmann number: Ha = B0 * a * sqrt(sigma / mu)
    """
    return B0 * a * np.sqrt(sigma / mu)

def calc_hartmann_boundary_layer(a, Ha):
    """
    Calculate Hartmann boundary layer thickness: delta_H = a / Ha [m]
    """
    return a / Ha

def calc_interaction_parameter(Ha, Re):
    """
    Calculate Interaction Parameter (Stuart number): N = Ha^2 / Re
    """
    return (Ha ** 2) / Re

def calc_mhd_pressure_drop(sigma, u0, B0, c_w=0.01):
    """
    Calculate MHD pressure drop in a rectangular channel [MPa/m]:
    dP/dz = - sigma * u0 * B0^2 * (c_w / (1 + c_w))
    """
    dp_dz_pascals = sigma * u0 * (B0 ** 2) * (c_w / (1.0 + c_w))
    return dp_dz_pascals / 1e6  # MPa/m
