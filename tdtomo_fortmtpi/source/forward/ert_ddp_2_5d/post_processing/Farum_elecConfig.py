#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# Program under the GNU General Public License <http://www.gnu.org/licenses/>.

"""
    FARUM ELECTRODE CONFIGURATIONS

    Split the x and y coordinates of the electrode ordered by quadrupole. This script is adapted to a specific 
    data-set acquired with time-lapse cross-borehole ERT during aquifer remediation.

    Last modification: July 2024 - Léa Lelimouzin

    :input:
        data_xx_yy: [array] coordinates of the electrodes per quadrupole
    :output: 
        [elec_configx, elec_configy]: [array] coordinates of the electrodes per quadrupole separated between x and y
        coordinates
"""


import numpy as np


def Farum_elecConfig(data_xx_yy):

    elecPos = data_xx_yy
    nrows = np.shape(elecPos)[0]
    elec_configx = np.zeros([nrows, 4])
    elec_configy = np.zeros([nrows, 4])

    for j in range(nrows):
        elec_configx[j, 0: 4] = elecPos[j, 0: 4] + 3  # Array of x coordinates: xA, xB, xM, xN. The +3 is a consequence of the C++ codes requirements to not have electrode at 0m.
        elec_configy[j, 0: 4] = -elecPos[j, 4: 8]     # Array of y coordinates: yA, yB, yM, yN

    return [elec_configx, elec_configy]
