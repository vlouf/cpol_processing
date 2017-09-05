"""
Codes for correcting and unfolding PHIDP as well as estimating KDP.

@title: phase_codes
@author: Valentin Louf <valentin.louf@monash.edu>
@institutions: Monash University and the Australian Bureau of Meteorology
@date: 15/08/2017

.. autosummary::
    :toctree: generated/

    phidp_giangrande
"""
# Other Libraries
import pyart
import numpy as np


def phidp_giangrande(radar, refl_field='DBZ', ncp_field='NCP',
                     rhv_field='RHOHV', phidp_field='PHIDP'):
    """
    Phase processing using the LP method in Py-ART. A LP solver is required,

    Parameters:
    ===========
        radar:
            Py-ART radar structure.
        refl_field: str
            Reflectivity field name.
        ncp_field: str
            Normalised coherent power field name.
        rhv_field: str
            Cross correlation ration field name.
        phidp_field: str
            Phase field name.

    Returns:
    ========
        phidp_gg: dict
            Field dictionary containing processed differential phase shifts.
        kdp_gg: dict
            Field dictionary containing recalculated differential phases.
    """
    phi = radar.fields[phidp_field]['data'].copy()
    # Unfolding phidp
    phi = np.ma.masked_where(phi > 0, phi)
    radar.add_field_like('PHIDP', "PHI_CORR", 360 + phi, replace_existing=True)

    # Processing PHIDP
    phidp_gg, kdp_gg = pyart.correct.phase_proc_lp(radar, 0.0,
                                                   LP_solver='cylp',
                                                   refl_field=refl_field,
                                                   ncp_field=ncp_field,
                                                   rhv_field=rhv_field,
                                                   phidp_field="PHI_CORR")

    # Removing tmp field
    radar.fields.pop("PHI_CORR")

    return phidp_gg, kdp_gg
