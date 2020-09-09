"""
Codes for correcting Doppler velocity.

@title: velocity
@author: Valentin Louf <valentin.louf@bom.gov.au>
@institutions: Monash University and the Australian Bureau of Meteorology
@creation: 11/12/2017
@date: 26/08/2020

.. autosummary::
    :toctree: generated/

    _check_nyquist_velocity
    unravel
"""
import pyart
import numpy as np


def _check_nyquist_velocity(radar, vel_name="VEL"):
    """
    Check if Nyquist velocity is present in the instrument parameters. If not,
    then it is created.
    """
    try:
        vnyq = radar.instrument_parameters["nyquist_velocity"]["data"]
        if vnyq is None:
            raise KeyError("Nyquist velocity does not exists.")
    except KeyError:
        vnyq = np.nanmax(radar.fields[vel_name]["data"])
        nray = len(radar.azimuth["data"])
        vnyq_array = np.array([vnyq] * nray, dtype=np.float32)
        nyquist_velocity = pyart.config.get_metadata("nyquist_velocity")
        nyquist_velocity["data"] = vnyq_array
        nyquist_velocity["_Least_significant_digit"] = 2
        radar.instrument_parameters["nyquist_velocity"] = nyquist_velocity

    return vnyq


def unravel(radar, gatefilter, vel_name="VEL", dbz_name="DBZ", nyquist=None):
    """
    Unfold Doppler velocity using Py-ART region based algorithm. Automatically
    searches for a folding-corrected velocity field.

    Parameters:
    ===========
    radar:
        Py-ART radar structure.
    gatefilter:
        Filter excluding non meteorological echoes.
    vel_name: str
        Name of the (original) Doppler velocity field.
    dbz_name: str
        Name of the reflecitivity field.
    nyquist: float
        Nyquist velocity co-interval.

    Returns:
    ========
    vel_meta: dict
        Unfolded Doppler velocity.
    """
    import unravel

    vnyq = _check_nyquist_velocity(radar, vel_name)
    if nyquist is None:
        if np.isscalar(vnyq):
            nyquist = vnyq

    unfvel = unravel.unravel_3D_pyart(
        radar, vel_name, dbz_name, gatefilter=gatefilter, alpha=0.8, nyquist_velocity=nyquist, strategy="long_range"
    )

    vel_meta = pyart.config.get_metadata("velocity")
    vel_meta["data"] = np.ma.masked_where(gatefilter.gate_excluded, unfvel).astype(np.float32)
    vel_meta["_Least_significant_digit"] = 2
    vel_meta["_FillValue"] = np.NaN
    vel_meta["comment"] = "UNRAVEL algorithm."
    vel_meta["long_name"] = "Doppler radial velocity of scatterers away from instrument"
    vel_meta["units"] = "m s-1"

    return vel_meta
