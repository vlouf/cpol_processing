CPOL Level 1b
=============

Code for producing the CPOL level 1b.

# Goals:
- A unique convention (CF/Radial)
- A unique file format (netcdf)

# Output parameters:
- name: DBZ
- name: VEL      # Unfolded velocity
- name: WIDTH
- name: ZDR
- name: PHIDP
- name: RHOHV
- name: sounding_temperature
- name: SNR
- name: KDP
- name: PHIDP_BRINGI
- name: KDP_BRINGI
- name: VEL_RAW      # Folded velocity
- name: specific_attenuation_zh
- name: specific_attenuation_zdr
- name: LWC       # Liquid water content
- name: IWC        # Ice water content

# Dependencies:
## Libraries
- [Py-ART][1]
- [Numpy][2]
- [Panda][3]

## Other requirements:

[1]: http://github.com/ARM-DOE/pyart
[2]: http://www.scipy.org/
[3]: http://pandas.pydata.org/
