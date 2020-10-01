#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from scipy.interpolate import RegularGridInterpolator

import uproot


rfile = uproot.open("sphenix3dmaprhophiz.root")
df = rfile['fieldmap'].pandas.df()

df = df.assign(**{
    'x':  df.rho*np.cos(df.phi),
    'y':  df.rho*np.sin(df.phi),
    'bx': df.brho*np.cos(df.phi),
    'by': df.brho*np.sin(df.phi),
    'b':  np.sqrt(df.brho**2 + df.bphi**2 + df.bz**2),
})

dfc = df[['x', 'y', 'z', 'bx', 'by', 'bz', 'b', 'brho', 'bphi', 'rho', 'phi']].copy()
dfc = dfc.round(3)
dfc.sort_values(['x', 'y', 'z'], inplace=True)

interp_bx = RegularGridInterpolator((x, y, z), dfc.bx.values.reshape(81, 81, 101))
interp_by = RegularGridInterpolator((x, y, z), dfc.by.values.reshape(81, 81, 101))
interp_bz = RegularGridInterpolator((x, y, z), dfc.bz.values.reshape(81, 81, 101))

p = np.linspace(0, 360, 37)
z = np.linspace(-100, 100, 101)
r = np.linspace(0, 80, 81)

pp, zz, rr = np.meshgrid(p, z, r)

coords = np.dstack(( rr.reshape(-1), zz.reshape(-1), pp.reshape(-1) )).reshape(-1,3)

df_new = pd.DataFrame(coords, columns=['r', 'z', 'phi'])

dfc_new = df_new.assign(**{
    'x': df_new.r*np.cos(df_new.phi/180*np.pi),
    'y': df_new.r*np.sin(df_new.phi/180*np.pi),
})

dfc_new = dfc_new.assign(**{
    'bx': interp_bx(np.dstack((dfc_new.x, dfc_new.y, dfc_new.z)).reshape(-1,3)),
    'by': interp_by(np.dstack((dfc_new.x, dfc_new.y, dfc_new.z)).reshape(-1,3)),
    'bz': interp_bz(np.dstack((dfc_new.x, dfc_new.y, dfc_new.z)).reshape(-1,3))
})

#dfc_new
