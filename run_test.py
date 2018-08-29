import sys

import pygsti
from pygsti.construction import std1Q_XYI
from create_data import generate_model_data, generate_model_data_rotation_error

import matplotlib.pyplot as plt

import NoiseSignal as _ns

maxLengths = [1,2,4,8,16,32]
gs_target = std1Q_XYI.gs_target
prep_fiducials = std1Q_XYI.prepStrs
meas_fiducials = std1Q_XYI.effectStrs

germ_strs = \
            [ 'Gi'
            , 'Gx'
            , 'Gy'
            , 'GxGy'
            , 'GxGyGi'
            , 'GxGiGy'
            , 'GxGiGi'
            , 'GyGiGi'
            , 'GxGxGiGy'
            , 'GxGyGyGi'
            , 'GxGxGyGxGyGy'
            ]
germs = pygsti.construction.gatestring_list(germ_strs)
listOfExperiments = pygsti.construction.make_lsgst_experiment_list(gs_target, prep_fiducials, meas_fiducials, germs, maxLengths)

"""
ds_filename = "amp_error_1.5em3_s20.txt"
ds = generate_model_data_rotation_error(listOfExperiments, 1000, 20, 20, 32, 'amplitude', 'fixed', 1.5e-3, 5e-4, 0.1, [1234])
results = pygsti.do_stdpractice_gst(ds, gs_target, prep_fiducials, meas_fiducials, germs, maxLengths, modes="TP,CPTP,Target")
pygsti.report.create_standard_report(results, filename="./L32-amp-error-1.5em3-s20", title="L=32 (amplitude error-1.5em3-s20)", verbosity=2)
pygsti.io.write_dataset(ds_filename, ds)

ds_filename = "amp_error_1.5em3_s1.txt"
ds = generate_model_data_rotation_error(listOfExperiments, 1000, 1, 20, 32, 'amplitude', 'fixed', 1.5e-3, 5e-4, 0.1, [1234])
results = pygsti.do_stdpractice_gst(ds, gs_target, prep_fiducials, meas_fiducials, germs, maxLengths, modes="TP,CPTP,Target")
pygsti.report.create_standard_report(results, filename="./L32-amp-error-1.5em3-s1", title="L=32 (amplitude error-1.5em3-s1)", verbosity=2)
pygsti.io.write_dataset(ds_filename, ds)
"""

ds_filename = "phase_damp_1em2_s20.txt"
ds = generate_model_data_rotation_error(listOfExperiments, 1000, 20, 20, 32, 'phase_damping', 'fixed', 1.0e-2, 5e-4, 0.1, [1234])
results = pygsti.do_stdpractice_gst(ds, gs_target, prep_fiducials, meas_fiducials, germs, maxLengths, modes="TP,CPTP,Target")
pygsti.report.create_standard_report(results, filename="./L32-phase-damp-error-1em2-s20", title="L=32 (phase damping error-1em2-s20)", verbosity=2)
pygsti.io.write_dataset(ds_filename, ds)
