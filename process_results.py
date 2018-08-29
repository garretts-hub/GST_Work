import sys
sys.path.append("../gst")
import pygsti
from pygsti.construction import std1Q_XYI
#from gst import *
from create_data import generate_model_data
import pickle

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

# Experimental data L=32
#
ds = pygsti.io.load_dataset("l32-data.dat")
results = pygsti.do_stdpractice_gst(ds, gs_target, prep_fiducials, meas_fiducials, germs, maxLengths, modes="TP,CPTP,Target")
estimated_gate_set = results.estimates['TP.Robust+'].gatesets['final iteration estimate']
Gi_SO = estimated_gate_set['Gi'].get_data()
print("Gi =")
print(Gi_SO)
Gx_SO = estimated_gate_set['Gx'].get_data()
print("Gx =")
print(Gx_SO)
pickle.dump(Gi_SO, open('GiL32_TP.pkl','wb'))
pickle.dump(Gx_SO, open('GxL32_TP.pkl','wb'))
pygsti.report.create_standard_report(results, filename="./L32-ion", title="L=32 (ion-data)", verbosity=2)

# Experimental data L=16
#
#ds = pygsti.io.load_dataset("l16-data.dat")
#results = pygsti.do_stdpractice_gst(ds, gs_target, prep_fiducials, meas_fiducials, germs, maxLengths, modes="TP,CPTP,Target")
#pygsti.report.create_standard_report(results, filename="./L16-ion", title="L=16 (ion-data)", verbosity=5)

# Experimental data L=16 (day one)
#
#ds = pygsti.io.load_dataset("l16-day1-data.dat")
#results = pygsti.do_stdpractice_gst(ds, gs_target, prep_fiducials, meas_fiducials, germs, maxLengths, modes="TP,CPTP,Target")
#pygsti.report.create_standard_report(results, filename="./L16-day1-ion", title="L=16 (day1-ion-data)", verbosity=2)

# Experimental data L=16 (day two)
#
#ds = pygsti.io.load_dataset("l16-day2-data.dat")
#results = pygsti.do_stdpractice_gst(ds, gs_target, prep_fiducials, meas_fiducials, germs, maxLengths, modes="TP,CPTP,Target")
#pygsti.report.create_standard_report(results, filename="./L16-day2-ion", title="L=16 (day2-ion-data)", verbosity=2)

