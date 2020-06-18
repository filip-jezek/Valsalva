# test FMPy optim
import fmpy

fmu = 'ADAN_0main_SystemicTree_Tilt_OlufsenTriSeg_0tiltable.fmu'

unzipdir = fmpy.extract(fmu)
# read the model description
model_description = fmpy.read_model_description(unzipdir)
# instantiate the FMU
fmu_instance = fmpy.instantiate_fmu(unzipdir, model_description, 'CoSimulation', fmi_call_logger=lambda s: print(s), debug_logging=True)

# the python breaks at fmu.terminate()
result = fmpy.simulate_fmu(unzipdir, fmu_instance=fmu_instance, stop_time=1, relative_tolerance=1e-6, output_interval=0.02)

# it wont come to this point without any trace
print(result)
fmpy.write_csv('test', result)
pass