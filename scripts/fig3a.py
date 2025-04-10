import numpy as np
import matplotlib.pyplot as plt
import tdgl, time, os, h5py
from scipy.signal import lombscargle
from create_device import create_device

film_radius = 3e3; T2 = Tmax = 0.98; d = 40
y = np.linspace(1.01 * film_radius, 0.79 * film_radius, 201)
x = np.zeros_like(y)
cross_section = np.array([x, y]).T
device, tau0 = create_device(film_radius, Tmax, d)
freqs = np.linspace(0.1, 40 * np.pi, 1000)
temps = np.linspace(0.70, 0.97, 55)
start = time.time()
for T1 in temps:
    tmax = round((4. / tau0)); dt = 1e-3; save = 50
    options = tdgl.SolverOptions(
        output_file = "ring-temp_sweep.h5",
        solve_time = tmax,
        skip_time = 2 * tmax,
        field_units = "uT",
        current_units = "uA",
        dt_init = 1e-4,
        dt_max = dt,
        save_every = save,
        adaptive = True,
    )

    def epsilon(r):
        x, y = r
        if y < 0:
            return 1. / T1 - 1.
        else:
            return 1. / T2 - 1.

    # ~ os.system("rm ring-temp_sweep.h5") # linux
    os.system("del ring-temp_sweep.h5") # windows

    solution = tdgl.solve(
        device,
        options,
        disorder_epsilon = epsilon,
    )

    times, jn  = tdgl.get_current_through_paths("./ring-temp_sweep.h5", cross_section, dataset = "normal_current", interp_method = 'linear', units = "uA", with_units = False, progress_bar = True)
    times, js  = tdgl.get_current_through_paths("./ring-temp_sweep.h5", cross_section, dataset = "supercurrent", interp_method = 'linear', units = "uA", with_units = False, progress_bar = True)
    total = jn + js
    times = times * tau0
    jn_max = np.max(jn)
    js_max = np.max(js)
    total_max = np.max(total)
    jn_min = np.min(jn)
    js_min = np.min(js)
    total_min = np.min(total)

    np.savetxt("./currents/js_" + f"{T1:2.2f}" , js)
    np.savetxt("./currents/jn_" + f"{T1:2.2f}" , jn)

    
    psd = lombscargle(times, jn + js, freqs, True, True)

    np.savetxt("./psd/psd_" + f"{T1:2.3f}" , psd)

xi, lam, eta, gam, tau0 = pars()
freq = np.linspace(0.01, 20, len(psd))
style = ":"
dT = np.linspace(0.03, 0.28, 55)
for i, t in enumerate(range(930, 700, -5)): 
    psd = np.loadtxt("./psd/psd_" + str(t))
    pic[:, i] = psd

ax.imshow(pic, origin = 'lower', extent = [0.0275, 0.2825, 0.0075, 20.25], aspect = "auto", interpolation = 'none', vmax = 0.2)
ax.plot(dT, eta * dT / (2. * np.pi * tau0), c = "w", linestyle = style)
ax.plot(dT, 2 * eta * dT / (2. * np.pi * tau0), c = (0.8, 0.8, 0.8), linestyle = style)
ax.plot(dT, 3 * eta * dT / (2. * np.pi * tau0), c = (0.6, 0.6, 0.6), linestyle = style)
ax.plot(dT, 4 * eta * dT / (2. * np.pi * tau0), c = (0.4, 0.4, 0.4), linestyle = style)
ax.plot(dT, 5 * eta * dT / (2. * np.pi * tau0), c = (0.2, 0.2, 0.2), linestyle = style)
plt.savefig("fig2a.svg")

end = time.time()
print((end - start) / 60.)

