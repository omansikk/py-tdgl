import numpy as np
import matplotlib.pyplot as plt
import tdgl, time, os, h5py
from scipy.signal import argrelextrema, lombscargle
from create_device import pars, create_device

film_radius = 3e3; T2 = Tmax = 0.98; d = 40
y = np.linspace(1.01 * film_radius, 0.79 * film_radius, 201)
x = np.zeros_like(y)
cross_section = np.array([x, y]).T
xi, lam, eta, gam, tau0 = pars()
device, tau0 = create_device(film_radius, Tmax, d)
device.layer.gamma = 0
nT, ndT = 24, 37
temps = np.linspace(0.75, 0.98, nT)
deltas = np.linspace(0.01, 0.1, ndT)
amplitudes = np.zeros((nT, ndT))
amplitudes2 = np.zeros((nT, ndT))
frequencies = np.zeros((nT, ndT))
freqs = np.linspace(0.001, 8 * np.pi, 500)
dt = 5e-2
save = 5
for i, Tmax in enumerate(temps):
    for j, dT in enumerate(deltas):
        freq = (eta ** 2 * (Tmax - T1) ** 2 - (1. / Tmax - 1.) * (1. / T1 - 1.)) / (eta * (Tmax - T1))
        tmax = round(10 / freq);
        def epsilon(r):
            x, y = r
            if y < 0:
                return 1. / (Tmax - dT) - 1.
            else:
                return 1. / Tmax - 1.

        # ~ os.system("rm ring.h5") # linux
        os.system("del ring.h5") # window
        options = tdgl.SolverOptions(
            output_file = "ring.h5",
            solve_time = tmax,
            skip_time = round(0.5 * tmax),
            field_units = "uT",
            current_units = "uA",
            dt_init = 1e-4,
            dt_max = dt,
            save_every = save,
            adaptive = True,
        )

        solution = tdgl.solve(
            device,
            options,
            disorder_epsilon = epsilon,
        )


        times, total  = tdgl.get_current_through_paths("./ring.h5", cross_section, interp_method = 'linear', units = "uA", with_units = False, progress_bar = False)
        
        psd = lombscargle(times * tau0, total, freqs, True, True)
        
        amplitudes2[i, j] = 0.5 * (np.max(total) - np.min(total))
        psd = np.clip(psd, 5e-2, 1.)
        peaks = argrelextrema(psd, np.greater, order = 10)[0]
        if len(peaks) > 0:
            amplitudes[i, j] = psd[peaks[0]]
            frequencies[i, j] = freqs[peaks[0]] / (2. * np.pi)
        
    itime = round((time.time() - start) / (i + 1))
    print(i, itime, (nT - i - 1) * itime)

np.savetxt("amps", amplitudes)
np.savetxt("amps_2", amplitudes2)
np.savetxt("freq", frequencies)

x = np.linspace(0.01, 0.1, 37)
y = np.linspace(0.75, 0.98, 24)
edge_x = 0.5 * (x[1] - x[0])
edge_y = 0.5 * (y[1] - y[0])
zf = zero_freq(eta, x)
fig, ax = plt.subplots(figsize = (18 / 2.54, 6 / 2.54))
c = ax.imshow(frequencies, origin = "lower", extent = [0.01 - edge_x, 0.1 + edge_x, 0.75 - edge_y, 0.98 + edge_y], aspect = "auto", interpolation = "none")
ax.set_ylim(0.75 - edge_y, 0.98 + edge_y)
ax.set_yticks([0.75, 0.8, 0.85, 0.9, 0.95])
ax.set_xlabel("Temperature difference $\Delta T/T_c$")
ax.set_ylabel("Temperature $T_2$")
fig.colorbar(c, ax = ax, label = "$\omega$ (Ghz)", ticks = [0, 3])
ax.plot(x, zf, c = "w")

end = time.time()
print((end - start) / 60.)
