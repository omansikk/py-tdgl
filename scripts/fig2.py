import numpy as np
import matplotlib.pyplot as plt
import tdgl, time, os, h5py
from tdgl.geometry import circle, box
from create_device import create_device

film_radius = 3e3; T2 = Tmax = 0.98; d = 40, T1 = Tmax - 0.05
y = np.linspace(1.01 * film_radius, 0.79 * film_radius, 201)
x = np.zeros_like(y)
cross_section = np.array([x, y]).T
device, tau0 = create_device(film_radius, Tmax, d)

options = tdgl.SolverOptions(
    output_file = "fig2.h5",
    solve_time = tmax,
    field_units = "uT",
    current_units = "nA",
    dt_init = 1e-5,
    dt_max = dt_max,
    save_every = 100,
    adaptive = True,
)

times, jn  = tdgl.get_current_through_paths("./fig2.h5", cross_section, dataset = "normal_current", interp_method = 'linear', units = "uA", with_units = False, progress_bar = True)
times, js  = tdgl.get_current_through_paths("./fig2.h5", cross_section, dataset = "supercurrent", interp_method = 'linear', units = "uA", with_units = False, progress_bar = True)

ns_up = np.zeros(len(times))
ns_down = np.zeros(len(times))
ns_left = np.zeros(len(times))
ns_right = np.zeros(len(times))

phase_up = np.zeros(len(times))
phase_down = np.zeros(len(times))

for i in range(0, len(times)):
    solution.solve_step = i
    psi = solution.interp_order_parameter(np.array([0.81 * film_radius, 0.]))[0]
    ns_right[i] = np.real(psi.conjugate() * psi)

    psi = solution.interp_order_parameter(np.array([-0.833 * film_radius, 0.]))[0]
    ns_left[i] = np.real(psi.conjugate() * psi)

    psi = solution.interp_order_parameter(np.array([0, 0.9 * film_radius]))[0]
    ns_up[i] = np.real(psi.conjugate() * psi)
    phase_up[i] = np.angle(psi)

    psi = solution.interp_order_parameter(np.array([0, -0.9 * film_radius]))[0]
    ns_down[i] = np.real(psi.conjugate() * psi)
    phase_down[i] = np.angle(psi)

font = {'family' : 'Latin Modern Roman',
        'size'   : 18}

plt.rc('font', **font)
plt.rcParams['lines.linewidth'] = 1.5
plt.rcParams['mathtext.fontset'] = "cm"

times = times * tau0

fig, ax = plt.subplots(4, figsize = (18 * cm, 20 * cm))
t = dynamics.time[indices] * tau0
x = np.array([t[0], t[-1]])
y = np.array([1, 1])
eps1 = 1. / T1 - 1.
eps2 = 1. / T2 - 1.
voltage = dynamics.voltage()[indices] * V0
mean_voltage = np.mean(voltage)
phase = dynamics.phase_difference()[indices]
unwrapped_phase = np.unwrap(phase)

# omit the last step because of a quick fix to occasional issues with the original py-tdgl and getting the time and current arrays to be same length
ax[0].plot(times[:-1], jn[:-1], "C0", label = "$j_q$")
ax[0].plot(times[:-1], js[:-1], "C1", label = "$j_s$")
ax[0].plot(times[:-1], js[:-1] + jn[:-1], "k--", label = "$j_q + j_s$")
ax[0].set_ylabel("$j$ (µA)")

ax[1].plot(t, voltage, "C0-")
ax[1].set_ylabel("$\\Delta\\mu$ (µV)")
ax[1].set_yticks([-10, -15, -5])

ax[2].plot(t, unwrapped_phase / (2 * np.pi), "C0")
ax[2].set_ylabel("$\\Delta\\theta/2\\pi$")
ax[2].set_ylim(-0.2, 6.2)
ax[2].set_yticks([0, 2, 4, 6])

ax[3].plot(times, ns_up)
ax[3].plot(times, ns_down)
ax[3].plot(times, ns_left)
ax[3].plot(times, ns_right)
ax[3].plot([times[0], times[-1]], [eps1, eps1], "k--")
ax[3].plot([times[0], times[-1]], [eps2, eps2], "k--")

ax[3].set_ylabel("$n_s$")
ax[3].set_xlabel("$t$ (ns)")

for i in range(4):
    ax[i].grid(axis = "both")

plt.savefig("fig2.svg")
