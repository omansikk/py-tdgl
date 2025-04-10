import numpy as np
import matplotlib.pyplot as plt
import tdgl, time, os, h5py
from create_device import pars, create_device

film_radius = 3e3; T2 = Tmax = 0.98; d = 40
y = np.linspace(1.01 * film_radius, 0.79 * film_radius, 201)
x = np.zeros_like(y)
cross_section = np.array([x, y]).T
device, tau0 = create_device(film_radius, Tmax, d)
dt = 1e-3; save = 50
T1, T2 = Tmax - 0.2, Tmax
tmax = round(1 / tau0)
options = tdgl.SolverOptions(
    output_file = "ring-field.h5",
    solve_time = tmax,
    skip_time = 100 * tmax,
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


phi0 = 90.24
flux_density =  np.linspace(-2 * phi0, 2 * phi0, 201)
start = time.time()
for i, fd in enumerate(flux_density):
    # ~ os.system("rm ring-field.h5") # linux
    os.system("del ring-field.h5") # windows

    solution = tdgl.solve(
        device,
        options,
        disorder_epsilon = epsilon,
        applied_vector_potential = fd,
        seed_solution = None
    )


    times, jn  = tdgl.get_current_through_paths("./ring-field.h5", cross_section, dataset = "normal_current", interp_method = 'linear', units = "uA", with_units = False, progress_bar = True)
    times, js  = tdgl.get_current_through_paths("./ring-field.h5", cross_section, dataset = "supercurrent", interp_method = 'linear', units = "uA", with_units = False, progress_bar = True)
    total = jn + js

    jn_max = np.max(jn)
    js_max = np.max(js)
    total_max = np.max(total)
    jn_min = np.min(jn)
    js_min = np.min(js)
    total_min = np.min(total)

    with open("field_output", "a") as file:
        file.write(f"{fd:2.8f} " + f"{0.5 * (jn_max - jn_min):2.8f} "+ f"{0.5 * (js_max - js_min):2.8f} " + f"{0.5 * (total_max - total_min):2.8f}\n")


output = np.readtxt(field_output)

fig, ax = plt.subplots(2, figsize = (18 / 2.54, 12 / 2.54), constrained_layout = True)
ax[0].plot(output[:, 0], output[:, 1])
ax[0].plot(output[:, 0], output[:, 2])
ax[0].grid(axis = "both")
ax[0].set_yticks([0., 1, 2, 3])
ax[0].set_ylabel("Amplitude $A$ (µA)")

ax[1].plot(output[:, 0], output[:, 3])
ax[1].grid(axis = "both")
ax[1].set_xlabel("Magnetic flux $\Phi / \Phi_0$")
ax[1].set_ylabel("Amplitude $A$ (µA)")
plt.savefig("fig4.svg")

end = time.time()
print((end - start) / 60.)
