import numpy as np
import matplotlib.pyplot as plt

y_0 = 1
v_0 = 0
input_data = [[5, 1], [3, 0.2],[5,0.1]]  # each element is a list with two parameters w_m and zeta


def undampedy_f(t, w_n):
    return (v_0 / w_n) * np.sin(w_n * t) + y_0 * np.cos((w_n * t))


def underdamped_f(t, w_n, zeta):
    w_d = w_n * np.sqrt(1 - zeta ** 2)
    c = y_0
    d = (v_0 + y_0 * zeta * w_n) / w_d
    return np.exp(-zeta * w_n * t) * (c * np.cos(w_d * t) + d * np.sin(w_d * t))


def critdampedy_f(t, w_n):
    return (y_0 + ((v_0 + y_0 * w_n) * t)) * np.exp(-w_n * t)


def overdamped_f(t, w_n, zeta):
    return ((v_0 - y_0 * w_n * (-zeta - np.sqrt(zeta ** 2 - 1))) / (2 * w_n * np.sqrt(zeta ** 2 - 1))) * np.exp(
        t * w_n * (-zeta + np.sqrt(zeta ** 2 - 1))) \
        + ((y_0 * w_n * (-zeta + np.sqrt(zeta ** 2 - 1)) - v_0) / (2 * w_n * np.sqrt(zeta ** 2 - 1))) * np.exp(
            t * w_n * (-zeta - np.sqrt(zeta ** 2 - 1)))


input_data = np.array(input_data)
maxt = 3 * (2 * np.pi / (np.max(input_data[:, 0])))
for i, j in input_data:
    y_data = [0]
    amplitude = np.sqrt(y_0 ** 2 + (v_0 / i) ** 2)
    t_data = np.arange(0, maxt, (2 * np.pi / i) * 0.01)
    if j == 0:
        y_data = undampedy_f(t_data, i)
    elif 0 < j < 1:
        y_data = underdamped_f(t_data, i, j)
    elif j == 1:
        y_data = critdampedy_f(t_data, i)
    elif 1 < j:
        y_data = overdamped_f(t_data, i, j)
    else:
        print('LAUDE negative zeta tera BAAP ne banaya hai kya')
    plt.hlines([amplitude, -amplitude], 0, maxt, color='black', lw=1, linestyles='dashed')
    plt.axhspan(amplitude, -amplitude, facecolor='#7FFFD4', alpha=0.1)
    plt.plot(t_data, y_data, label=f'\u03B6 = {j}')

#######################################################################################################################
plt.xlabel('Time (s)')
plt.ylabel('Displacement (m)')
plt.legend(loc=0, shadow=True)
plt.axhline(0, color='black')
plt.axvline(0, color='black')
plt.show()
