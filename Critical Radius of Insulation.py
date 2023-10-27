import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as spo

h_o = 10  # W/m^2k
L = 1  # m
K = 4  # W/mk
t_1 = 640  # k (t_i>t_1>t_2>t_o)
t_o = 298  # k (t_i>t_1>t_2>t_o)
r_1, r_2 = (0.11, 0.12)  # m , r_1 < r_2


def f1(x):  # f is equal 0
    if r_1 < r_2:
        return ((t_1 - t_o) * 2 * L * np.pi) / ((np.log(x / r_1) / K) + (1 / (x * h_o)))
    else:
        print('Pipes dimension is unrealistic')


def f2(x):  # f is equal f at r1
    return f1(x) - f1(r_1)


r_crit = K / h_o
r_cs = r_1
if r_crit > r_1:
    i = 1
    while r_cs == r_1:
        r_cs = spo.fsolve(f2, r_crit + 0.001 * i)[0]
        i += 1
t_2 = (K * t_1 + t_o * r_2 * h_o * np.log(r_2 / r_1)) / (K + r_2 * h_o * np.log(r_2 / r_1))
H_1 = f1(r_1)
H_2 = f1(r_2)
H_cs = f1(r_cs)
print(f"Heat loss Rate = {f1(r_2):.2f}W per unit Length"
      f" \nr\u2081 = {r_1}m"
      f"\nr\u2082 = {r_2}m"
      f"\nr_crit = {r_crit:.3f}m (Critical Radius of Insulation)"
      f"\nr_cs   = {r_cs:.3f}m (Radius of Crossover)"
      f"\nt\u2082 = {t_2:.2f}k")
if r_1 >= r_crit or r_2 > r_cs:
    print('The given Parameters makes a GOOD DESIGN '
          f'Heat Loss is {(f1(r_1) - f1(r_2)) * 100 / f1(r_1):.2f}% less than zero thickness Wall')
else:
    print('Change the material properties or increase the Outer Radius so that (r\u2082 > r_cs)\n'
          f'Heat Loss is {(f1(r_2) - f1(r_1)) * 100 / f1(r_1):.2f}% more than zero thickness Wall')

########################## Plot ##########
plt.style.use('bmh')
plt.axis([0, max(r_1, r_2, r_crit, r_cs) * 1.1, 0, f1(r_crit) * 1.05])
plt.title('Critical Radius of Insulation ''\n Heat Loss Rate v/s Outer Radius', fontdict={'fontsize': 15})
plt.xlabel("Outer Radius (Meter)")
plt.ylabel("Heat Loss Rate (Watt)")
# # # # # # # # # # # # # # # # # # # # #
x_list = np.linspace(0.001, max(r_1, r_2, r_crit, r_cs) * 1.05, 501)
y_list = f1(x_list)
plt.plot(x_list, y_list)
plt.plot(np.ones(2) * r_crit, np.linspace(0, f1(r_crit), 2), label='Critical Radius', c='red')
plt.plot(np.ones(2) * r_1, np.linspace(0, f1(r_1), 2), label='r1', c='green')
plt.plot(np.ones(2) * r_2, np.linspace(0, f1(r_2), 2), label='r2', c='green')
plt.legend()
plt.show()
