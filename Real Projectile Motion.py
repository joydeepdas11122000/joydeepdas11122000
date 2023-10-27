from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt


def dybdt(t, s, b, m, g):
    x, vx, y, vy = s
    return [vx,
            -(b / m) * np.sqrt(vx ** 2 + vy ** 2) * vx,
            vy,
            -g - (b / m) * np.sqrt(vx ** 2 + vy ** 2) * vy]


def impval(theta, u, b, m, g):  # Returns important Values
    maxt = (2 * u) / g
    n = int(np.log10(maxt)) + 1
    s = solve_ivp(dybdt, [0, maxt], y0=[0, u * np.cos(theta), 0, u * np.sin(theta)],
                  t_eval=np.arange(0, maxt, 0.1 ** (n + 5)), args=(b, m, g))
    # i_before_landing = np.where(np.diff(np.sign(s.y[2])) < 0)[0][0] #this also works
    i_before_landing = np.argmin(np.diff(np.sign(s.y[2])))
    i_before_maxheight = np.where(np.diff(np.sign(s.y[3])) < 0)[0][0]
    return [(s.t[i_before_landing] + s.t[i_before_landing + 1]) / 2,  # Total Time of Flight
            (s.y[2][i_before_maxheight] + s.y[2][i_before_maxheight + 1]) / 2,  # Maximum Height Reached
            (s.y[0][i_before_landing] + s.y[0][i_before_landing + 1]) / 2]  # Horizontal Range


def hrange(theta, u, b, m, g):  # Returns range
    maxt = (2 * u) / g
    n = int(np.log10(maxt)) + 1
    theta = theta * (np.pi / 180)
    s = solve_ivp(dybdt, [0, maxt], y0=[0, u * np.cos(theta), 0, u * np.sin(theta)],
                  t_eval=np.arange(0, maxt, 0.1 ** (n + 5)), args=(b, m, g))
    # i_before_landing = np.where(np.diff(np.sign(s.y[2])) < 0)[0][0] #this also works
    i_before_landing = np.argmin(np.diff(np.sign(s.y[2])))
    return (s.y[0][i_before_landing] + s.y[0][i_before_landing + 1]) / 2  # Horizontal Range


angles = [30,40, 45, 50]  # Initial Launch Angles in degrees
m = 0.25  # mass in kg
u = 1  # initial velocity in m/s
b = 50  # friction Coefficient
g = 9.806  # m/s^2

s = [0] * len(angles)
for i in range(len(angles)):
    theta = angles[i] * (np.pi / 180)
    T, H, R = impval(theta, u, b, m, g)
    s[i] = solve_ivp(dybdt, [0, T], y0=[0, u * np.cos(theta), 0, u * np.sin(theta)],
                     t_eval=np.linspace(0, T, 100),
                     args=(b, m, g))
    plt.plot(s[i].y[0], s[i].y[2], label=r'$\theta_0=' + str(angles[i]) + r'^{\circ}$')
    print(f'for launch angle = {angles[i]} degrees, Range = {R:.6f}m, Max. Height = {H:.3f}m, Total Time = {T:.3f}s')

plt.xlabel('$x (m)$', fontsize=20)
plt.ylabel('$y (m)$', fontsize=20)
plt.ylim(bottom=0)
plt.xlim(left=0)
plt.legend()
plt.show()


def plot_range_vs_theta(u, b, m, g=9.806):  # Returns important Values
    xlist = np.arange(0, 90, 1)
    rangelist = np.zeros(len(xlist))
    for j in range(len(xlist)):
        rangelist[j] = hrange(xlist[j], u, b, m, g)
    plt.plot(xlist, rangelist)
    plt.vlines((xlist[np.argmax(rangelist)]), 0, np.max(rangelist), color='red')
    plt.xlabel(r'$\theta_0$', fontsize=20)
    plt.ylabel('$Range (m)$', fontsize=20)
    plt.ylim(bottom=0)
    plt.xlim(left=0, right=90)
    plt.show()


plot_range_vs_theta(u, 50, m, g)
