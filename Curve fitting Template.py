import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

x_data = np.array([0., 0.15789474, 0.31578947, 0.47368421, 0.63157895,
                   0.78947368, 0.94736842, 1.10526316, 1.26315789, 1.42105263,
                   1.57894737, 1.73684211, 1.89473684, 2.05263158, 2.21052632,
                   2.36842105, 2.52631579, 2.68421053, 2.84210526, 3.])
y_data = np.array([2.95258285, 2.49719803, -2.1984975, -4.88744346,
                   -7.41326345, -8.44574157, -10.01878504, -13.83743553,
                   -12.91548145, -15.41149046, -14.93516299, -13.42514157,
                   -14.12110495, -17.6412464, -16.1275509, -16.11533771,
                   -15.66076021, -13.48938865, -11.33918701, -11.70467566])


def model_f(x, a, b, c):
    return a * ((x - b) ** 2) + c


p_opt, p_cov = curve_fit(model_f, x_data, y_data, p0=[3, 2, -16])
a, b, c = p_opt
# plotting the data
x_model = np.linspace(min(x_data), max(y_data), 100)
y_model = model_f(x_model, a, b, c)
plt.scatter(x_data, y_data)
plt.plot(x_model, y_model, c='green')
plt.show()