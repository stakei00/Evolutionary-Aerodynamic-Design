
# Implement the default Matplotlib key bindings.
from matplotlib.figure import Figure

import numpy as np

def plot_1():
    fig = Figure(figsize=(12, 4), dpi=100)
    t = np.arange(0, 3, .01)
    fig.add_subplot(111).plot(t, 2 * np.sin(2 * np.pi * t))
    return fig




