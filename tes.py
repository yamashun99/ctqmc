import numpy as np
import matplotlib.pyplot as plt
from numba import jit, f8
from numba.experimental import jitclass
import time
spec = [
    ('arr', f8[:, :, :]),
    ('new_arr', f8[:, :, :]),
]


@jitclass(spec)
class class_jit():
    def __init__(self, arr):
        self.arr = arr
        self.new_arr = np.zeros_like(arr, dtype=np.float64)

    def func(self):
        shape = self.arr.shape
        for x in range(shape[0]):
            for y in range(shape[1]):
                for z in range(shape[2]):
                    self.new_arr[x, y, z] = 1e6*x+1e3*y+z
        return self.new_arr

    def plot(self):
        plt.imshow(self.new_arr[0])


if __name__ == '__main__':
    n = 256
    arr = np.zeros([n]*3)
    st = time.time()
    obj1 = class_jit(arr)
    obj1.func()
    print(f'"class jit " elapsed_time {time.time()-st:2f}')
