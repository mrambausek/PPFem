import numpy as np


class Point(object):

    def __init__(self, *coords, index=None):
        self._coords = np.array(coords)
        self.index = index

    def __getitem__(self, item):
        return self._coords[item]

    def coords(self):
        return self._coords
