from matplotlib import pyplot as plt
import numpy
from model import GridCell, BoundingBox, hexagonal_earth_grid

GridCell.grid = hexagonal_earth_grid(
    BoundingBox(s=-1, n=1, w=-1, e=1), 450000000)
c = GridCell.random_cell()
