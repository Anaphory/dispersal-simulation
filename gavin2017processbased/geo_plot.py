def plot_hex_grid(function, all_gridcells):
    polygons = []
    values = []
    for cell in list(all_gridcells.values()):
        polygons.append(cell.polygon())
        values.append(function(cell))
    collection = matplotlib.collections.PolyCollection(
        polygons, facecolors=values)

    return collection

