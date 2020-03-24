import numpy

def TERMD2(land, **kwargs):
    if not land:
        return 0
    return 2000*356/2*1.5*10


from gavin2017processbased import PopulationCapModel, Point
p = PopulationCapModel()
p.alpha = 10 ** -8.07
p.beta = 2.64

def TERMD2(land, longitude, latitude):
    if not land:
        return 0
    else:
        return p.population_capacity(Point(longitude=longitude, latitude=latitude)) + numpy.random.random()
