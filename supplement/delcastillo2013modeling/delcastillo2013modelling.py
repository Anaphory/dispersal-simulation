"""delcastillo2013modeling

Del Castillo and colleagues (2013) build one of the few models I am aware of that model cultural traits as abstract numerical pattern associated to individuals, not connected to cultural practices, in a simulation of hunter-gatherer migration and population dynamics.

Castillo, F. del, J. A. Barceló, L. Mameli, F. Miguel & X. Vila. 2013. Modeling Mechanisms of Cultural Diversity and Ethnicity in Hunter–Gatherers. Journal of Archaeological Method and Theory 21(2). 364–384. doi:10.1007/s10816-013-9199-y.
"""

import numpy

class Household:
    maximum_steps = 1
    tech_efficiency = 1

    def __init__(self):
        self.birth = NOW
        cultural_identity = numpy.random.randint(size=(10, 10))

    @property
    def age(self):
        return NOW - self.birth

for time_step in range(2000):
    survive()
    move()
    identify_others()
    cooperate()
    if not cultural_consensus():
        work_alone()
    die()
