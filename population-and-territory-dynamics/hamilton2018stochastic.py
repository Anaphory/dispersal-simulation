"""hamilton2018stochastic

A replication of the model of hunter-gatherer population dynamics, with an
extinction threshold, described by Hamilton (2018).

Reference:
Hamilton, Marcus J & Walker, Robert S. 2018. A stochastic density-dependent
model of long-term population dynamics in hunter-gatherer populations.
Evolutionary Ecology Research 19(1). 85–102.
"""

from matplotlib import pyplot as plt
import numpy

TIME_STEPS_PER_YEAR = 1 #year
BAND_SIZE = 20
GENERATION_TIME = 28 * TIME_STEPS_PER_YEAR
demographic_variance = numpy.random.normal(
    0.24,
    scale = 0.02 / 1.96) # std in terms of the give 95% confidence interval
environmental_variance = numpy.random.uniform(0.005, 0.04)
catastrophe_rate = numpy.random.uniform(50, 400) * TIME_STEPS_PER_YEAR
carrying_capacity = 650

catastrophe_probability = 1/catastrophe_rate
intrinsic_growth_rate = 1/GENERATION_TIME

for _ in range(20):
    run = []
    log_pop_size = numpy.log(BAND_SIZE)
    while log_pop_size >= numpy.log(BAND_SIZE):
        run.append(log_pop_size)
        variance = (
            environmental_variance + demographic_variance /
            numpy.exp(log_pop_size)) / numpy.exp(
                2 * intrinsic_growth_rate) / TIME_STEPS_PER_YEAR
        log_pop_size = (
            # x(t)
            log_pop_size +
            # r(t)
            numpy.random.normal(
                loc=intrinsic_growth_rate - variance / 2, scale=variance ** 0.5) -
            # α
            (intrinsic_growth_rate / carrying_capacity) *
            # N(t)
            numpy.exp(log_pop_size) +
            # 1/f ln(1-s)
            (numpy.log(numpy.random.random())
            if numpy.random.random() < catastrophe_probability
            else 0)
        )
    run.append(0)
    plt.plot(numpy.arange(0, len(run)) / TIME_STEPS_PER_YEAR, run)
plt.ylim(0, 15)
plt.xlim(0, 1000)
plt.show()
