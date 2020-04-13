import csv
import enum
import numpy
import tifffile

import cartopy
import cartopy.crs as ccrs
import cartopy.geodesic as geodesic
from matplotlib import pyplot as plt

degree_celsius = float

def coordinates_to_index(points, resolution=2 * 60):
    """Convert long,lat coordinate pairs into indices in a TIF

    Convert a [..., 2] ndarray, or a pair of coordinates, into the matching
    grid indices of a Mercator projection pixel map with a given resolution in
    pixels per degree.

    Paramaters
    ==========
    points: ndarray-like, shape=(..., 2)
        An array of longitude. latitude pairs to be converted to grid indices

    resolution:
        The resolution of the grid in indices per degree

    Returns
    =======
    ndarray(int), shape=(..., 2)
        An integer array of grid indices

    """
    points = numpy.asarray(points)
    return numpy.stack(
        (numpy.round((-points[..., 1] + 90) * resolution).astype(int),
         numpy.round((points[..., 0] + 180) * resolution).astype(int)),
        -1)


def index_to_coordinates(indices, resolution=2 * 60):
    """Convert grid indices into long,lat coordinate pairs

    Convert a (ndarray of) grid indices of a Mercator projection pixel map with
    a given resolution in pixels per degree into geocoordinate pairs (longitude,
    latitude).

    Paramaters
    ==========
    indices: ndarray(int), shape=(..., 2)
        An integer array of grid indices

    resolution:
        The resolution of the grid in indices per degree

    Returns
    =======
    ndarray, shape=(..., 2)
        An array of longitude. latitude pairs

    """
    indices = numpy.asarray(indices)
    return numpy.stack(
        (indices[..., 1] / resolution - 180,
         90 - indices[..., 0] / resolution),
        -1)



def TERMD2(land, **kwargs):
    if not land:
        return 0
    return 2000*356/2*1.5*10


from gavin2017processbased import PopulationCapModel, Point
p = PopulationCapModel()
p.alpha = 4 * 10 ** -8.07
p.beta = 2.64

def TERMD2(land, longitude, latitude):
    if not land:
        return 0
    else:
        return p.population_capacity(Point(longitude=longitude, latitude=latitude)) + numpy.random.random()

class CLIM(enum.Enum):
    """Temperature ordination of climates CLIM."""
    # Table 4.02
    POLAR = 1
    BOREAL = 2
    COOL = 3
    WARM = 4
    SUBTROPICAL = 5
    TROPICAL = 6
    EQUATORIAL = 7

    @classmethod
    def from_effective_temperature(cls, effective_temperature: degree_celsius, **ignored):
        if effective_temperature < 10:
            return cls.POLAR
        elif effective_temperature < 12.50:
            return cls.BOREAL
        elif effective_temperature < 14.56:
            return cls.COOL
        elif effective_temperature < 16.62:
            return cls.WARM
        elif effective_temperature < 18.16:
            return cls.SUBTROPICAL
        elif effective_temperature < 25.58:
            return cls.TROPICAL
        else:
            return cls.EQUATORIAL

class SEASON(enum.Enum):
    """Which season is rain season?"""
    # p. 71, around Fig. 4.06
    SPRING = 1
    SUMMER = 2
    AUTUMN = 3
    WINTER = 4
    @classmethod
    def from_rain_heat_correlation(cls, rain_heat_correlation, **ignored):
        assert 0 <= rain_heat_correlation < 12, "rain_heat_correlation must be between 0 and 12!"
        if rain_heat_correlation < 3:
            return cls.SPRING
        elif rain_heat_correlation < 6:
            return cls.SUMMER
        elif rain_heat_correlation < 9:
            return cls.AUTUMN
        else: # rain_heat_correlation < 12:
            return cls.WINTER


precipitations = [
    tifffile.imread("/home/gereon/Public/settlement-of-americas/supplement/worldclim/wc2.1_5m_prec_{:02d}.tif".format(i + 1)).clip(0)
    for i in range(12)
    ]

temperature = [
    tifffile.imread("/home/gereon/Public/settlement-of-americas/supplement/worldclim/wc2.1_5m_tavg_{:02d}.tif".format(i + 1)).clip(0)
    for i in range(12)
    ]

class ClimateLocation:
    def __init__(self, longitude, latitude):
        self.longitude = longitude
        self.latitude = latitude
        self.index = tuple(coordinates_to_index((longitude, latitude), 60/5))

        self.mean_monthly_temperatures = [p[self.index] for p in temperature]
        self.mean_monthly_precipitations = [p[self.index] for p in precipitations]

        assert len(self.mean_monthly_temperatures) == 12
        assert len(self.mean_monthly_precipitations) == 12
    @property
    def MCM(self):
        """Compute mean warmest month temperature MWM.

        """
        try:
            return self.mean_coldest_month_temperature
        except AttributeError:
            self.mean_coldest_month_temperature = min(self.mean_monthly_temperatures)
            return self.mean_coldest_month_temperature

    @property
    def MWM(self):
        """Compute mean warmest month temperature MWM.

        """
        try:
            return self.mean_warmest_month_temperature
        except AttributeError:
            self.mean_warmest_month_temperature = max(self.mean_monthly_temperatures)
            return self.mean_warmest_month_temperature

    @property
    def ET(self):
        """Compute effective temperature ET.

        Effective temperature (ET) is a measure that was specifically designed by
        Bailey (1960) to examine the biological implications of ambient warmth,
        which is another way of referring to the amount of solar energy available
        at any given location. ET is calculated in terms of three empirically
        determined constants:

        1. The minimal mean temperature (18°C) of the coldest month of the year
        that will sustain tropical plant communities (those having 365-day
        growing seasons).
        2. The minimal mean temperature (10°C) expected at the beginning and end of
        the growing season along the zonal boundary between polar and boreal
        environments.
        3. The minimal mean temperature (8°C) at the beginning and end of the
        growing season, or—at the earth's poles—the warmest month.

        The choice of these constants permits the calculation of effective
        temperature values that relate directly to the major biological boundaries
        recognized empirically as marking transitions in biological activity. The
        relationship between ET and the earth’s biotic communities can be
        visualized as a scale along which ET values range. An ET value of 18°C or
        higher corresponds to those places in the world where there is no killing
        frost during the year. Locations that have an ET value of 10°C or less are
        characterized by fewer than thirty days of the year without a killing
        frost. Effective temperature provides more biologically relevant
        information than simple mean annual temperature, and, in addition, the
        magnitude of the value is directly indicative of the duration of the
        growing season.

        """
        # 4.01
        try:
            return self.effective_temperature
        except AttributeError:
            self.effective_temperature = (18 * self.MWM - 10 * self.MCM) / (self.MWM - self.MCM + 8)
            return self.effective_temperature

    @property
    def TEMP(self):
        """Compute temperateness TEMP.

        Bailey (1960:10) has developed a variable that he refers to as
        temperateness (TEMP), which combines aspects of evenness with a “comfort”
        judgment based on human experience. This measure tracks differences in
        temperature range between adjacent months at specific locations, with a
        positive bias in favor of locations where mean winter temperatures are
        above 0°C. When large month-to-month differences are noted, the location is
        said to lack temperateness, whereas the location is referred to as
        “temperate” when temperatures are very similar. Temperateness values may be
        calculated from the basic meteorological data available for each
        hunter-gatherer group in the sample, using the same constants used to
        calculate values of ET.

        """
        # 4.02
        return 161.7 - 41 * numpy.log((self.MWM - 10)**2 + (self.MCM - 18)**2)/numpy.log(10)

    @property
    def MTEMP(self):
        """Compute temperature extremeness measure MTEMP.

        In order to examine temperature constancy in greater detail, an additional
        measure, MTEMP, was calculated. In human terms, values of MTEMP provide a
        measure of the scale of temperature extremes with which an inhabitant of a
        specific locality would have to cope.

        """
        # 4.03
        try:
            return self.temperature_extremeness
        except AttributeError:
            self.temperature_extremeness = ((MCM + 45)/(MWM + 45)) * 100
            return self.temperature_extremeness

    @property
    def RHIGH(self):
        """Compute the maximum monthly rainfall RHIGH.

        """
        return max(self.mean_monthly_precipitations)

    @property
    def CRR(self):
        """Compute mean annual precipitation CRR.

        """
        return sum(self.mean_monthly_precipitations)

    @property
    def REVEN(self):
        """Compute rainfall evenness REVEN.

        REVEN is a straightforward measure resulting from a two-step calculation:
        first, total annual rainfall is divided by 12.0 to obtain the value
        expected for any one month, assuming that all months had identical
        rainfall. Then the value of the month observed to have the highest rainfall
        (RHIGH) is divided by the expected monthly value. REVEN is distributed so
        that a value of 1.0 would occur if all months received the same amount of
        rainfall, and it increases as the magnitude of the wettest month rises
        relative to the mean annual value. The highest possible value would be
        12.0, which indicates that all rain occurred in one month and no rain
        occurred in any other months.

        """
        # 4.04
        try:
            return self.rainfall_evenness
        except AttributeError:
            self.rainfall_evenness = self.RHIGH/(self.CRR/12)
            return self.rainfall_evenness

    @property
    def RCORR2(self):
        """Compute the rain/heat correlation RRCORR2.

        The number of months or parts thereof, positive or negative, that separate
        the wettest month from the warmest month was calculated and assigned the
        term RRCORR. All positive values of RRCORR refer to the number of months
        after the warmest month. Although these conventions were used in the
        initial tabulations of monthly weather data, in subsequent analysis the
        positive and negative scale became confusing, so a value of 4.5 was added
        to the recorded scale and any negative values remaining after this addition
        were subtracted from a value of 12.0. The result is a positive scale
        running from 0.0 to 12.0 and referred to as RRCORR2.

        """

        return (numpy.argmax(self.mean_monthly_precipitations) - numpy.argmax(self.mean_monthly_temperatures) + 4.5) % 12

    @property
    def SEASON(self):
        return SEASON.from_rain_heat_correlation(self.RCORR2)

    @property
    def CLIM(self):
        return CLIM.from_effective_temperature(self.ET)

def draw_and_plot_big_cities():
    """ Draw random cities with a population of at least 1000 inhabitants and plot them on a map.

    Returns the geocoordinates, as (longitude, latitude) pairs, of those cities.
    """
    random_big_cities = []
    with open("../geonames/cities1000.txt") as bigcities:
        for line in csv.reader(bigcities, delimiter="\t"):
            if numpy.random.random()<5e-3:
                random_big_cities.append((float(line[5]), float(line[4])))

    plt.figure()
    cmap = plt.get_cmap("viridis")
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines("50m")
    # ax.set_extent(gavin2017processbased.americas)
    plt.scatter(*zip(*random_big_cities))
    plt.show()
    return random_big_cities

def crr_by_lat(random_big_cities):
    """ Plot mean annual rainfall (CRR) vs. latitude"""
    # Figure 4.06
    plt.figure()
    data = [ClimateLocation(*l) for l in random_big_cities]
    for season in [SEASON.SPRING, SEASON.SUMMER, SEASON.AUTUMN, SEASON.WINTER]:
        y = [abs(l.latitude) for l in data
            if l.SEASON == season]
        x = [l.CRR for l in data
            if l.SEASON == season]
        s = "01234"[season.value]
        c = [None, "green", "orange", "red", "blue"][season.value]
        plt.scatter(x, y, label=season.name, marker=s, c=c)
    plt.legend()
    plt.xlim((-2000, 12000))
    plt.ylim((-10, 90))
    plt.xlabel("Mean Annual Rainfall (CRR) [mm]")
    plt.ylabel("Latitude (absolute value, North and South) [°]")
    plt.show()
    return data
