Check Sky Objects
=================

Jupyter notebook to inspect objects coordinates at thge Horizontal Coordinate System (Alt, Az) at a certain location and time. 

Dependencies
------------

In order to run this notebook, you will need python 3.9 or higher for the [zoneinfo](https://docs.python.org/3/library/zoneinfo.html). It also depends on [matpltotlib](https://matplotlib.org/) (v3.6 or higher) and [astropy](https://www.astropy.org/) (v5.1.1).

Usage
-----

For the correct calculation of the positions at a certain location, one should configure the correct site GPS coordinates at the function `get_myhome_location_and_time()`:

```python
def get_myhome_location_and_time():
    HOME_LAT = -27.480997*u.deg    # Latitude
    HOME_LON = -48.401418*u.deg    # Longitude
    HOME_HEI = 20*u.m              # Height
    HOME_TZ = 'America/Sao_Paulo'  # Timezone
    return get_earth_location_and_time(HOME_LAT, HOME_LON, HOME_HEI, HOME_TZ)
```

For different timezones, see [https://en.wikipedia.org/wiki/List_of_tz_database_time_zones].

A list of objects can be inserted at the dictionary `objects`. An example of this dictionary is on the notebook. The configuration is: `object name` as the key of the dictionary and a value-dictionary with the object equatorial system coordinates: right-ascencion (`ra`) and declination (`dec`). Example:

```python
objects = {
    'LMC': {'c': None, 'ra': 80.89417*u.deg, 'dec': -69.75611*u.deg},
    'SMC': {'c': None, 'ra': 13.15833*u.deg, 'dec': -72.80028*u.deg},
    'Antares': {'c': None, 'ra': 247.35191542*u.deg, 'dec': -26.4300261*u.deg},
    'Acrux': {'c': None, 'ra': 186.6495634*u.deg, 'dec': -63.09909286*u.deg},    
}
```

Contact
-------
	
Contact us: [dhubax@gmail.com](mailto:dhubax@gmail.com).