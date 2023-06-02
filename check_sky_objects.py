#!/usr/bin/env python3
import sys
import numpy as np
import configparser
import astropy.units as u
from os.path import basename
from zoneinfo import ZoneInfo
from astropy.time import Time
from matplotlib import colormaps
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
from datetime import datetime, timedelta
from matplotlib.dates import DateFormatter
from astropy.coordinates import AltAz, EarthLocation, SkyCoord, Angle, get_sun, get_moon, SkyCoord

__script_name__ = basename(sys.argv[0])
__defulat__ = 'deg'
__defulon__ = 'deg'
__t80slat__ = -30.1678639
__t80slon__ = -70.8056889
__t80shei__ = 2187
__t80stz__ = 'America/Santiago'

def get_earth_location_and_Time(latitude, longitude, height, location_dt):
    # SET EARTH LOCATION
    local_coordinates = EarthLocation(lat=latitude, lon=longitude, height=height)
    # SET TIME OF LOCATION
    local_now = Time(location_dt, location=local_coordinates)
    return local_coordinates, local_now

def create_location_AltAz_timeline(begin_time_ref, end_time_ref, location_time, location_coordinates, N=100):
    timeline = np.linspace(begin_time_ref, end_time_ref, N)*u.hour
    location_timeline = location_time + timeline
    location_timeline_AltAz = AltAz(obstime=location_timeline, location=location_coordinates)
    return location_timeline, location_timeline_AltAz, timeline

def create_sun_and_moon_location_AltAz_timeline(location_timeline, location_timeline_AltAz):
    sun_location_AltAz_timeline = get_sun(location_timeline).transform_to(location_timeline_AltAz)
    moon_location_AltAz_obstime = get_moon(location_timeline).transform_to(location_timeline_AltAz)
    return sun_location_AltAz_timeline, moon_location_AltAz_obstime

def get_location_moon_illumination(location_time):
    sun = get_sun(location_time)
    moon = get_moon(location_time)
    elongation = sun.separation(moon)
    moon_phase_angle = np.arctan2(
        sun.distance*np.sin(elongation), 
        moon.distance - sun.distance*np.cos(elongation)
    )
    moon_illumination = (1 + np.cos(moon_phase_angle))/2.0    
    return moon_illumination

def get_location_latlon(config):
    uconf = config['units']
    lconf = config['location']
    ulat = u.Unit(uconf.get('loc_lat_unit', __defulat__))
    ulon = u.Unit(uconf.get('loc_lat_unit', __defulon__))
    llat = lconf.getfloat('lat', __t80slat__)
    llon = lconf.getfloat('lon', __t80slon__)
    return Angle(llat, unit=ulat), Angle(llon, unit=ulon)

def get_config(filename):
    config = configparser.ConfigParser(
        # allow a variables be set without value
        allow_no_value=True,
        # allows duplicated keys in different sections
        strict=False,
        # deals with variables inside configuratio file
        interpolation=configparser.ExtendedInterpolation())
    config.optionxform = str
    config.read(filename)
    return config

def get_objects(config):
    uconf = config['units']
    ura = u.Unit(uconf.get('obj_ra_unit', __defulat__))
    udec = u.Unit(uconf.get('obj_dec_unit', __defulon__))
    oconf = config['objects']
    objects = {}
    for name, coords in oconf.items():
        print(name, coords)
        sRA, sDEC = coords.split(',')
        ra = Angle(sRA, unit=ura)
        dec = Angle(sDEC, unit=udec) 
        objects[name] = {'ra': ra, 'dec': dec}
    return objects

if __name__ == '__main__':
    try:
        config_file = sys.argv[1]
    except:
        print(f'usage: {__script_name__} CONFIG_FILENAME')
        sys.exit(1)        

    config = get_config(config_file)

    #### location
    llat, llon = get_location_latlon(config)
    lhei = config['location'].getfloat('hei', __t80shei__)*u.m
    ltz = ZoneInfo(config['location'].get('timezone', __t80stz__))

    #### datetime
    dtfmt = config['datetime'].get('datetime_format')
    ldt = config['datetime'].get('datetime')
    if 'now' in ldt:
         ldt = datetime.now(tz=ltz)
    else:
        ldt = datetime.strptime(ldt, dtfmt).astimezone(ltz)

    #### timeline
    begin_time_ref = config['datetime'].getint('begin_time_ref', 0)
    end_time_ref = config['datetime'].getint('end_time_ref', 24)

    #### objects
    objects = get_objects(config)

    lcoords, lTime = get_earth_location_and_Time(llat, llon, lhei, ldt)
    print('Location Earth Coordinates: ', lcoords)
    print('Location Time: ', lTime)  # Time displays UTC ?

    obstime, l_AltAz_obstime, _time_lapse = create_location_AltAz_timeline(
        begin_time_ref=begin_time_ref, 
        end_time_ref=end_time_ref,
        location_time=lTime, 
        location_coordinates=lcoords,
        N=100,
    )

    # Moon and sun positions at the location's timeline
    lsun_AltAz_obstime, lmoon_AltAz_obstime = create_sun_and_moon_location_AltAz_timeline(
        location_timeline=obstime, 
        location_timeline_AltAz=l_AltAz_obstime
    )

    lmoon_illumination = get_location_moon_illumination(lTime)

    # Objects positions at the location's timeline
    obj_AltAz_obstime = {}
    obj_moon_sep = {}
    for i, (k, v) in enumerate(objects.items()):
        c = SkyCoord(ra=v['ra'], dec=v['dec'], frame='icrs')
        c_AltAz_lTime = c.transform_to(AltAz(obstime=lTime, location=lcoords))
        objects[k]['AltAz_obstime'] = c_AltAz_lTime.transform_to(l_AltAz_obstime)
        _tmp = objects[k]['AltAz_obstime']
        objects[k]['moon_sep'] = _tmp.separation(lmoon_AltAz_obstime)

    mask_night = lsun_AltAz_obstime.alt < -0*u.deg
    mask_twilight = lsun_AltAz_obstime.alt < -18*u.deg

    time_lapse = np.asarray(
        [ldt + timedelta(hours=x.value) for x in _time_lapse]
    )

    # PLOT
    f, (ax1, ax2, ax3) = plt.subplots(1, 3)
    f.set_size_inches(15, 5)
    f.suptitle('Moon Illumination: {:.1f}%'.format(lmoon_illumination*100))

    ax1.plot(time_lapse, lsun_AltAz_obstime.alt, color='k', ls='--', lw=0.5)
    ax1.plot(time_lapse, lmoon_AltAz_obstime.alt, color='orange', ls='--', lw=0.5)
    ax1.plot(time_lapse[mask_twilight], lsun_AltAz_obstime.alt[mask_twilight], color='k', label='sun')
    ax1.plot(time_lapse[mask_twilight], lmoon_AltAz_obstime.alt[mask_twilight], color='orange', label='moon')
    ax1.fill_between(time_lapse, -90, 90, mask_night, color='0.8', zorder=0)
    ax1.fill_between(time_lapse, -90, 90, mask_twilight, color='0.6', zorder=0)
    ax1.axhline(y=0, ls='--', color='k')
    ax1.set_ylabel('Altitude [deg]')

    ax2.plot(time_lapse, lsun_AltAz_obstime.az, color='k', ls='--', lw=0.5)
    ax2.plot(time_lapse, lmoon_AltAz_obstime.az, color='orange', ls='--', lw=0.5)
    ax2.plot(time_lapse[mask_twilight], lsun_AltAz_obstime.az[mask_twilight], color='k', label='sun')
    ax2.plot(time_lapse[mask_twilight], lmoon_AltAz_obstime.az[mask_twilight], color='orange', label='moon')
    ax2.fill_between(time_lapse, 0, 360, mask_night, color='0.8', zorder=0)
    ax2.fill_between(time_lapse, 0, 360, mask_twilight, color='0.6', zorder=0)
    ax2.axhline(y=0, ls='--', color='k')
    ax2.legend(frameon=False, loc=1)

    ax3.fill_between(time_lapse, 0, 360, mask_night, color='0.8', zorder=0)
    ax3.fill_between(time_lapse, 0, 360, mask_twilight, color='0.6', zorder=0)
    ax3.set_title('Obj-Moon separation')
    ax3.set_ylabel('Separation [deg]')

    # PLOT OBJECTS
    n_obj = len(objects)
    cmap_key = config['plot'].get('objects_colormap', 'magma')
    cmap_norm = Normalize(vmin=0, vmax=n_obj-1)
    for i, (k, v) in enumerate(objects.items()):
        _c = colormaps[cmap_key](cmap_norm(i))
        _sep = v['moon_sep']
        _obj = v['AltAz_obstime']
        ax1.plot(time_lapse, _obj.alt, color=_c, ls='--', lw=0.5)
        ax1.plot(time_lapse[mask_twilight], _obj.alt[mask_twilight], color=_c, label=k)   
        ax2.plot(time_lapse, _obj.az, color=_c, ls='--', lw=0.5)
        ax2.plot(time_lapse[mask_twilight], _obj.az[mask_twilight], color=_c, label=k)
        ax3.plot(time_lapse, _sep, color=_c, label=k)
        if i:
            if max(_sep.value) > _range[1]:
                _range[1] = max(_sep.value)
            if min(_sep.value) < _range[0]:
                _range[0] = min(_sep.value)        
        else:
            _range = [min(_sep.value), max(_sep.value)]

    ax1.set_ylim(-90, 90)
    ax2.set_ylim(0, 360)
    ax3.set_ylim(_range[0]*0.95, _range[1]*1.05)

    for ax in [ax1, ax2, ax3]:
        nbins = (end_time_ref - begin_time_ref)
        formatter = DateFormatter('%H:%M:%S', tz=ltz)
        ax.xaxis.set_major_formatter(formatter)
        plt.setp(ax.get_xticklabels(), rotation=45)      
        ax.set_xlabel('Time')
        ax.legend(frameon=False, loc=1)
        ax.grid()

    _dpi = config["plot"].getint("image_dpi", 100)
    _imext = config["plot"].get("image_extention", "png")
    _plot_filename = f'location_sky.{_imext}'
    f.tight_layout()
    f.savefig(_plot_filename, dpi=_dpi)