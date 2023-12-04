import osmnx
import geopandas
import sys

date = sys.argv[1]
lon = sys.argv[2]
lat = sys.argv[3]
radius = sys.argv[4]
name=date+"_"+lon+"_"+lat+"_"+radius

osmnx.settings.overpass_settings="[out:json][timeout:1800][date:\""+str(date)+"-01-01T00:00:00Z\"]"
osmnx.settings.log_console = True

features = osmnx.features.features_from_point((float(lat),float(lon)),{'building': True}, float(radius))
features_reproj = features.to_crs(crs="EPSG:2154")

features_reproj.to_file('tmp/buildings_'+name+'.shp', engine='pyogrio')
