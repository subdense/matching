#!/bin/sh

# Test Toulouse center  lat=4.604 lon=1.444
# ./osmdata.sh 2020 1.444 4.604 100

year=$1
lon=$2
lat=$3
radius=$4
name=$year"_"$lon"_"$lat"_"$radius

rm tmp.overpassql
echo "" >> tmp.overpassql
echo "[timeout:1800][date:\"$year-01-01T00:00:00Z\"];" >> tmp.overpassql
echo "(" >> tmp.overpassql
echo "  way[building~\".\"](around:$radius,$lat,$lon);" >> tmp.overpassql
echo "  node(around:$radius,$lat,$lon);" >> tmp.overpassql
echo "  relation(around:$radius,$lat,$lon);" >> tmp.overpassql
echo ");" >> tmp.overpassql
echo "out;" >> tmp.overpassql
echo "" >> tmp.overpassql

wget -O "tmp/buildings_"$name".osm" --post-file=tmp.overpassql "https://overpass-api.de/api/interpreter"
ogr2ogr -skipFailures "tmp/buildings_"$name".shp" "tmp/buildings_"$name".osm" multipolygons -oo USE_CUSTOM_INDEXING=NO
