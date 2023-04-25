# Boiler plate stuff to start the module
import jpype
import jpype.imports
from jpype.types import *

# Launch the JVM
jpype.startJVM(classpath=['jars/geoxygene-matching-1.10-SNAPSHOT.jar'])

# import the Java modules
from fr.ign.cogit.geoxygene.contrib.appariement.surfaces import ParametresAppSurfaces, AppariementSurfaces
from fr.ign.cogit.geoxygene.util.conversion import ShapefileReader, ShapefileWriter

import sys

import shapely
from shapely import LineString
import geopandas

# Parameters of the matching algorithm
param = ParametresAppSurfaces()
param.surface_min_intersection = 1
param.pourcentage_min_intersection = 0.2
param.pourcentage_intersection_sur = 0.8
param.minimiseDistanceSurfacique = True
param.distSurfMaxFinal = 0.6
param.completudeExactitudeMinFinal = 0.3
param.regroupementOptimal = True
param.filtrageFinal = True
param.ajoutPetitesSurfaces = True
param.seuilPourcentageTaillePetitesSurfaces = 0.1
param.persistant = False
param.resolutionMin = 1
param.resolutionMax = 11

if len(sys.argv)!=3:
    layer1 = "./data/bati/bati_95430.shp"
    layer2 = "./data/bati/cadastre_bati_95430.shp"
else:
    layer1 = sys.argv[1]
    layer2 = sys.argv[2]

db1 = ShapefileReader.read(layer1, True)
db2 = ShapefileReader.read(layer2, True)
reader = ShapefileReader("./data/bati/bati_95430.shp", "bati_95430", None, False)
print(reader)
crs = reader.getCRS()
print(crs)
#print(reader.crs)
#print(reader.reader.getCRS())

# reduce multipolygons into single polygons
l1 = list()
for feat1 in db1:
    for i in range(0,feat1.getGeom().size()):
        single_feat = feat1.cloneGeom()
        single_feat.setGeom(feat1.getGeom().get(i))
        l1.append(single_feat)
db1.clear()
for f1 in l1:
    db1.add(f1)

l2 = list()
for feat2 in db2:
    for i in range(0,feat2.getGeom().size()):
        single_feat = feat2.cloneGeom()
        single_feat.setGeom(feat2.getGeom().get(i))
        l2.append(single_feat)
db2.clear()
for f2 in l2:
    db2.add(f2)

liensPoly = AppariementSurfaces.appariementSurfaces(db1, db2, param)

attrs = list()
geoms = list()
for f in liensPoly:
    for i in range(0,f.getGeom().size()):
        link = "LINESTRING ("+", ".join(list(map(lambda p: str(p.getX())+" "+str(p.getY()),f.getGeom().get(i).coord().getList())))+")"
        geoms.append(shapely.from_wkt(link))
        attrs.append(list(f.getSchema().getColonnes())) # no attributes?

#print(attrs)
links = geopandas.GeoDataFrame({'geometry':geoms}, crs = "EPSG:2154")
#print(links)
#links.crs = crs
#links.set_geometry()
#links.to_file('./data/bati/links.shp')

# geoxygen export fails
#AppariementSurfaces.writeShapefile(liensPoly, "./data/bati/appariement.shp")
#ShapefileWriter.write(liensPoly, "./data/bati/appariement.shp")

jpype.shutdownJVM()


