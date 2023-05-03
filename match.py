# Boiler plate stuff to start the module
import jpype
import jpype.imports
from jpype.types import *

# Launch the JVM
jpype.startJVM(classpath=['jars/geoxygene-matching-1.10-SNAPSHOT.jar'])

# import the Java modules
from fr.ign.cogit.geoxygene.contrib.appariement.surfaces import ParametresAppSurfaces, AppariementSurfaces
from fr.ign.cogit.geoxygene.util.conversion import ShapefileReader, ShapefileWriter

import sys,os

import shapely
from shapely.geometry import LineString
from shapely import from_wkt
from osgeo import ogr
import geopandas
import numpy

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

# param: index of IDs (same for both layers)
id_index = 1

if len(sys.argv)!=3:
    #layer1 = "./data/bati/bati_95430.shp"
    #layer2 = "./data/bati/cadastre_bati_95430.shp"
    #layer1 = "../../../Data/Test/Neudorf/valid_bati12.shp"
    #layer2 = "../../../Data/Test/Neudorf/valid_bati22.shp"
    layer1 = "./data/bati/neudorf_2012.shp"
    layer2 = "./data/bati/neudorf_2022.shp"
else:
    layer1 = sys.argv[1]
    layer2 = sys.argv[2]

path = layer1.split("/")
layer1name = os.path.splitext(path[len(path)-1])[0]
path.pop(len(path)-1)
path2 = layer2.split("/")
layer2name = os.path.splitext(path2[len(path2)-1])[0]

db1 = ShapefileReader.read(layer1, True)
db2 = ShapefileReader.read(layer2, True)

#reader = ShapefileReader("./data/bati/bati_95430.shp", "bati_95430", None, False)
#print(reader)
#crs = reader.getCRS()
#print(crs) # bug geoxygene: cant read CRS
#print(reader.crs)
#print(reader.reader.getCRS())

# get crs using python - assuming both layers have same CRS
# FIXME handle reprojections
#print(ogr.Open(layer1).GetLayer(0).GetSpatialRef().ExportToWkt()) # libgdal fails: ?
crs = geopandas.read_file(layer1).crs

# hashmap of all features
allfeatures = dict()

# reduce multipolygons into single polygons
l1 = list()
for feat1 in db1:
    for i in range(0,feat1.getGeom().size()):
        single_feat = feat1.cloneGeom()
        currentid = layer1name+"-"+str(single_feat.getAttributes()[id_index].toString())+"-"+str(i)
        single_feat.setGeom(feat1.getGeom().get(i))
        single_feat.setAttribute(0, currentid) # overwrite unused attr fid
        l1.append(single_feat)
        allfeatures[currentid] = single_feat
db1.clear()
for f1 in l1:
    db1.add(f1)

l2 = list()
for feat2 in db2:
    for i in range(0,feat2.getGeom().size()):
        single_feat = feat2.cloneGeom()
        currentid = layer2name+"-"+str(single_feat.getAttributes()[id_index].toString())+"-"+str(i)
        single_feat.setGeom(feat2.getGeom().get(i))
        single_feat.setAttribute(0, currentid)
        l2.append(single_feat)
        allfeatures[currentid] = single_feat
db2.clear()
for f2 in l2:
    db2.add(f2)

# call to geoxygene matching algorithm
liensPoly = AppariementSurfaces.appariementSurfaces(db1, db2, param)

# construct geopandas data frame
attrs = list()
geoms = list()

# to sort links by type
features_stable = list()
features_split = list()
features_merged = list()
features_aggregated = list()
all_link_targets = set()
all_link_sources = set()

for f in liensPoly:
    #print(f.getSchema().getColonnes()) # why do the matching links have no attributes - at least IDs should be kept
    #print(f.getNom()) # empty
    # FIXME getSchema.getColonnes() does not work: not initialised at shapefile reading? -> no attr names; here index hardcoded
    #for ref in f.getObjetsRef():
    #    print(ref.getAttributes()[id_index].toString())
    #print('--')
    #for comp in f.getObjetsComp():
    #    print(comp.getAttributes()[id_index].toString())
    #print('')

    # 1--1 : stability
    if len(f.getObjetsRef())==1 and len(f.getObjetsComp())==1:
        features_stable.append(f.getObjetsComp()[0])
        all_link_sources.add(f.getObjetsComp()[0].getAttribute(0))
        all_link_targets.add(f.getObjetsRef()[0].getAttribute(0))

    # FIXME add comparison geometries and semantic: ex height: important for densification

    # 1 -- n : split
    if len(f.getObjetsRef())==1 and len(f.getObjetsComp())>1:
        for comp in f.getObjetsComp():
            features_split.append(comp)
            all_link_sources.add(f.getObjetsRef()[0].getAttribute(0))
            all_link_targets.add(comp.getAttribute(0))

    # m -- 1 : merged
    if len(f.getObjetsRef())>1 and len(f.getObjetsComp())==1:
        for ref in f.getObjetsRef():
            features_merged.append(ref)
            all_link_sources.add(ref.getAttribute(0))
            all_link_targets.add(f.getObjetsComp()[0].getAttribute(0))

    # m -- n : aggregation
    if len(f.getObjetsRef())>1 and len(f.getObjetsComp())>1:
        for comp in f.getObjetsComp():
            features_aggregated.append(comp)
            all_link_targets.add(f.getObjetsComp()[0].getAttribute(0))
        for ref in f.getObjetsRef():
            all_link_sources.add(f.getObjetsRef()[0].getAttribute(0))

    # iterate over single links in the multiline
    for i in range(0,f.getGeom().size()):
        link = "LINESTRING ("+", ".join(list(map(lambda p: str(p.getX())+" "+str(p.getY()),f.getGeom().get(i).coord().getList())))+")"
        #print(link)
        geoms.append(from_wkt(link))
        attrs.append(list(f.getSchema().getColonnes())) # no attributes?

#print(attrs)
links = geopandas.GeoDataFrame({'geometry':geoms}, crs = crs)
#print(links)

# export links to shp
# geoxygen export fails
#AppariementSurfaces.writeShapefile(liensPoly, "./data/bati/appariement.shp")
#ShapefileWriter.write(liensPoly, "./data/bati/appariement.shp")

links.to_file('/'.join(path)+'/MATCHING-LINKS_'+layer1name+"_"+layer2name+'.shp')

# construct the evolution layer
# FIXME find a way to specify a generic interpretation of matching links
#print(features_stable)
#print(all_link_sources)
#print(features_split)

# construct appeared/disappeared layer

# 1 -- 0
features_disappeared = list()
for f1 in db1:
    if f1.getAttribute(0) not in all_link_sources:
        #print(f1.getAttribute(0))
        features_disappeared.append(f1)

# 0 -- 1
features_appeared = list()
for f2 in db2:
    if f2.getAttribute(0) not in all_link_targets:
         features_appeared.append(f2)

#print(features_appeared)

# export
evol_layer = features_appeared+features_disappeared+features_stable+features_split+features_merged+features_aggregated
#evol_attrs = numpy.repeat('appeared',len(features_appeared))+numpy.repeat('disappeared',len(features_disappeared))+numpy.repeat('stable',len(features_stable))+numpy.repeat('split',len(features_split))+numpy.repeat('merged',len(features_merged))+numpy.repeat('aggregated',len(features_aggregated))

# this does not work: issue with wkt import (reformatting of coordinates as sci notation)
#shapely.wkt.loads(["POLYGON ("+", ".join(list(map(lambda p: str(p.getX())+" "+str(p.getY()),x.getGeom().coord().getList())))+")" for x in evol_layer])
for x in evol_layer:
    wkt = "POLYGON ("+", ".join(list(map(lambda p: str(p.getX())+" "+str(p.getY()),x.getGeom().coord().getList())))+")"
    print(wkt)
    print(from_wkt(wkt))
#print([from_wkt() for x in evol_layer])

jpype.shutdownJVM()


