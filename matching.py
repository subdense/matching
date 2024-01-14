# import requires jpype JVM to already run
import json
from datetime import datetime

import jpype

from fr.ign.cogit.geoxygene.contrib.appariement.surfaces import ParametresAppSurfaces, AppariementSurfaces
from fr.ign.cogit.geoxygene.util.conversion import ShapefileReader
from fr.ign.cogit.geoxygene.feature import SchemaDefaultFeature, DefaultFeature, Population
from fr.ign.cogit.geoxygene.api.spatial.coordgeom import IPolygon
from fr.ign.cogit.geoxygene.schema.schemaConceptuelISOJeu import FeatureType, AttributeType
from java.util import HashMap
from fr.ign.cogit.geoxygene.util.index import Tiling, STRtreeJts

import sys,os
from os.path import exists

from shapely import from_wkt, Polygon
import geopandas
import numpy
import datetime

def default_params():
    params = dict()

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

    # test data by default: committed in data/bati
    layer1 = "./data/bati/neudorf_2012.shp"
    layer2 = "./data/bati/neudorf_2022.shp"

    params['algo_params']=param
    params['id_index']=id_index
    params['layer1']=layer1
    params['layer2']=layer2
    if exists('./data/parameters.json'):
        with open('./data/parameters.json') as parameter_file:
            params_from_file = json.load(parameter_file)
            print(params_from_file)
            params |= params_from_file
    
    return(params)

#' FIXME handle reprojections
#' FIXME paths must be the same
def get_data(params):
    if len(sys.argv)!=3:
        layer1 = params['layer1']
        layer2 = params['layer2']
    else:
        layer1 = sys.argv[1]
        layer2 = sys.argv[2]

    path1 = layer1.split("/")
    layer1name = os.path.splitext(path1[len(path1)-1])[0]
    path2 = layer2.split("/")
    layer2name = os.path.splitext(path2[len(path2)-1])[0]

    path1.pop() # dirty python mutables
    path = [".","output_data"]

    print("READ DB1")
    db1 = ShapefileReader.read(layer1, True)
    print("READ DB2")
    db2 = ShapefileReader.read(layer2, True)
    # get the list of feature attributes for db1
    print("DB1 attributes = " + str(db1.getFeatureType().getFeatureAttributes()))
    # get the list of feature attributes for db2
    print("DB2 attributes = " + str(db2.getFeatureType().getFeatureAttributes()))

    #crs = geopandas.read_file(layer1, engine="pyogrio").crs
    crs = geopandas.read_file(layer1).crs

    return(layer1name, layer2name, path, db1, db2, crs)

def get_data_and_preprocess(params, layer):
    layer = params[layer]

    path = layer.split("/")
    layername = os.path.splitext(path[len(path)-1])[0]

    path.pop() # dirty python mutables
    path = [".","output_data"]

    print(str(datetime.datetime.now())+ " - READ " + layer)
    db = ShapefileReader.read(layer, True)
    # get the list of feature attributes for db1
    print(str(datetime.datetime.now())+" - DB attributes = " + str(db.getFeatureType().getFeatureAttributes()))
    newdb = preprocess_layer(layername, db, 0)
    print(str(datetime.datetime.now())+" - preprocess done")
    crs = geopandas.read_file(layer,rows=1).crs
    return (layername, path, newdb, crs)

def preprocess_layer(layername, db, id_attribute):
    newFeatureType = FeatureType()
    newFeatureType.setTypeName("building")
    newFeatureType.setGeometryType(IPolygon.class_)
    id = AttributeType("id", "String")
    newFeatureType.addFeatureAttribute(id)
    schema = SchemaDefaultFeature()
    schema.setFeatureType(newFeatureType)
    newFeatureType.setSchema(schema)
    attLookup = HashMap()#<jpype.JInt, jpype.JString[:]>(0)
    attLookup.put(jpype.JInt(0), jpype.JString[:]@[id.getNomField(), id.getMemberName()])
    schema.setAttLookup(attLookup)
    newdb = Population(False, layername, DefaultFeature.class_, True)
    newdb.setFeatureType(newFeatureType)
    # reduce multipolygons into single polygons
    for feat in db:
        # print(str(feat1.getAttributes()))
        if feat.getAttributes()[id_attribute]:
            for i in range(0,feat.getGeom().size()):
                currentid = layername+"-"+str(feat.getAttributes()[id_attribute].toString())+"-"+str(i)
                n = newdb.nouvelElement(feat.getGeom().get(i))
                n.setSchema(schema)
                attributes = jpype.JObject[:]@[currentid]
                n.setAttributes(attributes)
        # else:
            # print("no attribute " + str(id_attribute) + " for " + str(feat.toString()))
            # newdb.initSpatialIndex(Tiling.class_, False)
            # print(str(datetime.datetime.now())+" - initSpatialIndex done")
            # return newdb
    index = STRtreeJts(newdb)
    newdb.setSpatialIndexToExisting(index)
    print(str(datetime.datetime.now())+" - initSpatialIndex done")
    # print("db preprocessed")
    return newdb

def preprocess_data(layer1name, layer2name, db1, db2, id_index):
    # hashmap of all features
    #allfeatures = dict()
    newFeatureType = FeatureType()
    newFeatureType.setTypeName("building")
    newFeatureType.setGeometryType(IPolygon.class_)
    id = AttributeType("id", "String")
    newFeatureType.addFeatureAttribute(id)
    schema = SchemaDefaultFeature()
    schema.setFeatureType(newFeatureType)
    newFeatureType.setSchema(schema)
                
    attLookup = HashMap()#<jpype.JInt, jpype.JString[:]>(0)
    attLookup.put(jpype.JInt(0), jpype.JString[:]@[id.getNomField(), id.getMemberName()])
    schema.setAttLookup(attLookup)
    newdb1 = Population(False, "db1", DefaultFeature.class_, True)
    newdb1.setFeatureType(newFeatureType)

    # reduce multipolygons into single polygons
    # l1 = list()
    for feat1 in db1:
        for i in range(0,feat1.getGeom().size()):
            single_feat = feat1#.cloneGeom()
            currentid = layer1name+"-"+str(single_feat.getAttributes()[id_index].toString())+"-"+str(i)
            # single_feat.setGeom(feat1.getGeom().get(i))
            # single_feat.setAttribute(0, currentid) # overwrite unused attr fid
            # l1.append(single_feat)

            n = newdb1.nouvelElement(feat1.getGeom().get(i))
            n.setSchema(schema)
            attributes = jpype.JObject[:]@[currentid]
            n.setAttributes(attributes)

            #allfeatures[currentid] = single_feat
    # db1.clear()
    # for f1 in l1:
    #     db1.add(f1)
    newdb1.initSpatialIndex(Tiling.class_, False)
    print("db1 preprocessed")
    newdb2 = Population(False, "db2", DefaultFeature.class_, True)#<DefaultFeature>
    newdb2.setFeatureType(newFeatureType)
    # l2 = list()
    for feat2 in db2:
        for i in range(0,feat2.getGeom().size()):
            single_feat = feat2#.cloneGeom()
            currentid = layer2name+"-"+str(single_feat.getAttributes()[id_index].toString())+"-"+str(i)
            # single_feat.setGeom(feat2.getGeom().get(i))
            # single_feat.setAttribute(0, currentid)
            # l2.append(single_feat)
            n = newdb2.nouvelElement(feat2.getGeom().get(i))
            n.setSchema(schema)
            attributes = jpype.JObject[:]@[currentid]
            n.setAttributes(attributes)

            #allfeatures[currentid] = single_feat
    # db2.clear()
    # for f2 in l2:
    #     db2.add(f2)
    newdb2.initSpatialIndex(Tiling.class_, False)
    print("db2 preprocessed")
    # force cleanup
    del db1
    del db2
    #return(db1, db2, allfeatures)
    return(newdb1, newdb2)

#'
#' FIXME find a way to specify a generic interpretation of matching links
#' FIXME add comparison geometries and semantic: ex height: important for densification
def post_process_links(lienspoly, db1, db2, crs):
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

    for f in lienspoly:
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
            all_link_sources.add(f.getObjetsRef()[0].getAttribute(0))
            all_link_targets.add(f.getObjetsComp()[0].getAttribute(0))

        # 1 -- n : split
        if len(f.getObjetsRef())==1 and len(f.getObjetsComp())>1:
            for comp in f.getObjetsComp():
                features_split.append(comp)
                all_link_sources.add(f.getObjetsRef()[0].getAttribute(0))
                all_link_targets.add(comp.getAttribute(0))

        # m -- 1 : merged
        if len(f.getObjetsRef())>1 and len(f.getObjetsComp())==1:
            for ref in f.getObjetsRef():
                #features_merged.append(ref)
                all_link_sources.add(ref.getAttribute(0))
            # for consistency with "merge ~ aggregation", target feature of the link only exported as evolution
            features_merged.append(f.getObjetsComp()[0])
            all_link_targets.add(f.getObjetsComp()[0].getAttribute(0))

        # m -- n : 20231130: new data model -> merge == agregation ~~aggregation~~
        if len(f.getObjetsRef())>1 and len(f.getObjetsComp())>1:
            for comp in f.getObjetsComp():
                #features_aggregated.append(comp)
                features_merged.append(comp)
                all_link_targets.add(f.getObjetsComp()[0].getAttribute(0))
            for ref in f.getObjetsRef():
                all_link_sources.add(f.getObjetsRef()[0].getAttribute(0))

        # iterate over single links in the multiline
        for i in range(0,f.getGeom().size()):
            link = "LINESTRING ("+", ".join(list(map(lambda p: str(p.getX())+" "+str(p.getY()),f.getGeom().get(i).coord().getList())))+")"
            geoms.append(from_wkt(link))
            #attrs.append(list(f.getSchema().getColonnes())) # no attributes?

    links = geopandas.GeoDataFrame({'geometry':geoms}, crs = crs)

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


    return(links, features_stable, features_split, features_merged, features_aggregated, all_link_targets, all_link_sources, features_disappeared, features_appeared)


def export_links(links, layer1name, layer2name, path, params):
    geojson_export(links, layer1name, layer2name, path, params)
    # export links to shp
    links.to_file('/'.join(path)+'/MATCHING-LINKS_'+layer1name+"_"+layer2name+'.shp')
    links.to_file('/'.join(path)+'/EVOLUTION_'+layer1name+"_"+layer2name+'.gpkg', layer='links', driver="GPKG")


def export(features_appeared, features_disappeared, features_stable, features_split, features_merged, features_aggregated, crs, layer1name, layer2name, path):
    # export

    evol_layer = features_appeared+features_disappeared+features_stable+features_split+features_merged+features_aggregated
    evol_attrs = list(numpy.repeat('appeared',len(features_appeared)))+list(numpy.repeat('disappeared',len(features_disappeared)))+list(numpy.repeat('stable',len(features_stable)))+list(numpy.repeat('split',len(features_split)))+list(numpy.repeat('merged',len(features_merged)))+list(numpy.repeat('aggregated',len(features_aggregated)))

    evol_polys = []
    for x in evol_layer:
        coordinates = []
        for p in x.getGeom().coord().getList():
            coordinates.append((p.getX(), p.getY()))
        polygon = Polygon(coordinates)
        evol_polys.append(polygon)

    evol_ids = [str(x.getAttribute(0)) for x in evol_layer]

    evol = geopandas.GeoDataFrame({'id':evol_ids, 'type':evol_attrs, 'geometry':evol_polys}, crs = crs)

    evol.to_file('/'.join(path)+'/EVOLUTION_'+layer1name+"_"+layer2name+'.shp')
    evol.to_file('/'.join(path)+'/EVOLUTION_'+layer1name+"_"+layer2name+'.gpkg', layer='evolution', driver="GPKG")


def geojson_export(links, layer1name, layer2name, path, params):
    links.to_file('/'.join(path)+'/MATCHING-LINKS_'+layer1name+"_"+layer2name+'.geojson', driver='GeoJSON')

    algoparams = params['algo_params']

    geojson_links = json.loads(links.to_json())
    geojson_links["@context"] = [
        "https://geojson.org/geojson-ld/geojson-context.jsonld",
        "http://www.w3.org/ns/anno.jsonld"]
    geojson_links["description"]="AppariementSurfaces (Atef Bel Hadj Ali)"
    current_dateTime = datetime.now()
    geojson_links["properties"] = {
        "created":str(current_dateTime),
        "Software":"https://github.com/subdense/matching",
        "linking":"https://hal.science/tel-03244834/document",
        "Text":"AppariementSurfaces (Atef Bel Hadj Ali), implemented in https://github.com/IGNF/geoxygene/blob/master/geoxygene-contrib/src/main/java/fr/ign/cogit/geoxygene/contrib/appariement/surfaces/AppariementSurfaces.java",
        "parameters": {
            "surface_min_intersection": algoparams.surface_min_intersection,
            "pourcentage_min_intersection": algoparams.pourcentage_min_intersection,
            "pourcentage_intersection_sur": algoparams.pourcentage_intersection_sur,
            "minimiseDistanceSurfacique": algoparams.minimiseDistanceSurfacique,
            "distSurfMaxFinal": algoparams.distSurfMaxFinal,
            "completudeExactitudeMinFinal": algoparams.completudeExactitudeMinFinal,
            "regroupementOptimal": algoparams.regroupementOptimal,
            "filtrageFinal": algoparams.filtrageFinal,
            "ajoutPetitesSurfaces": algoparams.ajoutPetitesSurfaces,
            "seuilPourcentageTaillePetitesSurfaces": algoparams.seuilPourcentageTaillePetitesSurfaces,
            "persistant": algoparams.persistant,
            "resolutionMin": algoparams.resolutionMin,
            "resolutionMax": algoparams.resolutionMax
        }
    }
    #print(json.dumps(geojson_links, indent=2))
    file_name = '/'.join(path)+'/MATCHING-LINKS_'+layer1name+"_"+layer2name+'-ld.geojson'
    with open(file_name, 'w', encoding='utf-8') as f:
        json.dump(geojson_links, f, ensure_ascii=False, indent=2)
