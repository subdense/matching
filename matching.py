# import requires jpype JVM to already run
import jpype.imports
from jpype.types import *

import json

import jpype

from fr.ign.cogit.geoxygene.contrib.appariement.surfaces import ParametresAppSurfaces, AppariementSurfaces
from fr.ign.cogit.geoxygene.util.conversion import ShapefileReader
from fr.ign.cogit.geoxygene.feature import SchemaDefaultFeature, DefaultFeature, Population
from fr.ign.cogit.geoxygene.api.spatial.coordgeom import IPolygon
from fr.ign.cogit.geoxygene.schema.schemaConceptuelISOJeu import FeatureType, AttributeType
from java.util import HashMap
from fr.ign.cogit.geoxygene.util.index import Tiling, STRtreeJts
from fr.ign.cogit.geoxygene.util.conversion import WktGeOxygene

import sys,os
from os.path import exists

from shapely import from_wkt, Polygon
from shapely.ops import transform
import geopandas
import pyproj
from pyproj import CRS, Transformer
import numpy
from datetime import datetime
from progress.bar import Bar

def get_params():
    params = dict()

    # Parameters of the matching algorithm
    param = ParametresAppSurfaces()
    # parameters for pre-matching
    param.surface_min_intersection = 1
    param.pourcentage_min_intersection = 0.2
    param.minimiseDistanceSurfacique = True
    # parameters for 'optimal' matching
    param.regroupementOptimal = True
    param.pourcentage_intersection_sur = 0.25
    # parameters for adding small areas
    param.ajoutPetitesSurfaces = True
    param.seuilPourcentageTaillePetitesSurfaces = 0.1
    # parameters for final filtering
    param.filtrageFinal = True
    param.distSurfMaxFinal = 0.6
    param.completudeExactitudeMinFinal = 0.3
    # parameter for db persistence
    param.persistant = False
    # parameters for 'robust' polygon operations (intersection, etc.)
    param.resolutionMin = 1
    param.resolutionMax = 11

    # param: index of IDs (same for both layers)
    id_index = 0

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
            params |= params_from_file

    if len(sys.argv)>=3:
        params['layer1'] = sys.argv[1]
        params['layer2'] = sys.argv[2]
    if len(sys.argv)>=4:
        params['crs'] = sys.argv[3]
    if len(sys.argv)>=5:
        print(sys.argv[4])
        params['attributes'] = json.loads(sys.argv[4])

    return(params)

#' FIXME handle reprojections
#' FIXME paths must be the same

def get_data_and_preprocess(params, layer):
    layer = params[layer]

    path = layer.split("/")
    layername = os.path.splitext(path[len(path)-1])[0]

    path.pop() # dirty python mutables
    path = [".","output_data"]

    attributes = []
    if "attributes" in params:
        attributes = params["attributes"]

    print(str(datetime.now())+ " - READ " + layer)
    db = geopandas.read_file(layer, engine="pyogrio", columns = attributes, fid_as_index=True)
    print(db.head())
    if "crs" in params:
        print(str(datetime.now())+ " - output CRS = " + str(params["crs"]))
        crs = CRS.from_user_input(params["crs"])
        if crs != db.crs:
            print(str(datetime.now())+" - reprojecting from CRS = " + str(db.crs))
            db = db.to_crs(crs)
    else:
        crs = db.crs
        params['crs'] = crs
        print(str(datetime.now())+" - output CRS from the first source = " + str(crs))
    # get the list of feature attributes for db1
    # print(str(datetime.datetime.now())+" - DB attributes = " + str(db.getFeatureType().getFeatureAttributes()))
    newdb = preprocess_layer(layername, db, attributes)
    print(str(datetime.now())+" - preprocess done")
    return (layername, path, newdb, crs)

def preprocess_layer(layername, layer, attributes):
    newFeatureType = FeatureType()
    newFeatureType.setTypeName("building")
    newFeatureType.setGeometryType(IPolygon.class_)
    id = AttributeType("id", "String")
    newFeatureType.addFeatureAttribute(id)
    schema = SchemaDefaultFeature()
    schema.setFeatureType(newFeatureType)
    newFeatureType.setSchema(schema)
    attLookup = HashMap()
    attLookup.put(jpype.JInt(0), jpype.JString[:]@[id.getNomField(), id.getMemberName()])
    for index, a in enumerate(attributes):
        att = AttributeType(a, "String")
        newFeatureType.addFeatureAttribute(att)
        attLookup.put(jpype.JInt(index+1), jpype.JString[:]@[att.getNomField(), att.getMemberName()])
    schema.setAttLookup(attLookup)
    newdb = Population(False, layername, DefaultFeature.class_, True)
    newdb.setFeatureType(newFeatureType)
    bar = Bar('Processing', max=len(layer))
    for index, feature in layer.iterrows():
        def addFeature(poly, attribute_values):
            n = newdb.nouvelElement(WktGeOxygene.makeGeOxygene(poly.wkt))
            n.setSchema(schema)
            attributes = jpype.JObject[:]@[*attribute_values]
            n.setAttributes(attributes)
        geom = feature["geometry"]
        feature_attributes = [index]
        for a in attributes:
            feature_attributes.append(feature[a])
        if geom.geom_type == 'MultiPolygon':
            polygons = list(geom.geoms)
            if len(polygons) == 1:
                addFeature(polygons[0],feature_attributes)
            else:
                for i in range(0,len(polygons)):
                    currentid = str(index)+"-"+str(i)#layername+"-"+
                    addFeature(polygons[i],currentid)
        else:
            addFeature(geom,feature_attributes)
        bar.next()
    bar.finish()
    index = STRtreeJts(newdb)
    newdb.setSpatialIndexToExisting(index)
    print(str(datetime.now())+" - initSpatialIndex done with " + str(newdb.size()))
    return newdb

#'
#' FIXME find a way to specify a generic interpretation of matching links
#' FIXME add comparison geometries and semantic: ex height: important for densification
def post_process_links(lienspoly, db1, db2, crs, layer1name, layer2name, id_index):
    # construct geopandas data frame
    attrs = list()
    geoms = list()
    source_ids = list()
    target_ids = list()

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
            all_link_sources.add(f.getObjetsRef()[0].getAttribute(id_index))
            all_link_targets.add(f.getObjetsComp()[0].getAttribute(id_index))
            # HACK - by JP: to get attributes for m-n links, use this internal Geoxygene structure to store objects
            f.getObjetsComp()[0].clearCorrespondants()
            f.getObjetsComp()[0].addCorrespondant(f.getObjetsComp()[0])

        # 1 -- n : split
        if len(f.getObjetsRef())==1 and len(f.getObjetsComp())>1:
            for comp in f.getObjetsComp():
                features_split.append(comp)
                all_link_sources.add(f.getObjetsRef()[0].getAttribute(id_index))
                all_link_targets.add(comp.getAttribute(id_index))
                # HACK
                comp.clearCorrespondants()
                comp.addCorrespondant(f.getObjetsRef()[0])

        # m -- 1 : merged
        if len(f.getObjetsRef())>1 and len(f.getObjetsComp())==1:
            for ref in f.getObjetsRef():
                #features_merged.append(ref)
                all_link_sources.add(ref.getAttribute(id_index))
            # for consistency with "merge ~ aggregation", target feature of the link only exported as evolution
            features_merged.append(f.getObjetsComp()[0])
            all_link_targets.add(f.getObjetsComp()[0].getAttribute(id_index))
            # HACK
            f.getObjetsComp()[0].clearCorrespondants()
            f.getObjetsComp()[0].addCorrespondant(f.getObjetsComp()[0])

        # m -- n : 20231130: new data model -> merge == agregation ~~aggregation~~
        if len(f.getObjetsRef())>1 and len(f.getObjetsComp())>1:
            # debug: neudorf: (2012 : BATIMENT0000000048824441 BATIMENT0000000048824481) -> (2022 : BATIMENT0000002223308372 BATIMENT0000002223308719 BATIMENT0000002223308720 BATIMENT0000002223308838)
            #if 'BATIMENT0000000048824441' in (ref.getAttribute(1) for ref in f.getObjetsRef()):
            #    print("Ref : ")
            #    print([ref.getAttribute(1) for ref in f.getObjetsRef()])
            #    print("Comp : ")
            #    print([ref.getAttribute(1) for ref in f.getObjetsComp()])
            for comp in f.getObjetsComp():
                # ! use features_merged but not features_aggregated (kept for reversibility)
                features_merged.append(comp)
                all_link_targets.add(comp.getAttribute(id_index))
                # HACK
                comp.clearCorrespondants()
                for ref in f.getObjetsRef():
                    comp.addCorrespondant(ref)
            for ref in f.getObjetsRef():
                all_link_sources.add(ref.getAttribute(id_index))

        # iterate over single links in the multiline
        for i in range(0,f.getGeom().size()):
            #if 'BATIMENT0000000048824441' in (ref.getAttribute(1) for ref in f.getObjetsRef()): # DEBUG
            #    print(f.getGeom().get(i))
            link = "LINESTRING ("+", ".join(list(map(lambda p: str(p.getX())+" "+str(p.getY()),f.getGeom().get(i).coord().getList())))+")"
            geoms.append(from_wkt(link))
            #attrs.append(list(f.getSchema().getColonnes())) # no attributes?

        # links in the multiline have same index than arcs -> construct origin/destination id lists
        for arc in f.getArcs():
            lien = jpype.JObject(arc.getCorrespondant(0), 'fr.ign.cogit.geoxygene.contrib.appariement.Lien')
            source_ids.append(layer1name+':'+str(lien.getObjetsRef()[0].getAttribute(id_index)))
            target_ids.append(layer2name+':'+str(lien.getObjetsComp()[0].getAttribute(id_index)))

    links = geopandas.GeoDataFrame({'geometry':geoms, 'source_id': source_ids, 'target_ids':target_ids}, crs = crs)

    # construct appeared/disappeared layer
    # 1 -- 0
    features_disappeared = list()
    for f1 in db1:
        if f1.getAttribute(id_index) not in all_link_sources:
            #print(f1.getAttribute(0))
            features_disappeared.append(f1)
            # HACK
            f1.clearCorrespondants()
    # 0 -- 1
    features_appeared = list()
    for f2 in db2:
        if f2.getAttribute(id_index) not in all_link_targets:
            features_appeared.append(f2)
            # HACK
            f2.clearCorrespondants()

    print(str(datetime.now())+" - " + str(len(links)) + " processed")
    return(links, features_stable, features_split, features_merged, features_aggregated, all_link_targets, all_link_sources, features_disappeared, features_appeared)


def export_links(links, layer1name, layer2name, path, params):
    output_dir = '/'.join(path)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    geojson_export(links, layer1name, layer2name, path, params)
    # export links to shp
    links.to_file('/'.join(path)+'/MATCHING-LINKS_'+layer1name+"_"+layer2name+'.shp')
    links.to_file('/'.join(path)+'/MATCHING-LINKS_'+layer1name+"_"+layer2name+'.gpkg', layer='links', driver="GPKG")

def export(features_appeared, features_disappeared, features_stable, features_split, features_merged, features_aggregated, crs, layer1name, layer2name, path, params):
    # export

    evol_layer = features_appeared+features_disappeared+features_stable+features_split+features_merged+features_aggregated
    evol_attrs = list(numpy.repeat('appeared',len(features_appeared)))+list(numpy.repeat('disappeared',len(features_disappeared)))+list(numpy.repeat('stable',len(features_stable)))+list(numpy.repeat('split',len(features_split)))+list(numpy.repeat('merged',len(features_merged)))+list(numpy.repeat('aggregated',len(features_aggregated)))

    attributes = dict()
    if "attributes" in params:
        for a in params["attributes"]:
            attributes[a+"_1"] = []
            attributes[a+"_2"] = []
    print("exporting attributes " + str(attributes.keys()))
    def to_str(list):
        return ",".join(map(str, list))

    #TODO clean than up
    if "attributes" in params:
        for x in features_appeared:
            for a in params["attributes"]:
                attributes[a+"_1"].append(to_str([]))
                attributes[a+"_2"].append(to_str([x.getAttribute(a)]))
        for x in features_disappeared:
            for a in params["attributes"]:
                attributes[a+"_1"].append(to_str([x.getAttribute(a)]))
                attributes[a+"_2"].append(to_str([]))
        for x in features_stable:
            for a in params["attributes"]:
                attributes[a+"_1"].append(to_str([x.getCorrespondant(0).getAttribute(a)]))
                attributes[a+"_2"].append(to_str([x.getAttribute(a)]))
        for x in features_split:
            for a in params["attributes"]:
                attributes[a+"_1"].append(to_str([x.getCorrespondant(0).getAttribute(a)]))
                attributes[a+"_2"].append(to_str([x.getAttribute(a)]))
        for x in features_merged:
            for a in params["attributes"]:
                attribute_values = []
                for c in range(0,x.getCorrespondants().size()):
                    attribute_values.append(x.getCorrespondant(c).getAttribute(a))
                attributes[a+"_1"].append(to_str(attribute_values))
                attributes[a+"_2"].append(to_str([x.getAttribute(a)]))
        #for x in features_aggregated:
        #    for a in params["attributes"]:
        #        attribute_values = []
        #        for c in range(0,x.getCorrespondants().size()):
        #            attribute_values.append(x.getCorrespondant(c).getAttribute(a))
        #        attributes[a+"_1"].append(to_str(attribute_values))
        #        attributes[a+"_2"].append(to_str([x.getAttribute(a)]))

    evol_polys = []
    for x in evol_layer:
        coordinates = []
        # TODO handle holes and not just the shell (exterior)
        for p in x.getGeom().getExterior().coord().getList():
            coordinates.append((p.getX(), p.getY()))
        holes = []
        for h in range(0,x.getGeom().sizeInterior()):
            ring = x.getGeom().getInterior(h)
            coords = []
            for p in ring.coord().getList():
                coords.append((p.getX(), p.getY()))
            holes.append(coords)
        polygon = Polygon(coordinates,holes)
        evol_polys.append(polygon)

    evol_ids = [str(x.getAttribute(0)) for x in evol_layer]

    evol = geopandas.GeoDataFrame({'id':evol_ids, 'type':evol_attrs, 'geometry':evol_polys}, crs = crs)
    for a, v in attributes.items():
        evol[a] = v

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
