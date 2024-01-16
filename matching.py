# import requires jpype JVM to already run
import jpype.imports
from jpype.types import *

import json
from datetime import datetime

from fr.ign.cogit.geoxygene.contrib.appariement.surfaces import ParametresAppSurfaces, AppariementSurfaces
from fr.ign.cogit.geoxygene.util.conversion import ShapefileReader


import sys,os

from shapely import from_wkt, Polygon
import geopandas
import numpy

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

    #layer1 = "./data/bati/bati_95430.shp"
    #layer2 = "./data/bati/cadastre_bati_95430.shp"
    #layer1 = "../../../Data/Test/Neudorf/valid_bati12.shp"
    #layer2 = "../../../Data/Test/Neudorf/valid_bati22.shp"
    #layer1 = "./data/bati/valid_strasbourg_2012.shp"
    #layer2 = "./data/bati/valid_strasbourg_2022.shp"
    # test data by default: committed in data/bati
    layer1 = "./data/bati/neudorf_2012.shp"
    layer2 = "./data/bati/neudorf_2022.shp"

    params['algo_params']=param
    params['id_index']=id_index
    params['layer1']=layer1
    params['layer2']=layer2

    return(params)


#'
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

    db1 = ShapefileReader.read(layer1, True)
    db2 = ShapefileReader.read(layer2, True)
    #print(db1.get(0).getAttributes())
    #print(db1.get(0).getSchema().getColonnes()) # empty
    #print(db1.getFeatureType().getFeatureAttributes())
    #print(db1.getSchema().getColonnes())

    #crs = geopandas.read_file(layer1, engine="pyogrio").crs
    crs = geopandas.read_file(layer1, rows = 1).crs

    return(layer1name, layer2name, path, db1, db2, crs)


def preprocess_data(layer1name, layer2name, db1, db2, id_index):
    # hashmap of all features
    #allfeatures = dict()

    # reduce multipolygons into single polygons
    l1 = list()
    for feat1 in db1:
        for i in range(0,feat1.getGeom().size()):
            single_feat = feat1.cloneGeom()
            currentid = layer1name+"-"+str(single_feat.getAttributes()[id_index].toString())+"-"+str(i)
            single_feat.setGeom(feat1.getGeom().get(i))
            single_feat.setAttribute(0, currentid) # overwrite unused attr fid
            l1.append(single_feat)
            #allfeatures[currentid] = single_feat
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
            #allfeatures[currentid] = single_feat
    db2.clear()
    for f2 in l2:
        db2.add(f2)

    #return(db1, db2, allfeatures)
    return(db1, db2)

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

        # 1 -- n : split
        if len(f.getObjetsRef())==1 and len(f.getObjetsComp())>1:
            for comp in f.getObjetsComp():
                features_split.append(comp)
                all_link_sources.add(f.getObjetsRef()[0].getAttribute(id_index))
                all_link_targets.add(comp.getAttribute(id_index))

        # m -- 1 : merged
        if len(f.getObjetsRef())>1 and len(f.getObjetsComp())==1:
            for ref in f.getObjetsRef():
                #features_merged.append(ref)
                all_link_sources.add(ref.getAttribute(id_index))
            # for consistency with "merge ~ aggregation", target feature of the link only exported as evolution
            features_merged.append(f.getObjetsComp()[0])
            all_link_targets.add(f.getObjetsComp()[0].getAttribute(id_index))

        # m -- n : 20231130: new data model -> merge == agregation ~~aggregation~~
        if len(f.getObjetsRef())>1 and len(f.getObjetsComp())>1:
            # debug: neudorf: (2012 : BATIMENT0000000048824441 BATIMENT0000000048824481) -> (2022 : BATIMENT0000002223308372 BATIMENT0000002223308719 BATIMENT0000002223308720 BATIMENT0000002223308838)
            #if 'BATIMENT0000000048824441' in (ref.getAttribute(1) for ref in f.getObjetsRef()):
            #    print("Ref : ")
            #    print([ref.getAttribute(1) for ref in f.getObjetsRef()])
            #    print("Comp : ")
            #    print([ref.getAttribute(1) for ref in f.getObjetsComp()])
            for comp in f.getObjetsComp():
                features_merged.append(comp)
                all_link_targets.add(comp.getAttribute(id_index))
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
    # 0 -- 1
    features_appeared = list()
    for f2 in db2:
        if f2.getAttribute(id_index) not in all_link_targets:
            features_appeared.append(f2)


    return(links, features_stable, features_split, features_merged, features_aggregated, all_link_targets, all_link_sources, features_disappeared, features_appeared)


def export_links(links, layer1name, layer2name, path, params):
    geojson_export(links, layer1name, layer2name, path, params)
    # export links to shp
    links.to_file('/'.join(path)+'/MATCHING-LINKS_'+layer1name+"_"+layer2name+'.shp')
    links.to_file('/'.join(path)+'/MATCHING-LINKS_'+layer1name+"_"+layer2name+'.gpkg', layer='links', driver="GPKG")


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
