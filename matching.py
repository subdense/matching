# import requires jpype JVM to already run
import jpype.imports
from jpype.types import *

import json

import jpype
from tqdm import tqdm

from fr.ign.cogit.geoxygene.contrib.appariement.surfaces import ParametresAppSurfaces, AppariementSurfaces
from fr.ign.cogit.geoxygene.feature import SchemaDefaultFeature, DefaultFeature, Population
from fr.ign.cogit.geoxygene.api.spatial.coordgeom import IPolygon
from fr.ign.cogit.geoxygene.schema.schemaConceptuelISOJeu import FeatureType, AttributeType
from java.util import HashMap
from fr.ign.cogit.geoxygene.util.index import STRtreeJts
from fr.ign.cogit.geoxygene.util.conversion import WktGeOxygene

import sys,os

from shapely import from_wkt, Polygon, LineString
from shapely.ops import transform
import geopandas
import pyproj
from pyproj import CRS, Transformer
import numpy
from datetime import datetime

def get_params(parameter_file = None, layer1 = None, layer2 = None, crs = None, attributes = None, output_prefix = None):
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

    params['algo_params']=param
    params['id_index']=id_index
    params['layer1']=layer1
    params['layer2']=layer2
    # if exists('./data/parameters.json'):
    #     with open('./data/parameters.json') as parameter_file:
    #         params_from_file = json.load(parameter_file)
    #         params |= params_from_file
    if parameter_file:
        params_from_file = json.load(parameter_file)
        params |= params_from_file
    if layer1:
        params['layer1'] = layer1
    if layer2:
        params['layer2'] = layer2
    if crs:
        params['crs'] = crs
    if attributes:
        params['attributes'] = json.loads(attributes)
    if output_prefix:
        params['output_prefix'] = output_prefix

        # test data by default: committed in data/bati
    if not layer1:
        layer1 = "./data/bati/neudorf_2012.shp"
    if not layer2:
        layer2 = "./data/bati/neudorf_2022.shp"

    return(params)

#' FIXME handle reprojections
#' FIXME paths must be the same

def get_data(params, layername):
    layer = params[layername]
    attributes = []
    if "attributes" in params: attributes = params["attributes"]
    print(f"{str(datetime.now())} - reading " + layer)
    db = geopandas.read_file(layer, engine="pyogrio", columns = attributes, fid_as_index=True)
    # print(db.head())
    if "crs" in params:
        # print(str(datetime.now())+ " - output CRS = " + str(params["crs"]))
        crs = CRS.from_user_input(params["crs"])
        if crs != db.crs:
            print(str(datetime.now())+" - reprojecting from CRS = " + str(db.crs))
            db = db.to_crs(crs)
    else:
        crs = db.crs
        params['crs'] = crs
        print(str(datetime.now())+" - output CRS from the first source = " + str(crs))
    return db, attributes

def get_data_and_preprocess(params, layername):
    db, attributes = get_data(params, layername)
    newdb = preprocess_layer(layername, db, attributes)
    print(str(datetime.now())+" - preprocess done")
    return newdb

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
    for index, feature in layer.iterrows():
        def addFeature(poly, feature, index):
            n = newdb.nouvelElement(WktGeOxygene.makeGeOxygene(poly.wkt))
            n.setSchema(schema)
            feature_attributes = [index]
            for a in attributes:
                feature_attributes.append(feature[a])
            new_attributes = jpype.JObject[:]@[*feature_attributes]
            n.setAttributes(new_attributes)
        geom = feature["geometry"]
        if geom.geom_type == 'MultiPolygon':
            polygons = list(geom.geoms)
            if len(polygons) == 1:
                addFeature(polygons[0],feature,index)
            else:
                for i in range(0,len(polygons)):
                    currentid = str(index)+"-"+str(i)#layername+"-"+
                    addFeature(polygons[i],feature,currentid)
        else:
            addFeature(geom,feature,index)
    index = STRtreeJts(newdb)
    newdb.setSpatialIndexToExisting(index)
    # print(str(datetime.now())+" - initSpatialIndex done with " + str(newdb.size()))
    return newdb

def match(idb1, idb2, attributes, params):
    def connected_components(idb1,idb2):
        print(str(datetime.now())+" - joining the data (sjoin)")
        join = geopandas.sjoin(idb1,idb2)
        join["index_left"] = join.index.map('L_{}'.format)
        join["index_right"] = join["index_right"].map('R_{}'.format)
        print(f"{str(datetime.now())} - computing the connected components with {len(idb1)} features on the left and {len(idb2)} features on the right")
        import networkx as nx
        G = nx.Graph()
        G.add_nodes_from(idb1.index.map('L_{}'.format))
        G.add_nodes_from(idb2.index.map('R_{}'.format))
        G.add_edges_from(list(map(tuple, join[["index_left","index_right"]].to_numpy())))
        return list(sorted(nx.connected_components(G), reverse=True))
    comp = connected_components(idb1, idb2)
    print(f"{str(datetime.now())} - found {len(comp)} connected components")
    feature_ids = list()
    feature_geoms = list()
    feature_evolution_types = list()
    feature_attributes = dict()
    for a in attributes:
        feature_attributes[a+"_1"] = []
        feature_attributes[a+"_2"] = []
    link_geoms = list()
    link_source_ids = list()
    link_target_ids = list()
    link_evaluation = list()
    for component in tqdm(comp, desc=f"{str(datetime.now())} - Processing connected components", position=0):
        db1_index = []
        db2_index = []
        for n in component:
            if n.startswith("L_"): db1_index.append(int(n[2:]))
            else: db2_index.append(int(n[2:]))
        # select corresponding input features
        c_db1 = preprocess_layer("layer1", idb1.loc[db1_index], attributes)
        # add them to the output data
        # db1.extend(c_db1)
        # select corresponding input features
        c_db2 = preprocess_layer("layer2", idb2.loc[db2_index], attributes)
        # add them to the output data
        # db2.extend(c_db2)
        if (not c_db1.isEmpty() and not c_db2.isEmpty()):
            c_links = AppariementSurfaces.appariementSurfaces(c_db1, c_db2, params["algo_params"])
            # links.extend(c_links)
        else:
            c_links = []
        c_feature_ids, c_feature_geoms, c_feature_evolution_types, c_feature_attributes, c_link_geoms, c_link_source_ids, c_link_target_ids, c_link_eval = process_links(c_links, c_db1, c_db2, attributes)
        feature_ids.extend(c_feature_ids)
        feature_geoms.extend(c_feature_geoms)
        feature_evolution_types.extend(c_feature_evolution_types)
        for a in attributes: 
            feature_attributes[a+"_1"].extend(c_feature_attributes[a+"_1"])
            feature_attributes[a+"_2"].extend(c_feature_attributes[a+"_2"])
        link_geoms.extend(c_link_geoms)
        link_source_ids.extend(c_link_source_ids)
        link_target_ids.extend(c_link_target_ids)
        link_evaluation.extend(c_link_eval)
    links = geopandas.GeoDataFrame({'geometry':link_geoms, 'source_id': link_source_ids, 'target_id':link_target_ids, 'evaluation': link_evaluation}, crs = params["crs"])
    evol = geopandas.GeoDataFrame({'id':feature_ids, 'type':feature_evolution_types, 'geometry':feature_geoms}, crs = params["crs"])
    for a, v in feature_attributes.items(): evol[a] = v
    return evol, links

def make_polygon(geom):
    coordinates = []
    for p in geom.getExterior().coord().getList():
        coordinates.append((p.getX(), p.getY()))
    holes = []
    for h in range(0,geom.sizeInterior()):
        ring = geom.getInterior(h)
        coords = []
        for p in ring.coord().getList():
            coords.append((p.getX(), p.getY()))
        holes.append(coords)
    return Polygon(coordinates,holes)
def make_linestring(geom):
    coordinates = []
    for p in geom.coord().getList():
        coordinates.append((p.getX(), p.getY()))
    return LineString(coordinates)
def to_str(list):
    return ",".join(map(str, list))

def process_links(input_links, db1, db2, attributes):
    # features
    feature_ids = list()
    feature_geoms = list()
    feature_evolution_types = list()
    feature_attributes = dict()
    for a in attributes:
        feature_attributes[a+"_1"] = []
        feature_attributes[a+"_2"] = []
    # links
    all_link_sources = set()
    all_link_targets = set()
    link_geoms = list()
    link_source_ids = list()
    link_target_ids = list()
    link_evaluation = list()
    for f in input_links:
        # 1--1 : stability
        if len(f.getObjetsRef())==1 and len(f.getObjetsComp())==1:
            ref = f.getObjetsRef()[0]
            comp = f.getObjetsComp()[0]
            feature_ids.append(comp.getAttribute(0))
            feature_geoms.append(make_polygon(comp.getGeom()))
            feature_evolution_types.append("stable")
            for a in attributes:
                feature_attributes[a+"_1"].append(to_str([ref.getAttribute(a)]))
                feature_attributes[a+"_2"].append(to_str([comp.getAttribute(a)]))
            all_link_sources.add(ref.getAttribute(0))
            all_link_targets.add(comp.getAttribute(0))
        # 1 -- n : split
        if len(f.getObjetsRef())==1 and len(f.getObjetsComp())>1:
            ref = f.getObjetsRef()[0]
            all_link_sources.add(ref.getAttribute(0))
            for comp in f.getObjetsComp():
                feature_ids.append(comp.getAttribute(0))
                feature_geoms.append(make_polygon(comp.getGeom()))
                feature_evolution_types.append("split")
                for a in attributes:
                    feature_attributes[a+"_1"].append(to_str([ref.getAttribute(a)]))
                    feature_attributes[a+"_2"].append(to_str([comp.getAttribute(a)]))
                all_link_targets.add(comp.getAttribute(0))
        # m -- 1 : merged
        if len(f.getObjetsRef())>1 and len(f.getObjetsComp())==1:
            comp = f.getObjetsComp()[0]
            for ref in f.getObjetsRef():
                all_link_sources.add(ref.getAttribute(0))
            ref_attribute_dict = dict()
            for a in attributes:
                ref_attribute_dict[a] = [ref.getAttribute(a) for ref in f.getObjetsRef()]
            # for consistency with "merge ~ aggregation", target feature of the link only exported as evolution
            feature_ids.append(comp.getAttribute(0))
            feature_geoms.append(make_polygon(comp.getGeom()))
            feature_evolution_types.append("merged")
            for a in attributes:
                feature_attributes[a+"_1"].append(to_str(ref_attribute_dict[a]))
                feature_attributes[a+"_2"].append(to_str([comp.getAttribute(a)]))
            all_link_targets.add(comp.getAttribute(0))
        # m -- n : 20231130: new data model -> merge == agregation ~~aggregation~~
        if len(f.getObjetsRef())>1 and len(f.getObjetsComp())>1:
            for ref in f.getObjetsRef():
                all_link_sources.add(ref.getAttribute(0))
            ref_attribute_dict = dict()
            for a in attributes:
                ref_attribute_dict[a] = [ref.getAttribute(a) for ref in f.getObjetsRef()]
            for comp in f.getObjetsComp():
                # ! use features_merged but not features_aggregated (kept for reversibility)
                all_link_targets.add(comp.getAttribute(0))
                feature_ids.append(comp.getAttribute(0))
                feature_geoms.append(make_polygon(comp.getGeom()))
                feature_evolution_types.append("merged")
                for a in attributes:
                    feature_attributes[a+"_1"].append(to_str(ref_attribute_dict[a]))
                    feature_attributes[a+"_2"].append(to_str([comp.getAttribute(a)]))
        # iterate over single links in the multiline
        for i in range(0,f.getGeom().size()):
            link_geoms.append(make_linestring(f.getGeom().get(i)))
        # links in the multiline have same index than arcs -> construct origin/destination id lists
        for arc in f.getArcs():
            lien = jpype.JObject(arc.getCorrespondant(0), 'fr.ign.cogit.geoxygene.contrib.appariement.Lien')
            link_source_ids.append(str(lien.getObjetsRef()[0].getAttribute(0)))
            link_target_ids.append(str(lien.getObjetsComp()[0].getAttribute(0)))
            link_evaluation.append(lien.getEvaluation())
    for f1 in db1:
        if f1.getAttribute(0) not in all_link_sources:
            feature_ids.append(f1.getAttribute(0))
            feature_geoms.append(make_polygon(f1.getGeom()))
            feature_evolution_types.append("disappeared")
            for a in attributes:
                feature_attributes[a+"_1"].append(to_str([f1.getAttribute(a)]))
                feature_attributes[a+"_2"].append(to_str([]))
    for f2 in db2:
        if f2.getAttribute(0) not in all_link_targets:
            feature_ids.append(f2.getAttribute(0))
            feature_geoms.append(make_polygon(f2.getGeom()))
            feature_evolution_types.append("appeared")
            for a in attributes:
                feature_attributes[a+"_1"].append(to_str([]))
                feature_attributes[a+"_2"].append(to_str([f2.getAttribute(a)]))
    return map(str, feature_ids), feature_geoms, feature_evolution_types, feature_attributes, link_geoms, link_source_ids, link_target_ids, link_evaluation


#' FIXME find a way to specify a generic interpretation of matching links
#' FIXME add comparison geometries and semantic: ex height: important for densification
def post_process_links(lienspoly, db1, db2, crs, id_index):
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

    for f in tqdm(lienspoly,desc=f"{str(datetime.now())} - Post-processing links", position=0):
        # 1--1 : stability
        if len(f.getObjetsRef())==1 and len(f.getObjetsComp())==1:
            ref = f.getObjetsRef()[0]
            comp = f.getObjetsComp()[0]
            features_stable.append(comp)
            all_link_sources.add(ref.getAttribute(id_index))
            all_link_targets.add(comp.getAttribute(id_index))
            # HACK - by JP: to get attributes for m-n links, use this internal Geoxygene structure to store objects
            comp.clearCorrespondants()
            comp.addCorrespondant(ref)

        # 1 -- n : split
        if len(f.getObjetsRef())==1 and len(f.getObjetsComp())>1:
            ref = f.getObjetsRef()[0]
            for comp in f.getObjetsComp():
                features_split.append(comp)
                all_link_sources.add(ref.getAttribute(id_index))
                all_link_targets.add(comp.getAttribute(id_index))
                # HACK
                comp.clearCorrespondants()
                comp.addCorrespondant(ref)

        # m -- 1 : merged
        if len(f.getObjetsRef())>1 and len(f.getObjetsComp())==1:
            comp = f.getObjetsComp()[0]
            for ref in f.getObjetsRef():
                #features_merged.append(ref)
                all_link_sources.add(ref.getAttribute(id_index))
            # for consistency with "merge ~ aggregation", target feature of the link only exported as evolution
            features_merged.append(comp)
            all_link_targets.add(comp.getAttribute(id_index))
            # HACK
            comp.clearCorrespondants()
            for ref in f.getObjetsRef():
                comp.addCorrespondant(ref)

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
            source_ids.append(str(lien.getObjetsRef()[0].getAttribute(id_index)))
            target_ids.append(str(lien.getObjetsComp()[0].getAttribute(id_index)))

    links = geopandas.GeoDataFrame({'geometry':geoms, 'source_id': source_ids, 'target_id':target_ids}, crs = crs)

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

    return(links, features_stable, features_split, features_merged, features_aggregated, features_disappeared, features_appeared)


def export_links(layer1, layer2, links, path, params, arguments):
    output_dir = '/'.join(path)
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    geojson_export(links, path, params)
    # export links to shp
    # links.to_file('/'.join(path)+'/MATCHING-LINKS_'+layer1name+"_"+layer2name+'.shp')
    # export links to GPKG
    # WARNING: does not work yet https://github.com/geopandas/geopandas/pull/2850
    algoparams = params['algo_params']
    metadata = {
        "title": "Matching links",
        "description": "AppariementSurfaces (Atef Bel Hadj Ali)",
        "created": str(datetime.now()),
        "Software":"https://github.com/subdense/matching",
        "linking":"https://hal.science/tel-03244834/document",
        "Abstract":"AppariementSurfaces (Atef Bel Hadj Ali), implemented in https://github.com/IGNF/geoxygene/blob/master/geoxygene-contrib/src/main/java/fr/ign/cogit/geoxygene/contrib/appariement/surfaces/AppariementSurfaces.java",
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
        },
        "Command line arguments": ' '.join(arguments)
    }
    links.to_file('/'.join(path)+f'/{params["output_prefix"]}_MATCHING-LINKS.gpkg', layer='links', driver="GPKG", metadata=metadata)
    layer1.to_file('/'.join(path)+f'/{params["output_prefix"]}_MATCHING-LINKS.gpkg', layer='layer1', driver="GPKG")
    layer2.to_file('/'.join(path)+f'/{params["output_prefix"]}_MATCHING-LINKS.gpkg', layer='layer2', driver="GPKG")
    with open('/'.join(path)+f'/{params["output_prefix"]}_MATCHING-LINKS_metadata.json', 'w', encoding='utf-8') as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2)

def export(features_appeared, features_disappeared, features_stable, features_split, features_merged, features_aggregated, crs, path, params):
    # export

    evol_layer = features_appeared+features_disappeared+features_stable+features_split+features_merged+features_aggregated
    evol_attrs = list(numpy.repeat('appeared',len(features_appeared)))+list(numpy.repeat('disappeared',len(features_disappeared)))+list(numpy.repeat('stable',len(features_stable)))+list(numpy.repeat('split',len(features_split)))+list(numpy.repeat('merged',len(features_merged)))+list(numpy.repeat('aggregated',len(features_aggregated)))

    attributes = dict()
    if "attributes" in params:
        for a in params["attributes"]:
            attributes[a+"_1"] = []
            attributes[a+"_2"] = []
    # print("exporting attributes " + str(attributes.keys()))
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
    for x in tqdm(evol_layer,desc=f"{str(datetime.now())} - Exporting features",position=0):
        evol_polys.append(make_polygon(x.getGeom()))

    evol_ids = [str(x.getAttribute(0)) for x in evol_layer]

    evol = geopandas.GeoDataFrame({'id':evol_ids, 'type':evol_attrs, 'geometry':evol_polys}, crs = crs)
    for a, v in attributes.items():
        evol[a] = v

    prefix = params["output_prefix"]
    # evol.to_file('/'.join(path)+f'/{prefix}_EVOLUTION.shp')
    evol.to_file('/'.join(path)+f'/{prefix}_EVOLUTION.gpkg', layer='evolution', driver="GPKG")


def geojson_export(links, path, params):
    links.to_file('/'.join(path)+f'/{params["output_prefix"]}_MATCHING-LINKS.geojson', driver='GeoJSON')

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
    file_name = '/'.join(path)+f'/{params["output_prefix"]}_MATCHING-LINKS-ld.geojson'
    with open(file_name, 'w', encoding='utf-8') as f:
        json.dump(geojson_links, f, ensure_ascii=False, indent=2)
