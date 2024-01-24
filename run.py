import jpype.imports
from jpype.types import *
import os
import datetime
import argparse
from tqdm import tqdm

from pyproj import CRS

def print_memory():
    # Total number of processors or cores available to the JVM
    print("Available processors (cores): " + str(Runtime.getRuntime().availableProcessors()))

    # Total amount of free memory available to the JVM
    print("Free memory (bytes): " + str(Runtime.getRuntime().freeMemory()))

    # This will return Long.MAX_VALUE if there is no preset limit
    maxMemory = Runtime.getRuntime().maxMemory()
    # Maximum amount of memory the JVM will attempt to use
    print("Maximum memory (bytes): " + ("no limit" if maxMemory == Long.MAX_VALUE else str(maxMemory)))

    # Total memory currently in use by the JVM
    print("Total memory (bytes): " + str(Runtime.getRuntime().totalMemory()))

parser = argparse.ArgumentParser(description='Process the matching of two datasets.')
parser.add_argument('-parameters', type = argparse.FileType('r'), required=False, help='a json parameter file')
parser.add_argument('-layer1', type = argparse.FileType('r'), required=True, help='layer1')
parser.add_argument('-layer2', type = argparse.FileType('r'), required=True, help='layer2')
parser.add_argument('-crs', type=str, required=False, help='crs')
parser.add_argument('-attributes', type=str, required=False, help='attributes as json array string')
parser.add_argument('-output_prefix', type=str, required=True, help='prefix for output files')
parser.add_argument('-java_memory', type=str, required=False, default = "16g", help='java memory to allocate')

args = parser.parse_args()

if jpype.isJVMStarted():
    jpype.shutdownJVM()

#jpype.startJVM(jpype.getDefaultJVMPath(),f"-Djava.class.path={os.path.dirname(__file__)}/jars/geoxygene-matching-1.10-SNAPSHOT.jar")
jpype.startJVM(jpype.getDefaultJVMPath(),f"-Djava.class.path={os.path.dirname(__file__)}/jars/geoxygene-matching-1.10-SNAPSHOT.jar",f"-Xmx{args.java_memory}")

from fr.ign.cogit.geoxygene.contrib.appariement.surfaces import AppariementSurfaces
from java.lang import Runtime, Long

import matching

params = matching.get_params(parameter_file=args.parameters.name if args.parameters else None, 
                             layer1=args.layer1.name,
                             layer2=args.layer2.name,
                             crs=args.crs if args.crs else None,
                             attributes=args.attributes if args.attributes else None,
                             output_prefix=args.output_prefix if args.output_prefix else None
                             )

# print_memory()
print(str(datetime.datetime.now())+" - reading the input data")
# db1 = matching.get_data_and_preprocess(params, "layer1")
# db2 = matching.get_data_and_preprocess(params, "layer2")
idb1,attributes = matching.get_data(params, "layer1")
idb2,_ = matching.get_data(params, "layer2")
import geopandas
print(str(datetime.datetime.now())+" - joining the data (sjoin)")
join = geopandas.sjoin(idb1,idb2)
join["index_left"] = join.index.map('L_{}'.format)
join["index_right"] = join["index_right"].map('R_{}'.format)
print(f"{str(datetime.datetime.now())} - computing the connected components with {len(idb1)} features on the left and {len(idb2)} features on the right")
import networkx as nx
G = nx.Graph()
G.add_nodes_from(idb1.index.map('L_{}'.format))
G.add_nodes_from(idb2.index.map('R_{}'.format))
G.add_edges_from(list(map(tuple, join[["index_left","index_right"]].to_numpy())))
comp = list(nx.connected_components(G))
print(f"{str(datetime.datetime.now())} - found {len(comp)} connected components")
db1 = []
db2 = []
liensPoly = []
for component in tqdm(comp, desc=f"{str(datetime.datetime.now())} - Processing connected components", position=0):
    db1_index = []
    db2_index = []
    for n in component:
        if n.startswith("L_"):
            db1_index.append(int(n[2:]))
        else:
            db2_index.append(int(n[2:]))
    # print("L = "+",".join(map(str,db1_index)))
    c_db1 = idb1.loc[db1_index]
    # print("R = "+",".join(map(str,db2_index)))
    c_db2 = idb2.loc[db2_index]
    newdb1 = matching.preprocess_layer("layer1", c_db1, attributes)
    newdb2 = matching.preprocess_layer("layer2", c_db2, attributes)
    db1.extend(newdb1)
    db2.extend(newdb2)
    if (not newdb1.isEmpty() and not newdb2.isEmpty()):
        c_links = AppariementSurfaces.appariementSurfaces(newdb1, newdb2, params['algo_params'])
        liensPoly.extend(c_links)
path = [".","output_data"]
crs = CRS.from_user_input(params["crs"])

# print(str(datetime.datetime.now())+" - preprocess done")
# print_memory()
# call to geoxygene matching algorithm
# from org.apache.logging.log4j import Level
# AppariementSurfaces.LOGGER.setLevel(Level.DEBUG)
# liensPoly = AppariementSurfaces.appariementSurfaces(db1, db2, params['algo_params'])
# print(str(datetime.datetime.now())+" - matching done")
# print_memory()
links, features_stable, features_split, features_merged, features_aggregated, all_link_targets, all_link_sources, features_disappeared, features_appeared = matching.post_process_links(liensPoly, db1, db2, crs, params['id_index'])
# print(str(datetime.datetime.now())+" - post_process_links done")
# print_memory()

matching.export_links(links, path, params)
# print(str(datetime.datetime.now())+" - export_links done")
# print_memory()

matching.export(features_appeared, features_disappeared, features_stable, features_split, features_merged, features_aggregated, crs, path, params)
print(str(datetime.datetime.now())+" - all done")
# print_memory()

jpype.shutdownJVM()
