import jpype.imports
from jpype.types import *
import os
import datetime

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

if jpype.isJVMStarted():
    jpype.shutdownJVM()

#jpype.startJVM(jpype.getDefaultJVMPath(),f"-Djava.class.path={os.path.dirname(__file__)}/jars/geoxygene-matching-1.10-SNAPSHOT.jar")
jpype.startJVM(jpype.getDefaultJVMPath(),f"-Djava.class.path={os.path.dirname(__file__)}/jars/geoxygene-matching-1.10-SNAPSHOT.jar","-Xmx16g")

from fr.ign.cogit.geoxygene.contrib.appariement.surfaces import AppariementSurfaces
from java.lang import Runtime, Long

import matching

params = matching.default_params()

print_memory()
print(str(datetime.datetime.now())+" - getdata")
layer1name, path, db1, crs = matching.get_data_and_preprocess(params, "layer1")
layer2name, path_, db2, crs_ = matching.get_data_and_preprocess(params, "layer2")
# layer1name, layer2name, path, db1, db2, crs = matching.get_data(params)

# print(str(datetime.datetime.now())+" - preprocess")
# db1 = matching.preprocess_data(layer1name, db1, id_index = 1)
# print(str(datetime.datetime.now())+" - preprocess 1 done")
# db2 = matching.preprocess_data(layer2name, db2, id_index = 1)
# print(str(datetime.datetime.now())+" - preprocess 2 done")
#db1, db2, allfeatures = matching.preprocess_data(layer1name, layer2name, db1, db2)
# db1, db2 = matching.preprocess_data(layer1name, layer2name, db1, db2, id_index = params['id_index'])

print(str(datetime.datetime.now())+" - preprocess done")
print_memory()
# call to geoxygene matching algorithm
liensPoly = AppariementSurfaces.appariementSurfaces(db1, db2, params['algo_params'])
print(str(datetime.datetime.now())+" - matching done")
print_memory()
links, features_stable, features_split, features_merged, features_aggregated, all_link_targets, all_link_sources, features_disappeared, features_appeared = matching.post_process_links(liensPoly, db1, db2, crs, layer1name, layer2name, params['id_index'])
print(str(datetime.datetime.now())+" - post_process_links done")
print_memory()

matching.export_links(links, layer1name, layer2name, path, params)
print(str(datetime.datetime.now())+" - export_links done")
print_memory()

matching.export(features_appeared, features_disappeared, features_stable, features_split, features_merged, features_aggregated, crs, layer1name, layer2name, path)
print(str(datetime.datetime.now())+" - export done")
print_memory()

jpype.shutdownJVM()
