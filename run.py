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

jpype.startJVM(jpype.getDefaultJVMPath(),f"-Djava.class.path={os.path.dirname(__file__)}/jars/geoxygene-matching-1.10-SNAPSHOT.jar",f"-Xmx{args.java_memory}")

from java.lang import Runtime, Long

import matching

params = matching.get_params(parameter_file=args.parameters.name if args.parameters else None, 
                             layer1=args.layer1.name,
                             layer2=args.layer2.name,
                             crs=args.crs if args.crs else None,
                             attributes=args.attributes if args.attributes else None,
                             output_prefix=args.output_prefix if args.output_prefix else None
                             )

print(str(datetime.datetime.now())+" - reading the input data")
idb1,attributes = matching.get_data(params, "layer1")
idb2,_ = matching.get_data(params, "layer2")
path = [".","output_data"]
crs = CRS.from_user_input(params["crs"])

evol, links = matching.match(idb1,idb2,attributes,params)

# links, features_stable, features_split, features_merged, features_aggregated, features_disappeared, features_appeared = matching.post_process_links(liensPoly, db1, db2, crs, params['id_index'])

matching.export_links(links, path, params)

# matching.export(features_appeared, features_disappeared, features_stable, features_split, features_merged, features_aggregated, crs, path, params)
prefix = params["output_prefix"]
evol.to_file('/'.join(path)+f'/{prefix}_EVOLUTION.gpkg', layer='evolution', driver="GPKG")
print(str(datetime.datetime.now())+" - all done")

jpype.shutdownJVM()
