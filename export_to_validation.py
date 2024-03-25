import os
from datetime import datetime
import argparse
from tqdm import tqdm
from pyproj import CRS
import geopandas
import networkx as nx
import json
import shapely
import random
from pyproj import CRS

parser = argparse.ArgumentParser(description='Export matches for validation.')
parser.add_argument('-layer1', type = argparse.FileType('r'), required=True, help='layer1')
parser.add_argument('-layer2', type = argparse.FileType('r'), required=True, help='layer2')
parser.add_argument('-date1', type = str, required=True, help='date1')
parser.add_argument('-date2', type = str, required=True, help='date2')
parser.add_argument('-wmts1', type = str, required=True, help='wmts1')
parser.add_argument('-wmts2', type = str, required=True, help='wmts2')
parser.add_argument('-output_directory', type=str, required=True, help='directory for output files')

args = parser.parse_args()

print(str(datetime.now())+" - reading the input data")
crs = CRS.from_epsg(4326)
idb1 = geopandas.read_file(args.layer1.name, engine="pyogrio", columns=[], fid_as_index=True).to_crs(crs)
idb2 = geopandas.read_file(args.layer2.name, engine="pyogrio", columns=[], fid_as_index=True).to_crs(crs)
path = [".","output_data"]
directory = args.output_directory
date1 = args.date1
date2 = args.date2
wmts1 = args.wmts1
wmts2 = args.wmts2
print(f"{str(datetime.now())} - computing connected components")
def connected_components(idb1,idb2):
    print(str(datetime.now())+" - joining the data (sjoin)")
    join = geopandas.sjoin(idb1,idb2)
    G = nx.Graph()
    # add prefix to node ids to be able to separate them later
    join["index_left"] = join.index.map('L_{}'.format)
    join["index_right"] = join["index_right"].map('R_{}'.format)
    G.add_nodes_from(idb1.index.map('L_{}'.format))
    G.add_nodes_from(idb2.index.map('R_{}'.format))
    G.add_edges_from(list(map(tuple, join[["index_left","index_right"]].to_numpy())))
    return list(sorted(nx.connected_components(G), reverse=True))
comp = connected_components(idb1, idb2)
# select 1000 samples randomly
comp = random.sample(comp, 1000)
print(f"{str(datetime.now())} - found {len(comp)} connected components")
os.makedirs(f'./output_data/{directory}', exist_ok=True)
fileList = []
for index, component in enumerate(tqdm(comp, desc=f"{str(datetime.now())} - Processing connected components", position=0)):
    features = []
    db1_index = []
    db2_index = []
    # use prefix to separate the nodes into db1 and db2
    for n in component:
        if n.startswith("L_"): db1_index.append(int(n[2:]))
        else: db2_index.append(int(n[2:]))
    def createFeature(f,d):
        return {
            "type":"Feature",
            "properties":{"date":d},
            "geometry": json.loads(shapely.to_geojson(f["geometry"]))
        }
    # select corresponding input features and add them to the output data
    for _, feature in idb1.loc[db1_index].iterrows(): features.append(createFeature(feature,date1))
    for _,feature in idb2.loc[db2_index].iterrows(): features.append(createFeature(feature,date2))
    geojson = {
        "type":"FeatureCollection",
        "features": features
    }
    with open('/'.join(path)+f'/{directory}/data_{index}.geojson', 'w', encoding='utf-8') as f:
        json.dump(geojson, f, ensure_ascii=False, indent=2)
    fileList.append(f'data_{index}.geojson')
jsonList = {
    "dates":[date1,date2],
    "wmts":[wmts1,wmts2],
     "files": fileList
}
with open('/'.join(path)+f'/{directory}/list.json', 'w', encoding='utf-8') as f:
        json.dump(jsonList, f, ensure_ascii=False, indent=2)
print(str(datetime.now())+" - all done")
