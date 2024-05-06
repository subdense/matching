from os.path import exists
import os
from datetime import datetime
import argparse
from tqdm import tqdm
from pyproj import CRS
import geopandas
import networkx as nx
import json
import shapely
from pyproj import CRS, Transformer
from shapely import STRtree
from shapely.ops import transform

parser = argparse.ArgumentParser(description='Export matches for validation.')
parser.add_argument('-layer1', type=argparse.FileType('r'),
                    required=True, help='layer1')
parser.add_argument('-layer2', type=argparse.FileType('r'),
                    required=True, help='layer2')
parser.add_argument('-date1', type=str, required=True, help='date1')
parser.add_argument('-date2', type=str, required=True, help='date2')
parser.add_argument('-wmts1', type=str, required=True, help='wmts1')
parser.add_argument('-wmts2', type=str, required=True, help='wmts2')
parser.add_argument('-select', type=argparse.FileType('r'),
                    required=True, help='selection data (points used to create tasks)')
parser.add_argument('-select_layer', type=str, required=True,
                    help='name of the layer to extract from the selection data')
# parser.add_argument('-attribute', type = str, required=True, help='attribute for the objects in selection layer used to give names to the output directories')
parser.add_argument('-radius', type=float, required=True,
                    help='radius to select features around the selection layer')
parser.add_argument('-output_directory', type=str,
                    required=True, help='directory for output files')

args = parser.parse_args()

print(str(datetime.now())+" - reading the input data")
# crs = CRS.from_epsg(4326)
wgs84 = CRS.from_epsg(4326)
idb1 = geopandas.read_file(
    args.layer1.name, engine="pyogrio", columns=[], fid_as_index=True)
crs = idb1.crs
idb2 = geopandas.read_file(args.layer2.name, engine="pyogrio", columns=[
], fid_as_index=True).to_crs(crs)
select_layer = args.select_layer
select = geopandas.read_file(args.select.name, layer=select_layer,
                             engine="pyogrio", columns=[], fid_as_index=True).to_crs(crs)
radius = args.radius
project = Transformer.from_crs(crs, wgs84, always_xy=True).transform
path = [".", "output_data"]
directory = args.output_directory
date1 = args.date1
date2 = args.date2
wmts1 = args.wmts1
wmts2 = args.wmts2
print(f"{str(datetime.now())} - computing connected components")


def connected_components(idb1, idb2):
    print(str(datetime.now())+" - joining the data (sjoin)")
    join = geopandas.sjoin(idb1, idb2)
    G = nx.Graph()
    # add prefix to node ids to be able to separate them later
    join["index_left"] = join.index.map('L_{}'.format)
    join["index_right"] = join["index_right"].map('R_{}'.format)
    G.add_nodes_from(idb1.index.map('L_{}'.format))
    G.add_nodes_from(idb2.index.map('R_{}'.format))
    G.add_edges_from(
        list(map(tuple, join[["index_left", "index_right"]].to_numpy())))
    return list(sorted(nx.connected_components(G), reverse=True))


comp = connected_components(idb1, idb2)
print(f"{str(datetime.now())} - found {len(comp)} connected components")
comp_geoms = []
comp_db1 = []
comp_db2 = []
for index, component in enumerate(tqdm(comp, desc=f"{str(datetime.now())} - Processing connected components", position=0)):
    geoms = []
    db1_index = []
    db2_index = []
    db1_features = []
    db2_features = []
    # use prefix to separate the nodes into db1 and db2
    for n in component:
        if n.startswith("L_"):
            db1_index.append(int(n[2:]))
        else:
            db2_index.append(int(n[2:]))
    # select corresponding input features and add them to the output data
    for _, feature in idb1.loc[db1_index].iterrows():
        geoms.append(feature["geometry"])
        db1_features.append(feature)
    for _, feature in idb2.loc[db2_index].iterrows():
        geoms.append(feature["geometry"])
        db2_features.append(feature)
    comp_geoms.append(shapely.union_all(geoms))
    comp_db1.append(db1_features)
    comp_db2.append(db2_features)
    # print(f'{index}: ({len(comp_geoms)},{len(comp_db1)},{len(comp_db2)}) with ({len(db1_features)},{len(db2_features)})')
# create the STRtree for querying
tree = STRtree(comp_geoms)
os.makedirs(f'./output_data/{directory}/{select_layer}', exist_ok=True)


def createFeature(f, d):
    reprojected = transform(project, f["geometry"])
    return {
        "type": "Feature",
        "properties": {"date": d},
        "geometry": json.loads(shapely.to_geojson(reprojected))
    }


# for each sample, select the corresponding components
task_list = []
for index, sample in tqdm(select.iterrows(), desc=f"{str(datetime.now())} - Processing sampling points", position=0):
    sample_geom = sample.geometry
    list = tree.query(sample_geom, predicate="dwithin",
                      distance=radius).tolist()
    # print(f"found {len(list)} components")
    if (len(list) > 0):
        os.makedirs(
            f'./output_data/{directory}/{select_layer}/{index}', exist_ok=True)
        fileList = []
        for comp_index in list:
            features = []
            # select corresponding input features and add them to the output data
            for feature in comp_db1[comp_index]:
                features.append(createFeature(feature, date1))
            for feature in comp_db2[comp_index]:
                features.append(createFeature(feature, date2))
            geojson = {
                "type": "FeatureCollection",
                "features": features
            }
            with open('/'.join(path)+f'/{directory}/{select_layer}/{index}/data_{comp_index}.geojson', 'w', encoding='utf-8') as f:
                json.dump(geojson, f, ensure_ascii=False, indent=2)
            fileList.append(
                f'{select_layer}/{index}/data_{comp_index}.geojson')
        jsonList = {
            # "dates":[date1,date2],
            # "wmts":[wmts1,wmts2],
            "tasks": fileList
        }
        with open('/'.join(path)+f'/{directory}/{select_layer}/{index}/tasks.json', 'w', encoding='utf-8') as f:
            json.dump(jsonList, f, ensure_ascii=False, indent=2)
        task_list.append(f'{select_layer}/{index}/tasks.json')

sampleList = {
    "dates": [date1, date2],
    "wmts": [wmts1, wmts2],
    "samples": task_list
}
with open('/'.join(path)+f'/{directory}/{select_layer}/samples.json', 'w', encoding='utf-8') as f:
    json.dump(sampleList, f, ensure_ascii=False, indent=2)


datasets = '/'.join(path)+f'/{directory}/datasets.json'
if exists(datasets):
    with open(datasets, "r", encoding='utf-8') as f:
        data = json.load(f)
else:
    data = {"datasets": []}
data["datasets"].append(f'{select_layer}/samples.json')
with open(datasets, 'w', encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=2)

print(str(datetime.now())+" - all done")
