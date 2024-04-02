import argparse
import pandas
import geopandas
from datetime import datetime
import numpy as np
from tqdm import tqdm
from shapely.ops import triangulate,unary_union
from libpysal.cg import voronoi_frames
from libpysal import weights
import momepy
import networkx as nx
from shapely import centroid, LineString, Polygon, MultiPolygon, union_all
geopandas.options.io_engine = "pyogrio"
import operator
from functools import reduce
from itertools import chain

parser = argparse.ArgumentParser(description='Create building blocks.')
parser.add_argument('-output', type = argparse.FileType('w'), required=True, help='output file')
parser.add_argument('-parcels', type = argparse.FileType('r'), required=True, help='input parcel file')
parser.add_argument('-buildings', type = argparse.FileType('r'), required=True, help="input building file")
parser.add_argument('files', nargs='+', metavar='input_file', type = argparse.FileType('r'), help="input linear files")
parser.add_argument('-threshold', type = float, required=False, default = 70.0, help="threshold for the overlay")
args = parser.parse_args()
print(f"{str(datetime.now())} - Create building blocks with {len(args.files)} files")
buildings = geopandas.read_file(args.buildings.name, engine="pyogrio", columns=[], fid_as_index=True)
building_points = buildings["geometry"].map(centroid)
common_crs = buildings.crs
buildings_centroids = geopandas.GeoDataFrame({'id': list(buildings.index), 'geometry':building_points}, crs = common_crs)
buildings_centroids.to_file(args.output.name, layer='building_centroids', driver="GPKG")
print(f"{str(datetime.now())} - Selecting buildings with threshold {args.threshold}")

delaunay_triangles_geometries = triangulate(unary_union(buildings_centroids.geometry.to_list()),edges=True)
delaunay_triangles_index = list(range(0,len(delaunay_triangles_geometries)))
delaunay_triangles = geopandas.GeoDataFrame({'id': list(delaunay_triangles_index), 'geometry':delaunay_triangles_geometries}, crs = common_crs)
delaunay_triangles.to_file(args.output.name, layer='delaunay_triangles', driver="GPKG")

to_drop = delaunay_triangles.index[delaunay_triangles.geometry.length > args.threshold].tolist()
print(f"{str(datetime.now())} - found {len(to_drop)} links longer than {args.threshold}")
delaunay_triangles.drop(to_drop, inplace=True)
delaunay_triangles.to_file(args.output.name, layer='delaunay_triangles_filtered', driver="GPKG")
 
lines = pandas.concat([geopandas.read_file(file.name, engine="pyogrio", columns=[]).to_crs(common_crs) for file in tqdm(args.files, desc=f"{str(datetime.now())} - Importing linear files")])
join = geopandas.sjoin(delaunay_triangles,lines,lsuffix="left_lines",rsuffix="right_lines")
to_drop = list(set(join.index))
print(f"{str(datetime.now())} - found {len(to_drop)} links intersecting lines")
delaunay_triangles.drop(to_drop, inplace=True)
delaunay_triangles.to_file(args.output.name, layer='delaunay_triangles_final', driver="GPKG")

print(f"{str(datetime.now())} - computing connected components")
G_primal = momepy.gdf_to_nx(delaunay_triangles, approach="primal", multigraph=False)
connected_components = list(sorted(nx.connected_components(G_primal), reverse=True))
node_dictionary = dict()
for index, comp in enumerate(tqdm(connected_components)):
    for n in comp: node_dictionary.update({n:index})
print(f"{str(datetime.now())} - found {len(connected_components)} connected components")
edge_dictionary = {(e[0], e[1]): {"component": node_dictionary[e[0]]} for e in G_primal.edges}
nx.set_edge_attributes(G_primal, edge_dictionary)
nodes, edges = momepy.nx_to_gdf(G_primal)
edges = edges.drop(columns=['mm_len','node_start','node_end'])
edges.to_file(args.output.name, layer='connected', driver="GPKG")
components = edges.dissolve(by='component')
components.to_file(args.output.name, layer='components', driver="GPKG")
parcels = geopandas.read_file(args.parcels.name, engine="pyogrio", columns=[], fid_as_index=True).to_crs(common_crs)
print(f"{str(datetime.now())} - join parcels and connected components")
join = geopandas.sjoin(parcels, components, lsuffix="left_parcel",rsuffix="right_component")
print(join.head())
G = nx.Graph()
join["index_left"] = join.index.map('L_{}'.format)
join["index_right"] = join["index_right_component"].map('R_{}'.format)
node_list = list(set(join["index_left"].to_list()+join["index_right"].to_list()))
G.add_nodes_from(node_list)
G.add_edges_from(list(map(tuple, join[["index_left","index_right"]].to_numpy())))
comp = list(sorted(nx.connected_components(G), reverse=True))
comp_index = []
comp_geom = []
for index, component in enumerate(tqdm(comp, desc=f"{str(datetime.now())} - Processing connected components")):
    parcel_ids = [int(n[2:]) for n in component if n.startswith("L_")]
    comp_index.append(index)
    comp_geom.append(union_all(parcels.loc[parcel_ids].geometry))
components_parcels = geopandas.GeoDataFrame({'id': comp_index, 'geometry':comp_geom}, crs = common_crs)
# components_parcels = join.dissolve(by='index_right_component')
def removeHoles(geom):
    if geom.geom_type == 'Polygon': return Polygon(shell=geom.exterior.coords)
    return MultiPolygon(list(map(lambda p: Polygon(shell=p.exterior.coords), geom.geoms)))
components_parcels["geometry"] = components_parcels["geometry"].map(removeHoles)
# components_parcels = geopandas.overlay(components_parcels, components_parcels, how='union')
print(components_parcels.head())
components_parcels.to_file(args.output.name, layer='components_parcels', driver="GPKG")

print(f"{str(datetime.now())} - Create building blocks done")
