import argparse
import pandas
import geopandas
from datetime import datetime
import pandas
from tqdm import tqdm
geopandas.options.io_engine = "pyogrio"

parser = argparse.ArgumentParser(description='Merging input datasets.')
parser.add_argument('-output', type = argparse.FileType('w'), required=True, help='output file')
parser.add_argument('-layer', type = str, required=True, help='input layer name')
parser.add_argument('-where', type = str, required=False, help="input filter, for instance: theme='Buildings'")
parser.add_argument('files', nargs='+', metavar='input_file', type = argparse.FileType('r'))
args = parser.parse_args()
print(f"{str(datetime.now())} - start merge with {len(args.files)} files")
gdf = pandas.concat([geopandas.read_file(file.name, engine="pyogrio", layer=args.layer, where=args.where) for file in tqdm(args.files, desc="Importing files")])
gdf.drop_duplicates(inplace=True)
gdf.rename(columns = {'fid':'ID'}, inplace = True) 
print(f"{str(datetime.now())} - saving to {args.output.name}")
gdf.to_file(args.output.name, layer='buildings', driver="GPKG")
print(f"{str(datetime.now())} - merge done")
