# matching

## Installation

Tested with python3.10 (does not work with anaconda: issues with importing shapely functions)

pip install -r requirements.txt

## Run

Pour Neudorf :
python run.py -layer1 ./data/bati/neudorf_2012.shp -layer2 ./data/bati/neudorf_2022.shp -attributes '["HAUTEUR","ID"]' -output_prefix FR_NEU -java_memory 4G

Pour Strasbourg :
python run.py -layer1 ./data/fr/strasbourg_building_2011.gpkg -layer2 ./data/fr/strasbourg_building_2021.gpkg -attributes '["HAUTEUR","ID"]' -output_prefix FR_STR -java_memory 16G

Pour Dortmund :
J'ai du faire un prétraitement pour avoir le même nom d'attribut (bon, c'est le nom de la zone donc ça ne sert à rien hein donc vous pouvez aussi ignorer les attributs)
ogrinfo ./data/de/dor_hu_2011_raw.gpkg -sql "ALTER TABLE dor_hu_2011_raw RENAME COLUMN 'AGS' TO 'ags'"
python run.py -layer1 ./data/de/dor_hu_2011_raw.gpkg -layer2 ./data/de/dor_hu_2021_raw.gpkg -attributes '["ags"]' -output_prefix DE_DOR -java_memory 16G

Pour Liverpool:
J'ai ajouté la fusion des fichiers de livraison (comme ça c'est reproductible et on a pas de doublons) et un petit script (merge_input_data.py):
python merge_input_data.py -output ./data/uk/LIV_2011.gpkg -layer Topographicarea -where "theme='Buildings'" ./data/uk/Liverpool\ 2011/*.gpkg
python merge_input_data.py -output ./data/uk/LIV_2021.gpkg -layer Topographicarea -where "theme='Buildings'" ./data/uk/Liverpool\ 2021/*.gpkg
python run.py -layer1 ./data/uk/LIV_2011.gpkg -layer2 ./data/uk/LIV_2021.gpkg -attributes '["ags"]' -output_prefix DE_DOR -java_memory 16G
