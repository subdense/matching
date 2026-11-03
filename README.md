# Matching

## Installation Options

- [UV - An ultra-fast python package manager](#uv) (Recommended)
- [Simple Python](#simple-python)
- [Docker](./README.Docker.md)

## Prerequisities

To run the script, ensure you have the following prerequisites installed:

- Java: Download [here](https://www.java.com/fr/download/manual.jsp)
- Python: Version 3.11.4
  - or alternately refer to the uv package manager section below

<h3 id="uv">UV - An Ultra-fast python package manager (Recommended)</h3>

1. Install uv. Follow the official installation guide [here](https://docs.astral.sh/uv/getting-started/installation/)

2. Set and locks the default python version for this current project

```py
uv python pin 3.11.4
```

3. Run the script directly with uv with the appropriate command for your terminal:

**Linux / macOS (Bash)**

```py
uv run run.py -layer1 ./data/bati/neudorf_2012.shp -layer2 ./data/bati/neudorf_2022.shp -attributes '["HAUTEUR","ID"]' -output_prefix FR_NEU -java_memory 4G
```

**Windows (PowerShell)**

_Note: Inner double quotes must be escaped using backticks._

```py
uv run run.py -layer1 ./data/bati/neudorf_2012.shp -layer2 ./data/bati/neudorf_2022.shp -attributes '[\"HAUTEUR",\"ID\"]' -output_prefix FR_NEU -java_memory 4G
```

### Simple Python

```shell
pip install -r requirements.txt
```

### Run

Pour Neudorf :

```shell
python run.py -layer1 ./data/bati/neudorf_2012.shp -layer2 ./data/bati/neudorf_2022.shp -attributes '["HAUTEUR","ID"]' -output_prefix FR_NEU -java_memory 4G
```

Pour Strasbourg :

```shell
python run.py -layer1 ./data/fr/strasbourg_building_2011.gpkg -layer2 ./data/fr/strasbourg_building_2021.gpkg -attributes '["HAUTEUR","ID"]' -output_prefix FR_STR -java_memory 16G
```

Pour Dortmund :
J'ai du faire un prétraitement pour avoir le même nom d'attribut (bon, c'est le nom de la zone donc ça ne sert à rien hein donc vous pouvez aussi ignorer les attributs)

```shell
ogrinfo ./data/de/dor_hu_2011_raw.gpkg -sql "ALTER TABLE dor_hu_2011_raw RENAME COLUMN 'AGS' TO 'ags'"
python run.py -layer1 ./data/de/dor_hu_2011_raw.gpkg -layer2 ./data/de/dor_hu_2021_raw.gpkg -attributes '["ags"]' -output_prefix DE_DOR -java_memory 16G
```

Pour Liverpool:
J'ai ajouté la fusion des fichiers de livraison (comme ça c'est reproductible et on a pas de doublons) et un petit script (merge_input_data.py):

```shell
python merge_input_data.py -output ./data/uk/LIV_2011.gpkg -layer Topographicarea -where "theme='Buildings'" ./data/uk/Liverpool\ 2011/*.gpkg
python merge_input_data.py -output ./data/uk/LIV_2021.gpkg -layer Topographicarea -where "theme='Buildings'" ./data/uk/Liverpool\ 2021/*.gpkg
python run.py -layer1 ./data/uk/LIV_2011.gpkg -layer2 ./data/uk/LIV_2021.gpkg -attributes '[]' -output_prefix UK_LIV -java_memory 16G
```

## Export for validation

```shell
python export_to_validation.py -layer1 data/fr/strasbourg_building_2011.gpkg -layer2 data/fr/strasbourg_building_2021.gpkg -date1 2011 -date2 2021 -wmts1 ORTHOIMAGERY.ORTHOPHOTOS2011 -wmts2 ORTHOIMAGERY.ORTHOPHOTOS -select data/fr/strasbourg_sample.gpkg -select_layer sample -radius 100.0 -output_directory VAL
```
