import jpype.imports
from jpype.types import *
jpype.startJVM(classpath=['jars/geoxygene-matching-1.10-SNAPSHOT.jar'])

from fr.ign.cogit.geoxygene.contrib.appariement.surfaces import AppariementSurfaces

import matching

params = matching.default_params()

layer1name, layer2name, path, db1, db2, crs = matching.get_data(params)

#db1, db2, allfeatures = matching.preprocess_data(layer1name, layer2name, db1, db2)
db1, db2 = matching.preprocess_data(layer1name, layer2name, db1, db2, id_index = params['id_index'])

# call to geoxygene matching algorithm
liensPoly = AppariementSurfaces.appariementSurfaces(db1, db2, params['algo_params'])

links, features_stable, features_split, features_merged, features_aggregated, all_link_targets, all_link_sources, features_disappeared, features_appeared = matching.post_process_links(liensPoly, db1, db2, crs, layer1name, layer2name, params['id_index'])

matching.export_links(links, layer1name, layer2name, path, params)

matching.export(features_appeared, features_disappeared, features_stable, features_split, features_merged, features_aggregated, crs, layer1name, layer2name, path)


jpype.shutdownJVM()
