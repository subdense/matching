import matching

# Launch the JVM
jpype.startJVM(classpath=['jars/geoxygene-matching-1.10-SNAPSHOT.jar'])

params = default_params()

layer1name, layer2name, db1, db2, crs = get_data(params)

#db1, db2, allfeatures = preprocess_data(layer1name, layer2name, db1, db2)
db1, db2 = preprocess_data(layer1name, layer2name, db1, db2)

# call to geoxygene matching algorithm
liensPoly = AppariementSurfaces.appariementSurfaces(db1, db2, params['algo_params'])

links, features_stable, features_split, features_merged, features_aggregated, all_link_targets, all_link_sources, features_disappeared, features_appeared = post_process_links(liensPoly)


# export links to shp
links.to_file('/'.join(path)+'/MATCHING-LINKS_'+layer1name+"_"+layer2name+'.shp')

geojson_export(links, layer1name, layer2name)

export(features_appeared, features_disappeared, features_stable, features_split, features_merged, features_aggregated)


jpype.shutdownJVM()


