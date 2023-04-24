# Boiler plate stuff to start the module
import jpype
import jpype.imports
from jpype.types import *

# Launch the JVM
jpype.startJVM(classpath=['jars/geoxygene-matching-1.10-SNAPSHOT.jar'])

# import the Java modules
from fr.ign.cogit.geoxygene.contrib.appariement.surfaces import ParametresAppSurfaces, AppariementSurfaces
from fr.ign.cogit.geoxygene.util.conversion import ShapefileReader, ShapefileWriter

# Copy in the patterns from the guide to replace the example code
param = ParametresAppSurfaces()
param.surface_min_intersection = 1
param.pourcentage_min_intersection = 0.2
param.pourcentage_intersection_sur = 0.8
param.minimiseDistanceSurfacique = True
param.distSurfMaxFinal = 0.6
param.completudeExactitudeMinFinal = 0.3
param.regroupementOptimal = True
param.filtrageFinal = True
param.ajoutPetitesSurfaces = True
param.seuilPourcentageTaillePetitesSurfaces = 0.1
param.persistant = False
param.resolutionMin = 1
param.resolutionMax = 11

db1 = ShapefileReader.read("./data/bati/bati_95430.shp", True)
db2 = ShapefileReader.read("./data/bati/cadastre_bati_95430.shp", True)
#print("db1 = ",db1.size())
for f in db1:
    f.setGeom(f.getGeom().get(0))
#print("db2 = ",db2.size())
for f in db2:
    f.setGeom(f.getGeom().get(0))
liensPoly = AppariementSurfaces.appariementSurfaces(db1, db2, param)
for f in liensPoly:
    if f.getGeom().size()>1:
        for i in range(0,f.getGeom().size()):
            print(f.getGeom().get(i))
    else:
        print(f.getGeom().get(0))
# AppariementSurfaces.writeShapefile(liensPoly, "./data/bati/appariement.shp")
#ShapefileWriter.write(liensPoly, "./data/bati/appariement.shp")
