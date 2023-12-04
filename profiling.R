
setwd(paste0(Sys.getenv('CS_HOME'),'/SuburbanDensification/Models/Matching/matching'))

library(peakRAM)

# Test
path = '../../../Data/DonneesBDTopoHist/'
layer1 = 'BATI_20.shp'
layer2 = 'BATI_22.shp'

peakRAM(system('python run.py'))
peakRAM(system(paste0('python run.py ',path,layer1,' ',path,layer2)))
# -> system call memory not taken into account by peakRAM

Rprof(tf <- "rprof.log", memory.profiling=TRUE)
system(paste0('python run.py ',path,layer1,' ',path,layer2))
Rprof(NULL)
summaryRprof(tf, memory="both")
# idem

# -> need to profile 'by hand'

# with JPype, JVM and python are running in the same process:
#  see https://jpype.readthedocs.io/en/latest/userguide.html ("Alternatives")

system.time(system('python run.py'))


# test mem profiling
# https://stackoverflow.com/questions/22261452/finding-memory-usage-of-a-process-in-linux
# https://stackoverflow.com/questions/7880784/what-is-rss-and-vsz-in-linux-memory-management
#system('./profiling.sh "python run.py"', wait=F)
system('./profiling.sh', wait=F)
system.time(system('python run.py'))
system('./kill.sh')
system('cat profiling.log; rm profiling.log')


# systematic test for 6 cities with an increasing number of buildings
# with OSM data (note: less "m:n" than with BDTopo, so perfs may be different? also test with bdtopo later)

# coords (lon,lat)
cities = list('Toulouse'=c(1.444,43.604),
           'Strasbourg'=c(7.750,48.584),
           'Dortmund'=c(7.465,51.514),
           'Frankfurt'=c(8.682,50.110),
           'Bristol'=c(-2.597,51.453),
           'Liverpool'=c(-2.991,53.407)
           )
radiuses = seq(from=100,to=2000,by=100)
t0=2015
t1=2022

res = data.frame()
for(city in cities){
  for(radius in radiuses){
    lon = cities[[city]][1]
    lat = cities[[city]][2]
    name=paste0(lon,'_',lat,'_',radius)
    name0 = paste0(t0,'_',name)
    name1 = paste0(t1,'_',name)
    
    # get data
    system(paste0('./osmdata.sh ',t0,' ',lon,' ',lat,' ',radius))
    system(paste0('./osmdata.sh ',t1,' ',lon,' ',lat,' ',radius))
    
    # run the algorithm with profiling
    system('./profiling.sh', wait=F)
    t = system.time(system(paste0('python run.py tmp/buildings_',name0,'.shp tmp/buildings_',name1,'.shp')))
    system('./kill.sh')
    profiling <- read.csv('profiling.log',sep=";",header=F)
    system('cat profiling.log; rm profiling.log')
  }
}

