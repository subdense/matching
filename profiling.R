
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





