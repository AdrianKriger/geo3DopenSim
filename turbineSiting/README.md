<<<<<<< HEAD
### geo3D OpenFOAM configuration files
=======
### *geo3D* OpenFOAM configuration files
>>>>>>> 5394befa792ee1b4cfde118623e335a01fb8d13a

![turbineSiting visualisation](../docs/turbineSiting.jpg)

LoD1 building models (.obj) for renewable wind facility optimisation 

<<<<<<< HEAD
openfoam commands are:
```
blockMesh

surfaceFeatures

surfaceCheck constant/geometry/terrain.obj | tee surfaceCheck_$(date +%Y%m%d_%H%M).log

snappyHexMesh | tee surfaceCheck_$(date +%Y%m%d_%H%M).log

checkMesh | tee checkMesh_$(date +%Y%m%d_%H%M).log

topoSet

checkMesh | tee checkMesh_$(date +%Y%m%d_%H%M).log

foamRun -solver incompressibleFluid | tee foamRun_$(date +%Y%m%d_%H%M).log

foamPostProcess -time 1000:
```
A taste of the results (wake field, power curve and power loss) is available in the [turbineSiting.ipynb](https://github.com/AdrianKriger/geo3DopenSim/blob/main/turbineSiting/turbineSiting.ipynb)
=======
OpenFOAM simulation results available [here](https://drive.google.com/file/d/1SCjUU2tgAF9MOSMQKba05FmNaqGnL6ew/view?usp=share_link)

openfoam commands are:
```
blockMesh
surfaceFeatures
surfaceCheck constant/geometry/terrain.obj | tee surfaceCheck_$(date +%Y%m%d_%H%M).log
snappyHexMesh
checkMesh | tee checkMesh_$(date +%Y%m%d_%H%M).log
topoSet
checkMesh | tee checkMesh_$(date +%Y%m%d_%H%M).log
foamRun -solver incompressibleFluid | tee foamRun_$(date +%Y%m%d_%H%M).log
```

A taste of the results (wake field, power curve and power loss) is available in the [turbineSiting.ipynb](https://github.com/AdrianKriger/geo3DopenSim/blob/main/turbineSiting/turbineSiting.ipynb)
>>>>>>> 5394befa792ee1b4cfde118623e335a01fb8d13a
