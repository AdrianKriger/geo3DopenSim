### geo3D OpenFOAM configuration files

![turbineSiting visualisation](../docs/turbineSiting.jpg)

LoD1 building models (.obj) for renewable wind facility optimisation

OpenFOAM simulation results available [here](https://drive.google.com/file/d/1SCjUU2tgAF9MOSMQKba05FmNaqGnL6ew/view?usp=share_link)

openfoam commands are:
```
blockMesh

surfaceFeatures

surfaceCheck constant/geometry/terrain.obj | tee surfaceCheck_$(date +%Y%m%d_%H%M).log

snappyHexMesh 

checkMesh| tee checkMesh_$(date +%Y%m%d_%H%M).log

topoSet

checkMesh | tee checkMesh_$(date +%Y%m%d_%H%M).log

foamRun -solver incompressibleFluid | tee foamRun_$(date +%Y%m%d_%H%M).log
```
a taste of the results (wake influence, power curve and power loss) are in the [turbineSiting.ipynb](https://github.com/AdrianKriger/geo3DopenSim/blob/main/turbineSiting/turbineSiting.ipynb)
