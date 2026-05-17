{\rtf1\ansi\ansicpg1252\cocoartf2869
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;\f1\fnil\fcharset0 AmericanTypewriter;}
{\colortbl;\red255\green255\blue255;\red1\green22\blue40;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c11373\c20784;\csgray\c0;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 ### geo3D OpenFOAM configuration files\
\
![turbineSiting visualisation](../docs/turbineSiting.jpg)\
\
LoD1 building models (.obj) for renewable wind facility optimisation \
\
OpenFOAM simulation results available [here](https://drive.google.com/file/d/1SCjUU2tgAF9MOSMQKba05FmNaqGnL6ew/view?usp=share_link)\
\
openfoam commands are:\
```
\f1 \
\pard\pardeftab720\partightenfactor0
\cf0 \expnd0\expndtw0\kerning0
blockMesh\
\
surfaceFeatures\
\
\cf2 surfaceCheck constant/geometry/terrain.obj \cf0 | tee surfaceCheck_$(date +%Y%m%d_%H%M).log\cf2 \
\cf0 \
snappyHexMesh \
\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardeftab720\pardirnatural\partightenfactor0
\cf3 \kerning1\expnd0\expndtw0 \CocoaLigature0 checkMesh \cf0 \expnd0\expndtw0\kerning0
\CocoaLigature1 | tee checkMesh_$(date +%Y%m%d_%H%M).log\
\pard\pardeftab720\partightenfactor0
\cf0 \
topoSet\
\pard\pardeftab720\partightenfactor0
\cf0 \kerning1\expnd0\expndtw0 \
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardeftab720\pardirnatural\partightenfactor0
\cf3 \CocoaLigature0 checkMesh \cf0 \expnd0\expndtw0\kerning0
\CocoaLigature1 | tee checkMesh_$(date +%Y%m%d_%H%M).log\kerning1\expnd0\expndtw0 \
\pard\pardeftab720\partightenfactor0
\cf0 \
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardeftab720\pardirnatural\partightenfactor0
\cf3 \CocoaLigature0 foamRun -solver incompressibleFluid | \cf0 \expnd0\expndtw0\kerning0
\CocoaLigature1 tee foamRun_$(date +%Y%m%d_%H%M).log\
```}