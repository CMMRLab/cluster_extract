# cluster_extract (version 16July2024)
cluster_extract is primarily meant to work with ReaxFF implemented in LAMMPS to analyze ReaxFF data files and bond order files. The code is meant to be able to "extract" clusters from a ReaxFF simulation via a mass or number of atoms cutoff. This allows for the "removal" of small products. After extracting cluster of atoms, cluster_extract performs in depth analysis of ring connectivity, hybridization state, ring alignment, and can spatially bin quantities to look at for example how the ring content and hybridization vary compared to a graphene sheet in a composite model.

cluster_extract can be run in a stand alone mode, can be run to iterate through a directory of data files and bond order files to perform batch analysis, or run within a PyLAMMPS script to switch between a LAMMPS simulation and an cluster_extract operation and analysis. The cluster_extract code base is the precursor to what become LUNAR (https://github.com/CMMRLab/LUNAR). Unfortunately, cluster_extract has not been properly documented like LUNAR and users of cluster_extract will have to use the examples and the comments to learn the code. Additionally, cluster_extract does not planned to be maintain as LUNAR does. 

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

## Terms of use
You agree to acknowledge the use of the LUNAR in any reports or publications of results obtained with the LUNAR by citing "Evolution of Glassy Carbon Derived from Pyrolysis of Furan Resin" at https://doi.org/10.1021/acsaenm.3c00360