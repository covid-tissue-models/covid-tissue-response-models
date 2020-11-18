The COVID-19 simulation’s source code is available in the GitHub repository 
https://github.com/covid-tissue-models/covid-tissue-response-models/tree/master/CC3D/Models/BiocIU/SARSCoV2MultiscaleVTM. 
The simulation is a model specification which runs in the CompuCell3D virtual-tissue simulation environment. To run the simulation requires installation of 
CompuCell3D version 4.1.1 or later. CompuCell3D is open-source and runs on Windows, Mac and Linux operating systems. It can be downloaded from 
https://compucell3d.org/SrcBin. Installers are available for Windows operating systems and Mac installation also does not require compilation. CompuCell3D’s manuals are 
available at https://compucell3d.org/Manuals.  The COVID-19 simulation can also be run online without requiring any installations or downloads on the nanoHUB servers 
at https://nanohub.org/tools/cc3dcovid19/. Use of nanoHUB is free but requires user registration. The simulation may take a few moments to load in its nanoHUB 
deployment; during load the simulation area will be blue.

*Twedit++* is a specialized text editor for CompuCell3D simulations which comes packaged with CompuCell3D. Tweedit++ can open cc3d file extensions which contain the 
simulation file structure for CompuCell3D simulations. To open the COVID-19 simulation click “Open CC3D Project” (Figure S12) and select *ViralInfectionVTM.cc3d* in
``<repository-folder>/covid-tissue-response-models/CC3D/Models/BiocIU/SARSCoV2MultiscaleVTM/Model. Once opened, ViralInfectionVTM.cc3d`` will appear in the left-hand panel 
“CC3D Simulation” (Figure S13).  Double click on it to open all simulation files. The main simulation files are: ``ViralInfectionVTM.xml``, which defines certain simulation 
properties (e.g., cell types, lattice size, energy-constraint plugins, diffusive fields); ``ViralInfectionVTMSteppables.py``, which defines the simulation’s initial 
conditions, main interactions and dynamics (e.g., cell initialization, immune-cell recruitment, secretion by cells into fields, uptake by cells from fields); 
``ViralInfectionVTMSteppableBasePy.py``, where the viral infection Antimony submodel is declared; ``ViralInfectionVTMModelInputs.py``, which sets the parameters. The submodels 
in ``ViralInfectionVTMSteppables.py`` are organized as python classes, making them easy to modify. Tweedit++ can also copy and rename the simulation project to a new 
directory by using CC3D Project; Save Project As. However the save as does not copy the folder ``<...>/Model/nCoVToolkit``, which must be copied into the new simulation 
directory separately.

CompuCell3D *Player* is a GUI tool which executes CompuCell3D simulations during desktop execution (or on nanoHUB). In order to run the simulation either right click *ViralInfectionVTM.cc3d* in the left hand panel and select “Open In Player” or open *CompuCell3D.exe* to open the CompuCell3D Player and select ``File; Open Simulation File`` and open *ViralInfectionVTM.cc3d*. Once the simulation is open in Player click play (on the nanoHub deployment the simulation should start automatically). Player will display windows with the cell lattice rendered (Figure S14), set the z-plane to 0 to visualize the epithelial cells and the z-plane to 1 to visualize the immune cells. More windows can be created (menu ``Window; New graphics Window``, Figure S15) as needed to visualize the virus, cytokine, oxidative agent fields. Each window has a drop-down menu containing the selection of fields that can be rendered (i.e., the chemical fields and the cell field, Figure S16). 

.. image::  https://github.com/covid-tissue-models/covid-tissue-response-models/blob/master/CC3D/Models/BiocIU/SARSCoV2MultiscaleVTM/Tutorials/Downloading%20and%20Running/Figure%20S12.png 
**Figure S12. Opening a CompuCell3D project in Tweedit++.**


.. image::  https://github.com/covid-tissue-models/covid-tissue-response-models/blob/master/CC3D/Models/BiocIU/SARSCoV2MultiscaleVTM/Tutorials/Downloading%20and%20Running/Figure%20S13.png 
**Figure S13. Tweedit++’s left hand panel with simulation project files openned.**


.. image::  https://github.com/covid-tissue-models/covid-tissue-response-models/blob/master/CC3D/Models/BiocIU/SARSCoV2MultiscaleVTM/Tutorials/Downloading%20and%20Running/Figure%20S14.png 
**Figure S14. Example of CompuCell3D’s player open with the COVID-19 simulation loaded.**


.. image::  https://github.com/covid-tissue-models/covid-tissue-response-models/blob/master/CC3D/Models/BiocIU/SARSCoV2MultiscaleVTM/Tutorials/Downloading%20and%20Running/Figure%20S15.png 
**Figure S15. How to open a new simulation render window in CompuCell3D Player.**


.. image::  https://github.com/covid-tissue-models/covid-tissue-response-models/blob/master/CC3D/Models/BiocIU/SARSCoV2MultiscaleVTM/Tutorials/Downloading%20and%20Running/Figure%20S16.png 
**Figure S16. Drop-down menu in simulation render window to select which field to render.**
