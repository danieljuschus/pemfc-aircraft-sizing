# Preliminary Propulsion System Sizing Methods for PEM Fuel Cell Aircraft - MSc thesis

This repository contains the main part of the code that I developed for my [Master thesis at the TU Delft in 2020/21](http://dx.doi.org/10.13140/RG.2.2.34701.72165).

For the required Python packages, see [requirements.txt](requirements.txt) (or [requirements_gui.txt](requirements_gui.txt), if you also want to use the GUI).

The basic way to use the models is shown in [main.py](app/main.py). The function [size_system](https://github.com/danieljuschus/pemfc-aircraft-sizing/blob/dee57f9b6d7d06745a7f22774348ed0542ef6bd1/app/main.py#L12) 
sizes a fuel cell system for an aircraft. There is also a GUI that you can use (experimental!).

For the methodology applied in the models used in the code please take a look at chapter 3 of my thesis report: 
http://dx.doi.org/10.13140/RG.2.2.34701.72165. The 
report also contains an introduction into FC aircraft and a summary of the results that were obtained. The directory [old](old) contains some old
scripts that I used to generate the results shown in my thesis. For the exact version of the code at the moment of my defence see https://github.com/danieljuschus/pemfc-aircraft-sizing/releases/tag/v1.0.

The models contained in the code are not being further developed in this repository. But if anyone finds a bug, feel free to open an issue and I will take a look at it. 
Any recent commits have just been due to me playing around with the GUI. 