# LeMoS URANS/LES model including shielding function (SLH)

### License

This file is part of OpenFOAM.  
  
OpenFOAM is free software: you can redistribute it and/or modify it  
under the terms of the GNU General Public License as published by  
the Free Software Foundation, either version 3 of the License, or  
(at your option) any later version.  
  
OpenFOAM is distributed in the hope that it will be useful, but WITHOUT  
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or  
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License  
for more details.  
  
You should have received a copy of the GNU General Public License  
along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.  

### Copyright Information
    
Copyright (C) 1991-2009 OpenCFD Ltd.  
Copyright (C) 2011-2024 Chair of Modeling and Simulation (LeMoS), University of Rostock, Germany

### Author
    
Ivan Shevchuk <ivan.shevchuk@uni-rostock.de>  
Website: https://www.lemos.uni-rostock.de

Changes implemented by Luise Draheim <luise.draheim@posteo.de>  

### Notes
    
Compatibility: OpenFOAM-5.x, OpenFOAM-6, OpenFOAM-7, v16x, v18x, v20x, v21x, v22x

### Description

Implementation of the hybrid turbulence modelling approach combining Menters  
k-Omega-SST model and Dynamic Smagorinsky SGS model for incompressible flows.  
  
The Hybrid approach used is described in:  
	
    @verbatim
		Kornev, N., Taranov, A., Shchukin, E. & Kleinsorge, L. 
		Development of hybrid URANS-LES methods for flow simulations in the ship stern area. 
		Ocean Engineering, 
		Vol. 38, 1831-1838. 
	@endverbatim

The computational domain is dynamically divided into URANS and LES regions, depending  
on the integral length scale. If the length scale is larger than the LES filter size, the  
turbulent viscosity is calculated using LES SGS model. Otherwise, the cell is in URANS  
region and turbulent viscosity is calculated using k-Omega SST model.  
  
At the moment the implementation does not support arbitrary switching between different  
URANS/LES models. This means, that only the combination of k-Omega-SST and Dynamic Smagorinsky  
model is available.  

### Instructions

1.) Enter the directory where the source code has been extracted, and compile it by typing: 

````
wmake libso  
````

2.) Add the following line to the controlDict of your case:

````
libs ( "libhybKOmegaSST.so" );
````

3.) Specify

````
RAS
{
    RASModel	hybKOmegaSST;

    turbulence	on;

    printCoeffs	on;

    filter	simple;

    delta	cubeRootVol;

    cubeRootVolCoeffs
    {
    	deltaCoeff      1;
    }

}
````
in turbulenceProperties.

### Disclaimer

This offering is not approved or endorsed by OpenCFD Limited, the producer  
of the OpenFOAM software and owner of the OPENFOAM®  and OpenCFD®  trade marks.  
  
Detailed information on the OpenFOAM trademark can be found at  
  
- http://www.openfoam.com/legal/trademark-policy.php  
- http://www.openfoam.com/legal/trademark-guidelines.php  

For further information on OpenCFD and OpenFOAM, please refer to  
  
- http://www.openfoam.com
