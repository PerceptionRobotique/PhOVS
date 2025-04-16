PhOVS: Photometric Omnidirectional Visual Servoing

Library developed first at UPJV's MIS lab, then at CNRS-AIST JRL from fall 2019; put on github.com/PerceptionRobotique in 2022, where it should have been from the beginning
Author: G. Caron
Contact: guillaume.caron@u-picardie.fr
Main associated research article (https://inria.hal.science/hal-00829822v1): 
Guillaume Caron, Eric Marchand, El Mustapha Mouaddib. Photometric visual servoing for omnidirectional cameras. Autonomous Robots, 2013, 35 (2-3), pp.177-193. ⟨10.1007/s10514-013-9342-3⟩.  

Prerequisities
0. CMake (version 3.28.1 tested)
2. ViSP (version 3.4.1 tested)

Configure and prepare PhOVS to build
1. Create a build directory in PhOVS
2. Use cmake to fill the build directory (select the ViSP build or install directory)
3. Open the project in build or use make in the latter directory to build the library

Finalize the library setting by pasting the content of the copyContentInLibDirectory directory in the build/lib.

####################################################################################
Copyright (C) 2017-2025 by MIS lab (UPJV). All rights reserved.

See http://mis.u-picardie.fr/~g-caron/fr/index.php?page=7 for more information.

This software was developed at:
MIS - UPJV
33 rue Saint-Leu
80039 AMIENS CEDEX
France

and at
CNRS - AIST JRL (Joint Robotics Laboratory)
1-1-1 Umezono, Tsukuba, Ibaraki
Japan

This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

Description:
Insight about how to set the project and build the program
Authors:
Guillaume CARON
