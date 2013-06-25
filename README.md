simplicial-nonlinear-principal-component-analysis
===============================================================================
Overview:
Simplicial nonlinear principal component analysis (SNPCA) is a manifold reconstruction algorithm that takes a point cloud that falls on the surface of a two-dimensional manifold as input, and outputs a triangulation that fits the manifold and point cloud. 
A description of the algorithm along with theoretical justification is available at:
http://arxiv.org/abs/1207.3374
===============================================================================
Driver files:
The easiest way to run SNPCA is with a driver file. A driver file typically consists of three sections:
A) Set algorithm parameters in SNPCA_params structure
B) Generate surface data point cloud
C) Run SNPCA algorithm by calling SNPCA_interleaved_main
A list of example driver files are presented in the next section.
===============================================================================
List of driver files:
SNPCA_cone.m
SNPCA_creased_sheet.m
SNPCA_Ford_LIDAR.m
SNPCA_sphere.m
SNPCA_swiss_roll.m
SNPCA_torus.m
===============================================================================
The file SCAN1000.mat contains point cloud data discussed in [1].
It is available from: http://robots.engin.umich.edu/SoftwareData/Ford.

[1] Gaurav Pandey, James R. McBride and Ryan M. Eustice, Ford campus vision and lidar data set. International Journal of Robotics Research, 30(13):1543-1552, November 2011.
===============================================================================
License:
Copyright (c) 2013, Thomas Hunt and Arthur J. Krener
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.
