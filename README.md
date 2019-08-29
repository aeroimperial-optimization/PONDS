# PONDS (Polynomial Optimisation of Non-linear systems)
An open-source MATLAB toolbox for finding lower bounds on long-time averaged polynomial magnitudes in hydrodynamic-type systems. Below is a quick guide to PONDS, for a more deatiled description and examples refer to the user manual, which you can find inside the `./docs/` folder.

Note: This code is under very active development and will be subject to changes and improvements, so watch this space!

* **Latest release:** 1.0
* **Release date:** 29 August 2019
* **Release notes:**
	- Toolbox [UODESys](https://github.com/aeroimperial-optimization/UODESys) has been included. UODESys interface has been modified to be more user friendly.
	- Functions in BoundSys folder are design to find bounds on long-time averaged polynomial magnitudes in hydrodynamic-type systems.
  - Examples have been included for Lorenz Attractor and N-dimensional truncations of Kuramoto-Sivashinsky equation.

## Contents
- [System requirements](#Requirements)
- [License](#License)
- [Installation](#Install)
- [Getting started](#GettingStarted)


## System requirements<a name="Requirements"></a>

In order to use QUINOPT, you will need:

1. A working version of [YALMIP](https://yalmip.github.io/), the MATLAB optimization modelling software by J. L&ouml;fberg.
2. Alternatively to YALMIP [SPOTLess](https://github.com/mmt/spotless/tree/master) can be used. However, YALMIP is still needed by UODESys. 
3. A suitable SDP solver. If you are using YALMIP, choices include [SeDuMi](https://github.com/sqlp/sedumi), [SDPT3](http://www.math.nus.edu.sg/~mattohkc/sdpt3.html), [SDPA](http://sdpa.sourceforge.net/), [Mosek](https://www.mosek.com/) (free for
    users in academia). If you are using SPOTLess you can chose between [Mosek](https://www.mosek.com/) and [SeDuMi](https://github.com/sqlp/sedumi).

## License<a name="License"></a>

POND is distributed under the [Apache 2.0 licence](http://www.apache.org/licenses/LICENSE-2.0):

Copyright 2019, M. Lino.

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at [http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

## Installation<a name="Install"></a>

To install PONDS:

1. [Download](https://yalmip.github.io/download/) and [install](https://yalmip.github.io/tutorial/installation/) YALMIP.
2. You can also [Download](https://github.com/mmt/spotless/tree/master) SPOTLess and install it by running spotinstall.m.
3. Install a semidefinite programming (SDP) solver compatible with YALMIP and/or SPOTLess, such as [Mosek](https://www.mosek.com/) or [SeDuMi](http://sedumi.ie.lehigh.edu/).
4. Add folder `lib/` and its subfolders to MATLAB's path by running:
```Matlab
>> addpath(genpath(’./lib’))
```
at MATLAB's command prompt.

## Getting started<a name="GettingStarted"></a>

To get started with PONDS, please look at the sample scripts on the folder './examples/'. A description of the examples can be found in the manual available in the './doc/' folder.

For more information on the main functions in BoundSys() and BoundUsys(), please consult the documentation or type

```Matlab
>> help BoundSys
>> help BoundUsys
```

at MATLAB's command prompt.

