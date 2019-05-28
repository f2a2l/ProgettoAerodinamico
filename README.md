# ProgettoAerodinamico

Credits:

- HessSmith solver originally developed by Federico Messanelli (federico.messanelli@polimi.it).

deleted:    FWcorrections.m

## Multi-element geometry

To deal with multiple geometries, the following conventions have been adopted:

- the reference aerodynamic length is the chord of the first (upstream) profile, and it’s unitary (= 1).
- a reference system attached to the first (upstream) profile is defined:
  - origin on the leading edge of the first profile;
  - x axis passing through the trailing edge of the first profile, pointing towards it;
    - this means that the point (1,0) corresponds to the trailing edge of the first profile;
  - y axis consequently.

Function `multiGeometry` can be found in folder `Geometry`; it can be used as such:

```MATLAB
[x, y, totLength] = multiGeometry(npoint, arflPar, alpha, dist, crel, pltFlag);
```

Input arguments:

- _npoint_ is half the number of nodes (actually, the number of nodes on the camber line, which thus define nodes on the edge of the profile).
- _arflPar_ is a matrix, each row corresponding to an airfoil (no. rows = no. airfoils); the 8 columns represent the 8 parameters needed by IGP parametrisation.
- _alpha_ is a row vector, containing angles of attack of the airfoils (all angles are expressed in an absolute frame).
- _dist_ is a matrix, each row corresponding to an airfoil; however, the first one is omitted, so that the number of rows is the the number of airfoils -1. Each row contains the distance vector of the corresponding profile from the trailing edge of the first one; columns correspond to the x and y components (expressed in the above mentioned frame).
- _crel_ is a column vector, each row/component corresponding to an airfoil; however, the first one is omitted, so that the number of rows/components is the the number of airfoils -1. Each component represents the chord of the corresponding profile adimensionalised over the chord of the first profile. Since the latter has unit chord, each component effectively contains the chord of the corresponding profile.
- _pltFlag_ is a boolean __optional__ input argument; if `true`, plots the geometry with nodes.

Outputs:

- _x_ and _y_ are matrices containing x and y coordinates for the nodes of each profile (using the above mentioned frame). Each column corresponds to a profile, each row to a node.
- _totLength_ is the maximum value of _x_ found in the above mentioned matrix; in other words, it is the total length of the multi-element profile, measured along the x axis (defined above). Keep in mind that $x = 1$ is the trailing edge of the first profile.

Example with two airfoils:

```MATLAB
arflPar = [0.3, 0.6, 0, 0, 0.3, 0.12, 0.3, 1.5;
					 0.3, 0.6, 0, 0, 0.3, 0.12, 0.3, 1.5];
[x, y, totLength] = multiGeometry(80, arflPar, [3 10], [-0.05 -0.05], 0.3, true);
```





## Hess smith solver

A Hess Smith solver is included in folder `HessSmith`. To start an UI version of this program, just run `StartHS` (after running `addPaths`). Alternatively, you can access the solver as a function with the instructions below.



### Using the solver

Function `solverHS` takes different inputs depending whether it is used for a single or multiple airfoils:



__Single airfoil case__: `[...] = solverHS(npoint, arflPar, alpha)`, where:
- _npoint_ is half the number of points used to generate geometry.
- _arflPar_ is a row vector, each component corresponding to the parameters describing airfoil geometry.
- _alpha_ is the angle of attack.

Example:
```MATLAB
arflPar = [0.3, 0.6, 0, 0, 0.3, 0.12, 0.3, 1.5];
[Cl, Cd] = solverHS(80, arflPar, 3);
```



__Multiple airfoil case__: `[...] = solverHS(npoint, aname, alpha, dist, crel)`; the inputs are exactly the same as the ones used for `multiGeometry` (see above).

Example:

```MATLAB
arflPar = [0.3, 0.6, 0, 0, 0.3, 0.12, 0.3, 1.5;
           0.3, 0.6, 0, 0, 0.3, 0.12, 0.3, 1.5];
[Cl, Cd] = solverHS(80, arflPar, [3, 6], [-0.05, -0.05], 0.3);
```



### Outputs

In most cases, it is sufficient to use (just like the examples above): 

```MATLAB	
[Cl, Cd] = solverHS(...)
```

_Cl_ and _Cd_ are column vectors containing the coefficients of lift and drag of each profile.

However, the full list of output is the following:

```MATLAB
[Cl, Cd, totLength, Cp, maxdCp, x, y, p, p1, SOL] = solverHS(...)
```

Any of them can be arbitrarily omitted, either by truncating the list of output arguments or by inserting `~` instead of them. Here’s a brief explaination of each output:

- _totLength_, _x_ and _y_ are the same as in `multiGeometry`.
- _Cp_ contains the chord distribution of pressure coefficient (possibly deprecated).
- _maxdCp_ is a matrix containing values of maximum $\Delta C_P$; each row corresponds to an airfoil, first column corresponts to lower part, second column to upper part.
- _p_, _p1_ are classes containing information about the panelisation; _SOL_ is the solution of the linear system associated to each geometry. These outputs are not meant for extarnal usage; however, they are used by the UI programs `StartHS`.





## Straight perfomance

The straight selected for the project is the straight between turn 2 and turn 3 of Baku City Circuit.

The evaluation of time needed to run this straight is performed by ```[T_sector] = sector(CL, CD, fig)```

- _CL_ is a vector containing the rear wing CL when DRS is closed and when DRS is open.
- _CD_ is a vector containing the rear wing CD when DRS is closed and when DRS is open.
- _fig_ creates figures of straight performances if ```true```.

Example: 

```matlab
CD = [1.169, 0.969];
CL = [4.846, 4.346];

[T_sector] = sector(CL, CD, true);
```





## Endplates correction: documentation
[t,s] = corr2dto3d(lambda), where:
- _lambda_ is the wing aspect ratio.
This code interploates a digitized graph of t and s correction facors, provided by Benzing, Ali/Wings p.63





## Xfoil-Matlab Interface: documentation.

The following is a function that prints out the outputs in Matlab from an XFoil analysis. 

`[polar, foil] = xfoil(coord,alpha,Re,Mach,varargin)`, where:


- _coord_ is the Naca airfoil type. ex: 'NACA 0012'. The four parameters can be modified to obtain the desired shape.
- _alpha_ is the angle of attack.
- _Re_ is the Reynolds number. If Re>0, the viscous analysis is automatically initiated.
- _Mach_ is the Mach number.
- _varargin_ allows the addition of other inputs such as 'Npanels' or 'Niter', corresponding to the number of panels for the airfoil, and the number of maximum iterations for the analysis.

Outputs:

- A plot of the airfoil countour with the correponding pressure distribution along its chord. 'Cp vs c/x' plot.
- Cl and Cd (lift and drag coefficients) at the specified AoA.
- Cdp - parasite drag coefficient (skin friction drag) , if viscous analysis is performed.
- Cm - moment coefficient at the specified AoA.
- Top_xtr and Bot_xtr: the top and bottom flow transition locations as a fraction of the airfoil's chord.
- Cf (Coefficient of friction) values and graph with respect to the airfoil's chord. 

