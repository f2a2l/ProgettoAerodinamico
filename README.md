# ProgettoAerodinamico


## Hess Smith: documentation

__Single airfoil case__: `solverHS(AirfNumb, aname, alpha)`, where:
- _AirfNumb_ is the number of airfoils (thus 1).
- _aname_ is a row vector, each component corresponding to the parameters describing airfoil geometry.
- _alpha_ is the angle of attack.

Example:
```MATLAB
arflPar = [0.3, 0.6, 0, 0, 0.3, 0.12, 0.3, 1.5];
solverHS(1, arflPar, 3);
```

__Multiple airfoil case__: `solverHS(AirfNumb, aname, alpha, dist, crel)`, where:
- _AirfNumb_ is the number of airfoils.
- _aname_ is a matrix; each row corresponds to an airfoil (no. rows = AirfNumb) and contains the parameters needed to define the airfoil.
- _alpha_ is a row vector, containing angles of attack of the airfoils (all angles are expressed in an absolute frame).
- _dist_ is a matrix, each row corresponding to an airfoil; columns correspond to the x and y components of the position of the leading edge, expressed in a body frame attached to the first airfoil (x axis parallel to chord, centered on the leading edge). Notice that the first airfoil is omitted (as its row would always be [0,0]), thus effectively reducing the number of rows to AirfNumb-1.
- _crel_ is a column vector; each component corresponds to the chord of the corresponding airfoil, adimensionalised on the chord of the first airfoil. Again, first airfoil is omitted (since its row would always be [1]), and the length of the vector is AirfNumb-1.

Keep in mind that the chord of the first airfoil is always 1 (thus the "relative chord" of the other airfoils is effectively ther chord).

Example:
```MATLAB
arflPar = [0.3, 0.6, 0, 0, 0.3, 0.12, 0.3, 1.5;
           0.3, 0.6, 0, 0, 0.3, 0.12, 0.3, 1.5];
solverHS(2, arflPar, [3, 6], [1.05, -0.05], 0.3);
```


## Coding guidelines
- __never__ push code to 'origin master'
- always run code from `ProgettoAerodinamico` folder; do not open subfolders in MATLAB
    - to do this, run the script `addPaths`; this will make all files and functions in subfolders available from the folder `ProgettoAerodinamico`
    - it might be useful to run `addPaths` even before programming, so that MATLAB linter can correctly access all files in subfolders
- keep in mind that commands such `checkout` and `pull` __overwrite your files__
- commit often and use meaningful commit messages
- __test your code__ before committing

## Useful stuff
- easy git guide: http://rogerdudler.github.io/git-guide/
- git cheatsheet: https://github.github.com/training-kit/
- how to ignore files with `.gitignore`: https://git-scm.com/docs/gitignore
