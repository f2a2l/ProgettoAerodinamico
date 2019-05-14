# ProgettoAerodinamico

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
