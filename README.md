# A Modified Phase-retrieval Algorithm
  A modified phase-retrieval algorithm, incorporating the π-half variant of charge flipping algorithm and direct-methods tangent formula into the relaxed alternating averaged reflection phase-retrieval (RAAR) algorithm, has been developed for the determination of anomalously scattering substructures from single-wavelength anomalous diffraction (SAD) data of macromolecules. 

## Requirement
- Fortran90
- Linux operating system
- [fgsl/gsl](http://www.gnu.org/software/gsl/) for random number generation
- [FFTW3](http://www.fftw.org) library for the fast Fourier transform
- [CCP4 subroutine libraries](https://www.ccp4.ac.uk/html/index.html) for basic crystallographic operations

## testting
1. Compiling  
Edit the [configure](configure) file and update the paths for  fgsl/gsl, FFTW3 and CCP4 subroutine libraries in the file. Then, run:
`./congfigure`
```$ npm install gitbook-plugin-copy-code-button```

To use it in your book, add this to book.json:
```
{
    "plugins": ["copy-code-button"]
}
```

## The dual-space iterative framework of the  Modified Phase-retrieval Algorithm
![image](fig/the_iterative_framework.png)
