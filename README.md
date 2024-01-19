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
Edit the [configure](configure) file and update the paths for  **fgsl/gsl**, **FFTW3** and **CCP4 subroutine libraries** in the file. Then, run:
```bash
./conafigure
```
2. A test case, 5OQ2  
There are two inputing parameters for the modified phase-retrieval algorithm:
	- **mtz file**. This should contain the conventional (CCP4) asymmetric unit of data with anomalous information (including F+/F- or I+/I-);  
	- **1 or 0**. "1" indicates refining phases in the dual-space iterative cycle; "0" indicate no refinement./<br>

The following is a case for protein with PDB entry 5OQ2:    
```bash
cd exp-test/7_5oq2
./run.sh 5oq2-sf.mtz 1
```
3. 
## The dual-space iterative framework of the  Modified Phase-retrieval Algorithm
<p align="center">
<img align="middle" src="fig/the_iterative_framework.png" width="500" alt="trg"/>
</p>
