!=====================================================================================
! (c) IHEP 2023. All rights reserved.
! Author: Geng Zhi
! Institute of High Energy Physics, Chinese Academy of Science (IHEP, CAS). 
! If you have any problem with the program, please contact author with the
! following 
! email: gengz@ihep.ac.cn
!======================================================================================

program main
   use functions
   implicit none

   real,allocatable :: x(:), y(:), z(:), asf(:,:)
   character(4),allocatable :: atms(:)

   call ccpfyp
   call xyzinit

   !read in atom coordinates from pdb file
   call pdbin(x,y,z,atms)

   !atom symmetry expansion
   call expand_atoms(x,y,z,atms,asf)
   
   !calculate structure factor from small molecule coordinates file cif
   call sfall(x,y,z,asf)

   stop 'Normally Terminated!'
end program main
