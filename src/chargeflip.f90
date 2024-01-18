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

   call ccpfyp

   !read in structure factor
   call mtzin


   !normalize structure factor
   !parameter: 1st : number of bins used in wilson statistics -- based
   !on sin(theta)/lambda
   !           2ed : the lowest resolution permitted in wilson
   !           statistics
   !           3th : average number of reflections in one bin in ecal
!   call ecal



   !main function for charge flipping
   !parameter: 1st : number of bins for wilson plot
   !           2ed : the lowest resolution truncated for wilson plot
   !           3th : fraction of density zone
   !           4th : number of charge-flipping iterations
   !           5th : try number of different initial random phases
   call charge_flipping

   stop 'Normally Terminated!'
end program main

