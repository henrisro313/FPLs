!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Metropolis-Hastings Monte Carlo sampling of FPL Hamiltonian using single spin flips
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
use MyDefs
implicit none
  character(len=150), dimension(:), allocatable :: args
  integer :: m, n, ix, num_args, bonds, NMC, nthermal, ncycles, i, j, no, neffective
  integer :: in, im, seedno
  integer, allocatable :: hor(:,:), ver(:,:)
  real(dbl) :: temp
  real(dbl) :: Energy, J0, u(3), a
  character(len=150) :: outfilename
  integer, allocatable :: seed(:)
  character(len=20) :: ba1, ba2, runno

  !! Read command line arguments:
  num_args = command_argument_count()
  allocate(args(num_args))
  do ix = 1, 4
    call get_command_argument(ix,args(ix))
  end do

  !! Dimensions of the grid:
  read(args(1),*) m
  read(args(2),*) n
  read(args(3),*) temp
  read(args(4),*) no
  write (ba1,'(I2.2)') m
  write (ba2,'(I2.2)') n
  write (runno,'(I2.2)') no
  !print *, "m, n: ", m, n

  !!!!!!!! PROBLEM PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!
  J0 = 1.0d0
  NMC = 1E6
  bonds = 8*m*n-2*m-2*n
  nthermal = NMC/10
  ncycles = NMC - nthermal + 1
  neffective = FLOOR(ncycles/FLOAT(bonds))
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(hor(2*m, 2*n-1), ver(2*m-1, 2*n))

  seedno=2
  call random_seed(size=seedno)
  allocate(seed(seedno))
  call random_seed(get=seed)
  !write (*, *) seed
  call random_seed(put=seed)

  !!!!!!! Initial state !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !hor(:,:) = +1 !! Initialize in the ferromagnetic all up state
  !ver(:,:) = +1 !! Initialize in the ferromagnetic all up state
  do im = 1, 2*m
    do in = 1, 2*n-1
      !! draw random number in [-1, 1] and set hor
      call random_number(a)
      if (a < 0.5d0) then
          hor(im,in) = +1
      else
          hor(im,in) = -1
      endif
    end do ! in
  end do ! im
  do im = 1, 2*m-1
    do in = 1, 2*n
      !! draw random number in [-1, 1] and set ver
      call random_number(a)
      if (a < 0.5d0) then
          ver(im,in) = +1
      else
          ver(im,in) = -1
      endif
    end do ! in
  end do ! im
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call MCsampling(m,n, NMC, temp, J0, ver, hor)

  deallocate(hor, ver, seed)

end program main


subroutine MCsampling(m, n, NMC, T, J0, ver_old, hor_old)
use MyDefs
implicit none
  integer, intent(in) :: m, n, NMC
  integer, intent(in) :: ver_old(2*m-1, 2*n), hor_old(2*m, 2*n-1)
  real(dbl), intent(in) :: T, J0
  integer :: ver_p(2*m-1, 2*n), hor_p(2*m, 2*n-1), ver_n(2*m-1, 2*n), hor_n(2*m, 2*n-1), ver_cur(2*m-1, 2*n), hor_cur(2*m, 2*n-1)
  integer :: sample, i_r, j_r, j
  real(dbl) :: u(4), Etest, deltaE, Energy, Mu, PBoltzmann, E_cur, E_p, SIG1, SIG2, eps
  character(len=150) :: output
  eps = 1E-8

  call random_seed()

  ver_p = ver_old
  hor_p = hor_old
  E_p = Energy(hor_old, ver_old, m, n, J0)
  SIG1 = 0.0d0
  SIG2 = 0.0d0

  do sample = 1, NMC
    hor_cur = hor_p
    ver_cur = ver_p
    E_cur = E_p

    call random_number(u)
    if(u(1) < 0.5d0) then
        !! Horizontal bond spin flips:
        i_r = 1+FLOOR(2*m*u(2)) ! [1,2m]
        j_r = 1+FLOOR((2*n-1)*u(3)) ! [1,2n-1]

        !!! IMPLEMENT MANUAL HORIZONTAL FLIP ENERGY COST:
        if ( (i_r .EQ. 1) .AND. (j_r .EQ. 1) ) then ! upper left corner
            SIG1 = hor_cur(i_r,j_r)+ver_cur(i_r,j_r) - 2.0d0
            SIG2 = hor_cur(i_r,j_r)+hor_cur(i_r,j_r+1)+ver_cur(i_r,j_r+1) - 1.0d0
        else if ( (i_r .EQ. 1) .AND. (j_r .EQ. 2*n-1) ) then ! upper right corner
            SIG1 = hor_cur(i_r,j_r)+ver_cur(i_r,j_r+1) - 2.0d0
            SIG2 = hor_cur(i_r,j_r)+hor_cur(i_r,j_r-1)+ver_cur(i_r,j_r) - 1.0d0
        else if ( (i_r .EQ. 2*m) .AND. (j_r .EQ. 1) ) then ! lower left corner
            SIG1 = hor_cur(i_r,j_r)+ver_cur(i_r-1,j_r) - 2.0d0
            SIG2 = hor_cur(i_r,j_r)+hor_cur(i_r,j_r+1)+ver_cur(i_r-1,j_r+1) - 1.0d0
        else if ( (i_r .EQ. 2*m) .AND. (j_r .EQ. 2*n-1) ) then ! lower right corner
            SIG1 = hor_cur(i_r,j_r)+ver_cur(i_r-1,j_r+1) - 2.0d0
            SIG2 = hor_cur(i_r,j_r)+hor_cur(i_r,j_r-1)+ver_cur(i_r-1,j_r) - 1.0d0

        else if ( (i_r .EQ. 1) .AND. (j_r < 2*n-1) .AND. (j_r > 1) ) then ! upper edge
            SIG1 = hor_cur(i_r,j_r)+hor_cur(i_r,j_r-1)+ver_cur(i_r,j_r) - 1.0d0
            SIG2 = hor_cur(i_r,j_r)+hor_cur(i_r,j_r+1)+ver_cur(i_r,j_r+1) - 1.0d0
        else if ( (i_r .EQ. 2*m) .AND. (j_r < 2*n-1) .AND. (j_r > 1) ) then ! lower edge
            SIG1 = hor_cur(i_r,j_r)+hor_cur(i_r,j_r-1)+ver_cur(i_r-1,j_r) - 1.0d0
            SIG2 = hor_cur(i_r,j_r)+hor_cur(i_r,j_r+1)+ver_cur(i_r-1,j_r+1) - 1.0d0
        else if ( (i_r > 1) .AND. (i_r < 2*m) .AND. (j_r .EQ. 1) ) then ! left edge
            SIG1 = hor_cur(i_r,j_r)+ver_cur(i_r,j_r)+ver_cur(i_r-1,j_r) - 1.0d0
            SIG2 = hor_cur(i_r,j_r)+hor_cur(i_r,j_r+1)+ver_cur(i_r,j_r+1)+ver_cur(i_r-1,j_r+1)
        else if ( (i_r > 1) .AND. (i_r < 2*m) .AND. (j_r .EQ. 2*n-1) ) then ! right edge
            SIG1 = hor_cur(i_r,j_r)+ver_cur(i_r,j_r+1)+ver_cur(i_r-1,j_r+1) - 1.0d0
            SIG2 = hor_cur(i_r,j_r)+hor_cur(i_r,j_r-1)+ver_cur(i_r,j_r)+ver_cur(i_r-1,j_r)

        else ! bulk
            SIG1 = hor_cur(i_r,j_r)+hor_cur(i_r,j_r-1)+ver_cur(i_r,j_r)+ver_cur(i_r-1,j_r)
            SIG2 = hor_cur(i_r,j_r)+hor_cur(i_r,j_r+1)+ver_cur(i_r,j_r+1)+ver_cur(i_r-1,j_r+1)
        end if
        Etest = E_cur + 4.0d0*J0*( 2.0d0 - (SIG1 + SIG2)*hor_cur(i_r,j_r) )
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        DeltaE = Etest-E_cur
        PBoltzmann = MINVAL( (/ 1.0d0, EXP(-DeltaE/T) /) )
        if(u(4) < PBoltzmann) then
            hor_cur(i_r, j_r) = -1*hor_cur(i_r, j_r) ! flipping horizontal spin
            E_cur = Etest
        end if

    else
        !! Vertical bond spin flips:
        i_r = 1+FLOOR((2*m-1)*u(2)) ! [1,2m-1]
        j_r = 1+FLOOR(2*n*u(3)) ! [1,2n]

        !!! IMPLEMENT MANUAL VERTICAL FLIP ENERGY COST:
        if ( (i_r .EQ. 1) .AND. (j_r .EQ. 1) ) then ! upper left corner
            SIG1 = ver_cur(i_r,j_r)+hor_cur(i_r,j_r) - 2.0d0
            SIG2 = ver_cur(i_r,j_r)+ver_cur(i_r+1,j_r)+hor_cur(i_r+1,j_r) - 1.0d0
        else if ( (i_r .EQ. 1) .AND. (j_r .EQ. 2*n) ) then ! upper right corner
            SIG1 = ver_cur(i_r,j_r)+hor_cur(i_r,j_r-1) - 2.0d0
            SIG2 = ver_cur(i_r,j_r)+ver_cur(i_r+1,j_r)+hor_cur(i_r+1,j_r-1) - 1.0d0
        else if ( (i_r .EQ. 2*m-1) .AND. (j_r .EQ. 1) ) then ! lower left corner
            SIG1 = ver_cur(i_r,j_r)+hor_cur(i_r+1,j_r) - 2.0d0
            SIG2 = ver_cur(i_r,j_r)+ver_cur(i_r-1,j_r)+hor_cur(i_r,j_r) - 1.0d0
        else if ( (i_r .EQ. 2*m-1) .AND. (j_r .EQ. 2*n) ) then ! lower right corner
            SIG1 = ver_cur(i_r,j_r)+hor_cur(i_r+1,j_r-1) - 2.0d0
            SIG2 = ver_cur(i_r,j_r)+ver_cur(i_r-1,j_r)+hor_cur(i_r,j_r-1) - 1.0d0

        else if ( (i_r > 1) .AND. (i_r < 2*m-1) .AND. (j_r .EQ. 1) ) then ! left edge
            SIG1 = ver_cur(i_r,j_r)+ver_cur(i_r-1,j_r)+hor_cur(i_r,j_r) - 1.0d0
            SIG2 = ver_cur(i_r,j_r)+ver_cur(i_r+1,j_r)+hor_cur(i_r+1,j_r) - 1.0d0
        else if ( (i_r > 1) .AND. (i_r < 2*m-1) .AND. (j_r .EQ. 2*n) ) then ! right edge
            SIG1 = ver_cur(i_r,j_r)+ver_cur(i_r-1,j_r)+hor_cur(i_r,j_r-1) - 1.0d0
            SIG2 = ver_cur(i_r,j_r)+ver_cur(i_r+1,j_r)+hor_cur(i_r+1,j_r-1) - 1.0d0
        else if ( (i_r .EQ. 1) .AND. (j_r > 1) .AND. (j_r < 2*n) ) then ! upper edge
            SIG1 = ver_cur(i_r,j_r)+hor_cur(i_r,j_r-1)+hor_cur(i_r,j_r) - 1.0d0
            SIG2 = ver_cur(i_r,j_r)+ver_cur(i_r+1,j_r)+hor_cur(i_r+1,j_r-1)+hor_cur(i_r+1,j_r)
        else if ( (i_r .EQ. 2*m-1) .AND. (j_r > 1) .AND. (j_r < 2*n) ) then ! lower edge
            SIG1 = ver_cur(i_r,j_r)+hor_cur(i_r+1,j_r-1)+hor_cur(i_r+1,j_r) - 1.0d0
            SIG2 = ver_cur(i_r,j_r)+ver_cur(i_r-1,j_r)+hor_cur(i_r,j_r-1)+hor_cur(i_r,j_r)

        else ! bulk
            SIG1 = ver_cur(i_r,j_r)+ver_cur(i_r-1,j_r)+hor_cur(i_r,j_r)+hor_cur(i_r,j_r-1)
            SIG2 = ver_cur(i_r,j_r)+ver_cur(i_r+1,j_r)+hor_cur(i_r+1,j_r)+hor_cur(i_r+1,j_r-1)
        end if
        Etest = E_p + 4.0d0*J0*( 2.0d0 - (SIG1 + SIG2)*ver_cur(i_r,j_r) )
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        DeltaE = Etest-E_cur
        PBoltzmann = MINVAL( (/ 1.0d0, EXP(-DeltaE/T) /) )
        if(u(4) < PBoltzmann) then
           ver_cur(i_r, j_r) = -1*ver_cur(i_r, j_r) ! flipping vertical spin
           E_cur = Etest
        end if
    end if

    E_p = E_cur ! update E_p for next step
    hor_p = hor_cur ! update hor_p for next step
    ver_p = ver_cur ! update ver_p for next step

  end do ! sample

  if (E_p < eps) then
    !! Write final state to text file:
    output = 'FinalState_hort.txt'
    open(unit=16,file=trim(output))
    do j = 1, 2*n-1
      write(16,*) hor_cur(:,j)
    end do ! j
    close(16)
    output = 'FinalState_vert.txt'
    open(unit=17,file=trim(output))
    do j = 1, 2*n
      write(17,*) ver_cur(:,j)
    end do ! j
    close(17)
  end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine MCsampling

function Energy(hor, ver, m, n, J0)
use MyDefs
implicit none
 integer, intent(in) :: m, n
 integer, intent(in) :: hor(2*m, 2*n-1), ver(2*m-1, 2*n)
 integer :: i, j
 real(dbl) :: Energy, J0

 Energy = 0.0d0
 !! Corner term:
 Energy = Energy + J0*(hor(1,1)+ver(1,1) - 2.0d0)**2 + J0*(hor(1,2*n-1)+ver(1,2*n) - 2.0d0)**2 &
          + J0*(hor(2*m,1)+ver(2*m-1,1) - 2.0d0)**2 + J0*(hor(2*m,2*n-1)+ver(2*m-1,2*n) - 2.0d0)**2
 !! Edge term:
 do j = 2, 2*n-1
    Energy = Energy + J0*(hor(1,j-1)+hor(1,j)+ver(1,j) - 1.0d0 )**2 &
              + J0*(hor(2*m,j-1)+hor(2*m,j)+ver(2*m-1,j) - 1.0d0 )**2
 end do ! j
 do i = 2, 2*m-1
    Energy = Energy + J0*(ver(i-1,1)+ver(i,1)+hor(i,1) - 1.0d0 )**2 &
              + J0*(ver(i-1,2*n)+ver(i,2*n)+hor(i,2*n-1) - 1.0d0 )**2
 end do ! i
 !! Bulk term:
 do i = 2, 2*m-1
   do j = 2, 2*n-1
     Energy = Energy + J0*(hor(i,j)+hor(i,j-1)+ver(i-1,j)+ver(i,j))**2
   end do ! j
end do ! i

end function Energy

function unit_to_interval(x_bar, a, b)
use MyDefs
implicit none
  real(dbl), intent(IN) :: x_bar, a, b
  real(dbl) :: unit_to_interval
  unit_to_interval = a+(b-a)*x_bar
end function unit_to_interval

subroutine init_random_seed()
   integer :: i, n, clock
   integer, dimension(:), allocatable :: seed
   call random_seed(size = n)
   allocate(seed(n))
   call system_clock(count=clock)
   seed = clock + 37 * (/ (i - 1, i = 1, n) /)
   call random_seed(put = seed)
   deallocate(seed)
end subroutine init_random_seed ! period for random_number : 2^32 = 4294967296
