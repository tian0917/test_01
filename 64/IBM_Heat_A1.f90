! MONOLITHIC PROJECTION METHOD WITH STAGGERED TIME STEP (Staggered_MPM)~~~
! IMMERSED BOUNDARY PROJECTION METHOD 
! 2D FORCED CONVECTION
!   
!   BOUNDARY CONDITION:
!   (1) U, V, T: HOMOGENEOUS DIRICHLET BOUNDARY in INLET
!   (2) U, V, T: NEUMANN BOUNDARY CONDITION IN OUTLET AND Y-DIRECTION
!   (2) P: NEUMANN BOUNDARY in X,Y-direction
!   (3) T=1 AT THE IB SURFACE
!                                                                           TIANTIAN XU
!                                                           COMPUTAIONAL FLUID DYNAMICS
!                                   DEPARTMENT OF COMPUTATIONAL SCIENCE AND ENGINEERING
!                                                                     YONSEI UNIVERSITY
!
!                                                                             2020/11/12

    Module Global
    implicit none  

    ! CONSTANT PARAMTERS
    integer, parameter :: ND = 3
    integer, parameter :: UNIFORM1 = 1, UNIFORM2 = 1
    double precision, parameter :: GAMMA1 = 3.d0, GAMMA2 = 3.d0
    double precision, parameter :: PI = acos(-1.d0)
    double precision, parameter :: dPTOL = 1.e-13

    ! PHYSICAL PARAMETERS
    double precision, parameter :: Thetam = 0.
    double precision, parameter :: Re = 100.0, Pr = 0.71, Ra = 1.0e6
    
    ! FORCED CONVECTION
!     double precision, parameter :: Cb = 0.0
!     double precision, parameter :: Cm = 1.0/Re, Ce = 1.0/(Re*Pr)

    ! NATURAL CONVECTION
    double precision, parameter :: Cb = Ra*Pr
    double precision, parameter :: Cm = Pr, Ce = 1.0

    double precision :: cfl = 3.0, cflm
    double precision :: umax, vmax 

    integer :: nprn = 10000, imore 


    ! ITERATION STEPS
    integer, parameter :: Tmax = 50000, DpNum = 1000, RelaxNum = 7, Tinter = 5000
    integer, parameter :: NK = 10


    ! COMPUTATIONAL SIZE FOR SPACE AND TIME DISCRETIZATIONS
    integer, parameter :: MaxLevel = 3     
    integer, parameter :: N1 = 401, N2 = 401
    integer, parameter :: N1m = N1-1, N2m = N2-1
    integer, parameter :: N4 = 3300, N5 = 400 

    ! SETTINGS FOR TIME STEP
    double precision :: dt = 1.0e-06 
    double precision, parameter :: MaxCFL = 5.0 


    ! DOMAIN SIZE FOR THE PHYSICAL PROBLEMS
    double precision, parameter :: r = 0.02, ALX = 1.0, ALY = 1.0

    double precision, parameter :: XMIN = -0.5d0, XMAX = XMIN + ALX
    double precision, parameter :: YMIN = -0.5d0, YMAX = YMIN + ALY

    double precision, parameter :: VOLUME = ALX*ALY


    ! INDEX VECTORS
    integer, dimension(N1m) :: IM_U, IM_V, IM_E, IP_A, IP_E 
    integer, dimension(N2m) :: JM_U, JM_V, JM_E, JP_A, JP_E


    ! MESH VECTORS
    double precision, dimension(0:N1) :: X, DX, DMX
    double precision, dimension(0:N2) :: Y, DY, DMY
    double precision, dimension(0:N1,MaxLevel) :: MGDX
    double precision, dimension(0:N2,MaxLevel) :: MGDY


    ! BOUNDARY CONDITION VECTORS
    double precision, dimension(0:N2, 1:2) :: UBC1, UBC2              ! in and out
    double precision, dimension(0:N1, 1:2) :: UBC3, UBC4              ! top and bottom
    double precision, dimension(0:N2) :: THETABC1, THETABC2           ! in and out
    double precision, dimension(0:N1) :: THETABC3, THETABC4           ! top and bottom

    ! VARIABLE NOTATIONS FOR PREVIOUS TIME STEP
    double precision, dimension(0:N1, 0:N2) :: Up, Vp, Pp, THETAp


    ! LAGRANGIAN PART
    integer :: NL, Nc 
    integer :: NLP(0:N4, N5)
    integer :: NP1(0:N4), NP2(0:N4), NNP(0:N4)
 
    double precision :: DS, DV, FSUM
    double precision, dimension(8) :: X0, Y0  
    double precision, dimension(0:N4) :: FP1, FP2 
    double precision, dimension(0:N4) :: UD1, UD2, QD  
    double precision, dimension(0:N4) :: ULAG1, ULAG2, ELAG


    ! TIME 
    double precision :: time(8) = 0.0 

    character*256 :: dirname

    integer :: TimeStep

    end Module Global
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


    !<<<<<<<<<<<<<<< MAIN PROGRAM >>>>>>>>>>>>>>>
    program main
    use omp_lib
    use Global
    implicit none

    character (len=90) :: filename, dirname1, dirname2

    integer :: i, j, l, m 
    integer :: ip, jp

    double precision :: t, tlength
    double precision :: t1, t2, ABS_DivU  

    ! Definition of unknown variables
    double precision, allocatable, dimension(:, :) :: U, V, P, THETA
    double precision, allocatable, dimension(:, :) :: DivU

    double precision, allocatable, dimension(:, :) :: BIBE1, BIBE2, BIBET 


    allocate(U(0:N1, 0:N2), V(0:N1, 0:N2), P(0:N1, 0:N2), THETA(0:N1, 0:N2))
    allocate(DivU(1:N1m, 1:N2m))

    allocate(BIBE1(0:N4,0:N4), BIBE2(0:N4,0:N4), BIBET(0:N4,0:N4))


    call get_environment_variable('PWD', dirname1)  
    dirname2 = '/scratch'
    dirname = './result/'  !trim(dirname2)//dirname1(6:)


   
    
    write (filename, '("LAG_POINTS.PLT" )') 
    open  (unit=16,file=trim(dirname)//filename)

    write (filename, '("COEF.PLT" )') 
    open  (unit=17,file=trim(dirname)//filename)

    write (filename, '("TIME.PLT" )' )
    open  (unit=19,file=trim(dirname)//filename)

    write (filename, '("HORI_U.PLT")')
    open  (unit=20,file=trim(dirname)//filename)      

    write (filename, '("HORI_V.PLT")')
    open  (unit=21,file=trim(dirname)//filename)

    write (filename, '("HORI_T.PLT")')
    open  (unit=22,file=trim(dirname)//filename)    

    write (filename, '("HORI_P.PLT")')
    open  (unit=23,file=trim(dirname)//filename)    

    write (filename, '("UAV.PLT")')
    open  (unit=28,file=trim(dirname)//filename)

    write (filename, '("EAV.PLT")')
    open  (unit=29,file=trim(dirname)//filename)
    
    write (filename, '("UT_HISTORY1.PLT")')
    open  (unit=30,file=trim(dirname)//filename)

    write (filename, '("UT_HISTORY2.PLT")')
    open  (unit=31,file=trim(dirname)//filename)

    write (filename, '("UT_HISTORY3.PLT")')
    open  (unit=32,file=trim(dirname)//filename)

    write (filename, '("UT_HISTORY4.PLT")')
    open  (unit=33,file=trim(dirname)//filename)

    write (filename, '("FSUM_HISTORY.PLT")')
    open  (unit=34,file=trim(dirname)//filename)

    call CPU_TIME(t1)
!     t1 = omp_get_wtime()


    call MESH 
    call INDICES
    call INIUPT(U, V, P, THETA)
!     call CoePoisson

    ! [INITIALIZATION BEGINS----------------------------
    Up = U; Vp = V; Pp = P; THETAp = THETA

    dt = cfl*DS 

    ! call Cfl_Con 
    call Cfl_Check

    call PrePrinting 
    call BCOND(t)  
    ! ----------------------------INITIALIZATION ENDS]

    write(*,*)'------------------------------------------------------------------------------'
    write(*,*)'>>>>>>>>>The preparation for Staggered_MPM has been finished!'
    write(*,*)'------------------------------------------------------------------------------'
    write(*,*) 


    ! FIND LAGRANGIAN POINTS
    call LAG_POINTS 

    ! CONSTRUCT BIBE MATRIX FOR IBM 
    call BIBEMAT(BIBE1, BIBE2, BIBET)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! [TIME DO LOOP BEGINS----------------------------
    do TimeStep = 1, Tmax

        write(*,*) '==================================================='
        write(*,*) 'Timestep, dt ', Timestep, dt 
        write(*,*) '==================================================='
        write(*,*) 

        t = t + dt

        call GETUPT(THETA, U, V, P, BIBE1, BIBE2, BIBET, t)

        Up = U
        Vp = V
        Pp = P
        THETAp = THETA

        ! call CheckCFL(t)  
        ! call Cal_dt(t, CFL)
      
     
        call CheckDiv(DivU)
        call Cfl_Check 
        call Cfl_Con
    
        write(*,*)
        write(*,*) 'MAX DIVERGENCE = ', sqrt(sum(abs(DivU)**2)/dble(N1m)/dble(N2m))


        ! PRINT THE FILES AT CERTAIN TIME
        if(mod(Timestep, nprn) == 0) then 
            imore = Timestep/nprn 
            call writeup
        end if 


        !TIME HISTORY OF VELOCITY AND TEMPERATURE
        write(30,1024) t, Up(N1m/4*3, N2m/2), Vp(N1m/4*3, N2m/2), THETAp(N1m/4*3, N2m/2)
        write(31,1024) t, Up(N1m/2, N2m/4*3), Vp(N1m/2, N2m/4*3), THETAp(N1m/2, N2m/4*3) 
        write(32,1024) t, Up(N1m/4, N2m/2), Vp(N1m/4, N2m/2), THETAp(N1m/4, N2m/2)
        write(33,1024) t, Up(N1m/2, N2m/4), Vp(N1m/2, N2m/4), THETAp(N1m/2, N2m/4)
    


    end do
    ! ------------------------------TIME DO LOOP ENDS] 

1024 format(4(3X, F22.10))

    call CPU_TIME(t2) 
!     t2 = omp_get_wtime()

    close(30)
    close(31)
    close(32)
    close(33)

    ! y=0
    do j = 1, N2 

        if(abs(Y(j)) < 1.0E-10) then 
            do i = 1, N1 
                write(20,*) X(i), U(i,j)
                write(21,*) X(i), V(i,j)
                write(22,*) X(i), THETA(i,j)
                write(23,*) X(i), P(i,j)
            end do 
        end if 

        if(abs(Y(j) + Y(j+1)) < 1.0E-10) then 
            do i = 1, N1 
                write(20,*) X(i), 0.5*(U(i,j) + U(i,j+1))
                write(21,*) X(i), 0.5*(V(i,j) + V(i,j+1))
                write(22,*) X(i), 0.5*(THETA(i,j) + THETA(i,j+1))
                write(23,*) X(i), 0.5*(P(i,j) + P(i,j+1))
            end do 
        end if 

    end do 

    close(20)
    close(21)
    close(22)
    close(23)    



    write(19,*) t2 - t1 
    write(19,*) TIME(1)
    write(19,*) TIME(2)
    write(19,*) TIME(3)
    write(19,*) TIME(4)
    write(19,*) TIME(5)
    write(19,*) TIME(6)
    write(19,*) TIME(7)
    write(19,*) TIME(8)




2000 format(3e14.6)
2002 format(6e15.7)  
2003 format(2i4, 3(3X, F12.7))
!     call PrintDiv(DivU)


    call output(t)  

    call streamfunction  


    close(16); close(17); close(19)


    write(*,*) 
    write(*,*) 
    write(*,*) 'COMPUTATION FINISH!'
    write(*,*) '-------------------'
    write(*,*) 
    write(*,*) 

    write(*,*) 'TOTAL TIME     ', t2 - t1 
    write(*,*) 'FEXTRAPOL      ', TIME(1)
    write(*,*) 'ENERGY TIME    ', TIME(2)
    write(*,*) 'MOMENTUM TIME  ', TIME(3)
    write(*,*) 'POISSON TIME   ', TIME(7)
    write(*,*) 
    write(*,*) 



    deallocate(U, V, P, THETA)
    deallocate(DivU)

    deallocate(BIBE1, BIBE2, BIBET)

   

    end program main
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       



    !<<<<<<<<<<<<<<< MESH >>>>>>>>>>>>>>>
    subroutine MESH
    use Global
    implicit none

    integer :: i, j
    integer :: level
    double precision :: temp



    ! X-DIRECTION
    X(0) = XMIN
    if(UNIFORM1 == 1) then 
        temp = (XMAX-XMIN)/dble(N1m)
        DS = temp 
        do i = 1, N1
            X(i) = dble(i-1)*temp + XMIN
        end do

    else
        do i = 1, N1
            X(i) = (XMAX-XMIN)*0.5*(1. + tanh(0.5*GAMMA1*(2.*dble(i-1)/dble(N1m)-1.0))/tanh(GAMMA1*0.5) ) + XMIN
        end do

    end if
    


    ! Y-DIRECTION
    Y(0) = YMIN
    if(UNIFORM2 == 1) then
        temp = (YMAX-YMIN)/dble(N2m)
        do j = 1, N2
            Y(j) = dble(j-1)*temp + YMIN
        end do

    else
        do j = 1, N2
            Y(j) = (YMAX-YMIN)*0.5*(1. + tanh(0.5*GAMMA2*(2.*dble(j-1)/dble(N2m)-1.0))/tanh(GAMMA2*0.5)) + YMIN
        end do

    end if


    ! GENERATE GRID SPACINGS 
    ! x-direction
    DX(0) = 0.
    do i = 1, N1m
        DX(i) = X(i+1) - X(i)
    end do
    DX(N1) = 0.

    DMX(0) = 0.
    do i = 1, N1
        DMX(i) = 0.5*(DX(i-1)+DX(i))
    end do


    ! y-direction 
    DY(0) = 0.
    do j = 1, N2m
        DY(j) = Y(j+1) - Y(j)
    end do
    DY(N2) = 0.

    DMY(0) = 0.
    do j = 1, N2
        DMY(j) = 0.5*(DY(j-1)+DY(j))
    end do


    ! FOR THE MULTIGRID PPOISSON SOLVER
    MGDX(0:N1,1:MaxLevel) = 0.
    MGDY(0:N2,1:MaxLevel) = 0.

    MGDX(0:N1,1) = DX(0:N1)
    MGDY(0:N2,1) = DY(0:N2)

    !  GENERATE MGDX, MGDY, MGDZ
    do level = 2, MaxLevel
        
        do i = 1, N1m/(2**(level-1))
            MGDX(i, level) = MGDX(2*i-1, level-1) + MGDX(2*i, level-1)
        end do

        do j = 1, N2m/(2**(level-1))
            MGDY(j, level) = MGDY(2*j-1, level-1) + MGDY(2*j, level-1)
        end do    

    end do


    return
    end subroutine MESH 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


    !<<<<<<<<<<<<<<< INDICES >>>>>>>>>>>>>>>
    subroutine INDICES 
    use Global
    implicit none

    integer :: i, j

    ! For VELOCITY IN X-DIRECRION
    do i = 1, N1m
        IP_A(i) = i + 1
        IM_U(i) = i - 1
        IM_V(i) = i - 1
    end do
    IP_A(N1m) = N1m 
    IM_U(2  ) = 2 
    IM_V(1  ) = 1 

    ! For VELOCITY IN Y-DIRECRION
    do j = 1, N2m
        JP_A(j) = j + 1
        JM_U(j) = j - 1
        JM_V(j) = j - 1
    end do
    JP_A(N2m) = N2m
    JM_U(1   ) = 1
    JM_V(2   ) = 2
    

    ! For TEMPERATURE X-DIRECTION
    do i = 1, N1m 
        IP_E(i) = i + 1 
        IM_E(i) = i - 1
    end do 
    IP_E(N1m) = N1m 
    IM_E(1  ) = 1

    ! For TEMPERATURE Y-DIRECRION
    do j = 1, N2m
        JP_E(j) = j + 1
        JM_E(j) = j - 1
    end do
    JP_E(N2m) = N2m
    JM_E(1  ) = 1


    return    
    end subroutine INDICES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 



    !<<<<<<<<<<<<<<< INIUPT >>>>>>>>>>>>>>>
    subroutine INIUPT(U, V, P, THETA)
    use omp_lib
    use Global
    implicit none

    integer :: i, j

    double precision, dimension(0:N1, 0:N2) :: U, V, P, THETA
!     double precision, dimension(0:N1, 0:N2) :: RanNum


!     call random_seed()
!     call random_number(RanNum)

    ! INITIALIZATION OF PHYSICAL FIELD
    !$omp parallel do
    do j = 0, N2
    do i = 0, N1

        U(i,j) = 0.
        V(i,j) = 0.
        P(i,j) = 0.
        THETA(i,j) = 0.

    end do
    end do
    !$omp end parallel do    

!     open(101, file='T_CENTER.PLT') 
!     open(102, file='P_CENTER.PLT')
!     open(103, file='U_EDGE.PLT')
!     open(104, file='V_EDGE.PLT')
!     do j = 0, N2 
!         do i = 0, N1 
!             read(101,*) THETA(i,j)
!             read(102,*) P(i,j)
!             read(103,*) U(i,j)
!             read(104,*) V(i,j)
!         end do 
!     end do 
!     close(101)
!     close(102)
!     close(103)
!     close(104)



    return    
    end subroutine INIUPT
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  

!     !<<<<<<<<<<<<<<< CoePoisson >>>>>>>>>>>>>>>
!     subroutine CoePoisson
!     use Global
!     implicit none

!     integer :: ip, jp


!     do j = 1, N2m
!     jp = j + 1
!     do i = 1, N1m
!     ip = i + 1

!         alpha(i,j,1) = 1./DX(i)/DMX(i )
!         alpha(i,j,2) = 1./DX(i)/DMX(ip)
!         alpha(i,j,3) = 1./DY(j)/DMY(j )
!         alpha(i,j,4) = 1./DY(j)/DMY(jp)

!         ! xz-Periodic
!         if(j .eq. 1  )  alpha(i,j,3) = 0.
!         if(j .eq. N2m)  alpha(i,j,4) = 0.


!         alpha(i,j,5) = -( alpha(i,j,1) + alpha(i,j,2) + alpha(i,j,3) + alpha(i,j,4))

!     end do
!     end do

    
!     return
!     end  subroutine CoePoisson
! !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!   ===============================================================
!               -IMMERSED BOUNDARY METHOD-
!   ===============================================================

    !<<<<<<<<<<<<<<< LAG_POINTS >>>>>>>>>>>>>
    subroutine LAG_POINTS
    use global 
    implicit none 

    integer :: i, j 
    integer :: L, K 
    double precision :: dthe 
    double precision :: distance, distancex, distancey 
    double precision :: xl0, yl0 
    

    X0 = (/ -7.0/18.0, -5.0/18.0, -1.0/6.0, -1.0/18.0, 1.0/18.0, 1.0/6.0, 5.0/18.0, 7.0/18.0 /)
    Y0 = (/ -7.0/18.0, -5.0/18.0, -1.0/6.0, -1.0/18.0, 1.0/18.0, 1.0/6.0, 5.0/18.0, 7.0/18.0 /)
    Nc = 8
    
    ! GENERATE LAGRANGIAN POINTS 
    NL = int(PI*2.0*r/DS) + 1
    dthe = 2.0*PI/dble(NL)

    L = 0 
    do j = 1, Nc
        yl0 = Y0(j)
        do i = 1, Nc 
            xl0 = X0(i)

            do K = 0, NL-1
                FP1(L) = xl0 + r*cos(dble(L)*dthe)
                FP2(L) = yl0 + r*sin(dble(L)*dthe)
                L = L + 1 
            end do 
        end do 
    end do 
    NL = NL*64


    ! PRINT LAG POINTS 
    do L = 0, NL-1
        write(16,*) FP1(L), FP2(L)
    end do 
    

    write(*,*) 'NUMBER OF WHOLE LAG POINTS = ', NL 


    DV = DS**2.0 
    
    ! INITIALIZE VELICITIES AT LAG POINTS 
    UD1(0:NL-1) = 0.0 
    UD2(0:NL-1) = 0.0 


    ! INITIALIZE TEMPERATURE AT LAG POINTS
    QD(0:NL-1) = 1.0  

    ! FIND NEIGHBORHOOD EULERIAN POINTS FROM LAG POINTS 
    do L = 0, NL-1 
        
        ! X-DIRECTION
        distance = abs(FP1(L) - X(0))
        do i = 1, N1 
            distancex = abs(FP1(L) - X(i))
            if(distancex < distance) then 
                distance = distancex 
                NP1(L) = i 
            end if 
        end do 

        ! Y-DIRECTION
        distance = abs(FP2(L) - Y(0))
        do j = 1, N2 
            distancey = abs(FP2(L) - Y(j))
            if(distancey < distance) then 
                distance = distancey 
                NP2(L) = j 
            end if 
        end do 


        ! FIND POINTS WITHIN (-3*DS,3*DS) FOR LAG 
        distance = dble(ND)*DS 

        do K = 0, NL-1 
            distancex = abs(FP1(K) - FP1(L))
            distancey = abs(FP2(K) - FP2(L))

            if(distancey < distance) then 
                if(distancex < distance) then 
                    NNP(L) = NNP(L) + 1
                    NLP(L, NNP(L)) = K 
                end if 
            end if 

        end do

    end do 



    return 

    end subroutine  LAG_POINTS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    !<<<<<<<<<<<<<<< BIBEMAT >>>>>>>>>>>>>>>>
    subroutine BIBEMAT(BIBE1, BIBE2, BIBET)
    use omp_lib
    use global 
    implicit none 

    double precision, dimension(0:N4, 0:N4) :: BIBE1, BIBE2, BIBET 

    integer :: i, j 
    integer :: L, M, N 
    double precision :: distancex, distancey 
    double precision, dimension(ND*2+1) :: DHX1, DHX2, DHY1, DHY2 


    ! INITIALIZATION
    BIBE1 = 0.0 
    BIBE2 = 0.0 
    BIBET = 0.0


    ! BIBE1
    DHX1 = 0.0 
    DHX2 = 0.0
    DHY1 = 0.0
    DHY2 = 0.0 
    !$omp parallel do                                   &
    !$omp & private(L, N, M, i, j)                      &
    !$omp & private(distancex, distancey)               &
    !$omp & private(DHX1, DHX2, DHY1, DHY2)             &
    !$omp & shared(BIBE1)
    do L = 0, NL-1
        do N = 1, NNP(L)
            M = NLP(L,N)

            ! FOR L 
            do j = NP2(L)-ND, NP2(L)+ND
                distancey = abs(0.5*(Y(j) + Y(j+1)) - FP2(L))
                call DELTA3P(distancey, DHY1(j-NP2(L)+ND+1))
            end do 

            do i = NP1(L)-ND, NP1(L)+ND 
                distancex = abs(X(i) - FP1(L))
                call DELTA3P(distancex, DHX1(i-NP1(L)+ND+1))
            end do 

            !FOR M 
            do j = NP2(L)-ND, NP2(L)+ND
                distancey = abs(0.5*(Y(j) + Y(j+1)) - FP2(M))
                call DELTA3P(distancey, DHY2(j-NP2(L)+ND+1))
            end do 

            do i = NP1(L)-ND, NP1(L)+ND 
                distancex = abs(X(i) - FP1(M))
                call DELTA3P(distancex, DHX2(i-NP1(L)+ND+1))
            end do 

            do j = 1, ND*2+1
                do i = 1, ND*2+1 

                    BIBE1(L,M) = BIBE1(L,M) + DHX1(i)*DHX2(i)*DHY1(j)*DHY2(j) 

                end do 
            end do 
        end do 
    end do 
    !$omp end parallel do 


    ! BIBE2
    DHX1 = 0.0 
    DHX2 = 0.0
    DHY1 = 0.0
    DHY2 = 0.0 
    !$omp parallel do                                   &
    !$omp & private(L, N, M, i, j)                      &
    !$omp & private(distancex, distancey)               &
    !$omp & private(DHX1, DHX2, DHY1, DHY2)             &
    !$omp & shared(BIBE2)
    do L = 0, NL-1
        do N = 1, NNP(L)
            M = NLP(L,N)

            ! FOR L 
            do j = NP2(L)-ND, NP2(L)+ND
                distancey = abs(Y(j) - FP2(L))
                call DELTA3P(distancey, DHY1(j-NP2(L)+ND+1))
            end do 

            do i = NP1(L)-ND, NP1(L)+ND 
                distancex = abs(0.5*(X(i) + X(i+1)) - FP1(L))
                call DELTA3P(distancex, DHX1(i-NP1(L)+ND+1))
            end do 

            !FOR M 
            do j = NP2(L)-ND, NP2(L)+ND
                distancey = abs(Y(j) - FP2(M))
                call DELTA3P(distancey, DHY2(j-NP2(L)+ND+1))
            end do 

            do i = NP1(L)-ND, NP1(L)+ND 
                distancex = abs(0.5*(X(i) + X(i+1)) - FP1(M))
                call DELTA3P(distancex, DHX2(i-NP1(L)+ND+1))
            end do 

            do j = 1, ND*2+1
                do i = 1, ND*2+1 

                    BIBE2(L,M) = BIBE2(L,M) + DHX1(i)*DHX2(i)*DHY1(j)*DHY2(j) 

                end do 
            end do 
        end do 
    end do 
    !$omp end parallel do 


    ! BIBET
    DHX1 = 0.0 
    DHX2 = 0.0
    DHY1 = 0.0
    DHY2 = 0.0 
    !$omp parallel do                                   &
    !$omp & private(L, N, M, i, j)                      &
    !$omp & private(distancex, distancey)               &
    !$omp & private(DHX1, DHX2, DHY1, DHY2)             &
    !$omp & shared(BIBET)
    do L = 0, NL-1
        do N = 1, NNP(L)
            M = NLP(L,N)

            ! FOR L 
            do j = NP2(L)-ND, NP2(L)+ND
                distancey = abs(0.5*(Y(j) + Y(j+1)) - FP2(L))
                call DELTA3P(distancey, DHY1(j-NP2(L)+ND+1))
            end do 

            do i = NP1(L)-ND, NP1(L)+ND 
                distancex = abs(0.5*(X(i) + X(i+1)) - FP1(L))
                call DELTA3P(distancex, DHX1(i-NP1(L)+ND+1))
            end do 

            !FOR M 
            do j = NP2(L)-ND, NP2(L)+ND
                distancey = abs(0.5*(Y(j) + Y(j+1)) - FP2(M))
                call DELTA3P(distancey, DHY2(j-NP2(L)+ND+1))
            end do 

            do i = NP1(L)-ND, NP1(L)+ND 
                distancex = abs(0.5*(X(i) + X(i+1)) - FP1(M))
                call DELTA3P(distancex, DHX2(i-NP1(L)+ND+1))
            end do 

            do j = 1, ND*2+1
                do i = 1, ND*2+1 

                    BIBET(L,M) = BIBET(L,M) + DHX1(i)*DHX2(i)*DHY1(j)*DHY2(j)

                end do 
            end do 
        end do 
    end do 
    !$omp end parallel do 


    

    return 

    end subroutine BIBEMAT  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    !<<<<<<<<<<<<<<< DELTA1 >>>>>>>>>>>>>>>>>
    subroutine DELTA1(j, i, L, dhx, dhy)
    use global 
    implicit none 

    integer :: i, j, L 
    double precision :: dhx, dhy 
    double precision :: distance, distancex, distancey 

    dhx = 0.0 
    dhy = 0.0 

    distancex = abs(X(i) - FP1(L))
    distancey = abs(0.5*(Y(j) + Y(j+1)) - FP2(L))

    call DELTA3P(distancex, dhx)
    call DELTA3P(distancey, dhy)

    
    return 

    end subroutine DELTA1 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    !<<<<<<<<<<<<<<< DELTA2 >>>>>>>>>>>>>>>>>
    subroutine DELTA2(j, i, L, dhx, dhy)
    use global 
    implicit none 

    integer :: i, j, L 
    double precision :: dhx, dhy 
    double precision :: distance, distancex, distancey 

    dhx = 0.0 
    dhy = 0.0 

    distancex = abs(0.5*(X(i) + X(i+1)) - FP1(L))
    distancey = abs(Y(j) - FP2(L))

    call DELTA3P(distancex, dhx)
    call DELTA3P(distancey, dhy)

    
    return 

    end subroutine DELTA2 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    !<<<<<<<<<<<<<<< DELTA >>>>>>>>>>>>>>>>>
    subroutine DELTA(j, i, L, dhx, dhy)
    use global 
    implicit none 

    integer :: i, j, L 
    double precision :: dhx, dhy 
    double precision :: distance, distancex, distancey 

    dhx = 0.0 
    dhy = 0.0 

    distancex = abs(0.5*(X(i) + X(i+1)) - FP1(L))
    distancey = abs(0.5*(Y(j) + Y(j+1)) - FP2(L))

    call DELTA3P(distancex, dhx)
    call DELTA3P(distancey, dhy)

    
    return 

    end subroutine DELTA
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    !<<<<<<<<<<<<<<< DELTA3P >>>>>>>>>>>>>>>>
    subroutine DELTA3P(distance, DH)
    use global 
    implicit none 

    double precision :: distance, DH 
    double precision :: temp 

    temp = distance/DS

    if(temp < 0.5) then 
        DH = 1.0/3.0*(1.0 + sqrt(1.0 - 3.0*temp**2.0))
    else if(temp < 1.5) then 
        DH = 1.0/6.0*(5.0 - 3.0*temp - sqrt(1.0 - 3.0*(1.0-temp)**2.0))
    else 
        DH = 0.0
    end if 


!     if(distance < 1.5*DS) then 
!         if(distance < 0.5*DS) then 
!             temp = distance 
!             DH = 1.0/(3.0*DS)*(1.0 + sqrt(1.0 - 3.0*(temp/DS)**2.0))
!         else
!             temp = distance 
!             DH = 1.0/(6.0*DS)*(5.0 - 3.0*temp/DS - sqrt(1.0 - 3.0*(1.0-temp/DS)**2.0))
!         end if 
!     end if

    return 
    end subroutine DELTA3P
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    !<<<<<<<<<<<<<<< DELTA4P >>>>>>>>>>>>>>>>
    subroutine DELTA4P(distance, DH)
    use global 
    implicit none 

    double precision :: distance, DH 
    double precision :: temp 
    
    temp = distance/DS 

    if(temp < 1.0) then 
        DH = 1.0/8.0*(3.0 - 2.0*temp + sqrt(1.0 + 4.0*temp - 4.0*temp**2.0))
    else if(temp < 2.0) then 
        DH = 1.0/8.0*(5.0 - 2.0*temp - sqrt(-7.0 + 12.0*temp - 4.0*temp**2.0))
    else
        DH = 0.0 
    end if 

    DH =DH/DS 


    return 
    end subroutine DELTA4P
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    !<<<<<<<<<<<<<<< SM_DELTA3P >>>>>>>>>>>>>>>>
    subroutine SM_DELTA3P(distance, DH)
    use global 
    implicit none 

    double precision :: distance, DH
    double precision :: temp 

    temp = distance/DS 

    if(temp < 1.0) then 
        DH = 17.0/48.0 + sqrt(3.0)*PI/108.0 + temp/4.0                          &
         & + (1.0 - 2.0*temp)/16.0*sqrt(-12.0*temp**2.0 + 12.0*temp + 1.0)      &
         & - sqrt(3.0)/12.0*asin(sqrt(3.0)/2.0*(2.0*temp - 1.0))
    else if(temp < 2.0) then 
        DH = 55.0/48.0 - sqrt(3.0)*PI/108.0 - 13.0*temp/12.0 + temp**2.0/4.0    &
         & + (2.0*temp - 3.0)/48.0*sqrt(-12.0*temp**2.0 + 36.0*temp - 23.0)     &
         & + sqrt(3.0)/36.0*asin(sqrt(3.0)/2.0*(2.0*temp - 3.0))  
    else
        DH = 0.0
    end if 

    DH = DH/DS 


    return 
    end subroutine SM_DELTA3P 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 


    !<<<<<<<<<<<<<<< SM_DELTA4P >>>>>>>>>>>>>>>>
    subroutine SM_DELTA4P(distance, DH)
    use global 
    implicit none 

    double precision :: distance, DH 
    double precision :: temp 
    
    temp = distance/DS 

    if(temp < 0.5) then 
        DH = 3.0/8.0 + pi/32.0 - temp**2.0/4.0 
    else if(temp < 1.5) then 
        DH = 1.0/4.0 + (1.0 - temp)/8.0*sqrt(-2.0 + 8.0*temp - 4.0*temp**2.0) - 1.0/8.0*asin(sqrt(2.0)*(temp - 1.0))
    else if(temp < 2.5) then 
        DH = 17.0/16.0 - PI/64.0 - 3.0*temp/4.0 + temp**2.0/8.0               & 
         & + (temp - 2.0)/16.0*sqrt(-14.0 + 16.0*temp - 4.0*temp**2.0)        &
         & + 1.0/16.0*asin(sqrt(2.0)*(temp - 2.0))
    else
        DH = 0.0 
    end if 

    DH =DH/DS 


    return 
    end subroutine SM_DELTA4P
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    !<<<<<<<<<<<<<<< B-SPLINE DELTA4P >>>>>>>>>>>>>>>>
    subroutine B_DELTA4P(distance, DH)
    use global 
    implicit none 

    double precision :: distance, DH 
    double precision :: temp 
    
    temp = distance/DS

    if(temp < 1.0) then 
        DH = 0.5*temp**3.0 - temp**2.0 + 2.0/3.0 
    else if(temp < 2.0) then 
        DH = -temp**3.0/6.0 + temp**2.0 - 2.0*temp + 4.0/3.0
    else
        DH = 0.0 
    end if 

    DH = DH/DS 

    return 

    end subroutine B_DELTA4P 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    !<<<<<<<<<<<<<<< B_DELTA5P >>>>>>>>>>>>>>>>>
    subroutine B_DELTA5P(distance, DH)
    use global 
    implicit none 

    double precision :: distance, DH
    double precision :: temp 

    temp = distance/DS 

    if(temp < 0.5) then 
        DH = 0.25*temp**4.0 - 5.0/8.0*temp**2.0 + 115.0/192.0 
    else if(temp < 1.5) then 
        DH = -1.0/6.0*temp**4.0 + 5.0/6.0*temp**3.0 - 5.0/4.0*temp**2.0 + 5.0/24.0*temp + 55.0/96.0 
    else if(temp < 2.5) then 
        DH = 1.0/24.0*temp**4.0 - 5.0/12.0*temp**3.0 + 25.0/16.0*temp**2.0 - 125.0/48.0*temp + 625.0/384.0 
    else
        DH = 0.0 
    end if 

    DH = DH/DS


    return 
    end subroutine B_DELTA5P     
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    !<<<<<<<<<<<<<<< B_DELTA6P >>>>>>>>>>>>>>>>>
    subroutine B_DELTA6P(distance, DH)
    use global 
    implicit none 

    double precision :: distance, DH 
    double precision :: temp 


    temp = distance/DS 

    if(temp < 1.0) then 
        DH = -1.0/12.0*temp**5.0 + 0.25*temp**4.0 - 0.5*temp**2.0 + 11.0/20.0 
    else if(temp < 2.0) then 
        DH = 1.0/24.0*temp**5.0 - 3.0/8.0*temp**4.0 + 5.0/4.0*temp**3.0 - 7.0/4.0*temp**2.0 + 5.0/8.0*temp + 17.0/40.0 
    else if(temp <3.0) then 
        DH = -1.0/120.0*temp**5.0 + 1.0/8.0*temp**4.0 - 3.0/4.0*temp**3.0 + 9.0/4.0*temp**2.0 - 27.0/8.0*temp + 81.0/40.0 
    else
        DH = 0.0 
    end if 

    DH = DH/DS 


    return 
    end subroutine B_DELTA6P 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    !<<<<<<<<<<<<<<< GETUPT >>>>>>>>>>>>>>>>>
    subroutine GETUPT(THETA, U, V, P, BIBE1, BIBE2, BIBET, t)
    use omp_lib 
    use global 
    implicit none 

    double precision :: t 
    double precision, dimension(0:N1, 0:N2) :: U, V, P, THETA 
    double precision, dimension(0:N4, 0:N4) :: BIBE1, BIBE2, BIBET 

    double precision :: UAV, EAV 
    double precision :: UHAV, EHAV
    double precision, dimension(0:N4) :: FL1, FL2, FLE   
    double precision, dimension(0:N1, 0:N2) :: dTHETA 
    double precision, dimension(0:N1, 0:N2) :: dUm, dVm, dP 
    double precision, dimension(0:N1, 0:N2) :: F1, F2, FE  
    double precision, dimension(0:N1, 0:N2) :: RHSe, RHSx, RHSy 

    integer :: i, j, k, L 
    double precision :: t1, t2 


    ! INITIALIZATION
    F1 = 0.0 
    F2 = 0.0
    FE = 0.0 

    FL1 = 0.0
    FL2 = 0.0
    FLE = 0.0  

    dP = 0.0 

    dUm = 0.0 
    dVm = 0.0
    dTHETA = 0.0 


    call BCOND(t)

    !===================================================
    !      IMMERSED BOUNDARY PROJECTION METHOD
    !===================================================
    EAV = 1.0 
    UAV = 1.0 
    EHAV = 1.0 
    UHAV = 1.0 


    ! CALCULATE THE INTERMEDIATE TEMPERATURE
    call CPU_TIME(t1)
    call Solve_dEnergy(dTHETA)
    call CPU_TIME(t2)
    time(2) = time(2) + t2 - t1
    
    ! CALCULATE THE INTERMEDIATE VELOCITY
    call CPU_TIME(t1)
    call SolveMomentum(dTHETA, t, dUm, dVm)
    call CPU_TIME(t2)   
    time(3) = time(3) + t2 - t1
    

    ! GET dUm and dTHETA AT LAGRANGIAN POINTS 
    call CPU_TIME(t1) 
    call UINTERPOL(dUm, dVm, dTHETA, UHAV, EHAV)
    call CPU_TIME(t2)
    time(4) = time(4) + t2 - t1 

    ! GET FORCING AT LAGRANGIAN POINTS
    call CPU_TIME(t1)
    call FLAG(BIBE1, BIBE2, BIBET, FL1, FL2, FLE)
    call CPU_TIME(t2)
    time(5) = time(5) + t2 - t1 

    ! GET FORCING AT EULERIAN POINTS 
    call CPU_TIME(t1)
    call FEXTRAPOL(F1, F2, FE, FL1, FL2, FLE)
    call CPU_TIME(t2)
    time(1) = time(1) + t2- t1 

    ! UPDATE INTERMEDIATE VELOCITY

    ! =================================================
    !           ALGORITHM 1 
    call CPU_TIME(t1)
    call UPDATEdUVT(dUm, dVm, dTHETA, F1, F2, FE)
    call CPU_TIME(t2)
    time(6) = time(6) + t2 - t1 
    ! #################################################
 

    call UINTERPOL(dUm, dVm, dTHETA, UHAV, EHAV)

    !==================================================
    !         POISSON EQUATION
    !==================================================
    
    !CALCULATE DP 
    call CPU_TIME(t1)
    call DPCALC(dUm, dVm, dP, t)
    call CPU_TIME(t2)
    time(7) = time(7) + t2 - t1


    ! UPDATE THE N+1 TIME STEP VELOCITY AND PRESSURE
    call CPU_TIME(t1)
    call UPDATEUVP(U, V, P, THETA, dUm, dVm, dP, dTHETA)
    call CPU_TIME(t2)
    time(8) = time(8) + t2 - t1 

    ! CALCULATE U^(N+1) AT LAGRANGIAN POINTS
    call CPU_TIME(t1) 
    call UINTERPOL(U, V, THETA, UAV, EAV)
    call CPU_TIME(t2)
    time(4) = time(4) + t2 - t1 


    write(*,*) 'UAV = ', UAV/max(umax, vmax) 
    write(*,*) 'EAV = ', EAV

    write(28,*) UAV/max(umax, vmax) , UHAV/max(umax, vmax) 
    write(29,*) EAV, EHAV 


    call COEFCALC(FL1, FL2, FLE, t)


    return 
    end subroutine GETUPT
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    !<<<<<<<<<<<<<<< BCOND >>>>>>>>>>>>>>>
    subroutine BCOND(t)
    use omp_lib
    use Global
    implicit none
    
    integer :: i, j
    double precision :: t
    double precision :: UC, RATE  
    double precision :: Q_IN, Q_EX, Q_UP, Q_DOWN   


    ! IN
    !$omp parallel do 
    do j = 0, N2 
        UBC1(j,1) = 0.0 
        UBC1(j,2) = 0.0 

        THETABC1(j) = 0.0 

        UBC2(j,1) = 0.0 
        UBC2(j,2) = 0.0 

        THETABC2(j) = 0.0 
    end do 
    !$omp end parallel do 


    ! OUT
!     !$omp parallel do 
!     do j = 0, N2 
!         UBC2(j,1) = Up(N1m,j)
!         UBC2(j,2) = Vp(N1m,j)

!         THETABC2(j) = THETAp(N1m,j)
!     end do 
!     !$omp end parallel do 



    ! UP & DOWN
    !$omp parallel do 
    do i = 0, N1

        UBC3(i,1) = 0.0
        UBC3(i,2) = 0.0   !Vp(i,2)
        THETABC3(i) = 0.0
        
        
        UBC4(i,1) = 0.0
        UBC4(i,2) = 0.0   !Vp(i,N2m)
        THETABC4(i) = 0.0

    end do
    !$omp end parallel do     


    ! CONVECTIVE BOUNDARY CONDITION
!     UC = 0.0 
!     do j = 1, N2m 
!         UC = UC + Up(N1,j)*DY(j)
!     end do 
!     UC = UC/ALY 

!     do j = 1, N2m 
!         UBC2(j,1) = Up(N1,j) - dt/DX(N1m)*UC*(Up(N1,j) - Up(N1m,j))
!         UBC2(j,2) = Vp(N1,j) - dt/DMX(N1)*UC*(Vp(N1,j) - Vp(N1m,j)) 
!         THETABC2(j) = THETAp(N1,j) - dt/DMX(N1)*UC*(THETAp(N1,j) - THETAp(N1m,j))
!     end do 

!     Q_IN = 0.0
!     Q_EX = 0.0
!     Q_UP = 0.0
!     Q_DOWN = 0.0 

!     do j = 1, N2m 
!         Q_IN = Q_IN + UBC1(j,1)*DY(j)
!         Q_EX = Q_EX + UBC2(j,1)*DY(j)
!     end do 

!     do i = 1, N1m
!         Q_UP   = Q_UP   + UBC4(i,2)*DX(i)
!         Q_DOWN = Q_DOWN + UBC3(i,2)*DX(i)
!     end do 

!     write(*,*) 
!     write(*,100) Q_IN, Q_UP, Q_DOWN, Q_EX

! 100 FORMAT('MASS FLUX: Q_IN = ', E11.4,X, 'Q_UP = ', E11.4,X, 'Q_DOWN = ', E11.4,X, 'Q_EX = ', E11.4)  !PARK

!     RATE = (Q_IN + Q_DOWN - Q_UP)/Q_EX

!     do j =1, N2m 
!         UBC2(j,1) = RATE*UBC2(j,1)
!         UBC2(j,2) = RATE*UBC2(j,2)
!         THETABC2(j) = RATE*THETABC2(j) 
!     end do 


    return    
    end subroutine BCOND
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

   
    !<<<<<<<<<<<<<<< FEXTRAPOL >>>>>>>>>>>>>>
    subroutine FEXTRAPOL(F1, F2, FE, FL1, FL2, FLE)
    use global
    implicit none 
    
    double precision, dimension(0:N4) :: FL1, FL2, FLE 
    double precision, dimension(0:N1, 0:N2) :: F1, F2, FE 

    integer :: i, j, L 
    double precision :: distancex, distancey 
    double precision, dimension(ND*2+1) :: dhx, dhy 
    

    !INITIALIZATION
    F1 = 0.0 
    F2 = 0.0 
    FE = 0.0 

    !F1 
    dhx = 0.0
    dhy = 0.0
    do L = 0, NL-1 

        do j = NP2(L)-ND, NP2(L)+ND
            distancey = abs(0.5*(Y(j) + Y(j+1)) - FP2(L))
            call DELTA3P(distancey, dhy(j-NP2(L)+ND+1))
        end do 

        do i = NP1(L)-ND, NP1(L)+ND 
            distancex = abs(X(i) - FP1(L))
            call DELTA3P(distancex, dhx(i-NP1(L)+ND+1))
        end do 

        do j = NP2(L)-ND, NP2(L)+ND
            do i = NP1(L)-ND, NP1(L)+ND 
                F1(i,j) = F1(i,j) + dhx(i-NP1(L)+ND+1)*dhy(j-NP2(L)+ND+1)*FL1(L)
            end do 
        end do 

    end do 


    !F2
    dhx = 0.0
    dhy = 0.0
    do L = 0, NL-1 

        do j = NP2(L)-ND, NP2(L)+ND
            distancey = abs(Y(j) - FP2(L))
            call DELTA3P(distancey, dhy(j-NP2(L)+ND+1))
        end do 

        do i = NP1(L)-ND, NP1(L)+ND 
            distancex = abs(0.5*(X(i) + X(i+1)) - FP1(L))
            call DELTA3P(distancex, dhx(i-NP1(L)+ND+1))
        end do 

        do j = NP2(L)-ND, NP2(L)+ND
            do i = NP1(L)-ND, NP1(L)+ND 
                F2(i,j) = F2(i,j) + dhx(i-NP1(L)+ND+1)*dhy(j-NP2(L)+ND+1)*FL2(L)
            end do 
        end do 

    end do 


    !FE
    dhx = 0.0
    dhy = 0.0
    do L = 0, NL-1 

        do j = NP2(L)-ND, NP2(L)+ND
            distancey = abs(0.5*(Y(j) + Y(j+1)) - FP2(L))
            call DELTA3P(distancey, dhy(j-NP2(L)+ND+1))
        end do 

        do i = NP1(L)-ND, NP1(L)+ND 
            distancex = abs(0.5*(X(i) + X(i+1)) - FP1(L))
            call DELTA3P(distancex, dhx(i-NP1(L)+ND+1))
        end do 

        do j = NP2(L)-ND, NP2(L)+ND
            do i = NP1(L)-ND, NP1(L)+ND 
                FE(i,j) = FE(i,j) + dhx(i-NP1(L)+ND+1)*dhy(j-NP2(L)+ND+1)*FLE(L)
            end do 
        end do 

    end do 


    return 

    end subroutine FEXTRAPOL
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    
    !<<<<<<<<<<<<<<< UINTERPOL >>>>>>>>>>>>>>>>>
    SUBROUTINE UINTERPOL(dUm, dVm, dTHETA, UAV, EAV)
    use global 
    implicit none 

    integer :: i, j, L 
    double precision :: EAV, UAV 
    double precision :: distancex, distancey  
    double precision, dimension(ND*2+1) :: dhx, dhy 
    double precision, dimension(0:N1, 0:N2) :: dUm, dVm, dTHETA  
    

    !INITIALIZATION
    ULAG1 = 0.0
    ULAG2 = 0.0
    ELAG  = 0.0

    EAV  = 0.0
    UAV = 0.0 


    !ULAG1
    dhx = 0.0 
    dhy = 0.0 
    !$omp parallel do                           &
    !$omp & private(L, i, j)                    &
    !$omp & private(distancex, distancey)       &
    !$omp & private(dhx, dhy)                   &
    !$omp & shared(ULAG1)                       
    do L = 0, NL-1

        do j = NP2(L)-ND, NP2(L)+ND
            distancey = abs(0.5*(Y(j) + Y(j+1)) - FP2(L))
            call DELTA3P(distancey, dhy(j-NP2(L)+ND+1))
        end do 

        do i = NP1(L)-ND, NP1(L)+ND 
            distancex = abs(X(i) - FP1(L))
            call DELTA3P(distancex, dhx(i-NP1(L)+ND+1))
        end do 

        do j = NP2(L)-ND, NP2(L)+ND
            do i = NP1(L)-ND, NP1(L)+ND 
                ULAG1(L) = ULAG1(L) + dhx(i-NP1(L)+ND+1)*dhy(j-NP2(L)+ND+1)*dUm(i,j)
            end do 
        end do 

    end do 
    !$omp end parallel do 


    !ULAG2
    dhx = 0.0 
    dhy = 0.0 
    !$omp parallel do                           &
    !$omp & private(L, i, j)                    &
    !$omp & private(distancex, distancey)       &
    !$omp & private(dhx, dhy)                   &
    !$omp & shared(ULAG2) 
    do L = 0, NL-1

        do j = NP2(L)-ND, NP2(L)+ND
            distancey = abs(Y(j) - FP2(L))
            call DELTA3P(distancey, dhy(j-NP2(L)+ND+1))
        end do 

        do i = NP1(L)-ND, NP1(L)+ND 
            distancex = abs(0.5*(X(i) + X(i+1)) - FP1(L))
            call DELTA3P(distancex, dhx(i-NP1(L)+ND+1))
        end do 

        do j = NP2(L)-ND, NP2(L)+ND
            do i = NP1(L)-ND, NP1(L)+ND 
                ULAG2(L) = ULAG2(L) + dhx(i-NP1(L)+ND+1)*dhy(j-NP2(L)+ND+1)*dVm(i,j)
            end do 
        end do 

    end do
    !$omp end parallel do  



    !ELAG
    dhx = 0.0 
    dhy = 0.0 
    !$omp parallel do                           &
    !$omp & private(L, i, j)                    &
    !$omp & private(distancex, distancey)       &
    !$omp & private(dhx, dhy)                   &
    !$omp & shared(ELAG) 
    do L = 0, NL-1

        do j = NP2(L)-ND, NP2(L)+ND
            distancey = abs(0.5*(Y(j) + Y(j+1)) - FP2(L))
            call DELTA3P(distancey, dhy(j-NP2(L)+ND+1))
        end do 

        do i = NP1(L)-ND, NP1(L)+ND 
            distancex = abs(0.5*(X(i) + X(i+1)) - FP1(L))
            call DELTA3P(distancex, dhx(i-NP1(L)+ND+1))
        end do 

        do j = NP2(L)-ND, NP2(L)+ND
            do i = NP1(L)-ND, NP1(L)+ND 
                ELAG(L) = ELAG(L) + dhx(i-NP1(L)+ND+1)*dhy(j-NP2(L)+ND+1)*dTHETA(i,j)
            end do 
        end do 

    end do 
    !$omp end parallel do 



    UAV = sqrt(sum((ULAG1(0:NL-1) - UD1(0:NL-1))**2. + (ULAG2(0:NL-1) - UD2(0:NL-1))**2.) /dble(NL))
    EAV = sqrt(sum((ELAG(0:NL-1)  -  QD(0:NL-1))**2.) /dble(NL))



    return 

    end subroutine UINTERPOL
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
    !<<<<<<<<<<<<<<< FLAG >>>>>>>>>>>>>>>>>>>
    subroutine FLAG(BIBE1, BIBE2, BIBET, FL1, FL2, FLE)
    use global 
    implicit none 

    double precision, dimension(0:N4) :: FL1, FL2, FLE
    double precision, dimension(0:N4, 0:N4) :: BIBE1, BIBE2, BIBET 

    integer :: L 
    double precision, dimension(0:N4) :: H1, H2, HE 
    double precision, dimension(0:N4) :: d1, d2, de 


    ! CALCULATE RHS OF BIBE*DF=(UD-UL)/DT & BIBE*DF=(QD-QL)/DT
    !$omp parallel do 
    do L = 0, NL-1 
        H1(L) = (UD1(L) - ULAG1(L))/dt 
        H2(L) = (UD2(L) - ULAG2(L))/dt 
        HE(L) = (QD(L ) - ELAG(L ))/dt 
    end do 
    !$omp end parallel do 
    

    ! CALCULATE DELTA FORCING 
!     call CG(BIBE1, H1, dFL1)
!     call CG(BIBE2, H2, dFL2)
!     call CG(BIBET, HE, dFLE)

    !THE EXPLICIT METHOD
    do L = 0, NL-1 
        d1(L) = sum(BIBE1(0:NL-1,L))
        d2(L) = sum(BIBE2(0:NL-1,L))
        de(L) = sum(BIBET(0:NL-1,L))
    end do 

    do L = 0, NL-1
        FL1(L) = H1(L)/d1(L)
        FL2(L) = H2(L)/d2(L)
        FLE(L) = HE(L)/de(L)
    end do 


    FSUM = sum(FLE(0:NL-1))

    write(34,*) sum(FL1(0:NL-1)), sum(FL2(0:NL-1)), sum(FLE(0:NL-1))


    return

    end subroutine FLAG
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    !<<<<<<<<<<<<<<< CG >>>>>>>>>>>>>>>>>>>>>
    subroutine CG(BIBE, H, SOL)
    use global 
    implicit none 

    double precision :: BIBE(0:N4, 0:N4)
    double precision :: H(0:N4), SOL(0:N4) 

    integer :: J, K, L, M 
    integer :: ITERMAX  
    double precision :: TOL, RESI, DUMY  
    double precision :: ALPHA, BETA, ADR, RDR
    double precision, dimension(0:N4) :: RNEW, ROLD, P 

    TOL = 1.0E-12
    RESI = 1.0 
    ITERMAX = 1000

    ! INITIALIZATION
    SOL = 0.0 
    RNEW = H
    ROLD = H 
    P = H 

    ! START ITERATION 
    K = 0 

    
    do while (RESI > TOL .and. K < ITERMAX)
        
        ! GET ALPHA 
        ADR = 0.0 
        RDR = 0.0

        do L = 0, NL-1 
            
            ADR = ADR + ROLD(L)**2.0 
            
            DUMY = 0.0
            do J = 1, NNP(L)
                M = NLP(L,J)
                DUMY = DUMY + BIBE(L,M)*P(M)
            end do 
            RDR = RDR + P(L)*DUMY 
        end do 

        if(RDR == 0.0) then
            exit 
        end if 

        ALPHA = ADR/RDR 

        ! UPDATE SOLUTION & RESIDUAL 
        RESI = 0.0 
        do L = 0, NL-1 
                
            SOL(L) = SOL(L) + ALPHA*P(L)
            RNEW(L) = ROLD(L)
                
            DUMY = 0.0 
            do J = 1, NNP(L)
                M = NLP(L,J)
                DUMY = DUMY + BIBE(L,M)*P(M)
            end do 

            RNEW(L) = ROLD(L) - ALPHA*DUMY 
            RESI = RESI + abs(RNEW(L))

        end do 

        RESI = RESI/dble(NL)

        ! GET BETA 
        ADR = 0.0
        RDR = 0.0 

        do L = 0, NL-1
            ADR = ADR + RNEW(L)**2.0
            RDR = RDR + ROLD(L)**2.0 
        end do 

        BETA = ADR/RDR 

        do L = 0, NL-1 
            P(L) = RNEW(L) + BETA*P(L)
            ROLD(L) = RNEW(L)
        end do 

        K = K + 1 

    end do 


    write(*,*) 'RESIDUAL = ', RESI, 'CGITER = ', K 

    
    return 

    end subroutine CG 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    !<<<<<<<<<<<<<<< UPDATEdUVWm >>>>>>>>>>>>
    subroutine UPDATEdUVT(dUm, dVm, dTHETA, F1, F2, FE)
    use global 
    implicit none 
    
    double precision, dimension(0:N1, 0:N2) :: dUm, dVm, dTHETA
    double precision, dimension(0:N1, 0:N2) :: F1, F2, FE  

    integer :: i, j
    double precision :: fc, Thetay 
    
    ! TEMPERATURE
    !$omp parallel do 
    do j = 1, N2m 
        do i = 1, N1m 
            dTHETA(i,j) = dTHETA(i,j) + dt*FE(i,j)
        end do 
    end do 
    !$omp end parallel do 

    
    ! VELOCITY IN X-DIRECTION
    !$omp parallel do
    do j = 1, N2m
        do i = 2, N1m 
            dUm(i,j) = dUm(i,j) + dt*F1(i,j)
        end do 
    end do 
    !$omp end parallel do 

    ! VELOCITY IN Y-DIRECTION
    !$omp parallel do 
    do j = 2, N2m 
        do i = 1, N1m 
            dVm(i,j) = dVm(i,j) + dt*F2(i,j)
        end do 
    end do
    !$omp end parallel do 

    
    return 

    end subroutine UPDATEdUVT 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    !<<<<<<<<<<<<<< DPCALC >>>>>>>>>>>>>>>>> 
    subroutine DPCALC(dUm, dVm, dP, t)
    use global 
    implicit none 

    double precision :: t 
    double precision, dimension(0:N1, 0:N2) :: dUm, dVm
    double precision, dimension(0:N1, 0:N2) :: dP, RDP 

    call RHSDP(RDP, dUm, dVm)

    call PoissonSolver_FFTW1dy(RDP, dP)

!     call PoissonSolverMG(RDP(1:N1m,1:N2m), dP, t)


    end subroutine DPCALC
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    !<<<<<<<<<<<<<<< RHSDP >>>>>>>>>>>>>>>>>>
    subroutine RHSDP(RDP, dUm, dVm)
    use omp_lib
    use global 
    implicit none

    double precision, dimension(0:N1, 0:N2) :: RDP 
    double precision, dimension(0:N1, 0:N2) :: dUm, dVm
    
    integer :: i, j
    integer :: ip, ium, iup, jp, jvm, jvp

    double precision :: cbc, DivU 


    !$omp parallel do                                      
    do j = 1, N2m
    jp = j + 1
    jvp = JP_A(j) - j
    jvm = j - JM_U(j)    
    do i = 1, N1m
    ip = i + 1
    iup = IP_A(i) - i
    ium = i - IM_V(i)
  

        DivU = ( dble(iup)*dUm(ip,j) - dble(ium)*dUm(i,j))/DX(i)    &
           & + ( dble(jvp)*dVm(i,jp) - dble(jvm)*dVm(i,j))/DY(j) 

        cbc = dble(1. - ium)*UBC1(j,1)/DX(i)    &
          & - dble(1. - iup)*UBC2(j,1)/DX(i)    &
          & + dble(1. - jvm)*UBC3(i,2)/DY(j)    &
          & - dble(1. - jvp)*UBC4(i,2)/DY(j)


        RDP(i,j) = (DivU - cbc)/dt

    end do
    end do
    !$omp end parallel do


    return    

    end subroutine RHSDP
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    !<<<<<<<<<<<<<<< GetdUVW >>>>>>>>>>>>>>>
    subroutine UPDATEUVP(U, V, P, THETA, dUm, dVm, dP, dTHETA)
    use omp_lib
    use Global
    implicit none
    
    integer :: i, j
    integer :: ip, im, jp, jm

    double precision :: t

    double precision, dimension(0:N1, 0:N2) :: dUm, dVm, dP, dTHETA 
    double precision, dimension(0:N1, 0:N2) :: U, V, P, THETA  
    double precision :: tempx(2,2), tempy(2,2), temp(2,2)


    ! call omp_set_nested(.true.)

    ! TEMPERATURE
    !$omp parallel do 
    do j = 1, N2m 
        do i = 1, N1m 
            THETA(i,j) = dTHETA(i,j)
        end do 
    end do 
    !$omp end parallel do 

    ! BOUNDARY OF THETA
    ! y = ymin & ymax
    !$omp parallel do      
    do i = 0, N1
        THETA(i,0 ) = THETABC3(i)
        THETA(i,N2) = THETABC4(i)
    end do
    !$omp end parallel do

    ! x = xmin & xmax
    !$omp parallel do    
    do j = 0, N2
        THETA(0, j) = THETABC1(j)
        THETA(N1,j) = THETABC2(j)
    end do
    !$omp end parallel do

    ! UPDATE P
    !$omp parallel do    
    do j = 1, N2m
    do i = 1, N1m
        P(i,j) = Pp(i,j) + dP(i,j)
    end do
    end do
    !$omp end parallel do   


    ! BOUNDARY OF P
    tempx(1,1) = (DMX(1 )+DMX(2  ))**2./((DMX(1 )+DMX(2  ))**2.-DMX(1 )**2.)
    tempx(2,1) = DMX(1)**2.            /((DMX(1 )+DMX(2  ))**2.-DMX(1 )**2.)
    tempx(1,2) = (DMX(N1)+DMX(N1m))**2./((DMX(N1)+DMX(N1m))**2.-DMX(N1)**2.)
    tempx(2,2) =  DMX(N1)**2.          /((DMX(N1)+DMX(N1m))**2.-DMX(N1)**2.)  

    tempy(1,1) = (DMY(1 )+DMY(2  ))**2./((DMY(1 )+DMY(2  ))**2.-DMY(1 )**2.)
    tempy(2,1) = DMY(1 )**2.           /((DMY(1 )+DMY(2  ))**2.-DMY(1 )**2.)
    tempy(1,2) = (DMY(N2)+DMY(N2m))**2./((DMY(N2)+DMY(N2m))**2.-DMY(N2)**2.)
    tempy(2,2) = DMY(N2)**2.           /((DMY(N2)+DMY(N2m))**2.-DMY(N2)**2.)


    ! x = xmin & xmax
    P(0, 1:N2m) = tempx(1,1)*P(1,   1:N2m) - tempx(2,1)*P(2,   1:N2m)     
    P(N1,1:N2m) = tempx(1,2)*P(N1-1,1:N2m) - tempx(2,2)*P(N1-2,1:N2m)

    ! y = ymin & ymax
    P(1:N1m,0 ) = tempy(1,1)*P(1:N1m,1   ) - tempy(2,1)*P(1:N1m,2   )
    P(1:N1m,N2) = tempy(1,2)*P(1:N1m,N2-1) - tempy(2,2)*P(1:N1m,N2-2)
    

    ! POINTS
    P(0, 0) = 0.5*tempy(1,1)*P(0,   1) - 0.5*tempy(2,1)*P(0,   2)       &
          & + 0.5*tempx(1,1)*P(1,   0) - 0.5*tempx(2,1)*P(2,   0)
    P(N1,0) = 0.5*tempy(1,1)*P(N1,  1) - 0.5*tempy(2,1)*P(N1,  2)       &
          & + 0.5*tempx(1,2)*P(N1-1,0) - 0.5*tempx(2,2)*P(N1-2,0)
    
    P(0, N2) = 0.5*tempy(1,2)*P(0,   N2-1) - 0.5*tempy(2,2)*P(0,   N2-2)       &
           & + 0.5*tempx(1,1)*P(1,   N2  ) - 0.5*tempx(2,1)*P(2,   N2  )       
    P(N1,N2) = 0.5*tempy(1,2)*P(N1,  N2-1) - 0.5*tempy(2,2)*P(N1,  N2-2)       &
           & + 0.5*tempx(1,2)*P(N1-1,N2  ) - 0.5*tempx(2,2)*P(N1-2,N2  )


    
    ! UPDATE U
    !$omp parallel do private(j, i, im) shared(U, dUm)
    do j = 1, N2m
    do i = 2, N1m
        im = i - 1 
        U(i,j) = dUm(i,j) - dt*(dP(i,j) - dP(im,j))/DMX(i)
    end do
    end do
    !$omp end parallel do


    ! BOUNDARY OF U 
    ! y = ymin & ymax
    !$omp parallel do      
    do i = 1, N1
        U(i,0 ) = UBC3(i,1)
        U(i,N2) = UBC4(i,1)
    end do
    !$omp end parallel do

    ! x = xmin & xmax
    !$omp parallel do    
    do j = 0, N2
        U(1, j) = UBC1(j,1)
        U(N1,j) = UBC2(j,1)
    end do
    !$omp end parallel do


    ! UPDATE V
    !$omp parallel do private(j, jm, i) shared(V, dVm)    
    do j = 2, N2m
    jm = j - 1
    do i = 1, N1m
        V(i,j) = dVm(i,j) - dt*(dP(i,j) - dP(i,jm))/DMY(j)
    end do
    end do
    !$omp end parallel do
    

    ! BOUNDARY OF V
    ! y = ymin & ymax
    !$omp parallel do      
    do i = 0, N1        
        V(i,1 ) = UBC3(i,2)
        V(i,N2) = UBC4(i,2)
    end do
    !$omp end parallel do

    
    ! x = xmin & xmax
    !$omp parallel do      
    do j = 1, N2
        V(0, j) = UBC1(j,2)
        V(N1,j) = UBC2(j,2)
    end do
    !$omp end parallel do    


    
    return    
    end subroutine UPDATEUVP 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


    !<<<<<<<<<<<<<< COEFCALC >>>>>>>>>>>>>>
    subroutine COEFCALC(FL1, FL2, FLE, t)  
    use global 
    implicit none 

    double precision :: t 
    double precision, dimension(0:N4) :: FL1, FL2, FLE 
    
    integer :: L 
    double precision :: CD, CL, Nu  

    CD = 0.0 
    CL = 0.0
    Nu = 0.0  

    do L = 0, NL-1 
        CD = CD + FL1(L)*DV 
        CL = CL + FL2(L)*DV 
        Nu = Nu + FLE(L)*DV 
    end do 

    CD = -2.0*CD 
    CL = -2.0*CL 

    write(17, 107) t, CD, CL, Nu  

107 format(4(3X, F22.10))

    return 

    end subroutine COEFCALC
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 



    !<<<<<<<<<<<<<<< Cfl_Con >>>>>>>>>>>>>>>
    subroutine Cfl_Con
    use global 
    implicit none  

    integer :: i, j 
    double precision :: utmp, vtmp 


    umax = 0.0 
    vmax = 0.0 

    do j = 1, N2 
        do i = 1, N1 
            utmp = Up(i,j)
            vtmp = Vp(i,j)

            if(utmp > umax) umax = utmp 
            if(vtmp > vmax) vmax = vtmp 

        end do 
    end do 

    if(umax == 0.0 .and. vmax == 0.0) then 
        dt = dt 
    else 
        dt = cfl*DS/(umax + vmax)
    end if 

    end subroutine Cfl_Con 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


    !<<<<<<<<<<<<<<< Cal_dt >>>>>>>>>>>>>>>
    subroutine Cfl_Check
    use global 
    implicit none 

    integer :: i, j, ip, jp 
    double precision :: cfll 

    cflm = 0.0 
    do j = 1, N2m 
        jp = j +1 
        do i = 1, N1m 
            ip = i + 1
            cfll = abs(Up(i,j) + Up(ip,j))*(0.5/DX(i))    &
               & + abs(Vp(i,j) + Vp(i,jp))*(0.5/DY(j)) 
            
            cflm = dmax1(cflm, cfll)
        end do 
    end do 

    end subroutine Cfl_Check
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!     !<<<<<<<<<<<<<<< Cal_dt >>>>>>>>>>>>>>>
!     subroutine Cal_dt(t, CFL)
!     use Global
!     implicit none
    
!     integer :: i, j
!     integer :: ip, jp, ic, jc

!     double precision :: t, tmp, CFL_Idt, CFL
!     double precision :: xc, yc

!     xc = 0.
!     yc = 0.
!     CFL_Idt = 0.


!     do j = 1, N2m
!     jp = j + 1
!     do i = 1, N1m
!     ip = i + 1
!         tmp = ABS(Up(i,j) + Up(ip,j))*(0.5/DX(i))    &
!           & + ABS(Vp(i,j) + Vp(i,jp))*(0.5/DY(j)) 

!         if(CFL_Idt .le. tmp)then
!             CFL_Idt = tmp
!             xc = X(i); ic = i
!             yc = Y(j); jc = j
!         end if

!     end do
!     end do


!     ! Here, Pi is needed to be considered, since we are using the 
!     ! non-uniform Chebyshev nodes
!     ! dt = min(MaxCfl/(CFL_Idt*PI), 0.75*dt+0.25*MaxCfl/(CFL_Idt*PI))

!     dt = min(MaxCfl/(CFL_Idt), dt)


! !     write(*,*) dt

!     ! dt = dtstart
!     CFL = CFL_Idt*dt


! !     write(28, '(1e11.5,3e20.12)') t, CFL_Idt, CFL, dt
! !     write(29, 2001) dt, xc, ic, yc, jc, zc, kc

! 2001 format('DTmax=',E16.8, 3X, 'x_c = ', E11.5, i5, 3X, 'y_c = ', E11.5, i5)

!     return
!     end subroutine Cal_dt
! ! !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 



    !<<<<<<<<<<<<<<< CheckDiv >>>>>>>>>>>>>>>
    subroutine CheckDiv(DivU)
    use omp_lib
    use Global
    implicit none

    integer :: i, j
    integer :: ip, jp
    double precision, dimension(1:N1m, 1:N2m) :: DivU

!     call omp_set_nested(.true.)

    !$omp parallel do private(j, jp, i, ip) shared(DivU)
    do j = 1, N2m
    jp = j + 1
    do i = 1, N1m
    ip = i + 1

        DivU(i,j) = (Up(ip,j)-Up(i,j))/DX(i)      &
                & + (Vp(i,jp)-Vp(i,j))/DY(j)       

    end do
    end do
    !$omp end parallel do


    return
    end subroutine CheckDiv
! !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    !<<<<<<<<<<<<<<< PrintDiv >>>>>>>>>>>>>>>
    subroutine PrintDiv(DivU)
    use Global
    implicit none

    integer :: i, j
    integer :: ip, jp
    double precision :: xc, yc
    double precision, dimension(1:N1m, 1:N2m) :: DivU


    write(16,*) 'zone t="',1,'"','j=',N2m,'i=',N1m


    do j = 1, N2m
    jp = j + 1
    yc = 0.5*(Y(j) + Y(jp))
    do i = 1, N1m
    ip = i + 1
    xc = 0.5*(X(i) + X(ip))

        write(16,'(2(e11.5,2x),e20.12)') xc, yc, DivU(i,j)

    end do
    end do


    return
    end subroutine PrintDiv
! !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

 

    !**************************************************************************
    !               sunroutine for TDMA
    !**************************************************************************
    subroutine tdma(a,b,c,d,x,n)
    implicit none
    !        a - sub-diagonal (means it is the diagonal below the main
    !        diagonal)
    !        b - the main diagonal
    !        c - sup-diagonal (means it is the diagonal above the main
    !        diagonal)
    !        d - right part
    !        x - the answer
    !        n - number of equations

    integer,intent(in) :: n
    double precision,dimension(n),intent(in) :: a,b,c,d
    double precision,dimension(n),intent(out) :: x
    double precision,dimension(n) :: cp,dp
    double precision :: m
    integer i

    ! initialize c-prime and d-prime
    cp(1) = c(1)/b(1)
    dp(1) = d(1)/b(1)
    ! solve for vectors c-prime and d-prime
    do i = 2,n
        m = b(i)-cp(i-1)*a(i)
        cp(i) = c(i)/m
        dp(i) = (d(i)-dp(i-1)*a(i))/m
    enddo
    ! initialize x
    x(n) = dp(n)
    ! solve for x from the vectors c-prime and d-prime
    do i = n-1, 1, -1
        x(i) = dp(i)-cp(i)*x(i+1)
    enddo

    end subroutine tdma


    !<<<<<<<<<<<<<<< PrePrinting >>>>>>>>>>>>>>>
    subroutine PrePrinting
    use Global
    implicit none


    ! [PRE-PRITING BEGINS----------------------------
    write(*,*)'-------------------------------------------------'
    write(*,*) 'We are using'
    if(UNIFORM1==1 ) write(*,*)' Uniform mesh in x direction'
    if(UNIFORM2==1 ) write(*,*)' Uniform mesh in y direction'

    if(UNIFORM1==0 ) write(*,*)' Non-uniform in x with GAMMA1=', GAMMA1
    if(UNIFORM2==0 ) write(*,*)' Non-uniform in y with GAMMA2=', GAMMA2

    write(*,*) ' Conservative form for nonlinear term'
    write(*,*)'-------------------------------------------------'
    write(*,*)
    write(*,*)

    write(*,1995) Thetam
    write(*,1994) Pr, Ra
    write(*,1993) N1, N2, dt, maxcfl  
    write(*,1992) dPTOL          
    write(*,1991) XMIN, XMAX
    write(*,1990) YMIN, YMAX

    write(*,*)'-------------------------------------------------'
    write(*,*)
    write(*,*)    



    ! ----------------------------PRE-PRITING ENDS] 

1995 format('Thetam                             ->', 1e20.12)
1994 format('Pr, Ra                             ->', 2e20.12)
1993 format('N1, N2, dt, maxcfl                 ->', 2i8,2e20.12)
1992 format('Stopping criteria                  ->', 1e20.12)
1991 format('XMIN, XMAX                         ->', 2f20.12) 
1990 format('YMIN, YMAX                         ->', 2f20.12) 



    return 
    end subroutine PrePrinting





    !<<<<<<<<<<<<<<< Solve_dEnergy >>>>>>>>>>>>>>>        
    subroutine Solve_dEnergy(dTHETA)
    use omp_lib
    use Global
    implicit none

    integer :: i, j
    integer :: im, ip, jm, jp
    integer :: iep, iem, jem, jep

    double precision :: u1, u2, e1, e2
    double precision :: v3, v4, e3, e4

    double precision :: dedx1, dedx2
    double precision :: dedy3, dedy4

    double precision :: convect_e1, convect_e2, convect
    double precision :: viscous_e1, viscous_e2, viscous

    double precision :: ebc_in, ebc_out, ebc_down, ebc_up, ebc

    double precision :: eAMI, eACI, eAPI
    double precision :: eAMJ, eACJ, eAPJ
    double precision :: RHS_e


    double precision, dimension(0:N1, 0:N2) :: RHSe, dTHETAi
    double precision, dimension(0:N1, 0:N2) :: dTHETA

    double precision, dimension(1:N1m) :: ami, aci, api, adi, ri
    double precision, dimension(1:N2m) :: amj, acj, apj, adj, rj

    

    ! INITIALIZATION
    RHSe(0:N1, 0:N2) = 0.
    ! call omp_set_nested(.true.)


    ! For the RHS 

    !$omp parallel do                                                               &
    !$omp &  private(jp, jm, jem, jep, ip, im, iem, iep, j, i)                      & 
    !$omp &  private(u1, u2, v3, v4)                                                & 
    !$omp &  private(dedx1, dedx2)                                                  & 
    !$omp &  private(dedy3, dedy4)                                                  & 
    !$omp &  private(convect_e1, convect_e2, convect)                               & 
    !$omp &  private(viscous_e1, viscous_e2, viscous)                               & 
    !$omp &  private(eAMI, eACI, eAPI, eAMJ, eACJ, eAPJ)                            & 
    !$omp &  private(RHS_e, ebc_in, ebc_out, ebc_up, ebc_down, ebc)                 & 
    !$omp &   shared(RHSe)
    do j = 1, N2m
    jp = j + 1
    jm = j - 1
    jep = JP_E(j) - j
    jem = j - JM_E(j)

    do i = 1, N1m
    ip = i + 1
    im = i - 1
    iep = IP_E(i) - i 
    iem = i - IM_E(i)

    ! CONVECTION TERM
        u1 = Up(i ,j)
        u2 = Up(ip,j)
        dedx1 = (THETAp(i, j) - THETAp(im,j))/DMX(i )
        dedx2 = (THETAp(ip,j) - THETAp(i, j))/DMX(ip)  
       

        v3 = Vp(i,j )
        v4 = Vp(i,jp)
        dedy3 = (THETAp(i,j ) - THETAp(i,jm))/DMY(j )
        dedy4 = (THETAp(i,jp) - THETAp(i,j ))/DMY(jp)


    ! 1/2 ( Up(THETAp)x + Vp(THETAp)y + Wp(THETAp)z )
        convect_e1 = (u1*dedx1 + u2*dedx2)*0.5
        convect_e2 = (v3*dedy3 + v4*dedy4)*0.5  

        convect = 0.5*(convect_e1 + convect_e2)


    ! DIFFUSION TERM
        viscous_e1 = 1./DX(i)*(dedx2 - dedx1)
        viscous_e2 = 1./DY(j)*(dedy4 - dedy3)

        viscous = 0.5*Ce*(viscous_e1 + viscous_e2) 


    ! ebc
    ! From Convection Terms
        ! For X-direction
        ebc_in = 0.25*u1/DMX(i)*THETABC1(j)

        ebc_out = - 0.25*u2/DMX(ip)*THETABC2(j)

        ! For Y-direction
        ebc_down = 0.25*v3/DMY(j)*THETABC3(i)    

        ebc_up = -0.25*v4/DMY(jp)*THETABC4(i)   

    ! From Diffusion Terms
        ! For X-direction
        ebc_in = ebc_in + 0.5*Ce/DX(i)/DMX(i)*THETABC1(j)

        ebc_out = ebc_out + 0.5*Ce/DX(i)/DMX(ip)*THETABC2(j)

        ! For Y-direction
        ebc_down = ebc_down + 0.5*Ce/DY(j)/DMY(j)*THETABC3(i)        

        ebc_up = ebc_up + 0.5*Ce/DY(j)/DMY(jp)*THETABC4(i)        

        ebc = dble(1. - iem)*ebc_in         &
          & + dble(1. - iep)*ebc_out        &
          & + dble(1. - jem)*ebc_down       &
          & + dble(1. - jep)*ebc_up


    ! A THETAP                    

        ! X-DIRECTION
        eAPI = -0.5*Ce/DX(i)/DMX(ip)                                    &
           & + 0.25*u2/DMX(ip)  
        eAPI = eAPI*dble(iep)              

        eACI = 0.5*Ce/DX(i)*(1.0/DMX(ip) + 1.0/DMX(i))                  &
           & + (0.25*u1/DMX(i) - 0.25*u2/DMX(ip))          

        eAMI = -0.5*Ce/DX(i)/DMX(i)                                     &
           & - 0.25*u1/DMX(i)  
        eAMI = eAMI*dble(iem)    

        ! Y-DIRECTION
        eAPJ = -0.5*Ce/DY(j)/DMY(jp)                                    &
           & + 0.25*v4/DMY(jp)                
        eAPJ = eAPJ*dble(jep)

        eACJ = 0.5*Ce/DY(j)*(1.0/DMY(jp) + 1.0/DMY(j))                  &
           & + (0.25*v3/DMY(j ) - 0.25*v4/DMY(jp)) 

        eAMJ = -0.5*Ce/DY(j)/DMY(j)                                     &
            &  - 0.25*v3/DMY(j)                     
        eAMJ = eAMJ*dble(jem)

        RHS_e = eAPJ*THETAp(i,jp) + eACJ*THETAp(i,j) + eAMJ*THETAp(i,jm)      &
            & + eAPI*THETAp(ip,j) + eACI*THETAp(i,j) + eAMI*THETAp(im,j)

    ! RIGHT-HAND SIDE      
        RHSe(i,j) = THETAp(i,j)/dt - convect + viscous + ebc        &
                & - (THETAp(i,j)/dt + RHS_e)       


    end do
    end do
    !$omp end parallel do

    ! Compute the TDMA coefficients for x-, y-, and z- directions


    ! SOLVE IN Y-DIRECTION
     dTHETAi(0:N1, 0:N2) = 0.

    !$omp parallel do                                       &
    !$omp &  private(i, j, jp, jm, jep, jem)                &
    !$omp &  private(v3, v4)                                &
    !$omp &  private(apj, acj, amj, adj, rj)                & 
    !$omp &   shared(RHSe, dTHETAi)
    do i = 1, N1m

        do j = 1, N2m
        jp = j + 1
        jm = j - 1
        jep = JP_E(j) - j
        jem = j - JM_E(j)

        v3 = Vp(i,j )
        v4 = Vp(i,jp)

        apj(j) = -0.5*Ce/DY(j)/DMY(jp)                                &
             & + 0.25*v4/DMY(jp)      
        apj(j) = apj(j)*dble(jep)*dt

        acj(j) = 0.5*Ce/DY(j)*(1.0/DMY(jp) + 1.0/DMY(j)    )          &
             & + (0.25*v3/DMY(j ) - 0.25*v4/DMY(jp))  
        acj(j) = acj(j)*dt + 1.

        amj(j) = -0.5*Ce/DY(j)/DMY(j)                                 &
            &  - 0.25*v3/DMY(j)                    
        amj(j) = amj(j)*dble(jem)*dt

        adj(j) = RHSe(i,j)*dt 

        end do

        call TDMA(amj(1:N2m), acj(1:N2m), apj(1:N2m), adj(1:N2m), rj(1:N2m), N2m)
        
        do j = 1, N2m
            dTHETAi(i,j) = rj(j)
        end do

    end do
    !$omp end parallel do


    ! SOLVE IN X-DIRECTION
    dTHETA(0:N1, 0:N2) = 0.

    !$omp parallel do                                       &
    !$omp &  private(i, j, ip, im, iep, iem)                &
    !$omp &  private(u1, u2)                                &
    !$omp &  private(api, aci, ami, adi, ri)                & 
    !$omp &   shared(dTHETAi)
    do j = 1, N2m
        do i = 1, N1m
        ip = i + 1
        im = i - 1 
        iep = IP_E(i) - i 
        iem = i - IM_E(i)

        u1 = Up(i ,j)
        u2 = Up(ip,j)

        api(i) = -0.5*Ce/DX(i)/DMX(ip)                           &
             & + 0.25*u2/DMX(ip)                 
        api(i) = api(i)*dble(iep)*dt

        aci(i) = 0.5*Ce/DX(i)*(1.0/DMX(ip) + 1.0/DMX(i))         &
             & + (0.25*u1/DMX(i) - 0.25*u2/DMX(ip))          
        aci(i) = aci(i)*dt + 1.

        ami(i) = -0.5*Ce/DX(i)/DMX(i)                            &
             & - 0.25*u1/DMX(i)
        ami(i) = ami(i)*dble(iem)*dt

        adi(i) = dTHETAi(i,j) 

        end do

        call TDMA(ami(1:N1m), aci(1:N1m), api(1:N1m), adi(1:N1m), ri(1:N1m), N1m)

        do i = 1, N1m
            dTHETA(i,j) = ri(i)   
            dTHETA(i,j) = dTHETA(i,j) + THETAp(i,j)
        end do

    end do
    !$omp end parallel do


    return
    end subroutine Solve_dEnergy
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       



    !<<<<<<<<<<<<<<< SolveMomentum >>>>>>>>>>>>>>>        
    subroutine SolveMomentum(dTHETA, t, dUm, dVm)
    use omp_lib
    use Global
    implicit none

    double precision :: t

    double precision, dimension(0:N1, 0:N2) :: dTHETA
    double precision, dimension(0:N1, 0:N2) :: ddU, dUm, RHSx
    double precision, dimension(0:N1, 0:N2) :: ddV, dVm, RHSy

    integer :: i, j


    call RHS1(t, RHSx)
    call RHS2(dTHETA, t, RHSy)


    call GetddU(RHSx, t, ddU)
    call GetddV(RHSy, t, ddU, ddV)

    call GetdUVm(ddU, ddV, t, dUm, dVm)


    !X-DIRECTION
    !$omp parallel do 
    do j = 1, N2m 
        do i = 2, N1m 
            dUm(i,j) = dUm(i,j) + Up(i,j)
        end do 
    end do 
    !$omp end parallel do 

    !Y-DIRECTION
    !$omp parallel do 
    do j = 2, N2m 
        do i = 1, N1m 
            dVm(i,j) = dVm(i,j) + Vp(i,j)
        end do 
    end do 
    !$omp end parallel do 

    

    return

    end subroutine SolveMomentum
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 


    !<<<<<<<<<<<<<<< RHS1 >>>>>>>>>>>>>>>
    subroutine RHS1(t, RHSx)
    use omp_lib
    use Global
    implicit none

    integer :: i, j
    integer :: im, ip, ium, iup 
    integer :: jm, jp, jum, jup

    double precision :: t

    double precision :: u1, u2
    double precision :: v3, v4

    double precision :: dudx1, dudx2 
    double precision :: dudy3, dudy4

    double precision :: viscous_u1, viscous_u2, viscous_u12, viscous, dpdx
    double precision :: xbc_in, xbc_out, xbc_down, xbc_up, xbc
    double precision :: mAPI, mACI, mAMI, mAPJ, mACJ, mAMJ

    double precision :: M11Un, M12Vn

    double precision, dimension(0:N1, 0:N2) :: RHSx


    ! INITIALIZATION    
    RHSx(0:N1, 0:N2) = 0.
    
    ! call omp_set_nested(.true.)

    !$omp parallel do                                                               &
    !$omp &  private(jp, jm, jup, jum, ip, im, iup, ium, j, i)                      &
    !$omp &  private(u1, u2, v3, v4)                                                &
    !$omp &  private(dudx1, dudx2)                                                  &
    !$omp &  private(dudy3, dudy4)                                                  &
    !$omp &  private(viscous_u1, viscous_u2, viscous, dpdx)                         &
    !$omp &  private(xbc_in, xbc_out, xbc_down, xbc_up, xbc)                        &
    !$omp &  private(mAPI, mACI, mAMI, mAPJ, mACJ, mAMJ)                            &
    !$omp &  private(M11Un, M12Vn)                                                  &
    !$omp &  shared(RHSx) 
    do j = 1, N2m
    jp = j + 1
    jm = j - 1
    jup = JP_A(j) - j
    jum = j - JM_U(j)

    do i = 2, N1m
    ip = i + 1 
    im = i - 1 
    iup = IP_A(i) - i 
    ium = i - IM_U(i)

    ! DIFFUSION TERM
        dudx1 = (Up(i, j) - Up(im,j))/DX(im)
        dudx2 = (Up(ip,j) - Up(i, j))/DX(i )

        dudy3 = (Up(i,j ) - Up(i,jm))/DMY(j )
        dudy4 = (Up(i,jp) - Up(i,j ))/DMY(jp)


        viscous_u1 = 1./DMX(i)*(dudx2 - dudx1)
        viscous_u2 = 1./DY(j )*(dudy4 - dudy3)


        viscous = 0.5*Cm*(viscous_u1 + viscous_u2)

    ! PRESSURE TERM
        dpdx = (Pp(i,j) - Pp(im,j))/DMX(i)


    ! ubc
        ! FROM CONVECTION TERM
        ! X-direction
        u1 = 0.5*(Up(im,j) + Up(i ,j))
        u2 = 0.5*(UP(i ,j) + Up(ip,j))

        xbc_in = (0.25*u1/DX(im)*UBC1(j,1)                           &
             & -  0.25*dudx1*0.5*UBC1(j,1)) 
 
        xbc_out = (-0.25*u2/DX(i )*UBC2(j,1)                         &
              & -   0.25*dudx2*0.5*UBC2(j,1)) 


        ! Y-direction
        v3 = (DX(im)*Vp(i,j ) + DX(i)*Vp(im,j ))*(0.5/DMX(i))
        v4 = (DX(im)*Vp(i,jp) + DX(i)*Vp(im,jp))*(0.5/DMX(i))

        xbc_down = (0.25*v3/DMY(j)*UBC3(i,1)                                             &
               & -  0.25*dudy3*(DX(i)*UBC3(im,2) + DX(im)*UBC3(i,2))*(0.5/DMX(i)))  
 
        xbc_up = (-0.25*v4/DMY(jp)*UBC4(i,1)                                             &
             & -   0.25*dudy4*(DX(i)*UBC4(im,2) + DX(im)*UBC4(i,2))*(0.5/DMX(i)))  


        ! FROM DIFFUSION TERM
        ! X-direction         
        xbc_in  = xbc_in  + 0.5*Cm/DMX(i)/DX(im)*UBC1(j,1)
        xbc_out = xbc_out + 0.5*Cm/DMX(i)/DX(i )*UBC2(j,1)


        ! Y-direction                
        xbc_down = xbc_down + 0.5*Cm/DY(j)/DMY(j )*UBC3(i,1)
        xbc_up   = xbc_up   + 0.5*Cm/DY(j)/DMY(jp)*UBC4(i,1)

        xbc = dble(1. - ium)*xbc_in      &
          & + dble(1. - iup)*xbc_out     & 
          & + dble(1. - jum)*xbc_down    &
          & + dble(1. - jup)*xbc_up


    ! M11Un
        ! X-direction

        ! 'u1, u2, dudx1, dudx2' HAVE BEEN DEFINED

        mAPI = -0.5*Cm/DMX(i)/DX(i)                        &
           & + (0.25*u2/DX(i)  + 0.25*dudx2*0.5)  
        mAPI = mAPI*dble(iup) 

        mACI = 0.5*Cm/DMX(i)*(1.0/DX(i) + 1.0/DX(im))      &
           & + (0.25*dudx2*0.5 + 0.25*dudx1*0.5             &
           & -  0.25*u2/DX(i)  + 0.25*u1/DX(im) )    

        mAMI = -0.5*Cm/DMX(i)/DX(im)                       &
           & + (-0.25*u1/DX(im) + 0.25*dudx1*0.5) 
        mAMI = mAMI*dble(ium)  

        ! Y-direction
        
        ! 'v3, v4' HAVE BEEN DEFINED
        mAPJ = -0.5*Cm/DY(j)/DMY(jp)                       &
           & + 0.25*v4/DMY(jp)                      
        mAPJ = mAPJ*dble(jup)

        mACJ = 0.5*Cm/DY(j)*(1.0/DMY(jp) + 1.0/DMY(j))     &
           & + (-0.25*v4/DMY(jp) + 0.25*v3/DMY(j ) )   

        mAMJ = -0.5*Cm/DY(j)/DMY(j)                        &
           & - 0.25*v3/DMY(j)                       
        mAMJ = mAMJ*dble(jum)


        M11Un = mAPI*Up(ip,j ) + mACI*Up(i,j) + mAMI*Up(im,j )      &
            & + mAPJ*Up(i, jp) + mACJ*Up(i,j) + mAMJ*Up(i, jm)   

    ! M12Vn
    
    ! 'v3, v4, dudy3, dudy4' HAVE BEEN DEFINED

        M12Vn = 0.25*(v4*dudy4*dble(jup) + v3*dudy3*dble(jum))                          


    ! Final RHS for the momentum equation in x-direction
        RHSx(i,j) =  Up(i,j)/dt + viscous + xbc - dpdx    &
                & - (Up(i,j)/dt + M11Un + M12Vn)        

    end do
    end do
    !$omp end parallel do      

    return
    end subroutine RHS1
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    !<<<<<<<<<<<<<<< RHS2 >>>>>>>>>>>>>>>
    subroutine RHS2(dTHETA, t, RHSy)
    use omp_lib
    use Global
    implicit none

    integer :: i, j, ii, jj 
    integer :: im, ip, ivp, ivm 
    integer :: jm, jp, jvp, jvm

    double precision :: t
    double precision :: xl0, yl0 

    double precision :: u1, u2
    double precision :: v3, v4
    double precision :: ec

    double precision :: dvdx1, dvdx2
    double precision :: dvdy3, dvdy4


    double precision :: viscous_v1, viscous_v2, viscous, dpdy, THETAy
    double precision :: ybc_in, ybc_out, ybc_down, ybc_up, ybc
    double precision :: mAPI, mACI, mAMI, mAPJ, mACJ, mAMJ

    double precision :: M21Un, M22Vn

    double precision, dimension(0:N1, 0:N2) :: dTHETA
    double precision, dimension(0:N1, 0:N2) :: RHSy

    double precision :: xc, yc 

    ! INITIALIZATION    
    RHSy(0:N1, 0:N2) = 0.

    ! call omp_set_nested(.true.)

    !$omp parallel do                                                             &
    !$omp &  private(jp, jm, jvp, jvm, ip, im, ivp, ivm, j, i)                    &
    !$omp &  private(u1, u2, v3, v4, ec)                                          &
    !$omp &  private(dvdx1, dvdx2, dvdy3, dvdy4)                                  &
    !$omp &  private(viscous_v1, viscous_v2, viscous, dpdy, THETAy)               &
    !$omp &  private(ybc_in, ybc_out, ybc_down, ybc_up, ybc)                      &
    !$omp &  private(mAPI, mACI, mAMI, mAPJ, mACJ, mAMJ)                          &
    !$omp &  private(M21Un, M22Vn)                                                &
    !$omp &   shared(RHSy, dTHETA)   
    do j = 2, N2m
    jp = j + 1
    jm = j - 1
    jvp = JP_A(j) - j
    jvm = j - JM_V(j)
    do i = 1, N1m
    ip = i + 1 
    im = i - 1 
    ivp = IP_A(i) - i 
    ivm = i - IM_V(i)

    
    ! DIFFUSION TERM
        dvdx1 = (Vp(i, j) - Vp(im,j))/DMX(i )
        dvdx2 = (Vp(ip,j) - Vp(i, j))/DMX(ip)
    
        dvdy3 = (Vp(i,j ) - Vp(i,jm))/DY(jm)
        dvdy4 = (Vp(i,jp) - Vp(i,j ))/DY(j )

        viscous_v1 = 1./DX(i )*(dvdx2 - dvdx1)
        viscous_v2 = 1./DMY(j)*(dvdy4 - dvdy3)

        viscous = 0.5*Cm*(viscous_v1 + viscous_v2)

    ! PRESSURE TERM
        dpdy = (Pp(i,j) - Pp(i,jm))/DMY(j)

    ! BUOYANCY TERM
        ec = (DY(jm)*dTHETA(i,j) + DY(j)*dTHETA(i,jm))*(0.5/DMY(j))
        THETAy = Cb*(ec - Thetam)

        xc = 0.5*(X(i) + X(ip))
        yc = Y(j)

        do jj = 1, Nc
            yl0 = Y0(jj)
            do ii = 1, Nc 
                xl0 = X0(ii)

                if((xc - xl0)**2.0 + (yc - yl0)**2.0 < r**2.0) then 
                    THETAy = 0.0
                end if 
            end do 
        end do 

    ! vbc
        ! FROM CONVECTION TERM
        ! X-direction
        u1 = 0.5/DMY(j)*(DY(j)*Up(i ,jm) + DY(jm)*Up(i ,j))
        u2 = 0.5/DMY(j)*(DY(j)*Up(ip,jm) + DY(jm)*Up(ip,j))

        ybc_in  = 0.25*u1/DMX(i)*UBC1(j,2)                                          & 
              & - 0.25*dvdy3*0.5/DMY(j)*(DY(j)*UBC1(jm,1) + DY(jm)*UBC1(j,1))

        ybc_out = -0.25*u2/DMX(ip)*UBC2(j,2)                                        &
              & - 0.25*dvdy4*0.5/DMY(j)*(DY(j)*UBC2(jm,1) + DY(jm)*UBC2(j,1))

        ! Y-direction
        v3 = 0.5*(Vp(i,jm) + Vp(i,j ))
        v4 = 0.5*(Vp(i,j ) + Vp(i,jp))

        ybc_down = 0.25/DY(jm)*v3*UBC3(i,2)        &
               & - 0.25*0.5*dvdy3*UBC3(i,2)

        ybc_up = -0.25/DY(j )*v4*UBC4(i,2)         &
             & -  0.25*0.5*dvdy4*UBC4(i,2) 

        ! FROM DIFFUSION TERM
        ! X-direction
        ybc_in  = ybc_in  + 0.5*Cm/DX(i)/DMX(i )*UBC1(j,2)
        ybc_out = ybc_out + 0.5*Cm/DX(i)/DMX(ip)*UBC2(j,2)

        ! Y-direction
        ybc_down = ybc_down + 0.5*Cm/DMY(j)/DY(jm)*UBC3(i,2)
        ybc_up   = ybc_up   + 0.5*Cm/DMY(j)/DY(j )*UBC4(i,2)


        ybc = dble(1. - ivm)*ybc_in        &
          & + dble(1. - ivp)*ybc_out       &
          & + dble(1. - jvm)*ybc_down      &
          & + dble(1. - jvp)*ybc_up


    ! M21Un 
    ! 'u1, u2, dvdx1, dvdx2' HAVE BEEN DEFINED
        M21Un = 0.25*(u2*dvdx2*dble(ivp) + u1*dvdx1*dble(ivm))   

    ! M22Vn

        ! Y-DIRECTION

        ! 'v3, v4, dvdy3, dvdy4' HAVE BEEN DEFINED
        mAPJ = -0.5*Cm/DMY(j)/DY(j)                     &
           & + (0.25*v4/DY(j)  + 0.25*dvdy4*0.5)    
        mAPJ = mAPJ*dble(jvp)

        mACJ = 0.5*Cm/DMY(j)*(1.0/DY(j) + 1.0/DY(jm))   &
           & + (0.25*dvdy4*0.5 + 0.25*dvdy3*0.5          &
           & -  0.25*v4/DY(j)  + 0.25*v3/DY(jm) )      

        mAMJ = -0.5*Cm/DMY(j)/DY(jm)                    &
           & + (-0.25*v3/DY(jm) + 0.25*dvdy3*0.5)    
        mAMJ = mAMJ*dble(jvm)


        ! X-Direction
        ! 'u1, u2' HAVE BEEN DEFINED 
        mAPI = -0.5*Cm/DX(i)/DMX(ip)                   &
           & + 0.25*u2/DMX(ip)   
        mAPI = mAPI*dble(ivp)                  

        mACI = 0.5*Cm/DX(i)*(1.0/DMX(ip) + 1.0/DMX(i)) &
           & + (-0.25*u2/DMX(ip) + 0.25*u1/DMX(i ))                      

        mAMI = -0.5*Cm/DX(i)/DMX(i)                    &
           & - 0.25*u1/DMX(i)  
        mAMI = mAMI*dble(ivm)                      
  

        M22Vn = mAPI*Vp(ip,j ) + mACI*Vp(i,j) + mAMI*Vp(im,j )      &
            & + mAPJ*Vp(i, jp) + mACJ*Vp(i,j) + mAMJ*Vp(i, jm)  

    
    ! Final RHS for the momentum equation in y-direction
        RHSy(i,j) =   Vp(i,j)/dt + viscous + ybc - dpdy + THETAy        &
                & - ( Vp(i,j)/dt + M21Un + M22Vn)       

    end do
    end do
    !$omp end parallel do


    return
    end subroutine RHS2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


    !<<<<<<<<<<<<<<< GetddU >>>>>>>>>>>>>>>
    subroutine GetddU(RHSx, t, ddU)
    use omp_lib
    use Global
    implicit none

    integer :: i, j
    integer :: im, ip, ium, iup 
    integer :: jm, jp, jum, jup

    double precision :: t


    double precision :: u1, u2
    double precision :: v3, v4

    double precision :: dudx1, dudx2 

    double precision, dimension(1:N1m) :: ami, aci, api, adi, ri
    double precision, dimension(1:N2m) :: amj, acj, apj, adj, rj

    double precision, dimension(0:N1, 0:N2) :: RHSx, ddU, ddUj

    ! INITIALIZATION
    ddU(0:N1, 0:N2) = 0.
    ddUj(0:N1, 0:N2) = 0.

    ! call omp_set_nested(.true.)

    ! SOLVE IN X-DIRECTION

    !$omp parallel do                           &
    !$omp &  private(j, i, ip, im, iup, ium)    &
    !$omp &  private(u1, u2)                    &
    !$omp &  private(dudx1, dudx2)              &
    !$omp &  private(api, aci, ami, adi, ri)    &
    !$omp &   shared(RHSx)
    do j = 1, N2m
    do i = 2, N1m
    ip = i + 1
    im = i - 1 
    iup = IP_A(i) - i 
    ium = i - IM_U(i)

        u1 = 0.5*(Up(im,j) + Up(i,j))
        u2 = 0.5*(Up(ip,j) + Up(i,j))

        dudx1 = (Up(i, j) - Up(im,j))/DX(im)
        dudx2 = (Up(ip,j) - Up(i, j))/DX(i )

        api(i) = -0.5*Cm/DMX(i)/DX(i)                        &
             & + (0.25*u2/DX(i)  + 0.25*dudx2*0.5)    
        api(i) = api(i)*dble(iup)*dt           

        aci(i) = 0.5*Cm/DMX(i)*(1.0/DX(i) + 1.0/DX(im))      &
             & + (0.25*dudx2*0.5 + 0.25*dudx1*0.5             &
             & -  0.25*u2/DX(i)  + 0.25*u1/DX(im) )      
        aci(i) = aci(i)*dt + 1.

        ami(i) = -0.5*Cm/DMX(i)/DX(im)                       &
             & + (-0.25*u1/DX(im) + 0.25*dudx1*0.5)   
        ami(i) = ami(i)*dble(ium)*dt

        adi(i) = RHSx(i,j)*dt

        end do

        call TDMA(ami(2:N1m), aci(2:N1m), api(2:N1m), adi(2:N1m), ri(2:N1m), N1m-1)

        do i = 2, N1m
            ddUj(i,j) = ri(i)
        end do

    end do
    !$omp end parallel do



    ! SOLVE IN Y-DIRECTION
    !$omp parallel do                                       &
    !$omp &  private(i, im, j, jp, jm, jup, jum)            &
    !$omp &  private(v3, v4)                                &
    !$omp &  private(apj, acj, amj, adj, rj)                &
    !$omp &   shared(ddUj)   
    do i = 2, N1m
    im = i - 1
    ip = i + 1 
    do j = 1, N2m
    jp = j + 1
    jm = j - 1
    jup = JP_A(j) - j
    jum = j - JM_U(j)


        v3 = (DX(im)*Vp(i,j ) + DX(i)*Vp(im,j ))*(0.5/DMX(i))
        v4 = (DX(im)*Vp(i,jp) + DX(i)*Vp(im,jp))*(0.5/DMX(i))

        apj(j) = -0.5*Cm/DY(j)/DMY(jp)                         &
             & + 0.25*v4/DMY(jp)                   
        apj(j) = apj(j)*dble(jup)*dt

        acj(j) = 0.5*Cm/DY(j)*(1.0/DMY(jp) + 1.0/DMY(j))       &
             & + (-0.25*v4/DMY(jp) + 0.25*v3/DMY(j ) )  
        acj(j) = acj(j)*dt + 1.

        amj(j) = -0.5*Cm/DY(j)/DMY(j)                          &
             & - 0.25*v3/DMY(j)                     
        amj(j) = amj(j)*dble(jum)*dt

        adj(j) = ddUj(i,j)

    end do

        call TDMA(amj(1:N2m), acj(1:N2m), apj(1:N2m), adj(1:N2m), rj(1:N2m), N2m)

        do j = 1, N2m
            ddU(i,j) = rj(j)
        end do

    end do
    !$omp end parallel do


    return
    end subroutine GetddU
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    


    !<<<<<<<<<<<<<<< GetddV >>>>>>>>>>>>>>>
    subroutine GetddV(RHSy, t, ddU, ddV)
    use omp_lib
    use Global
    implicit none

    integer :: i, j
    integer :: im, ip, ivp, ivm 
    integer :: jm, jp, jvp, jvm

    double precision :: t

    double precision :: u1, u2
    double precision :: v3, v4

    double precision :: dvdx1, dvdx2
    double precision :: dvdy3, dvdy4

    double precision :: ddu1, ddu2
    
    double precision :: M21ddU

    double precision, dimension(1:N1m) :: ami, aci, api, adi, ri
    double precision, dimension(2:N2m) :: amj, acj, apj, adj, rj

    double precision, dimension(0:N1, 0:N2) :: RHSy, ddU, ddV, ddVi

    ! INITIALIZATION
    ddV(0:N1, 0:N2) = 0.
    ddVi(0:N1, 0:N2) = 0.

    ! call omp_set_nested(.true.)


    ! Y-direction
    !$omp parallel do                                               &
    !$omp &  private(i, ip, im, ivp, ivm, j, jp, jm, jvp, jvm)      &
    !$omp &  private(v3, v4)                                        &
    !$omp &  private(ddu1, ddu2)                                    &
    !$omp &  private(dvdx1, dvdx2)                                  &
    !$omp &  private(dvdy3, dvdy4)                                  &
    !$omp &  private(apj, acj, amj, adj, rj, M21ddU)                &
    !$omp &  shared(RHSy) 
    do i = 1, N1m
    ip = i + 1
    im = i - 1 
    ivp = IP_A(i) - i 
    ivm = i - IM_V(i)
    do j = 2, N2m
    jp = j + 1
    jm = j - 1
    jvp = JP_A(j) - j
    jvm = j - JM_V(j)

        v3 = 0.5*(Vp(i,jm) + Vp(i,j ))
        v4 = 0.5*(Vp(i,j ) + Vp(i,jp))

        dvdy3 = (Vp(i,j ) - Vp(i,jm))/DY(jm)
        dvdy4 = (Vp(i,jp) - Vp(i,j ))/DY(j )

        apj(j) = -0.5*Cm/DMY(j)/DY(j)                         &
             & + (0.25*v4/DY(j)  + 0.25*dvdy4*0.5)  
        apj(j) = apj(j)*dble(jvp)*dt

        acj(j) = 0.5*Cm/DMY(j)*(1.0/DY(j) + 1.0/DY(jm))       &
             & + (0.25*dvdy4*0.5 + 0.25*dvdy3*0.5                &
             & -  0.25*v4/DY(j)  + 0.25*v3/DY(jm) )     
        acj(j) = acj(j)*dt + 1.


        amj(j) = -0.5*Cm/DMY(j)/DY(jm)                        &
             & + (-0.25*v3/DY(jm) + 0.25*dvdy3*0.5)  
        amj(j) = amj(j)*dble(jvm)*dt

        ! M21 ddU
        ddu1 = (DY(jm)*ddU(i, j) + DY(j)*ddU(i, jm))*(0.5/DMY(j))
        ddu2 = (DY(jm)*ddU(ip,j) + DY(j)*ddU(ip,jm))*(0.5/DMY(j))

        dvdx1 = (Vp(i, j) - Vp(im,j))/DMX(i )
        dvdx2 = (Vp(ip,j) - Vp(i, j))/DMX(ip)

        M21ddU = 0.25*(ddu2*dvdx2*dble(ivp) + ddu1*dvdx1*dble(ivm))           

        adj(j) = RHSy(i,j)*dt - M21ddU*dt

    end do

        call TDMA(amj(2:N2m), acj(2:N2m), apj(2:N2m), adj(2:N2m), rj(2:N2m), N2m-1)

        do j = 2, N2m
          ddVi(i,j) = rj(j)
        end do

    end do
    !$omp end parallel do


    ! X-direction
    !$omp parallel do                                     &
    !$omp &  private(j, jm, i, ip, im, ivp, ivm)          &
    !$omp &  private(u1, u2)                              &
    !$omp &  private(api, aci, ami, adi, ri)              &
    !$omp &  shared(ddVi)
    do j = 2, N2m
    jm = j - 1
    jp = j + 1
    do i = 1, N1m
    ip = i + 1 
    im = i - 1 
    ivp = IP_A(i) - i 
    ivm = i - IM_V(i)

    
        u1 = (DY(jm)*Up(i, j) + DY(j)*Up(i, jm))*(0.5/DMY(j))
        u2 = (DY(jm)*Up(ip,j) + DY(j)*Up(ip,jm))*(0.5/DMY(j))

        api(i) = -0.5*Cm/DX(i)/DMX(ip)                         &
             & + 0.25*u2/DMX(ip)                   
        api(i) = api(i)*dble(ivp)*dt

        aci(i) = 0.5*Cm/DX(i)*(1.0/DMX(ip) + 1.0/DMX(i))       &
             & + (-0.25*u2/DMX(ip) + 0.25*u1/DMX(i ))                
        aci(i) = aci(i)*dt + 1.

        ami(i) = -0.5*Cm/DX(i)/DMX(i)                          &
             & - 0.25*u1/DMX(i)                  
        ami(i) = ami(i)*dble(ivm)*dt  

        adi(i) = ddVi(i,j)

    end do

        call TDMA(ami(1:N1m), aci(1:N1m), api(1:N1m), adi(1:N1m), ri(1:N1m), N1m)

        do i = 1, N1m
          ddV(i,j) = ri(i)
        end do

    end do
    !$omp end parallel do


    return
    end subroutine GetddV
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  




    !<<<<<<<<<<<<<<< GetdUVWm >>>>>>>>>>>>>>>
    subroutine GetdUVm(ddU, ddV, t, dUm, dVm)
    use omp_lib
    use Global
    implicit none

    double precision :: t

    integer :: i, j
    integer :: ip, im
    integer :: jp, jm, jum, jup

    double precision :: dvm3, dvm4
    double precision :: dudy3, dudy4
    double precision :: M12dVm

    double precision, dimension(0:N1, 0:N2) :: ddU, ddV, dUm, dVm


    ! call omp_set_nested(.true.)

    dUm = 0.0 
    dVm = 0.0  

    !$omp parallel do private(i, j)
    do j = 2, N2m
    do i = 1, N1m
        dVm(i,j) = ddV(i,j)
    end do
    end do
    !$omp end parallel do


    !$omp parallel do                                             &
    !$omp &  private(j, jp, jm, jup, jum, i, im)                  &
    !$omp &  private(dvm3, dvm4)                                  &
    !$omp &  private(dudy3, dudy4)                                &
    !$omp &  private(M12dVm)                                      &
    !$omp &  shared(ddU)
    do j = 1, N2m
    jp = j + 1
    jm = j - 1
    jup = JP_A(j) - j
    jum = j - JM_U(j)

    do i = 2, N1m
    im = i - 1

        ! M12 dVm
        dvm3 = (DX(i)*dVm(im,j ) + DX(im)*dVm(i,j ))*(0.5/DMX(i))
        dvm4 = (DX(i)*dVm(im,jp) + DX(im)*dVm(i,jp))*(0.5/DMX(i))
        
        dudy3 = (Up(i,j ) - Up(i,jm))/DMY(j )
        dudy4 = (Up(i,jp) - Up(i,j ))/DMY(jp)

        M12dVm = 0.25*(dvm4*dudy4*dble(jup) + dvm3*dudy3*dble(jum))                   
  

        dUm(i,j) = ddU(i,j) - dt*M12dVm


    end do
    end do
    !$omp end parallel do


    
    return
    end subroutine GetdUVm
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  



    !<<<<<<<<<<<<<<< PoissonSolverMG >>>>>>>>>>>>>>>
    subroutine PoissonSolverMG(f, dP, t)
    use Global
    implicit none

    integer :: i, j
    integer :: im, ip, jm, jp

    integer :: Vlevel
    integer :: iter, PoiMaxIter

    double precision :: t, dPNorm, tmp
    double precision :: OMEGA, tep, tempx(2,2), tempy(2,2), temp(2,2)

    double precision, dimension(1:N1m, 1:N2m) :: f
    double precision, dimension(0:N1,  0:N2 ) :: dP, dPold

    OMEGA = 2./(1.+pi/sqrt(2.0)/((dble(N1m*N2m))**(1./2.)))

    dPNorm = 0.
    dPold = dP

    tempx(1,1) = (dmx(1 )+dmx(2  ))**2./((dmx(1 )+dmx(2  ))**2.-dmx(1 )**2.)
    tempx(2,1) = dmx(1 )**2.           /((dmx(1 )+dmx(2  ))**2.-dmx(1 )**2.)
    tempx(1,2) = (dmx(N1)+dmx(N1m))**2./((dmx(N1)+dmx(N1m))**2.-dmx(N1)**2.)
    tempx(2,2) = dmx(N1)**2.           /((dmx(N1)+dmx(N1m))**2.-dmx(N1)**2.)

    tempy(1,1) = (dmy(1 )+dmy(2  ))**2./((dmy(1 )+dmy(2  ))**2.-dmy(1 )**2.)
    tempy(2,1) = dmy(1 )**2.           /((dmy(1 )+dmy(2  ))**2.-dmy(1 )**2.)
    tempy(1,2) = (dmy(N2)+dmy(N2m))**2./((dmy(N2)+dmy(N2m))**2.-dmy(N2)**2.)
    tempy(2,2) = dmy(N2)**2.           /((dmy(N2)+dmy(N2m))**2.-dmy(N2)**2.)

    
    do iter = 1, DpNum
       
        Vlevel = 1

        call MG_Vcycle(N1m, N2m, f, dP, DX, DY, Vlevel)


!         do j = 1, N2m
!         jp = j + 1
!         jm = j - 1
!         do i = 1, N1m
!         ip = i + 1
!         im = i - 1

!             tmp = alpha(i,j,1)*dP(im,j) + alpha(i,j,2)*dP(ip,j)     &
!               & + alpha(i,j,3)*dP(i,jm) + alpha(i,j,4)*dP(i,jp)   

!             dP(i,j) = OMEGA*(f(i,j) - tmp)/alpha(i,j,5) + (1.d0 - OMEGA)*dP(i,j)

!         end do
!         end do

        dPNorm = sqrt(sum((dP(1:N1m,1:N2m) - dpold(1:N1m,1:N2m))**2.)/dble(N1m*N2m))/   &
               & sqrt(sum((dP(1:N1m,1:N2m))**2.)/dble(N1m*N2m))

        dPold = dP

!         write(*,*) iter, dPNorm

        if(dPNorm .le. DPTOL) then
        ! UPDATE THE BOUNDARY INFORMATION
        tempx(1,1) = (DMX(1 )+DMX(2  ))**2./((DMX(1 )+DMX(2  ))**2.-DMX(1 )**2.)
        tempx(2,1) = DMX(1)**2.            /((DMX(1 )+DMX(2  ))**2.-DMX(1 )**2.)
        tempx(1,2) = (DMX(N1)+DMX(N1m))**2./((DMX(N1)+DMX(N1m))**2.-DMX(N1)**2.)
        tempx(2,2) =  DMX(N1)**2.          /((DMX(N1)+DMX(N1m))**2.-DMX(N1)**2.)  

        tempy(1,1) = (DMY(1 )+DMY(2  ))**2./((DMY(1 )+DMY(2  ))**2.-DMY(1 )**2.)
        tempy(2,1) = DMY(1 )**2.           /((DMY(1 )+DMY(2  ))**2.-DMY(1 )**2.)
        tempy(1,2) = (DMY(N2)+DMY(N2m))**2./((DMY(N2)+DMY(N2m))**2.-DMY(N2)**2.)
        tempy(2,2) = DMY(N2)**2.           /((DMY(N2)+DMY(N2m))**2.-DMY(N2)**2.)


        ! x = xmin & xmax
        dP(0, 1:N2m) = tempx(1,1)*dP(1,   1:N2m) - tempx(2,1)*dP(2,   1:N2m)     
        dP(N1,1:N2m) = tempx(1,2)*dP(N1-1,1:N2m) - tempx(2,2)*dP(N1-2,1:N2m)

        ! y = ymin & ymax
        dP(1:N1m,0 ) = tempy(1,1)*dP(1:N1m,1   ) - tempy(2,1)*dP(1:N1m,2   )
        dP(1:N1m,N2) = tempy(1,2)*dP(1:N1m,N2-1) - tempy(2,2)*dP(1:N1m,N2-2)
    

        ! POINTS
        dP(0, 0) = 0.5*tempy(1,1)*dP(0,   1) - 0.5*tempy(2,1)*dP(0,   2)       &
               & + 0.5*tempx(1,1)*dP(1,   0) - 0.5*tempx(2,1)*dP(2,   0)
        dP(N1,0) = 0.5*tempy(1,1)*dP(N1,  1) - 0.5*tempy(2,1)*dP(N1,  2)       &
               & + 0.5*tempx(1,2)*dP(N1-1,0) - 0.5*tempx(2,2)*dP(N1-2,0)
    
        dP(0, N2) = 0.5*tempy(1,2)*dP(0,   N2-1) - 0.5*tempy(2,2)*dP(0,   N2-2)       &
                & + 0.5*tempx(1,1)*dP(1,   N2  ) - 0.5*tempx(2,1)*dP(2,   N2  )       
        dP(N1,N2) = 0.5*tempy(1,2)*dP(N1,  N2-1) - 0.5*tempy(2,2)*dP(N1,  N2-2)       &
                & + 0.5*tempx(1,2)*dP(N1-1,N2  ) - 0.5*tempx(2,2)*dP(N1-2,N2  )


        PoiMaxIter = iter
        go to 1001

        end if

    end do

    write(*,*) 'Poisson equation is not converged'

1001 continue

    return
    end subroutine PoissonSolverMG
! !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    !<<<<<<<<<<<<<<< MG_Vcycle >>>>>>>>>>>>>>>
    !*************************************************************************  
    !                  recursive subroutine MG_Vcycle
    !
    !       get the V-Cycle recursive subroutine with by using 
    !                 Gauss-Seidel MG_Relaxation method.
    !*************************************************************************  
    recursive subroutine MG_Vcycle(ng1, ng2, mg_f, mg_p, mg_DX, mg_DY, Vlevel)
    use Global
    implicit none
    
    integer :: i, j
    integer :: ng1, ng2, Vlevel

    double precision, dimension(ng1, ng2) :: mg_f
    double precision, dimension(ng1, ng2, 5) :: mg_alpha
    double precision, dimension(0:ng1+1, 0:ng2+1) :: mg_p

    double precision, dimension(0:ng1+1) :: mg_DX 
    double precision, dimension(0:ng2+1) :: mg_DY

    double precision, allocatable, dimension(:,:) :: mg_rhs, mg_v, mg_vf

    allocate(mg_rhs(ng1/2,ng2/2))
    allocate(mg_v(0:ng1/2+1,0:ng2/2+1))
    allocate(mg_vf(0:ng1+1,0:ng2+1))

    ! Coefficients in the multi-grid spacing
    call MG_CoePoisson(ng1, ng2, mg_DX, mg_DY, mg_alpha)

    ! Pre-smoothing
    call MG_Relax(ng1, ng2, mg_DX, mg_DY, mg_p, mg_f, mg_alpha, Vlevel)

    if (Vlevel < MaxLevel) then

        ! Restriction
        call MG_Restriction(ng1, ng2, mg_DX, mg_DY, mg_p, mg_f, mg_rhs, mg_alpha, Vlevel)

        ! Recursive Procedure
        mg_v = 0.
        mg_DX = MGDX(0:ng1/2+1,Vlevel+1)
        mg_DY = MGDY(0:ng2/2+1,Vlevel+1)
        call MG_Vcycle(ng1/2, ng2/2, mg_rhs, mg_v, mg_DX, mg_DY, Vlevel+1)

        ! Interpolation or Prolongation
        mg_vf = 0.
        call MG_Prolongation(ng1/2, ng2/2, mg_v, mg_vf, Vlevel+1)

        do j = 1, ng2
        do i = 1, ng1
            mg_p(i,j) = mg_p(i,j) + mg_vf(i,j)
        end do
        end do

        ! Post-Smoothing
        mg_DX = MGDX(0:ng1+1, Vlevel)
        mg_DY = MGDY(0:ng2+1, Vlevel)
        call MG_Relax(ng1, ng2, mg_DX, mg_DY, mg_p, mg_f, mg_alpha, Vlevel)


    end if

    call MG_Uniq(mg_p, ng1, ng2)

    deallocate(mg_rhs, mg_v, mg_vf)


    end subroutine MG_Vcycle
! !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    !<<<<<<<<<<<<<<< MG_CoePoisson >>>>>>>>>>>>>>>
    subroutine MG_CoePoisson(ng1, ng2, mg_DX, mg_DY, mg_alpha)
    implicit none

    integer :: i,   j
    integer :: ip,  jp
    integer :: ng1, ng2

    double precision, dimension(0:ng1+1) :: mg_DX, mg_DMX
    double precision, dimension(0:ng2+1) :: mg_DY, mg_DMY

    double precision, dimension(ng1,ng2,5) :: mg_alpha

    mg_DMX(0) = 0.
    do i = 1, ng1+1
      mg_DMX(i) = 0.5*(mg_DX(i-1) + mg_DX(i))
    end do  
    mg_DMX(1) = mg_DMX(1) + mg_DMX(ng1+1)
    mg_DMX(ng1+1) = mg_DMX(1) 

    mg_DMY(0) = 0.
    do j = 1, ng2+1
      mg_DMY(j) = 0.5*(mg_DY(j-1) + mg_DY(j))
    end do 

    
    do j = 1, ng2
    jp = j + 1
    do i = 1, ng1
    ip = i + 1

        mg_alpha(i,j,1) = 1./mg_DX(i)/mg_DMX(i )
        mg_alpha(i,j,2) = 1./mg_DX(i)/mg_DMX(ip)
        mg_alpha(i,j,3) = 1./mg_DY(j)/mg_DMY(j )
        mg_alpha(i,j,4) = 1./mg_DY(j)/mg_DMY(jp)

        if (j==1  )    mg_alpha(i,j,3) = 0.
        if (j==ng2)    mg_alpha(i,j,4) = 0.

        mg_alpha(i,j,5) = - ( mg_alpha(i,j,1) + mg_alpha(i,j,2) + mg_alpha(i,j,3) + mg_alpha(i,j,4))


    end do
    end do        

    return
    end subroutine MG_CoePoisson
! !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    !<<<<<<<<<<<<<<< MG_Relax >>>>>>>>>>>>>>>
    subroutine MG_Relax(ng1, ng2, mg_DX, mg_DY, mg_p, mg_f, mg_alpha, level)
    use Global
    implicit none

    integer :: i, j
    integer :: ip, im, jp, jm
    integer :: ng1, ng2
    integer :: iter, level

    double precision :: tmp, OMEGA

    double precision, dimension(ng1, ng2) :: mg_f
    double precision, dimension(ng1, ng2, 5) :: mg_alpha
    double precision, dimension(0:ng1+1, 0:ng2+1) :: mg_p, mg_po
    double precision, dimension(0:ng1+1) :: mg_DX
    double precision, dimension(0:ng2+1) :: mg_DY

    OMEGA = 1.

    do iter = 1, RelaxNum
        mg_po = mg_p

    
        do j = 1, ng2
        jp = j + 1
        jm = j - 1
        do i = 1, ng1
        ip = i + 1  
        im = i - 1    

            tmp = mg_alpha(i,j,1)*mg_p(im,j) + mg_alpha(i,j,2)*mg_po(ip,j)       &
              & + mg_alpha(i,j,3)*mg_p(i,jm) + mg_alpha(i,j,4)*mg_po(i,jp)   

            mg_p(i,j) = OMEGA*(mg_f(i,j) - tmp)/mg_alpha(i,j,5) + (1.-OMEGA)*mg_po(i,j)

        end do
        end do

    end do

    return
    end subroutine MG_Relax
! !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    !<<<<<<<<<<<<<<< MG_Restriction >>>>>>>>>>>>>>>
    subroutine MG_Restriction(ng1, ng2, mg_DX, mg_DY, mg_p, mg_f, mg_fc, mg_alpha, level)
    use Global
    implicit none

    integer :: i, j
    integer :: ip, im, jp, jm
    integer :: i2, j2
    integer :: ng1, ng2, level

    double precision, dimension(5) :: VLM

    double precision, dimension(0:ng1+1) :: mg_DX
    double precision, dimension(0:ng2+1) :: mg_DY

    double precision, dimension(ng1, ng2, 5) :: mg_alpha
    double precision, dimension(ng1, ng2) :: mg_f, rhs_f
    double precision, dimension(ng1/2, ng2/2) :: mg_fc
    double precision, dimension(0:ng1+1, 0:ng2+1) :: mg_p

    rhs_f(1:ng1, 1:ng2) = 0.
    mg_fc(1:ng1/2, 1:ng2/2) = 0.

   
    do j = 1, ng2
    jp = j + 1
    jm = j - 1
    do i = 1, ng1
    ip = i + 1
    im = i - 1

        rhs_f(i,j) = mg_f(i,j)                                                              &
                &  - mg_alpha(i,j,1)*mg_p(im,j) - mg_alpha(i,j,2)*mg_p(ip,j)          &
                &  - mg_alpha(i,j,3)*mg_p(i,jm) - mg_alpha(i,j,4)*mg_p(i,jp)          &
                &  - mg_alpha(i,j,5)*mg_p(i,j)  
    end do
    end do


    do j = 1, ng2/2
    j2 = j*2
    do i = 1, ng1/2
    i2 = i*2
        ! M_20160830 (Misunderstand the relationship between 'Volume' and 'rhd_f')
        VLM(1) = MGDX(i2-1,level)*MGDY(j2-1,level)
        VLM(2) = MGDX(i2-1,level)*MGDY(j2  ,level)
        VLM(3) = MGDX(i2  ,level)*MGDY(j2-1,level)
        VLM(4) = MGDX(i2  ,level)*MGDY(j2  ,level)
        VLM(5) = sum(VLM(1:4))

        mg_fc(i,j) = (VLM(1)*rhs_f(i2,  j2  )                   &
                & +  VLM(2)*rhs_f(i2,   j2-1)                   & 
                & +  VLM(3)*rhs_f(i2-1, j2  )                   &
                & +  VLM(4)*rhs_f(i2-1, j2-1))/VLM(5)

    end do
    end do

    return 
    end subroutine MG_Restriction
! !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    !<<<<<<<<<<<<<<< MG_Prolongation >>>>>>>>>>>>>>>
    subroutine MG_Prolongation(ng1, ng2, mg_vc, mg_vf, level)
    use Global
    implicit none

    integer :: i, j
    integer :: ng1, ng2
    integer :: i2, j2
    integer :: level

    double precision :: tmp

    double precision, dimension(0:ng1+1, 0:ng2+1) :: mg_vc
    double precision, dimension(0:ng1*2+1, 0:ng2*2+1) :: mg_vf

    
    do j = 1, ng2
    j2 = j*2
    do i = 1, ng1
    i2 = i*2

        tmp = mg_vc(i,j)

        mg_vf(i2-1,j2-1) = tmp
        mg_vf(i2-1,j2-1) = tmp
        mg_vf(i2-1,j2  ) = tmp
        mg_vf(i2-1,j2  ) = tmp
        mg_vf(i2  ,j2-1) = tmp
        mg_vf(i2  ,j2-1) = tmp
        mg_vf(i2  ,j2  ) = tmp
        mg_vf(i2  ,j2  ) = tmp

    end do
    end do

    return
    end subroutine MG_Prolongation
 ! !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


    !<<<<<<<<<<<<<<< MG_Uniq >>>>>>>>>>>>>>>
    subroutine MG_Uniq(mg_p, ng1, ng2)
    implicit none

    integer :: i, j
    integer :: ng1, ng2

    double precision :: AVER

    double precision, dimension(0:ng1+1, 0:ng2+1) :: mg_p

    AVER = sum(mg_p(1:ng1,1:ng2))/dble(ng1)/dble(ng2)

    do j = 1, ng2
    do i = 1, ng1
        mg_p(i,j) = mg_p(i,j) - AVER
    end do
    end do

    
    return
    end subroutine MG_Uniq
! !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



!    <<<<<<<<<<<<<<<  PoissonSolver_FFTW1dy >>>>>>>>>>>>>>>    Xiaomin
     subroutine PoissonSolver_FFTW1dy(RDP, dP)
     use omp_lib
     use Global
     implicit none
     include 'fftw3.f'

     DOUBLE PRECISION, DIMENSION(0:N1,0:N2) :: dP, RDP  

     integer :: nn2 
     complex(8), parameter :: II = (0.d0, 1.d0)

     integer(8) :: plan
     integer :: ix, iy, iix, iiy, ip, jp, im, jm
     integer :: ivm, ivp


     double precision :: OMEGA, tempx(2,2), tempy(2,2)

     double precision :: ddy, temp
     double precision :: kx2, ky2
     double precision :: t_x, t_y


     double precision, dimension(1:N1m, 1:N2m) ::  fphi
     double precision, dimension(0:N1, 0:N2) :: phi

     double precision, dimension(N2m) :: tempf
     double precision, dimension(N1m) :: am,  ac,  ap , rR, rI, M_ac
     double precision, dimension(N1m) :: ad
     double precision, dimension(N1m) :: amR, acR, apR
     double precision, dimension(N1m) :: adR


     double precision, allocatable :: tempfk(:)

     double precision, allocatable  :: fk(:,:), phik(:,:)

     nn2 = N2M

     allocate(tempfk(nn2))
     allocate(fk(n1m,nn2), phik(n1m,nn2))

     fphi = RDP(1:N1M,1:N2M)


    tempx(1,1) = (dmx(1 )+dmx(2  ))**2./((dmx(1 )+dmx(2  ))**2.-dmx(1 )**2.)
    tempx(2,1) = dmx(1 )**2.           /((dmx(1 )+dmx(2  ))**2.-dmx(1 )**2.)
    tempx(1,2) = (dmx(N1)+dmx(N1m))**2./((dmx(N1)+dmx(N1m))**2.-dmx(N1)**2.)
    tempx(2,2) = dmx(N1)**2.           /((dmx(N1)+dmx(N1m))**2.-dmx(N1)**2.)

    tempy(1,1) = (dmy(1 )+dmy(2  ))**2./((dmy(1 )+dmy(2  ))**2.-dmy(1 )**2.)
    tempy(2,1) = dmy(1 )**2.           /((dmy(1 )+dmy(2  ))**2.-dmy(1 )**2.)
    tempy(1,2) = (dmy(N2)+dmy(N2m))**2./((dmy(N2)+dmy(N2m))**2.-dmy(N2)**2.)
    tempy(2,2) = dmy(N2)**2.           /((dmy(N2)+dmy(N2m))**2.-dmy(N2)**2.)

!      tempx(1,1) = 1.
!      tempx(2,1) = 0.
!      tempx(1,2) = 1.
!      tempx(2,2) = 0.

!      tempy(1,1) = 1.
!      tempy(2,1) = 0.
!      tempy(1,2) = 1.
!      tempy(2,2) = 0.



     ddy = dy(5)   ! Y-direction is uniform
     fk = 0.


                       
    ! Forward FFT
    ! compute the fourier transform of the right-hand side
     
     call dfftw_plan_r2r_1d(plan, N2m, tempf(1:N2m), tempfk(1:nn2), FFTW_REDFT10, FFTW_ESTIMATE)
   
     do ix = 1, N1m

        tempf(1:N2m) = fphi(ix, 1:N2m)
        
        call dfftw_execute_r2r(plan, tempf(1:N2m), tempfk(1:nn2))
        
        fk(ix, 1:nn2) = tempfk(1:nn2)

     end do
     call dfftw_destroy_plan(plan)


    !  STATEMENT OF THE TDMA COEFFICIENTS
     phik = 0.

     
     do ix = 1, N1m
       ip = ix + 1
       im = ix - 1
      ivm = ix - IM_V(ix)
      ivp = IP_A(ix) - ix

        am(ix) = 1./dx(ix)/dmx(ix)*dble(ivm)
        ap(ix) = 1./dx(ix)/dmx(ip)*dble(ivp)
        ac(ix) = - am(ix) - ap(ix)
     end do


     do iy = 1, nn2

        iiy = iy - 1 ! 20151005 BE CAREFUL!

        ky2 = 4./ddy/ddy*sin(0.5*(iy-1)*Pi/N2m)*sin(0.5*(iy-1)*Pi/N2m)


        do ix = 1, N1m
            M_ac(ix) = ac(ix) - ky2
        end do


        ad(1:N1m) = fk(1:N1m, iy)

        amR(1:N1m) =   am(1:N1m)
        acR(1:N1m) = M_ac(1:N1m)
        apR(1:N1m) =   ap(1:N1m)
        adR(1:N1m) =   ad(1:N1m)


        if(iiy == 0) then
            amR(1) = 0.
            acR(1) = 1.
            apR(1) = 0.
            adR(1) = 0.
        end if


        call TDMA(amR(1:N1m), acR(1:N1m), apR(1:N1m), dble(adR(1:N1m)), rR(1:N1m), N1m)

           
        phik(1:N1m, iy) = rR(1:N1m) 


     end do   


     phi = 0.

    ! Backward FFTW
    ! check the correction of fftw's calling

     call dfftw_plan_r2r_1d(plan, N2m, tempfk(1:nn2), tempf(1:N2m), FFTW_REDFT01, FFTW_ESTIMATE)

     do ix = 1, N1m

        tempfk(1:nn2) = phik(ix, 1:nn2)
        
        call dfftw_execute_r2r(plan,  tempfk(1:nn2), tempf(1:N2m) )        

        phi(ix,1:N2m) = tempf(1:N2m)/dble(N2m)/2.

     end do
     call dfftw_destroy_plan(plan)


   

     temp = sum(phi(1:N1m,1:N2m))/dble(N1m)/dble(N2m)
     phi(1:N1m,1:N2m) = phi(1:N1m,1:N2m) - temp

    ! UPDATE THE BOUNDARY OF phi 20151001
    ! x= xmin & xmax 
     phi(0, 1:N2m) = tempx(1,1)*phi(1,   1:N2m) - tempx(2,1)*phi(2,   1:N2m)
     phi(N1,1:N2m) = tempx(1,2)*phi(N1-1,1:N2m) - tempx(2,2)*phi(N1-2,1:N2m)
 
    ! y = ymin & ymax
     phi(1:N1m,0 ) = tempy(1,1)*phi(1:N1m,1   ) - tempy(2,1)*phi(1:N1m,2   )
     phi(1:N1m,N2) = tempy(1,2)*phi(1:N1m,N2-1) - tempy(2,2)*phi(1:N1m,N2-2)

    ! POINTS
     phi(0, 0) = 0.5*tempy(1,1)*phi(0,   1) - 0.5*tempy(2,1)*phi(0,   2)    &
             & + 0.5*tempx(1,1)*phi(1,   0) - 0.5*tempx(2,1)*phi(2,   0)
     phi(N1,0) = 0.5*tempx(1,2)*phi(N1-1,0) - 0.5*tempx(2,2)*phi(N1-2,0)    &
             & + 0.5*tempy(1,1)*phi(N1,  1) - 0.5*tempy(2,1)*phi(N1,  2)

     phi(0, N2) = 0.5*tempy(1,2)*phi(0,   N2-1) - 0.5*tempy(2,2)*phi(0,   N2-2)    &
              & + 0.5*tempx(1,1)*phi(1,   N2  ) - 0.5*tempx(2,1)*phi(2,   N2  )
     phi(N1,N2) = 0.5*tempx(1,2)*phi(N1-1,N2  ) - 0.5*tempx(2,2)*phi(N1-2,N2  )    &
              & + 0.5*tempy(1,2)*phi(N1,  N2-1) - 0.5*tempy(2,2)*phi(N1,  N2-2)

     DP(0:N1,0:N2) = phi(0:N1,0:N2)

     deallocate(tempfk)
     deallocate(fk, phik)

     return

     end subroutine PoissonSolver_FFTW1dy
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    !<<<<<<<<<<<<<<< streamfunction >>>>>>>>>>>>>>>
    subroutine streamfunction
    use global 
    implicit none 

    double precision :: t 

    integer :: i, j, im, jm  
    double precision :: psi(0:N1,0:N2)


    psi(0:N1,0:N2) = 0.

    open  (unit=1001,file=trim(dirname)//'STREAM.PLT')
    write(1001,*) 'zone t="',1,'"','j=',N2,'i=',N1
    !Calculate stream function
    do j = 1, N2
        jm = j - 1 
        do i = 1, N1

            psi(i,j) = psi(i,jm) + Up(i,jm)*DY(j)

            write(1001,'(2(e11.5,2x), e20.12)') X(i), Y(j), psi(i,j)

        end do 
    end do 
    close(1001)

    
    return
    end subroutine streamfunction
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    !<<<<<<<<<<<<<<< output >>>>>>>>>>>>>>>
    subroutine output(t)
    use global 
    implicit none 

    double precision :: t 

    integer :: i, j, im, jm 
    double precision :: Pup, Pdown, Eup, Edown 
    double precision :: Uu, Vv, P1, Ee, VOR 
    double precision :: DUDY, DVDX 



    open  (unit=120,file=trim(dirname)//'RESULT.PLT')  
    open  (unit=121,file=trim(dirname)//'THETA.PLT')
    

    write(120,*) 'VARIABLES=X, Y, U, V, P, THETA, VORTICITY'
    write(120,*) 'zone t="',1,'"','j=',N2,'i=',N1


    do j = 1, N2 
        jm = j - 1
        do i = 1, N1 
            im = i - 1 
                
            if(i == 1) then 
                Vv    = Vp(im,j)
                Pup   = Pp(i,j )
                Pdown = Pp(i,jm)
                Eup   = THETAp(i,j )
                Edown = THETAp(i,jm)
                DVDX = 2.0/DX(i)*(Vp(i,j) - Vp(im,j))
            else if(i == N1) then 
                Vv    = Vp(i,j)
                Pup   = Pp(i,j )
                Pdown = Pp(i,jm)
                Eup   = THETAp(i,j )
                Edown = THETAp(i,jm)
                DVDX = 2.0/DX(im)*(Vp(i,j) - Vp(im,j))
            else
                Vv    = 0.5/DMX(i)*(DX(im)*Vp(i,j ) + DX(i)*Vp(im,j ))
                Pup   = 0.5/DMX(i)*(DX(im)*Pp(i,j ) + DX(i)*Pp(im,j ))
                Pdown = 0.5/DMX(i)*(DX(im)*Pp(i,jm) + DX(i)*Pp(im,jm))
                Eup   = 0.5/DMX(i)*(DX(im)*THETAp(i,j ) + DX(i)*THETAp(im,j ))
                Edown = 0.5/DMX(i)*(DX(im)*THETAp(i,jm) + DX(i)*THETAp(im,jm))
                DVDX  = 1.0/DMX(i)*(Vp(i,j) - Vp(im,j))
            end if 

            if(j == 1) then 
                Uu = Up(i,jm)
                P1 = Pup 
                Ee = Eup 
                DUDY = 2.0/DY(j)*(Up(i,j) - Up(i,jm))
            else if(j == N2) then
                Uu = Up(i,j)
                P1 = Pdown 
                Ee = Edown 
                DUDY = 2.0/DY(jm)*(Up(i,j) - Up(i,jm))
            else
                Uu = 0.5/DMY(j)*(DY(jm)*Up(i,j) + DY(j)*Up(i,jm))
                P1 = 0.5/DMY(j)*(DY(jm)*Pup       + DY(j)*Pdown)
                Ee = 0.5/DMY(j)*(DY(jm)*Eup       + DY(j)*Edown)
                DUDY = 1.0/DMY(j)*(Up(i,j) - Up(i,jm))
            end if 

            VOR = DUDY - DVDX 

            write(120, 1200) X(i), Y(j), Uu, Vv, P1, Ee, VOR 
            write(121, *) Ee 

        end do 
    end do 

1200 FORMAT(7(3X, F22.10))


    close(120)
    close(121)


    return 


    return
    end subroutine output
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    subroutine writeup
    use global 
    implicit none 

    character*256 :: cfile 

    integer :: i, j, im, jm 
    double precision :: Pup, Pdown, Eup, Edown 
    double precision :: Uu, Vv, P1, Ee, VOR 
    double precision :: DUDY, DVDX 
    


   
    write (cfile,*) imore 
    open  (unit=112,file=trim(dirname)//'RESULT_'//trim(adjustl(cfile))//'.PLT')  

    write(112,*) 'VARIABLES=X, Y, U, V, P, THETA, VORTICITY'
    write(112,*) 'zone t="',1,'"','j=',N2,'i=',N1

    do j = 1, N2 
        jm = j - 1
        do i = 1, N1 
            im = i - 1
            
            if(i == 1) then 
                Vv    = Vp(im,j)
                Pup   = Pp(i,j )
                Pdown = Pp(i,jm)
                Eup   = THETAp(i,j )
                Edown = THETAp(i,jm)
                DVDX = 2.0/DX(i)*(Vp(i,j) - Vp(im,j))
            else if(i == N1) then 
                Vv    = Vp(i,j)
                Pup   = Pp(i,j )
                Pdown = Pp(i,jm)
                Eup   = THETAp(i,j )
                Edown = THETAp(i,jm)
                DVDX = 2.0/DX(im)*(Vp(i,j) - Vp(im,j))
            else
                Vv    = 0.5/DMX(i)*(DX(im)*Vp(i,j ) + DX(i)*Vp(im,j ))
                Pup   = 0.5/DMX(i)*(DX(im)*Pp(i,j ) + DX(i)*Pp(im,j ))
                Pdown = 0.5/DMX(i)*(DX(im)*Pp(i,jm) + DX(i)*Pp(im,jm))
                Eup   = 0.5/DMX(i)*(DX(im)*THETAp(i,j ) + DX(i)*THETAp(im,j ))
                Edown = 0.5/DMX(i)*(DX(im)*THETAp(i,jm) + DX(i)*THETAp(im,jm))
                DVDX  = 1.0/DMX(i)*(Vp(i,j) - Vp(im,j))
            end if 
                

            if(j == 1) then 
                Uu = Up(i,jm)
                P1 = Pup 
                Ee = Eup 
                DUDY = 2.0/DY(j)*(Up(i,j) - Up(i,jm))
            else if(j == N2) then
                Uu = Up(i,j)
                P1 = Pdown 
                Ee = Edown 
                DUDY = 2.0/DY(jm)*(Up(i,j) - Up(i,jm))
            else
                Uu = 0.5/DMY(j)*(DY(jm)*Up(i,j) + DY(j)*Up(i,jm))
                P1 = 0.5/DMY(j)*(DY(jm)*Pup       + DY(j)*Pdown)
                Ee = 0.5/DMY(j)*(DY(jm)*Eup       + DY(j)*Edown)
                DUDY = 1.0/DMY(j)*(Up(i,j) - Up(i,jm))
            end if 

            VOR = DUDY - DVDX 

            write(112, 1120) X(i), Y(j), Uu, Vv, P1, Ee, VOR 

        end do 
    end do 

1120 FORMAT(7(3X, F22.10))


    close(112)


    return 

    end subroutine writeup 
