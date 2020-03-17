module FVW_Tests 

   use NWTC_Library

   use FVW_Types
   use FVW_Subs
   use FVW_VortexTools
   use FVW_Wings
   use FVW_IO
   use FVW_BiotSavart

   implicit none

   public :: FVW_RunTests
   private

   character(len=255),save :: testname
   interface test_equal; module procedure &
         test_equal_i1, &
         test_equal_i0
   end interface
    interface test_almost_equal; module procedure &
          test_almost_equal_0, &
          test_almost_equal_1, &
          test_almost_equal_2
    end interface
contains
   ! --------------------------------------------------------------------------------
   ! --- Helper functions (should be part of NWTC library)
   ! --------------------------------------------------------------------------------
    subroutine test_success(info,bPrint_in)
        character(len=*), intent(in) :: info
        logical, intent(in), optional  ::  bPrint_in
        if(present(bPrint_in)) then
            if(bPrint_in) then
                write(*,'(A)')'[ OK ] '//trim(testname)//': '//trim(Info)
            endif
        else
            write(*,'(A)')'[ OK ] '//trim(testname)//': '//trim(Info)
        endif
    end subroutine

    subroutine test_fail(info,bPrint_in,bStop_in)
        character(len=*), intent(in) :: info
        logical, intent(in), optional  ::  bPrint_in
        logical, intent(in), optional  ::  bStop_in
        if(present(bPrint_in)) then
            if(bPrint_in) then
                write(*,'(A)')'[FAIL] '//trim(testname)//': '//trim(Info)
            endif
        else
            write(*,'(A)')'[FAIL] '//trim(testname)//': '//trim(Info)
        endif
        if(present(bStop_in)) then
            if(bStop_in) then
                STOP 
            endif
        else
            STOP
        endif
    end subroutine

    subroutine  test_equal_i0(Var,iTry,iRef)
        ! Arguments
        character(len=*), intent(in) :: Var
        integer, intent(in) :: iTry         !< 
        integer, intent(in) :: iRef         !< 
        ! Variables
        character(len=255) :: InfoAbs
        if(iRef/=iTry) then
            write(InfoAbs,'(A,I0,A,I0)') trim(Var),iRef,'/',iTry
            call test_fail(InfoAbs)
            STOP -1 !OTHER-COMPILER
            STOP ! COMPAQ-COMPILER
        else
            write(InfoAbs,'(A,A,I0)') trim(Var),' ok ',iRef
            call test_success(InfoAbs)
        endif
    end subroutine

    subroutine  test_equal_i1(Var,VecTry,VecRef,bTest,bPrintOnly,bPassed)
        ! Arguments
        character(len=*), intent(in) :: Var
        integer, dimension(:), intent(in) :: VecTry         !< 
        integer, dimension(:), intent(in) :: VecRef         !< 
        logical, intent(in) :: bTest
        logical, intent(in) :: bPrintOnly
        logical, intent(out),optional :: bPassed
        ! Variables
        character(len=255) :: InfoAbs
        integer :: i,cpt
        ! 
        cpt=0
        do i=1,size(VecRef)
            if(VecRef(i)/=VecTry(i)) then
                cpt=cpt+1
            endif
        enddo
        if(cpt>0) then
            write(InfoAbs,'(A,I0)') trim(Var)//' Elements different: ',cpt
            if(present(bPassed)) then
                bPassed=.false.
            endif
        else
            write(InfoAbs,'(A)') trim(Var)//' reproduced to identity'
            if(present(bPassed)) then
                bPassed=.true.
            endif
        endif
        if(bPrintOnly) then
            print'(A)',trim(InfoAbs)
        endif
        if(bTest) then
            if(cpt>0) then
                call test_fail(InfoAbs)
                STOP 
            else
                call test_success(InfoAbs)
            endif
        endif
    end subroutine

    subroutine  test_almost_equal_0(Var,Ref,Try,MINNORM,bStop,bPrint,bPassed)
        ! Arguments
        character(len=*), intent(in) :: Var
        real(ReKi), intent(in) :: Ref         !< 
        real(ReKi), intent(in) :: Try         !< 
        real(ReKi), intent(in) :: MINNORM
        logical, intent(in) :: bStop
        logical, intent(in) :: bPrint
        logical, intent(out),optional :: bPassed
        ! Variables
        character(len=255) :: InfoAbs
        real(ReKi) :: delta
        integer :: cpt
        ! 
        cpt=0
        delta=abs(Ref-Try)
        if(delta>MINNORM) then
            write(InfoAbs,'(A,ES8.1E2,A,ES8.1E2,A,I0)') trim(Var)//' tol: ',MINNORM,', mean: ',delta,' - Failed:',cpt
            call test_fail(InfoAbs,bPrint,bStop)
        else
            write(InfoAbs,'(A,ES8.1E2,A,ES8.1E2)') trim(Var)//' tol: ',MINNORM,', mean: ',delta
            call test_success(InfoAbs,bPrint)
        endif
        if(present(bPassed)) then
            bPassed=delta>MINNORM
        endif
    end subroutine
    subroutine  test_almost_equal_1(Var,VecRef,VecTry,MINNORM,bStop,bPrint,bPassed)
        ! Arguments
        character(len=*), intent(in) :: Var
        real(ReKi), dimension(:), intent(in) :: VecRef         !< 
        real(ReKi), dimension(:), intent(in) :: VecTry         !< 
        real(ReKi), intent(in) :: MINNORM
        logical, intent(in) :: bStop
        logical, intent(in) :: bPrint
        logical, intent(out),optional :: bPassed
        ! Variables
        character(len=255) :: InfoAbs
        integer :: i,cpt
        real(ReKi) :: delta
        real(ReKi) :: delta_cum
        ! 
        cpt=0
        delta_cum=0.0_ReKi
        do i=1,size(VecRef,1)
            delta=abs(VecRef(i)-VecTry(i))
            delta_cum=delta_cum+delta
            if(delta>MINNORM) then
                cpt=cpt+1
            endif
        enddo
        delta_cum=delta_cum/size(VecRef)

        if(cpt>0) then
            write(InfoAbs,'(A,ES8.1E2,A,ES8.1E2,A,I0)') trim(Var)//' tol: ',MINNORM,', mean: ',delta_cum,' - Failed:',cpt
            call test_fail(InfoAbs,bPrint,bStop)
        else
            write(InfoAbs,'(A,ES8.1E2,A,ES8.1E2)') trim(Var)//' tol: ',MINNORM,', mean: ',delta_cum
            call test_success(InfoAbs,bPrint)
        endif
        if(present(bPassed)) then
            bPassed=(cpt==0)
        endif
    end subroutine
    subroutine  test_almost_equal_2(Var,VecRef,VecTry,MINNORM,bStop,bPrint,bPassed)
        ! Arguments
        character(len=*), intent(in) :: Var
        real(ReKi), dimension(:,:), intent(in) :: VecRef         !< 
        real(ReKi), dimension(:,:), intent(in) :: VecTry         !< 
        real(ReKi), intent(in) :: MINNORM
        logical, intent(in) :: bStop
        logical, intent(in) :: bPrint
        logical, intent(out),optional :: bPassed
        ! Variables
        real(ReKi), dimension(:),allocatable :: VecRef2    !< 
        real(ReKi), dimension(:),allocatable :: VecTry2   !<
        integer :: p, i,j,n1,n2,nCPs
        ! 
        n1 = size(VecRef,1); n2 = size(VecRef,2); nCPs=n1*n2
        allocate ( VecRef2 (n1*n2)  ) ; allocate ( VecTry2 (n1*n2)  ) 
        p=0
        do j=1,n2; do i=1,n1
            p=p+1
            VecRef2(p)=VecRef(i,j)
            VecTry2(p)=VecTry(i,j)
        enddo; enddo;
        call  test_almost_equal(Var,VecRef2,VecTry2,MINNORM,bStop,bPrint,bPassed)
    end subroutine

   ! --------------------------------------------------------------------------------}
   ! --- Specific FVW tests 
   ! --------------------------------------------------------------------------------{
   !>
   subroutine Test_BiotSavart_Sgmt(ErrStat, ErrMsg)
      integer(IntKi)      , intent(out) :: ErrStat !< Error status of the operation
      character(ErrMsgLen), intent(out) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
      integer(IntKi)       :: ErrStat2
      character(ErrMsgLen) :: ErrMsg2
      ! 
      real(ReKi), dimension(3) :: P1,P2,P3,CP
      real(ReKi), dimension(3) :: U1,U2
      real(ReKi) :: SegGamma1 !< Circulation  [m^2/s]
      real(ReKi) :: RegParam1 !< 
      integer(IntKi) :: i1,i2
      integer(IntKi) :: RegFunction 
      integer(IntKi), parameter :: nSegTot  = 2
      integer(IntKi), parameter :: nSegPTot = 3
      integer(IntKi), parameter :: nCPsTot  = 1
      real(ReKi),     dimension(3,nCPsTot) :: CPs            !< Control points
      real(ReKi),     dimension(3,nSegPTot) :: SegPoints     !< Segment points
      integer(IntKi), dimension(2,nSegTot) :: SegConnct      !< Connectivity, indices of segments points iSeg1, iSeg2
      real(ReKi),     dimension(nSegTot)   :: SegGamma       !< Segment circulation
      real(ReKi),     dimension(nSegTot)   :: RegParam       !< Regularization parameter
      real(ReKi),     dimension(3,nCPsTot) :: Uind_out       !< Induced velocity vector - Side effects!!!
      real(ReKi),     dimension(3,4) :: CPs_test   !< 
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""
      ! --- Test that the two functions return the same values 
      P1=(/0.  ,0.,-1./)
      P2=(/0.  ,0., 1./)
      CPs_test(:,1) = (/ 0.0,  0., 0.0  /) ! Middle
      CPs_test(:,2) = P1                   ! Extremity
      CPs_test(:,3) = (/ 0.05, 0., -0.5 /) ! Close
      CPs_test(:,4) = (/ 10.,  0., 0.0  /) ! Far
      do i2 = 1, size(CPs_test,2)
         ! Segment param
         CP=CPs_test(:,i2)
         SegGamma1=1
         RegParam1=0.5
         ! One segment param
         SegConnct(:,1)=(/1,2/)
         SegPoints(:,1) = P1
         SegPoints(:,2) = P2
         SegGamma(:) = SegGamma1
         RegParam(:) = RegParam1
         CPs (:,1)   = CP
         do i1=1,5
            RegFunction = idRegVALID(i1)
            ! Method 1
            Uind_out =0.0_ReKi
            call ui_seg(1, 1, nCPsTot, CPs, &
                  1, 1, nSegTot, nSegPTot, SegPoints, SegConnct, SegGamma,  &
                  RegFunction, RegParam, Uind_out)
            ! Method 2
            call ui_seg_11(CP-P1, CP-P2, SegGamma1, RegFunction, RegParam1, U1)
            ! Test
            print*,'Reg function', RegFunction, 'CP',CP
            print*,'Uind_out',Uind_out
            print*,'U1      ',U1
            call test_almost_equal('Uind method1/2', U1, Uind_out(:,1), 1e-4, .true.,.true.)
            !call test_almost_equal('Uind method1/2', U1, Uind_out(:,1), 1e-4, .false.,.true.)
         enddo
      enddo

      ! --- Test that the two segments or one segment returns the same value
      P1=(/0.  ,0.,-1./)
      P2=(/0.  ,0., 1./)
      P3=(/0.  ,0., 0./)
      CPs_test(:,1) = (/ 0.0,  0., 0.0  /) ! Middle
      CPs_test(:,2) = P1                   ! Extremity
      CPs_test(:,3) = (/ 0.05, 0., -0.5 /) ! Close
      CPs_test(:,4) = (/ 100.,  0., -0.5  /) ! Far
      do i2 = 1,size(CPs_test,2)
         ! Segment param
         CP=CPs_test(:,i2)
         SegGamma1=1
         RegParam1=0.5
         ! One segment param
         SegConnct(:,1)=(/1,2/)
         SegConnct(:,2)=(/2,3/)
         SegPoints(:,1) = P1
         SegPoints(:,2) = P3
         SegPoints(:,3) = P2
         SegGamma(:) = SegGamma1
         RegParam(:) = RegParam1
         CPs (:,1)   = CP
         do i1=1,4 ! NOTE stopping at 4 since Offset is not linear
            RegFunction = idRegVALID(i1)
            ! Method 1
            Uind_out =0.0_ReKi
            call ui_seg(1, 1, nCPsTot, CPs, &
                  1, 2, nSegTot, nSegPTot, SegPoints, SegConnct, SegGamma,  &
                  RegFunction, RegParam, Uind_out)
            ! Method 2
            call ui_seg_11(CP-P1, CP-P2, SegGamma1, RegFunction, RegParam1, U1)
            !print*,'Reg function', RegFunction, 'CP',CP
            !print*,'Uind_out',Uind_out
            !print*,'U1      ',U1
            call test_almost_equal('Uind 1seg/2seg', U1, Uind_out(:,1), 1e-4, .true.,.true.)
         enddo
      enddo
   end subroutine

   !>
   subroutine Test_LatticeToSegment(iStat)
      integer(IntKi), intent(  out)  :: iStat !< Status for test
      ! Local
      integer(IntKi),dimension(:,:), allocatable :: SegConnct !< Segment connectivity
      real(ReKi),    dimension(:,:), allocatable :: SegPoints !< Segment Points
      real(ReKi),    dimension(:)  , allocatable :: SegGamma  !< Segment Circulation
      real(ReKi),    dimension(:),   allocatable :: SegEpsilon !< 
      !
      real(ReKi),    dimension(:,:,:), allocatable :: LatticePoints1 !< Lattice Points
      real(ReKi),    dimension(:,:,:), allocatable :: LatticePoints2 !< Lattice Points
      real(ReKi),    dimension(:,:),   allocatable :: LatticeGamma1  !< Lattice Circulation
      real(ReKi),    dimension(:,:),   allocatable :: LatticeGamma2  !< Lattice Circulation
      real(ReKi),    dimension(:,:),   allocatable :: CPs   !< ControlPoints
      real(ReKi),    dimension(:,:),   allocatable :: Uind  !< Induced velocity
      integer(IntKi) :: iHeadC
      integer(IntKi) :: iHeadP
      integer(IntKi) :: i,j,k
      integer(IntKi) :: nP
      integer(IntKi) :: nC
      integer(IntKi) :: nP1, nP2
      integer(IntKi) :: nC1, nC2
      integer(IntKi) :: nDepth, nSpan
      integer(IntKi) :: SmoothModel

      ! --- Creating two lattice
      allocate(LatticePoints1(3,2,2)) 
      allocate(LatticePoints2(3,3,4)) 
      allocate(LatticeGamma1(1,1)) ; 
      allocate(LatticeGamma2(2,3)) ; 
      LatticeGamma1=1
      ! Test shed vorticity
      LatticeGamma2(:,1)=1
      LatticeGamma2(:,2)=2
      LatticeGamma2(:,3)=3
      ! Test trailed vorticity
!       LatticeGamma2(1,:)=1
!       LatticeGamma2(2,:)=2
      CALL MeshMe(LatticePoints1,(/0.,0.,0./))
      CALL MeshMe(LatticePoints2,(/0.,0.,1./))

      CALL WrVTK_Lattice('Points1.vtk',LatticePoints1, LatticeGamma1)
      CALL WrVTK_Lattice('Points2.vtk',LatticePoints2, LatticeGamma2)

      ! --- Convert lattice 1 to segments
      nSpan  = size(LatticePoints1,2)
      nDepth = size(LatticePoints1,3)
      nP1 = nSpan*nDepth
      nC1 = 2*(nSpan*nDepth)-nSpan-nDepth
      allocate(SegConnct(1:2,1:nC1)); SegConnct=-1
      allocate(SegPoints(1:3,1:nP1)); SegPoints=-1
      allocate(SegGamma (1:nC1)    ); SegGamma=-999
      allocate(SegEpsilon(1:nC1)    ); SegEpsilon=0.0_ReKi

      iHeadP=1
      iHeadC=1
      CALL LatticeToSegments(LatticePoints1, LatticeGamma1, 1, SegPoints, SegConnct, SegGamma, iHeadP, iHeadC, .true. )
      CALL printall()
      CALL WrVTK_Segments('Points1_seg.vtk', SegPoints, SegConnct, SegGamma, SegEpsilon) 

      allocate(Uind(1:3,1) ); Uind=0.0_ReKi
      allocate(CPs (1:3,1) ); 
      CPs(1:3,1)=(/1.5,1.5,0./)
      SegEpsilon=100.0_ReKi
      SmoothModel=0 ! No smooth
      CALL ui_seg(1, 1, 1, CPs, &
      1, nC1, nC1, nP1, SegPoints, SegConnct, SegGamma,   &
      SmoothModel, SegEpsilon, Uind)
      print*,'Uind',Uind

      ! --- Convert lattice 2 to segments
      nSpan  = size(LatticePoints2,2)
      nDepth = size(LatticePoints2,3)
      nP2 = nSpan*nDepth
      nC2 = 2*(nSpan*nDepth)-nSpan-nDepth
      deallocate(SegConnct)
      deallocate(SegPoints)
      deallocate(SegGamma)
      allocate(SegConnct(1:2,1:nC2)); SegConnct=-1
      allocate(SegPoints(1:3,1:nP2)); SegPoints=-1
      allocate(SegGamma (1:nC2)    ); SegGamma=-9999
      iHeadP=1
      iHeadC=1
      CALL LatticeToSegments(LatticePoints2, LatticeGamma2, 1, SegPoints, SegConnct, SegGamma, iHeadP, iHeadC , .true.)
      CALL printall()
      CALL WrVTK_Segments('Points2_seg.vtk', SegPoints, SegConnct, SegGamma, SegEpsilon) 

      ! --- Concatenate both
      nP = nP1 + nP2
      nC = nC1 + nC2
      iHeadP=1
      iHeadC=1
      deallocate(SegConnct)
      deallocate(SegPoints)
      deallocate(SegGamma)
      allocate(SegConnct(1:2,1:nC)); SegConnct=-1
      allocate(SegPoints(1:3,1:nP)); SegPoints=-1
      allocate(SegGamma (1:nC)    ); SegGamma=-9999
      CALL LatticeToSegments(LatticePoints1, LatticeGamma1, 1, SegPoints, SegConnct, SegGamma, iHeadP, iHeadC, .true. )
      CALL LatticeToSegments(LatticePoints2, LatticeGamma2, 1, SegPoints, SegConnct, SegGamma, iHeadP, iHeadC, .true. )
      CALL printall()
      CALL WrVTK_Segments('PointsBoth_seg.vtk', SegPoints, SegConnct, SegGamma, SegEpsilon) 


   contains
      subroutine printall()
         print*,'Points'
         do i=1,size(SegPoints,2)
            print*,'i',i,'Coords:', SegPoints(1:3,i)
         enddo
         print*,'Connectivity'
         do i=1,size(SegConnct,2)
            print*,'i',i,'Conn:', SegConnct(1:2,i),'Gam:', SegGamma(i)
         enddo
         print*,'-----------------------------'
      endsubroutine

      subroutine MeshMe(M,offset)
         real(ReKi), dimension(:,:,:), intent(inout) :: M
         real(ReKi), dimension(3)    , intent(in   ):: offset
         do j=1,size(M,3)
            do i=1,size(M,2)
               M(1,i,j)=i + offset(1)
               M(2,i,j)=j + offset(2)
               M(3,i,j)=0 + offset(3)
            enddo
         enddo 
      endsubroutine 
   endsubroutine Test_LatticeToSegment

   !> Main test function 
   subroutine FVW_RunTests(ErrStat,ErrMsg)
      integer(IntKi)      , intent(out) :: ErrStat !< Error status of the operation
      character(ErrMsgLen), intent(out) :: ErrMsg  !< Error message if ErrStat /= ErrID_None
      integer(IntKi)       :: ErrStat2
      character(ErrMsgLen) :: ErrMsg2
      ! Initialize ErrStat
      ErrStat = ErrID_None
      ErrMsg  = ""
      call Test_BiotSavart_Sgmt(ErrStat2, ErrMsg2)
   end subroutine FVW_RunTests

end module FVW_Tests
