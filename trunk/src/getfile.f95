
MODULE mod_getfile
  USE mod_param
  USE netcdf
  
  INTEGER, DIMENSION(1)                      :: start1D  ,count1D
  INTEGER, DIMENSION(4)                      :: start2D  ,count2D ,map2D
  INTEGER, DIMENSION(4)                      :: start3D  ,count3D ,map3D
  INTEGER, DIMENSION(4)                      :: start4D  ,count4D ,map4D
  INTEGER                                    :: ncTpos

  INTEGER                                    :: ierr, varid,ncid
  
CONTAINS

  !===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

  function get1DfieldNC (fieldFile ,varName)
    CHARACTER (len=*)                       :: fieldFile ,varName 
    REAL, ALLOCATABLE,   DIMENSION(:)       :: get1dfieldNC
    INTEGER,             DIMENSION(1)       :: d    
    INTEGER                                 :: varid ,ncid
  
    d=count1d(1)+start1d(1)-1
    allocate ( get1DfieldNC(d(1)) )

    ierr=NF90_OPEN(trim(fieldFile) ,NF90_NOWRITE ,ncid)
    if(ierr.ne.0) call printReadError(1)    
    ierr=NF90_INQ_VARID(ncid ,varName ,varid)
    if(ierr.ne.0) call printReadError(2)
    ierr=NF90_GET_VAR(ncid ,varid ,get1DfieldNC ,start1d ,count1d)
    if(ierr.ne.0) call printReadError(3)
    ierr=NF90_CLOSE(ncid)
  end function get1dfieldNC

  !===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

  function get2DfieldNC (fieldFile ,varName)
    CHARACTER (len=*)                       :: fieldFile ,varName 
    REAL, ALLOCATABLE,   DIMENSION(:,:)     :: get2DfieldNC
    REAL, ALLOCATABLE,   DIMENSION(:,:)     :: field
    INTEGER,             DIMENSION(4)       :: d, s, c, dimids
    INTEGER                                 :: ncid
    
    s = start2d(map2d)
    c = count2d(map2d)
    d = start2d + count2d - 1    
    allocate ( field(d(1),d(2)), get2dfieldNC(imt+2,jmt) )
    
    ierr=NF90_OPEN(trim(fieldFile) ,NF90_NOWRITE ,ncid)
    if(ierr.ne.0) call printReadError(1)
    ierr=NF90_INQ_VARID(ncid ,varName ,varid)
    if(ierr.ne.0) call printReadError(2)
    ierr=NF90_GET_VAR(ncid ,varid , field, start2d, count2d)
    if(ierr.ne.0) call printREadError(3)
    r=NF90_inquire_variable(ncid, varid, dimids = dimids)   
    ierr=NF90_CLOSE(ncid)

    if ( map2d(2) > map2d(1) ) then
       get2DfieldNC(1:imt,:) = field
    else
       do i = 1,imt
          get2DfieldNC(i,:) = field(:,i)
       end do
    end if

  end function get2DfieldNC
  
  !===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===

  function get3DfieldNC (fieldFile ,varName)
    CHARACTER (len=*)                       :: fieldFile ,varName 
    REAL, ALLOCATABLE, DIMENSION(:,:,:)     :: field
    REAL, ALLOCATABLE, DIMENSION(:,:,:)     :: get3dfieldNC
    INTEGER,             DIMENSION(4)       :: d, s, c
    INTEGER                                 :: i,j,k

    start3d(1) = ncTpos
    s = start3d(map3d)
    c = count3d(map3d)
    d = c + s - 1
    allocate ( field(d(1),d(2),d(3)), get3dfieldNC(imt+2,jmt,km) )

    ierr = NF90_OPEN(trim(fieldFile) ,NF90_NOWRITE ,ncid)
    if(ierr.ne.0) call printReadError(1)
    ierr=NF90_INQ_VARID(ncid ,varName ,varid)
    if(ierr.ne.0) call printReadError(2)
    ierr=NF90_GET_VAR(ncid ,varid ,field, s, c)
    if(ierr.ne.0) call printReadError(3)

    if  (all(map3d == (/1,2,3,4/),DIM=1) ) then
       get3DfieldNC = field
    elseif (all(map3d == (/2,3,4,1/),DIM=1) )  then
       get3DfieldNC = field
    elseif (all(map3d == (/3,2,4,1/),DIM=1) )  then
       do k=1,km
          do i=1,imt
             get3DfieldNC(i,:,k) = field(:,i,k)
          end do
       end do
    elseif (all(map3d == (/3,2,4,1/),DIM=1) )  then
       do k=1,km
          do i=1,imt
             get3DfieldNC(i,:,k) = field(:,i,k)
          end do
       end do
    else
       print *,"ERROR!"
       print *,"==================================================="
       print *," This combination of dimensions in the indata file"
       print *," isn't implemented yet. Please contact"
       print *," Bror Jonsson (brorfred@gmail.com) for help."
       stop 999
    end if
    ierr=NF90_CLOSE(ncid)

  end function get3DfieldNC


  subroutine printReadError(sel)
    USE netcdf
    INTEGER                                 :: sel, ndims, v,dimln
    INTEGER, DIMENSION(4)                   :: dimvec
    CHARACTER(len=20)                       :: dimname
    select case (sel)
    case (1)
       print *,'Error when trying to open the file'
       print *,'   ' ,getfile
       print *,'    Error: ' , NF90_STRERROR(ierr)
       stop 
    case (2)
       print *,'Error when trying to open the field   ',ncvar
       print *,'    Error: ' , NF90_STRERROR(ierr)
       stop 
    case (3)

       print *,'Error when trying to read the field   ',ncvar
       print *, 'start1d =  ' ,start1d
       print *, 'count1d =  ' ,count1d
       print *, 'start2d =  ' ,start2d
       print *, 'count2d =  ' ,count2d
       print *, 'start3d =  ' ,start3d(map3d)
       print *, 'count3d =  ' ,count3d(map3d)
       print *, '-------------------------------------------------'
       print *,'Error:      ' ,NF90_STRERROR(ierr)
      
       ierr= nf90_inquire_variable(ncid,varid,ndims=ndims,dimids=dimvec)
       do v = 1,ndims
          ierr = nf90_inquire_dimension(ncid, dimvec(v), dimname, dimln)
          print *,dimname,dimln
       end do
       stop
       
       !r=NF90_inquire_variable(ncid, varid, dimids = dimids)
       !do i=1,4
       !   r=NF90_inquire_dimension(ncid, dimids(i), len=d(i))
       !end do
       !print * ,'Error when trying to read the field   ',varName
       !print * ,'start2d =  ' ,start2d
       !print * ,'count2d =  ' ,count2d
       !print * ,'Dimensions: ' ,d
       !print * ,'Error:      ' ,NF90_STRERROR(ierr)
       !stop

    end select
  end subroutine printReadError
end MODULE mod_getfile





!###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###     
!  function get4dfield ()
!    REAL, ALLOCATABLE,   DIMENSION(:,:,:,:) :: get4dfield
!    INTEGER,             DIMENSION(4)       :: d
!    INTEGER                                 :: varid ,ncid
!  !===   ===   ===   ===   ===   ===   ===   ===   ===   ===   ===   === 
!    d=count4d+start4d-1
!    allocate ( get4dfield(d(1),d(2),d(3),d(4)) )
!
!    ierr=NF90_OPEN(trim(getfile) ,NF90_NOWRITE ,ncid)
!    fileError: if(ierr.ne.0) then
!       print *,'Error when trying to open the file'
!       print *,'   ' ,getfile
!       print *,'    Error code: ' , ierr
!       stop 3001
!    end if fileError
!    ierr=NF90_INQ_VARID(ncid ,varName ,varid)
!    varError: if(ierr.ne.0) then
!       print *,'Error when trying to find the field   ',varName
!       print *,'Error code: ' , ierr
!       stop 3001
!    end if varError
!    ierr=NF90_GET_VAR(ncid ,varid ,get4dfield ,start1d ,count1d)
!    if(ierr.ne.0) stop 3100
!    ierr=NF90_CLOSE(ncid)
!  end function get4dfield