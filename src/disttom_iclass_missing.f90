subroutine realmin(ain,size_ain,icls,themin) 
  ! outputs 
  ! icls: the position of the vector containing the min
  ! the min: the min value

  implicit none
  integer :: i,size_ain,icls 
  double precision, dimension(size_ain) :: ain
  double precision :: themin
  icls=1
  themin=ain(1)
  if1: if (size_ain == 1) then
     themin=ain(1)
  else 
  do i=2,size_ain
     if2: if (themin > ain(i)) then
          themin=ain(i)
          icls=i
     end if if2
  end do
  end if if1
end subroutine realmin
!
!
subroutine get_index(index01,n,n1,out_index)
!
  implicit none
  ! Return the positions of the vector
  ! which contains 1 (1 = nonmissing)

  ! length of index01 and index
  integer :: n,n1 
  ! vector containing 0 or 1
  integer, dimension(n) :: index01 
  ! output
  integer,dimension(n1) :: out_index 
  ! dummies
  integer :: i,count 
  !
  count = 1
  do i = 1,n
     if1: if (index01(i)==1) then
        out_index(count) = i
        count = count + 1
     end if if1 
  end do
end subroutine get_index


subroutine disttom_iclass_missing(data,miss_data,nrow,ncol,mu,miss_mu,k,iclass,disttom,W,sumW)

  integer :: k,j,l,nrow,ncol,icls,n_miss_vec
  double precision :: themin, sumW, unadj_d
  double precision, dimension(k) :: dj
  double precision, dimension(nrow,ncol) :: data
  double precision, dimension(k,ncol) :: mu

  ! weights for scaling the distances 
  double precision, dimension(ncol) :: W,unsummed_d
  ! out_index records the positions of the distance vector
  ! which contains nonmissing values. 
  ! Its length depends on the distance vector of interest.
  integer, dimension(:), allocatable :: out_index
  integer, dimension(ncol) :: p_miss_vec, jth_miss_d, kth_miss_mu

  ! miss_data(i,j) = 1 if the (i,j) entry is non-missing in data 
  integer, dimension(nrow,ncol) :: miss_data
  ! miss_mu(i,j) = 1 if the (i,j) entry is non-missing in mu
  integer, dimension(k,ncol) :: miss_mu

  ! outputs of interests
  ! cluster labels
  integer, dimension(nrow) :: iclass
  ! The distances of cases to their closest cluster centers
  double precision, dimension(nrow) :: disttom

  doj: do j=1,nrow
     !
     jth_miss_d = miss_data(j,:)
     !
     dol1: do l =1,k
        kth_miss_mu = miss_mu(l,:)
        !
        ! p_miss_vec(i,j) = 1 iff jth_miss_d(i,j) and kth_miss_mu(i,j)
        ! are both non-missing
        p_miss_vec = jth_miss_d*kth_miss_mu
        ! n_miss_vec = # NON-missing entries in the distance vector
        n_miss_vec = sum(p_miss_vec)
        !
        allocate( out_index(n_miss_vec))
        call get_index(p_miss_vec,ncol,n_miss_vec,out_index)
        !
        unsummed_d = (data(j,:)-mu(l,:) )**2
        unadj_d = sum( unsummed_d(out_index) )
        dj(l)=unadj_d*sumW/sum(W(out_index))
        !
        deallocate(out_index)
     end do dol1
     !
     call realmin(dj,k,icls,themin)
     disttom(j)=themin
     iclass(j)=icls
  end do doj
end subroutine disttom_iclass_missing

subroutine disttom_iclass(data,nrow,ncol,mu,k,iclass,disttom)

  integer :: k,j,l,nrow,ncol,icls
  double precision :: themin
  double precision, dimension(k) :: dj
  double precision, dimension(nrow,ncol) ::data
  double precision, dimension(k,nrow) :: mu

  ! outputs of interests

  integer, dimension(nrow) :: iclass
  double precision, dimension(nrow) :: disttom

  doj: do j=1,nrow

     dj=(/(0,i=1,k)/)

     dol1: do l =1,k
        dj(l)=sum( (data(j,:)-mu(l,:) )**2)
     end do dol1

     call realmin(dj,k,icls,themin)
     disttom(j)=themin
     iclass(j)=icls

  end do doj

end subroutine disttom_iclass
