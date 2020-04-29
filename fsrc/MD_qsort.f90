MODULE srtmod

   TYPE ipair

      INTEGER :: v     ! sort key

      INTEGER :: a     ! associated data

   END TYPE

END MODULE srtmod
subroutine Merge(A,NA,B,NB,C,NC)

use srtmod

integer, intent(in) :: NA,NB,NC      ! Normal usage: NA+NB = NC

type(ipair), intent(in out) :: A(NA) ! B overlays C(NA+1:NC)

type(ipair), intent(in) :: B(NB)

type(ipair), intent(in out) :: C(NC)
integer :: I,J,K
I = 1; J = 1; K = 1;

do while(I <= NA .and. J <= NB)

   if (A(I)%v <= B(J)%v) then

      C(K) = A(I)

      I = I+1

   else

      C(K) = B(J)

      J = J+1

   endif

   K = K + 1

enddo
do while (I <= NA)

   C(K) = A(I)

   I = I + 1

   K = K + 1

enddo
return

end subroutine merge
recursive subroutine MergeSort(A,N,T)

use srtmod

integer, intent(in) :: N

type(ipair), dimension(N), intent(in out) :: A

type(ipair), dimension((N+1)/2), intent (out) :: T
integer :: NA,NB

type(ipair) :: VT
if (N < 2) return

if (N == 2) then

   if (A(1)%v > A(2)%v) then

      VT  = A(1)

      A(1) = A(2)

      A(2) = VT

   endif

   return

endif
NA=(N+1)/2

NB=N-NA

call MergeSort(A,NA,T)

call MergeSort(A(NA+1),NB,T)
if (A(NA)%v > A(NA+1)%v) then

   T(1:NA)=A(1:NA)

   call Merge(T,NA,A(NA+1),NB,A,N)

endif
return

end subroutine MergeSort
program TestMergeSort

use srtmod

integer :: i,N

type(ipair), dimension(:), allocatable :: A , T
read (*,*) N

allocate (A(N),T((N+1)/2))

read(*,*)(A(i)%v,A(i)%a,i=1,N)
call MergeSort(A,N,T)
write(*,'(A)')' Sorted Data:'

write(*,'(/,(2I2))')(A(i)%v,A(i)%a,i=1,N)

end program TestMergeSort
