! Module containing an example subroutine that calculates the error between 2 vectors

! Copyright (C) 2002 Jon May

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.
!===========================================================================

module sims

implicit none

CONTAINS


subroutine Example(Models,ValMax,Misfits,RptMod,Debug)

implicit none

real*8,dimension(:,:),intent(in) :: Models
real*8,dimension(:),intent(in) :: ValMax
real*8,dimension(:),intent(out) :: Misfits
integer,intent(in),optional :: RptMod
logical,intent(in),optional :: Debug

real*8,dimension(:),allocatable,save :: Aim
logical,save :: called=.False.
integer :: i


! Assign values to Aim
if(.Not. called) then
  allocate(Aim(size(Models,2)))

! Specific
  ! Aim = (/8.5,7.5/)
  Aim = (/8.5,7.5,7.5/)
  ! Aim = (/9.0,3.0,3.0/)

! Random
  ! call Random_number(Aim)
  ! do i =1,size(Aim) ! Adjust Aim as Models are range 0-ValMax(i)
    ! Aim(i) = Aim(i) * ValMax(i)
  ! end do

  if(present(Debug) .and. (debug)) print*,'Aim:',Aim
called = .True.
end if

do i = 1,size(Models,1)
  if(present(RptMod) .and. i==RptMod) then
    if(Debug) print'(A,I0,A,I0)',' Misfit subroutine skipping model ',RptMod,' / ',size(Models,1)
    cycle
  end if
  Misfits(i) = (sum((Models(i,:)-Aim(:))**2))**0.5
end do

end subroutine


end module
