! Module containing Grid Search subroutines

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

module subs

implicit none

CONTAINS


subroutine CalcGrid(VarValsMin,VarValsMax,NumMods,Grid,Step)

real*8, dimension(:),intent(in) :: VarValsMin
real*8, dimension(:),intent(in) :: VarValsMax
integer,dimension(:),intent(in) :: NumMods

real*8,dimension(:,:),allocatable,intent(out) :: Grid ! 2D grid of vectors, max(NumMods)*NumVariables [each column is a variable to keep column access]
real*8,dimension(:),allocatable,intent(out) :: Step

integer :: i,j


allocate(Grid(maxval(NumMods),size(NumMods)),Step(size(NumMods)))

do i = 1,size(NumMods)
  Step(i) = (VarValsMax(i) - VarValsMin(i) )/(NumMods(i))
  Grid(1,i) = VarValsMin(i) + (Step(i)/2.0)
  Grid(NumMods(i),i) = VarValsMax(i) - (Step(i)/2.0)
    do j = 2,NumMods(i)-1
      Grid(j,i) = (VarValsMin(i)+(Step(i)/2.0)) + (Step(i)*(j-1))
    end do
end do

end subroutine


subroutine GenerateModels(NumMods,TotNumMods,Grid,Models)

integer,dimension(:),intent(in)   :: NumMods
integer,intent(in) :: TotNumMods
real*8, dimension(:,:),intent(in) :: Grid

real*8,dimension(:,:),intent(out) :: Models

integer :: i,j,k,ModsBefore,ModsAfter,jump


do k = size(Grid,2),1,-1 ! k controls accessed variable (starting from final)
  ModsBefore = 1

  do i = 1,size(NumMods) - (1 + size(Grid,2) - k)
    ModsBefore = ModsBefore*NumMods(i)
  end do

  ModsAfter = (TotNumMods/ModsBefore)/NumMods(i) ! ModsAfter is repeat number

  do i = 1,ModsBefore
    jump = (i-1)*NumMods(k)*ModsAfter
    do j = 1,NumMods(k)
      Models( ((j-1)*ModsAfter)+1+jump:(j*ModsAfter)+jump ,k) = Grid(j,k) ! create all models in Models
    end do
  end do
end do

end subroutine


subroutine MisfitAnalysis(misfits,Cells,Keep)

real*8,dimension(:),intent(in) :: misfits
integer,intent(in) :: Cells

integer,dimension(:),allocatable,intent(out) :: Keep

logical,dimension(:),allocatable :: mask
integer :: i


allocate(Keep(Cells),mask(size(misfits,1)))

mask = .True.
do i = 1,Cells
  Keep(i) = minloc(misfits,dim=1,mask=mask)
  mask(minloc(misfits,mask)) = .False.
end do

end subroutine


subroutine ReadCLA(Debug,Writeall,Kill,TargMisM,TargMisK)

logical,intent(out) :: Debug
logical,intent(out) :: Writeall,Kill
real*8,intent(out) :: TargMisM,TargMisK

character(len=10) :: arg
integer :: i
logical :: exists


Writeall = .True.
Kill     = .False.
Debug    = .False.
TargMisM = 0.0
TargMisK = 0.0

print'(A,I0,A)','Program called with ',command_argument_count(),' CLA(s)'

do i = 1,command_argument_count()
if(arg == '-m') then
  call get_command_argument(i,arg)
  read(arg,'(F9.5)') TargMisM
  print'(A,F9.5)','User defined Misfit Storage Target ',TargMisM
  cycle
elseif(arg == '-k') then
  call get_command_argument(i,arg)
  read(arg,'(F9.5)') TargMisK
  print'(A,F9.5)','User defined Misfit Kill Target ',TargMisK
  cycle
end if

  call get_command_argument(i,arg)

  select case (arg)

    case ('-h')
      call print_help()
      stop "STOP as -h flag prevents program running"

    case ('-i')
      call print_inst()
      stop "STOP as -i flag prevents program running"

    case ('-c')
      call print_cite()
      stop "STOP as -c flag prevents program running"

    case ('-d')
      print*,'debug option set to TRUE'
      Debug = .True.

    case ('-m')
	  print*,'Program will write only best models to "Misfits-m.txt"'
	  inquire(file="Misfits-m.txt", exist=exists)
	  if(exists) call execute_command_line('rm Misfits-m.txt')

    case ('-k')
	  Kill = .True.

    case ('-nw')
  	  Writeall = .False.

    case default
      print '(2a, /)', 'Unrecognised command-line option: ', arg
      call print_help()
      stop

  end select
end do

if(command_argument_count() == 0) print"(A,/)",'The default program will run and produce "AllModels.txt"'

end subroutine


subroutine ReadInput(VarNames,VarValsMin,VarValsMax,NumMods,Cells,Levels,TotNumMods,ModRpt,RptMod)

character(len=20),dimension(:),allocatable,intent(out) :: VarNames
real*8, dimension(:),allocatable,intent(out) :: VarValsMin
real*8, dimension(:),allocatable,intent(out) :: VarValsMax
integer,dimension(:),allocatable,intent(out) :: NumMods
integer,intent(out) :: Cells, Levels
integer,intent(out) :: TotNumMods
logical,intent(out) :: ModRpt
integer,intent(out) :: RptMod

integer :: i,NumVars


TotNumMods = 1
ModRpt = .True.

open(unit=1,file='Input.txt')
  read(1,*) NumVars,Cells,Levels
  allocate(VarNames(NumVars),VarValsMin(NumVars),VarValsMax(NumVars),NumMods(NumVars))
  read(1,*)

  do i = 1,NumVars
    read(1,*) VarNames(i)
	read(1,*) VarValsMin(i)
    read(1,*) VarValsMax(i)
    read(1,*) NumMods(i)
	TotNumMods = TotNumMods*NumMods(i)
	if(mod(NumMods(i),2)==0) ModRpt = .False.
	read(1,*)
  end do
close(1)

! Find repeated model - will always be central model
if(ModRpt) then
  RptMod = (TotNumMods+1)/2
else
  RptMod = 0
end if

end subroutine


subroutine UpdateSetup(Models,Step,VarValsMin,VarValsMax)

real*8,dimension(:),intent(in) :: Models
real*8,dimension(:),intent(in) :: Step

real*8,dimension(:),intent(out) :: VarValsMin
real*8,dimension(:),intent(out) :: VarValsMax

integer :: i


do i = 1,size(Models,1)
  VarValsMin(i) = Models(i) - Step(i)/2.0
  VarValsMax(i) = Models(i) + Step(i)/2.0
end do

end subroutine


subroutine WriteMisfitTarget(AllModels,TargMisM,LvlNumMods,Levels,Step,iLevel)

real*8,dimension(:,:),intent(in) :: AllModels
real*8,intent(in) :: TargMisM
integer,intent(in):: LvlNumMods,Levels
real*8,dimension(:),intent(in) :: Step
integer,intent(in),optional :: iLevel

integer :: i,j
character(len=100) :: frmt
logical,save :: Called = .False.


if(.not. Called) then
  open(1,file='Misfits-m.txt')
  Called = .True.
else
  open(1,file='Misfits-m.txt',access='append',status='old')
end if

write(frmt,'(A,I0,A)') "(A,I0,A,I0,A,I0,A,",size(Step),"F11.5)"
write(1,frmt) 'Level ', iLevel ,' of ',Levels,' with ',LvlNumMods,' models and Cell sizes +-:',Step/2.0

write(frmt,'(A,I0,A)') "(",size(AllModels,2),"F11.5)"
do i = 1,size(AllModels,1)
  if(AllModels(i,size(AllModels,2)) <= TargMisM) write(1,frmt) (AllModels(i,j), j=1,size(AllModels,2))
end do
write(1,*)

close(1)

end subroutine


subroutine WriteOut(Part,LvlNumMods,Levels,AllModels,Step,iLevel)

character(len=3),intent(in) :: Part
integer,intent(in):: LvlNumMods,Levels ! LvlNumMods = Cells*TotNumMods
real*8,dimension(:,:),intent(in) :: AllModels
real*8,dimension(:),intent(in) :: Step
integer,intent(in),optional :: iLevel

integer :: i,j
character(len=100) :: frmt


select case(Part)
  case('1st')! Write first Level output to file
    open(unit=1,file='AllModels.txt')
	  write(frmt,'(A,I0,A)') "(A,I0,A,I0,A,",size(Step),"F11.5)"
      write(1,frmt) 'Initial Level of ',Levels,' with ',LvlNumMods,' models and Cell sizes +-:',Step/2.0
	  write(frmt,'(A,I0,A)') "(",size(AllModels,2),"F11.5)"
      do i = 1,size(AllModels,1)
        write(1,frmt) (AllModels(i,j), j=1,size(AllModels,2))
      end do
	  write(1,*)
    close(1)

  case('Nth')! Write Nth level output to file
    open(unit=1,file='AllModels.txt',access='append',status='old')
	  write(frmt,'(A,I0,A)') "(A,I0,A,I0,A,I0,A,",size(Step),"F11.5)"
      write(1,frmt) 'Level ', iLevel ,' of ',Levels,' with ',LvlNumMods,' models and Cell sizes +-:',Step/2.0
      write(frmt,'(A,I0,A)') "(",size(AllModels,2),"F11.5)"
	  do i = 1,size(AllModels,1)
        write(1,frmt) (AllModels(i,j), j=1,size(AllModels,2))
      end do
	  write(1,*)
    close(1)

end select

end subroutine


!-------------------- Non-running CLA ---------------------------------------------

subroutine print_help()

  print '(A,/)', 'command-line options:'
  print '(A)'  , '  -h          print command line argument information and stop'
  print '(A)'  , '  -i          print instructions for program use and stop'
  print '(A)'  , '  -c          print citation information and stop'
  print '(A)'  , '  -d          set debug option to TRUE - print extra output to terminal'
  print '(A)'  , '  -k          set misfit value at which search will exit'
  print '(A)'  , '  -m          set misfit value below which all models are stored in Misfits-m.txt'
  print '(A)'  , '  -nw         prevents creation of AllModels.txt output file'
  print '(A,/)', '  no option   normal execution of program - run search to completion and create AllModels.txt'

  print '(A)', 'Copyright (C) 2002 Jon May'
  print '(A)', 'This program comes with ABSOLUTELY NO WARRANTY.'
  print '(A)', 'This is free software, you are welcome to redistribute it under conditions'
  print '(A)', 'set out in GNU General Public License version 3, or any later version.' 
  print '(A,/)', 'See the included license file, or https://www.gnu.org/licenses/, for details.'

end subroutine print_help


subroutine print_inst()

  print '(A)'  ,'The user must complete the input file "Input.txt" with the information required for the'
  print '(A,/)','Grid Search. Below is a description of the file with 2 variables:'

  print '(A)'  ,'Number of variables, number of best Cells to keep, number of search Levels (incl. initial)'
  print '(A)'  ,'----------'
  print '(A)'  ,'Var 1 Name'
  print '(A)'  ,'Var 1 Min Value'
  print '(A)'  ,'Var 1 Max Value'
  print '(A)'  ,'Var 1 Number of models'
  print '(A)'  ,'----------'
  print '(A)'  ,'Var 2 Name'
  print '(A)'  ,'Var 2 Min Value'
  print '(A)'  ,'Var 2 Max Value'
  print '(A)'  ,'Var 2 Number of models'
  print '(A,/)','----------'

  print '(A)'  ,'Variable names are not currently used in the search algorithm but are present incase a'
  print '(A,/)','linked forward model requires them as input.'

  print '(A)'  ,'The user should update the subroutine calls to "Example" inside "main.f90" at line 74 & 133'
  print '(A)'  ,'to link to their desired forward model or misfit calculator. This new subroutine can take'
  print '(A)'  ,'any input but must return the misfit values to the program (in a vector) where they are'
  print '(A)'  ,'stored and analysied in order to select the best cell(s) for the next level. The newly'
  print '(A)'  ,'linked subroutine need not be written in Fortran or contained within the files provided'
  print '(A,/)','but it is up to the user to ensure any mixed language programs compile and link properly.'

  print '(A)'  ,'Currently the first call to "Example" accepts 3 inputs (Models,VarValsMax,Misfits) but the'
  print '(A)'  ,'second call takes 5 (Models,VarValsMax,Misfits,RptMod,Debug), both "RptMod" & "Debug" are'
  print '(A)'  ,'optional and only needed for the scoring subroutine call within the main program loops.'
  print '(A)'  ,'This is because "RptMod" allows the subroutine to skip a single model if it is a repeated'
  print '(A)'  ,'model from the previous level, the "Debug" variable is simply to print information from'
  print '(A)'  ,'the subroutine if the "-d" flag is used and isn''t strictly required.'
  print '(A)'  ,'Therefore any new subroutine linked should also have the option of "RptMod" in order to'
  print '(A,/)','prevent repeated tests unless the user programs another method.'

  print '(A)'  ,'The default befaviour of the program is to generate an output file "AllModels.txt", this'
  print '(A)'  ,'file stores all models run in order, seperated into levels, with their respective misfit'
  print '(A)'  ,'values. Should the user run the program with the "-nw" flag then the program will not'
  print '(A)'  ,'write "AllModels.txt", it is therefore important to pair this with "-m" in order to'
  print '(A,/)','generate the "Misfits-m.txt" output file'

  print '(A)', 'Copyright (C) 2002 Jon May'
  print '(A)', 'This program comes with ABSOLUTELY NO WARRANTY.'
  print '(A)', 'This is free software, you are welcome to redistribute it under conditions'
  print '(A)', 'set out in GNU General Public License version 3, or any later version.' 
  print '(A,/)', 'See the included license file, or https://www.gnu.org/licenses/, for details.'

end subroutine


subroutine print_cite()

  print '(A)'  ,'The program should be cited if used in any work.'
  print '(A,/)','Below is a bibtex format citation example:'

  print '(A)'  ,'@misc{May2022GridSearch,'
  print '(A)'  ,'author       = {Jon May},'
  print '(A)'  ,'howpublished = {Online, Zenodo},'
  print '(A)'  ,'title        = {Fortran Grid Search (1.0)},'
  print '(A)'  ,'month        = {June},'
  print '(A)'  ,'year         = {2022},'
  print '(A)'  ,'url          = {url},'
  ! print '(A)'  ,'Git url      = {Giturl},'
  print '(A)'  ,'doi          = {doi},'
  print '(A,/)','}'

  print '(A)', 'Copyright (C) 2002 Jon May'
  print '(A)', 'This program comes with ABSOLUTELY NO WARRANTY.'
  print '(A)', 'This is free software, you are welcome to redistribute it under conditions'
  print '(A)', 'set out in GNU General Public License version 3, or any later version.' 
  print '(A,/)', 'See the included license file, or https://www.gnu.org/licenses/, for details.'

end subroutine


end module
