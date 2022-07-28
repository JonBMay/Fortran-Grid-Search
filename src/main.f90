! Grid Search Program

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

program GridSearch

use subs
use sims

implicit none

character(len=20),dimension(:),allocatable :: VarNames
real*8, dimension(:),  allocatable :: VarValsMin,VarValsMax
real*8, dimension(:,:),allocatable :: Grid
real*8, dimension(:,:),allocatable :: Models
real*8, dimension(:,:),allocatable :: AllModels
real*8, dimension(:),  allocatable :: Step,StepOld
real*8, dimension(:),  allocatable :: Misfits
integer,dimension(:),  allocatable :: NumMods
integer,dimension(:),  allocatable :: Keep
real*8,dimension(:,:), allocatable :: KeepCells

integer :: TotNumMods,Cells,Levels
integer :: i,j,iLevel,iCell,RptMod
logical :: Debug,Writeall,Kill,MisExit,ModRpt
real*8 :: TargMisM,TargMisK


!----------------------------------------------------------------------------------------------------
!----------------------------------Fixed Program (edit at own risk)----------------------------------
! Read Command Line Arguments
call ReadCLA(Debug,Writeall,Kill,TargMisM,TargMisK)
! Read initial user input & allocate
call ReadInput(VarNames,VarValsMin,VarValsMax,NumMods,Cells,Levels,TotNumMods,ModRpt,RptMod)
allocate(Models(TotNumMods,size(NumMods)))
allocate(Misfits(TotNumMods))
allocate(AllModels(TotNumMods,1+size(NumMods)))
allocate(StepOld(size(NumMods)))
allocate(KeepCells(Cells,1+size(NumMods)))

if(Debug) then
  print*,'The number of models in each dimenison is:',NumMods
  if(ModRpt) then
    print'(X,A)','Since all dimensions have an odd number of models the central model will be a repeat of the'
    print'(X,A)','model from the previous level, this model will therefore be skipped and the previous result used.'
    print'(X,A,I0,A,I0,/)','Number of repeated model inside next level: ',RptMod,' / ',size(Models,1)
  else
    print*,'Since at least 1 dimension has an even number of models there will be no repeated models'
  end if
end if

! Calculate the Grid of values
call CalcGrid(VarValsMin,VarValsMax,NumMods,Grid,Step)
! Generate Models from Grid
call GenerateModels(NumMods,TotNumMods,Grid,Models)

!----------------------------------------------------------------------------------------------------
!---------------------------------------------User Edit----------------------------------------------
! User call to misfit calculation
call Example(Models,VarValsMax,Misfits,RptMod=TotNumMods+1,Debug=Debug)

!----------------------------------------------------------------------------------------------------
!----------------------------------Fixed Program (edit at own risk)----------------------------------
AllModels(1 : TotNumMods,1:size(Models,2)) = Models
AllModels(1 : TotNumMods,1+size(Models,2)) = Misfits

! Decide which cells to keep
call MisfitAnalysis(AllModels(:,size(AllModels,2)),Cells,Keep)
do i = 1,size(Keep)
  KeepCells(i,1:size(KeepCells,2)-1) = AllModels(Keep(i),1:size(AllModels,2)-1)
  KeepCells(i,  size(KeepCells,2))   = AllModels(Keep(i),  size(AllModels,2))
end do

! Write first Level output to file
if(Writeall) call WriteOut('1st',TotNumMods,Levels,AllModels,Step)
if(minval(KeepCells(:,size(KeepCells,2)))<TargMisM) call WriteMisfitTarget(AllModels,TargMisM,TotNumMods,Levels,Step,1)

if(minval(KeepCells(:,size(KeepCells,2))) > TargMisK) then
  do iLevel = 2,Levels
    if(iLevel==2) then
	  deallocate(AllModels)
	  allocate(AllModels(Cells*TotNumMods,1+size(Models,2)))
    end if
	StepOld = Step

    if(Debug) then
	  if(iLevel==2) then
	    print*,'Analysed models:',1,':',(1+((iLevel-2)*Cells))*TotNumMods
	  else
	    print*,'Analysed models:',1+(1+((iLevel-3)*Cells))*TotNumMods,':',(1+((iLevel-2)*Cells))*TotNumMods
      end if
	  do i = 1,size(Keep)
        print*,'Keeping model:',Keep(i), &
             & 'Keeping model:',KeepCells(i,1:size(KeepCells,2)-1), &
		     & ' With misfit:' ,KeepCells(i,size(KeepCells,2))
      end do
	  print*,''
    end if

! Loop through, create models and misfits for all kept cells
    do iCell = 1,Cells
! Update setup values using selected cell values
      call UpdateSetup(KeepCells(iCell,1:size(KeepCells,2)-1),StepOld,VarValsMin,VarValsMax)
      if(Debug) then
        print*,'Cell:',Keep(iCell),' = ',KeepCells(iCell,1:size(KeepCells,2)-1)
        print*,'VarValsMin:',VarValsMin
        print*,'VarValsMax:',VarValsMax
		print*,''
      end if

! Calculate a new Grid
      call CalcGrid(VarValsMin,VarValsMax,NumMods,Grid,Step)
! Generate all models from Grid
      call GenerateModels(NumMods,TotNumMods,Grid,Models)

!----------------------------------------------------------------------------------------------------
!---------------------------------------------User Edit----------------------------------------------
! User call to misfit calculation
      call Example(Models,VarValsMax,Misfits,RptMod=RptMod,Debug=Debug)

!----------------------------------------------------------------------------------------------------
!----------------------------------Fixed Program (edit at own risk)----------------------------------
      AllModels(1+(iCell-1)*TotNumMods : iCell*TotNumMods,1:size(AllModels,2)-1) = Models
  	  AllModels(1+(iCell-1)*TotNumMods : iCell*TotNumMods,  size(AllModels,2))   = Misfits
! Copy repeated model data
      if(ModRpt) then
        AllModels(RptMod+(iCell-1)*TotNumMods,1:size(AllModels,2)-1) = KeepCells(iCell,1:size(KeepCells,2)-1)
  	    AllModels(RptMod+(iCell-1)*TotNumMods,  size(AllModels,2))   = KeepCells(iCell,  size(KeepCells,2))
      end if
    end do ! loop over Cells

! Write level output to file
    if(Writeall) call WriteOut('Nth',Cells*TotNumMods,Levels,AllModels,Step,iLevel)
    if(minval(AllModels(:,size(AllModels,2)))<TargMisM) call WriteMisfitTarget(AllModels,TargMisM,TotNumMods,Levels,Step,iLevel)

    call MisfitAnalysis(AllModels(:,size(AllModels,2)),Cells,Keep)
	do i = 1,size(Keep)
      KeepCells(i,1:size(KeepCells,2)-1) = AllModels(Keep(i),1:size(AllModels,2)-1)
      KeepCells(i,  size(KeepCells,2))   = AllModels(Keep(i),  size(AllModels,2))
	end do

! Check Kill misfit value
    if(Kill .and. minval(KeepCells(:,size(KeepCells,2))) <= TargMisK) then
	  MisExit = .True.
	  exit
	end if
  end do ! end loop over Levels
else ! minval(KeepCells(:,size(KeepCells,2))) <= TargMisK at first level
  MisExit = .True.
end if

if(MisExit) then
  if(Writeall) then
    open(unit=1,file='AllModels.txt',access='append',status='old')
	write(1,'(/,A,F10.7)') 'Program was halted as a misfit was <= ',TargMisK
  end if

  if(TargMisM /= 0.0) then
    open(1,file='Misfits-m.txt',access='append',status='old')
	write(1,'(/,A,F10.7)') 'Program was halted as a misfit was <= ',TargMisK
  end if
end if

end program
