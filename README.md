# Fortran Grid Search
This program is a simple serial grid search algorithm written in Fortran.

### Program explanation

File "main.f90" contains the main working of the program.
File "subs.f90" contains the subroutines necessary for general running of the program.
File "sims.f90" contains an Example subroutine and acts as a potential location for future user
	defined Fortran forward models or misfit calculators. It is called on lines 58 & 119 of "main.f90".

The user is under no obligation to connect Fortran subroutines; indeed Fortran has the ability to interface  
with routines written in other languages, however the operability of any mixed language connection(s) is up  
to the user.  
  
Currently the first call to "Example" accepts 4 inputs (Models,VarValsMax,Misfits,Debug) but the second call  
takes 5 (Models,VarValsMax,Misfits,RptMod,Debug), both "RptMod" & "Debug" are optional. The "Debug" variable  
is simply to print information from the subroutine if the "-d" flag is used and isn't strictly required.  
"RptMod" allows the subroutine to skip a single model if it is a repeated model from a previous level and is  
therefore only required in the level loop. Therefore any new subroutine linked should also have the option  
of "RptMod" in order to prevent repeated tests unless the user programs another method.  
"Example" has programmed two options to generate a 'target' model, either a user defined vector or  
a loop to generate a random number for each variable, the user should comment/uncomment as they wish.  
  
  
The program can accept numerous Command Line Arguments (CLAs), the complete list of all CLAs is:  
no flag -> run the program in a default mode (outlined below)  
"-h"    -> print information about CLA options & STOP  
"-i"    -> print program instructions & STOP  
"-c"    -> print citation information and STOP  
"-d"    -> run the program with debug mode enabled (extra output to terminal) & generate "AllModels.txt"  
"-k"    -> user defined target misfit value that will exit the search & generate "AllModels.txt"  
"-m"    -> user defined target misfit value under which all models are stored, generate "Misfits-m.txt" & "AllModels.txt"  
"-nw"   -> prevent writing of "AllModels.txt"  
  
Since "-h", "-i" and "-c" call a STOP they may not be combined, however the "-d","-k", "-m" and "-nw" flags  
may be combined on any single program invocation. If "-m" and/or "-k" is selected the next CLA should be a  
real number, this will define the user's target misfit value for storage and stopping respectively.  
E.G: "./gridsearch -m 2.5 -k 1.5 -nw" will store all models with misfit < 2.5 inside "Misfits-m.txt",  
automatically STOP when any misfit is <= 1.5 and will not generate "AllModels.txt".  
Whereas "./gridsearch -m 2.5 -d" will store all models with misfit < 2.5 inside "Misfits-m.txt", print all   
debugging information and will generate "AllModels.txt".  
For obvious reasons any "-k" option should be smaller than any "-m" option if paired.  
  
  
The default behaviour of the program is to generate an output file "AllModels.txt", this file stores all  
models run in finishing order, separated into levels, with their respective misfit values. By default the  
program will stop if any misfit=0.0.  
  
  
To compile: make gridsearch  
To clean: make clean OR make cleanall   
  
  
### Used in
This program was modified for use in the creation of the ShellSet program (DOI:)
