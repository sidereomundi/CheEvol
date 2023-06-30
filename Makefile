F77 = gfortran

FLAGS  =  -O2 -Wall 

FILES2  =  IMFsub.f chemFin6omo.f

.f.o:
	$(F77) $(FLAGS) -c $<

Fede:  gradfede.o $(FILES2) Makefile
	$(F77) $(FLAGS) -o fede.x  gradfede.o $(FILES2) 
