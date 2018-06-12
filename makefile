	SHELL	=	/bin/sh
	EXEC	=	spic1d
	FC	=	g++
objects=diag.o init.o main.o maxw.o move.o pic.o pois.o setv.o solvfiel.o 
$ (EXEC) :$(objects)
#echo objects
	$(FC) $(CFLAGS) -o $(EXEC) $(objects) $(lib)

$(objects): %.o: %.cpp
	$(FC) -c $(CFLAGS)  $(INCS) $< -o $@

/PHONY :clean
clean:
	rm *.o $(EXEC)	
