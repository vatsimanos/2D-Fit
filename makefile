#makefile fuer Patch Clamp Analysis



SRCDIR = ~/2DFit/src/
OBJDIR = ~/2DFit/obj/
BINDIR = ~/2DFit/obj/
BINFILE= 2DFit64

CC = mpic++	#gcc c/c++ #mpiicpc #mpic++

VPATH = %.o $(OBJDIR) 
VPATH = $(SRCDIR)


INCDIRS = -I ~/2DFit/galib247
 
LIBDIRS = -L ~/2DFit/galib247/ga 
DEFS = -Dcplusplus #-Ddebug_simulat #-Ddebug_mem_matrix2 #-Ddebug_mem_patchio #-Ddebug_mem_datarray #-Ddebug_mem_setfile #-Ddebug_mem_rausch #-Ddebug_mem_simulat  #-Ddebug_datarray #-Ddebug_simulat  #-DFORTIFY 

CFLAGS =  -O2 -march=native $(DEFS) $(INCDIRS)  #-O3 -march=native -pipe #`--cflags --pkg-config libs gtk+-3.0`

 

#-O3  	#optimization high lvl
#-march=native  #native processor optimazion   
#-pedantic 
# -g3          debug infos high lvl        
# -Wall       alle Warnings an
# -Ddef	      define def 1
# -mwindows   ohne DOS Fenster
# -mno-cygwin mingw mit cygwin 
# -mpentium   
# -fnative-struct

LDFLAGS = $(LIBDIRS) -lga #`pkg-config --cflags --libs gtk+-3.0 atk`

O = \
	$(OBJDIR)\gtksig.o\
	$(OBJDIR)\gtkmain.o\
	$(OBJDIR)\gtkfunc.o\
	$(OBJDIR)\declare.o\
	$(OBJDIR)\error.o\
	$(OBJDIR)\daten.o\
	$(OBJDIR)\daten2.o\
	$(OBJDIR)\detector.o\
	$(OBJDIR)\fit.o\
	$(OBJDIR)\y_shut.o\
	$(OBJDIR)\steady.o\
	$(OBJDIR)\matu.o\
	$(OBJDIR)\steuerung.o\
	$(OBJDIR)\dwell2d.o\
	$(OBJDIR)\testout.o\
	$(OBJDIR)\sim_fit.o\
	$(OBJDIR)\err.o\
	$(OBJDIR)\patchio.o\
	$(OBJDIR)\rausch.o\
	$(OBJDIR)\setfile.o\
	$(OBJDIR)\timeseries.o\
	$(OBJDIR)\filter.o\
	$(OBJDIR)\channel.o\
	$(OBJDIR)\matrix2.o\
	$(OBJDIR)\datarray.o\
	$(OBJDIR)\bessel.o\
	$(OBJDIR)\MultiLevel2D.o\
	$(OBJDIR)\fft.o\
	

main:	$(O) 
	$(CC)  -o $(BINDIR)\$(BINFILE) $(O) $(LDFLAGS)

$(OBJDIR)\\%.o: $(SRCDIR)/%.c
	$(CC) -c $(CFLAGS) -o $@ $<  


include makefile.dep


.PHONY: clean cleanbin cleanall start dep

cleanbin:       
	del $(BINDIR)\$(BINFILE)

clean:  
	del $(OBJDIR)\*.o

cleanall:	clean cleanbin

start:
	$(BINDIR)\$(BINFILE)

dep:
	$(CC) -MM -E $(INCDIRS) src/*.c > makefile.dep

