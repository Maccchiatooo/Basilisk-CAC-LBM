CC = qcc
SRC = multiphase.c
LIB = /home/chunheng/basilisk/src/gl
CFLAGS = -Wall -O3
LDLIBS = -lm -L${LIB} -lglutils -lfb_glx -lGLU -lGLEW -lGL -lX11 
sessil:	${SRC}
	${CC} ${CFLAGS} ${SRC} -o czhao29 ${LDLIBS}
clean:
	@echo "Cleaning files..."
	rm -rf *.mp4 *.dat *.exe *.html dump *.vtk
