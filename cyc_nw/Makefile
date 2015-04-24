ifdef SystemRoot
	RM = del /Q
	EXT = .exe
else
	ifeq ($(shell uname), Linux)
		RM = rm -f
		EXT =
	endif
endif

all:
	g++ -std=c++11 -D_USE_64 -msse3 -O3 -fomit-frame-pointer -funroll-loops -DNDEBUG -lm -I . cyc_nw.cpp -o cyc_nw$(EXT)

clean:
	$(RM) cyc_nw.o cyc_nw
