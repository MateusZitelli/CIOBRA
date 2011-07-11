all:
	gcc simula.c -o simula.bin `sdl-config --cflags` `sdl-config --libs` -lSDL_gfx
