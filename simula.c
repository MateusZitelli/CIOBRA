//CIOBRA - C implementation of Barnes-Hut Algorithm
//Copyright (C) 2011 Mateus Zitelli (zitellimateus@gmail.com)

//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.

//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.

//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "SDL.h"
#include "SDL/SDL_gfxPrimitives.h"
#define WIDTH 600
#define HEIGHT 600
#define BPP 4
#define DEPTH 32
#define MAXBODIES 50
#define PRECISION 0.01

long last;

struct Vect {
	double x;
	double y;
};

struct Body {
	struct Vect P;
	struct Vect V;
	struct Vect F;
	struct Quad *Q;
	double M;
	double C;
};

struct Quad {
	struct Vect S;
	struct Vect E;
	struct Body **B;
	struct Quad *QA;
	struct Quad *QB;
	struct Quad *QC;
	struct Quad *QD;
	struct Quad *QUP;
	struct Body *COM;
	long N;
	int I;
};

void setpixel(SDL_Surface * screen, int x, int y, Uint8 r, Uint8 g, Uint8 b)
{
	if(x > WIDTH - 1 || x < 0 || y > HEIGHT - 1 || y < 0)
		return; 
	Uint32 *pixmem32;
	Uint32 colour;
	colour = SDL_MapRGB(screen->format, r, g, b);
	pixmem32 = (Uint32 *) screen->pixels + WIDTH * y + x;
	*pixmem32 = colour;
}

void set_quad(struct Quad *Q, double sx, double sy, double ex, double ey,
	      struct Quad *QUP)
{
	Q->S.x = sx;
	Q->S.y = sy;
	Q->E.x = ex;
	Q->E.y = ey;
	//free(Q->B);
	//Q->B = (struct Body **)malloc(sizeof(struct Body *) * MAXBODIES);
	Q->QA = NULL;
	Q->QB = NULL;
	Q->QC = NULL;
	Q->QD = NULL;
	Q->QUP = QUP;
	Q->COM = NULL;
	Q->N = 0;
}

void set_body(struct Body *B, double x, double y, double charge)
{
	B->P.x = x;
	B->P.y = y;
	B->V.x = 0;
	B->V.y = 0;
	B->F.x = 0;
	B->F.y = 0;
	B->M = 1;
}

void divide(struct Quad *Q, struct Quad *LQ)
{
	int i;
	int num = 0;
	if (Q == NULL)
		return;
	if (Q->N <= 1)
		return;
	double sx, sy, mx, my, ex, ey;
	sx = Q->S.x;
	sy = Q->S.y;
	ex = Q->E.x;
	ey = Q->E.y;
	mx = (ex + sx) / 2.0;
	my = (ey + sy) / 2.0;
	for (i = 0; i < Q->N; i++) {
		if(Q->B[i]->P.x > ex || Q->B[i]->P.y > ey || Q->B[i]->P.x < sx || Q->B[i]->P.y < sy)
			continue;
		if (Q->B[i]->P.x < mx && Q->B[i]->P.y < my) {
			if (!(num & 1)) {
				set_quad(&LQ[last], sx, sy, mx, my, Q);
				Q->QA = &LQ[last++];
				num |= 1;
			}
			Q->QA->B[Q->QA->N++] = Q->B[i];
		} else if (Q->B[i]->P.x > mx && Q->B[i]->P.y < my) {
			if (!(num & 2)) {
				set_quad(&LQ[last], mx, sy, ex, my, Q);
				Q->QB = &LQ[last++];
				num |= 2;
			}
			Q->QB->B[Q->QB->N++] = Q->B[i];
		} else if (Q->B[i]->P.x > mx && Q->B[i]->P.y > my) {
			if (!(num & 4)) {
				set_quad(&LQ[last], mx, my, ex, ey, Q);
				Q->QC = &LQ[last++];
				num |= 4;
			}
			Q->QC->B[Q->QC->N++] = Q->B[i];
		} else if (Q->B[i]->P.x < mx && Q->B[i]->P.y > my) {
			if (!(num & 8)) {
				set_quad(&LQ[last], sx, my, mx, ey, Q);
				Q->QD = &LQ[last++];
				num |= 8;
			}
			Q->QD->B[Q->QD->N++] = Q->B[i];
		}
	}
	//printf("%li\n", last);
	divide(Q->QA, LQ);
	divide(Q->QB, LQ);
	divide(Q->QC, LQ);
	divide(Q->QD, LQ);
}

void get_com(struct Quad *Q)
{
	//if (Q->COM != NULL)
	//	return *Q->COM;
	long i;
	struct Body BO;
	BO.P.x = 0;
	BO.P.y = 0;
	BO.M = 0;
	for (i = 0; i < Q->N; i++){
		BO.P.x += Q->B[i]->P.x * Q->B[i]->M;
		BO.P.y += Q->B[i]->P.y * Q->B[i]->M;
		BO.M += Q->B[i]->M;
	}
	if(BO.C){
		BO.P.x /= BO.M;
		BO.P.y /= BO.M;
	}
	Q->COM = &BO;
	//return (BO);
}

void apply_force(struct Body *b1, struct Body b2)
{
	double F, D, CX, CY;
	CX = (b1->P.x - b2.P.x) * 10E3;
	CY = (b1->P.y - b2.P.y) * 10E3;
	D = CX * CX + CY * CY;	//Disntace squared
	if(D > 10)
		F = 9 * 10E5 * b1->M * b2.C / D;
	b1->F.x -= F * CX / D;
	b1->F.y -= F * CY / D;
}

void solve(struct Body *B, struct Quad *Q, int ind, int mode)
{
	if(Q->N > 0){
		struct Body COM;
		double DX, DY;
		if (Q->QA != NULL && Q->QA->I != ind) {
			COM = *(Q->QA->COM);
			DX = COM.P.x - B->P.x;
			DY = COM.P.y - B->P.y;
			printf("%f\n", COM.M);
			if((Q->E.x - Q->S.x) / sqrt(DX * DX + DY * DY) < PRECISION && Q->N > 1){
				solve(B, Q->QA, ind, 0);
			}else{
				apply_force(B,COM);
			}
		}
		if (Q->QB != NULL && Q->QB->I != ind) {
			COM = *(Q->QB->COM);
			DX = COM.P.x - B->P.x;
			DY = COM.P.y - B->P.y;
			if((Q->E.x - Q->S.x) / sqrt(DX * DX + DY * DY) < PRECISION && Q->N > 1){
				solve(B, Q->QB, ind, 0);
			}else{
				apply_force(B,COM);
			}
		}
		if (Q->QC != NULL && Q->QC->I != ind) {
			COM = *(Q->QC->COM);
			DX = COM.P.x - B->P.x;
			DY = COM.P.y - B->P.y;
			if((Q->E.x - Q->S.x) / sqrt(DX * DX + DY * DY) < PRECISION && Q->N > 1){
				solve(B, Q->QC, ind, 0);
			}else{
				apply_force(B,COM);
			}
		}
		if (Q->QD != NULL && Q->QD->I != ind) {
			COM = *(Q->QD->COM);
			DX = COM.P.x - B->P.x;
			DY = COM.P.y - B->P.y;
			if((Q->E.x - Q->S.x) / sqrt(DX * DX + DY * DY) < PRECISION && Q->N > 1){
				solve(B, Q->QD, ind, 0);
			}else{
				apply_force(B,COM);
			}
		}
	}
	if (Q->QUP == NULL)
		return;
	if(mode)
		solve(B, Q->QUP, ind, 1);
}

void apply_forces(struct Body *B)
{
	//printf("%f, %f\n", B->P.x, B->P.y);
	B->V.x += B->F.x / B->M;
	B->V.y += B->F.y / B->M;
	B->P.x += B->V.x;
	B->P.y += B->V.y;
	B->F.x = 0;
	B->F.y = 0;
}

int main(void)
{
	srand(time(0));
	long i = 0;
	SDL_Surface *screen;
	SDL_Event event;
	int keypress = 0;
	int h = 0;
	if (SDL_Init(SDL_INIT_VIDEO) < 0)
		return 1;
	if (!(screen = SDL_SetVideoMode(WIDTH, HEIGHT, DEPTH, SDL_HWSURFACE))) {
		SDL_Quit();
		return 1;
	}
	struct Quad *quads =
	    (struct Quad *)malloc(sizeof(struct Quad) * MAXBODIES * 10);
	for (i = 0; i < MAXBODIES * 10; i++) {
		quads[i].I = i;
		quads[i].B = (struct Body **)malloc(sizeof(struct Body *) * MAXBODIES);
	}
	struct Body *bodies =
	    (struct Body *)malloc(sizeof(struct Body) * MAXBODIES);
	set_quad(&quads[0], 0, 0, WIDTH, HEIGHT, NULL);
	for (i = 0; i < MAXBODIES; i++) {
		set_body(&bodies[i], rand() % 10000 / 10000.0 * WIDTH,
			 rand() % 10000 / 10000.0 * WIDTH, 1);
		quads[0].B[i] = &bodies[i];
	}
	quads[0].N = MAXBODIES - 1;
	while (!keypress) {
		last = 1;
		divide(&quads[0], quads);
		for(i = 0; i < last; i++){
			get_com(&quads[i]);
		}
		for (i = 0; i < last; i++) {
			#if 0
			rectangleRGBA(screen, quads[i].S.x, quads[i].S.y,
				      quads[i].E.x, quads[i].E.y, 0, 255, 0,
				      255);
			#endif
			#if 1
			if(quads[i].COM != NULL)
				setpixel(screen, quads[i].COM->P.x, quads[i].COM->P.y, 255, 255, 255);
			#endif
			if (quads[i].N == 1) {
				solve(quads[i].B[0], quads[i].QUP, quads[i].I, 1);
			}
		}
		SDL_Flip(screen);
		SDL_FillRect(screen, NULL, 0x000000);
		while (SDL_PollEvent(&event)) {
			switch (event.type) {
			case SDL_QUIT:
				keypress = 1;
				break;
			case SDL_KEYDOWN:
				keypress = 1;
				break;
			}
		}
		for (i = 0; i < MAXBODIES -1; i++) {
			apply_forces(&bodies[i]);
			setpixel(screen, bodies[i].P.x, bodies[i].P.y, 255, 0, 0);
		}
		for(i = 1; i < last; i++)
			quads[i].N = 0;
		//printf("%li\n", last);
	}
	free(quads);
	free(bodies);
	SDL_Quit();
	return (0);
}

