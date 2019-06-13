/*
UDF for droplet detection and localisation of phase areas in vof cases.


Settings:
	_PHASE_IDX : Phase index phase which should be detected
	_MIN_VOL_FRAC : Lower Limit for detection


WARNING: THIS IS AN EARLY VERSION THERE MAY BE INEXPECTED BUGS!

Ver: 0.01

Current Limitations & Todos:
- Only work in seriell (but u could save data and run this function via 
						scheme script on saved *.dat files)
- Only finds cell neighbors which are face neighbors (maybe you may also want 
														edge (node) neighbors)
- Only double precision solver supported, due to file i/o
- By now only single domain fluid/case supported 
									(due to loop over all threads in fluid domain (1))

Copyright 2019 Christian Schubert

License (MIT):

Permission is hereby granted, free of charge, to any person obtaining a copy of 
this software and associated documentation files (the "Software"), to deal in the 
Software without restriction, including without limitation the rights to use, copy, 
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, 
and to permit persons to whom the Software is furnished to do so, subject to the 
following conditions:

The above copyright notice and this permission notice shall be included in all 
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "udf.h"
#include "stdlib.h"
#include "mem.h"
#include "sg_mphase.h"
 

/* Settings */
#define _PHASE_IDX 0 /* for phase to detect droplets of */
#define _MIN_VOL_FRAC 0.01 /* Lower Limit for phase detection */


/* ------------------------------------------------------------------------- */

#define _MYDEBUG 0
#define STATE_OK 1
#define STATE_ERROR 0

#define CELL_CHECKED 1

/* ------------------------------------------------------------------------- */

struct Droplet
{
	cell_t *cell_list;
	int    no_cells;
	real   com[ND_ND];						/* center of mass (mass weighted average) */
	real   sumed_com_cell_weights[ND_ND];	/* intermediate values for com determination*/
	real   mass;
	real   vol;
};


int initDroplet(struct Droplet *d, cell_t c)
{
	int state = STATE_OK;
	int i;

	(*d).no_cells = 1;
	(*d).mass = 0;
	(*d).vol  = 0;
	(*d).cell_list = NULL;
	(*d).cell_list = (cell_t*) calloc((*d).no_cells, sizeof(cell_t));
	(*d).cell_list[0] = c;

	for(i = 0; i < ND_ND; ++i)
	{
		(*d).com[i] = 0;
		(*d).sumed_com_cell_weights[i] = 0;
	}
	
	if((*d).cell_list != NULL) 
	{}
	else 
	{
		DebugMessage("Error (calloc): No free memory for D.cells "
						"available!\n");
		state = STATE_ERROR;
	}

	return state;
}

int dropletCellsAppend(struct Droplet *d, cell_t val)
{
	int state = STATE_OK;

	(*d).no_cells ++;
	(*d).cell_list = (cell_t*) realloc((*d).cell_list, 
											(*d).no_cells*sizeof(cell_t));

	if((*d).cell_list != NULL) 
	{
		(*d).cell_list[(*d).no_cells-1] = val; /* append to end of array */
	}
	else 
	{
		DebugMessage("Error (calloc): No free memory for D.cells "
					"available!\n");
		state = STATE_ERROR;
	}	

	return state;
}


void updateDropletWeights(struct Droplet *d, real c_mass, real x_c[ND_ND])
{
	int i = 0;

	for(i = 0; i<ND_ND; ++i)
	{
		(*d).sumed_com_cell_weights[i] += x_c[i] * c_mass;
 	}
}


void setFinalDropletCOM(struct Droplet *d)
{
	int i = 0;

	for(i = 0; i<ND_ND; ++i)
	{
		(*d).com[i] = (*d).sumed_com_cell_weights[i] / (*d).mass;
	}
}


int printDList(char filename[], struct Droplet *DList, int sizeDList)
{
	FILE *f = NULL;
	int state = STATE_OK;
	int i,j;
	
	f = fopen(filename, "w");

	if(f == NULL)
	{
		Message("Error (printDList()): Unable to open file %s "
		        "for writing!\n", filename);
		state = STATE_ERROR;
	}
	else
	{
		for(i = 0; i < sizeDList; ++i)
		{
			fprintf(f, "Droplet: %i\n", i);
			fprintf(f, "Droplet cells: %i\n", DList[i].no_cells);
			fprintf(f, "Droplet mass: %lE\n", DList[i].mass);
			fprintf(f, "Droplet volume: %lE\n", DList[i].vol);
			fprintf(f, "Droplet com:");

			for(j = 0; j <ND_ND; ++j)
			{
				fprintf(f, " %lf", DList[i].com[j]);
			}
			fprintf(f, "\n");

			#if _MYDEBUG
			fprintf(f, "\nList of cells in droplet:\n");
			for(j = 0; j < DList[i].no_cells; ++j)
			{
				fprintf(f,"%i\n", DList[i].cell_list[j]);
			} 
			#endif

			fprintf(f,"\n----------------------------------------------\n");
		}
		
		fclose(f);
	}

	return state;
}


/* ------------------------------------------------------------------------- */

void DebugMessage(char Msg[])
{
	#if _MYDEBUG
	Message(Msg);
	#endif
}

/* ------------------------------------------------------------------------- */

int setCellAsChecked(int *cell_checked_arr, int cell_count, int i)
{
	int state = STATE_OK;

	if(i >= 0 && i < cell_count)
	{
		cell_checked_arr[i] = CELL_CHECKED;
	}
	else
	{
			Message("Warning index checked cells out of bounds");
			state = STATE_ERROR;
	}

	return state;
}

/* ------------------------------------------------------------------------- */

int copyCell_tArray(cell_t **toarr, cell_t *fromarr, int arrsize)
{
	int i;
	int state = STATE_OK;

	*toarr = (cell_t *) calloc(arrsize, sizeof(cell_t));

	if(*toarr != NULL) 
	{
		for(i = 0; i < arrsize; ++i)
		{
			(*toarr)[i] = fromarr[i];
		}
	}else 
	{
		DebugMessage("\nError: Kein freier Speicher für Copy Calloc vorhanden.\n");
		state = STATE_ERROR;
	}

	return state;
}


int copyRealND_ND_Array(real (*toarr)[ND_ND], real *fromarr, int arrsize)
{
	int i;
	int state = STATE_OK;

	for(i = 0; i < arrsize; ++i)
	{
		(*toarr)[i] = fromarr[i];
	}

	return state;
}


/* ------------------------------------------------------------------------- */


DEFINE_ON_DEMAND(DropletDetermination)
 {
	/*
	
	*/
	#if !RP_HOST
	cell_t c, cx, c0x, c1x;					/*Cell index*/
	Domain *domain =Get_Domain(1);			/*Fluid domain*/
	 
	int state = STATE_OK;
	Thread *ct;								/*Cell Thread pointer*/
	Thread *tf;								/*Face Thread pointer*/
	Thread **pt;
	face_t f;								/*Face index*/
	real x_c[ND_ND];
	real c_mass = 0;
	int n;									/*Local face index number*/
	int i;
	struct Droplet D;						/*Single droplet*/
	int dci;								/*Droplet cell index*/		
	int cell_count = 0;
	int *cell_checked_arr = NULL;			/*0=unchecked, 1=checked*/
	struct Droplet *DList = NULL;			/*Droplet List (dynamic allocation)*/
	int i_DList = 0;						/* Current droplet index in droplet list */
	int DList_length = 0;					/*Number droplets in droplet list */

	thread_loop_c(ct, domain)	
	{
	 	cell_count += THREAD_N_ELEMENTS_INT(ct);
	} 

	cell_checked_arr = (int*) calloc(cell_count,  sizeof(int));
	
	if(cell_checked_arr != NULL) 
	{}
	else
	{
		DebugMessage("Error (calloc): No free memory for cell_checked_arr available!\n");
		state = STATE_ERROR;
	}
	
	if(state == STATE_OK)
	{
		thread_loop_c(ct, domain)	
		{
			begin_c_loop(c, ct)
			{	
				pt = THREAD_SUB_THREADS(ct);
				/* find initial droplet cell */
				if (
					cell_checked_arr[c] != CELL_CHECKED 
					&& C_VOF(c, pt[_PHASE_IDX]) > _MIN_VOL_FRAC
					&& state == STATE_OK
				   )
				{
					state = initDroplet(&D, c);
					DList_length++;

					if (state == STATE_OK)
					{
						/* find connected droplet cells */
						for (dci = 0;  dci < D.no_cells; dci++)
						{
							cx = D.cell_list[dci];

							/* Update droplet properties */
							c_mass = C_VOF(cx, pt[_PHASE_IDX]) * C_VOLUME(cx, ct) 
								 * C_R(cx, pt[_PHASE_IDX]);
							D.mass += c_mass;
							D.vol += C_VOF(cx, pt[_PHASE_IDX]) * C_VOLUME(cx, ct);
							C_CENTROID(x_c, cx, ct);
							updateDropletWeights(&D, c_mass, x_c);
						
							if (cell_checked_arr[cx] != CELL_CHECKED)
							{	
								setCellAsChecked(cell_checked_arr, cell_count, cx);

								c_face_loop(cx, ct, n)
								{
									f = C_FACE(cx,ct,n);					
									tf = C_FACE_THREAD(cx,ct,n);

									if(!BOUNDARY_FACE_THREAD_P(tf))
									{
										c0x = F_C0(f, tf);
										c1x = F_C1(f, tf);

										if (
											c0x == cx 
											&& c1x >= 0 
											&& c1x < cell_count
										   )
										{
											pt = THREAD_SUB_THREADS(THREAD_T1(tf));
											if(C_VOF(c1x, pt[_PHASE_IDX]) > _MIN_VOL_FRAC)	
											{	
												state = dropletCellsAppend(&D, c1x);		
											}
											else
											{
												setCellAsChecked(cell_checked_arr, cell_count, c1x);
											}
										}
										else if(
												c1x == cx 
												&& c0x >= 0 
												&& c0x < cell_count
												)
										{
											pt = THREAD_SUB_THREADS(THREAD_T0(tf));
											if(C_VOF(c0x, pt[_PHASE_IDX]) > _MIN_VOL_FRAC)	
											{
												state = dropletCellsAppend(&D, c0x);					
											}
											else
											{
												setCellAsChecked(cell_checked_arr, cell_count, c0x);
											}
										}	
									}
								}
							}
						}		

						setFinalDropletCOM(&D);

						if(DList_length == 0)
						{
								DList = (struct Droplet*) calloc(DList_length,
														 sizeof(struct Droplet));
						}
						else
						{
							DList = ( struct Droplet*) realloc(DList, 
										(DList_length) * sizeof(struct Droplet));
						}
							
						if(DList != NULL) 
						{
							state = copyCell_tArray(&(DList[i_DList].cell_list),
														 D.cell_list, D.no_cells);
							state = copyRealND_ND_Array(&(DList[i_DList].com),
																 D.com, ND_ND);
							DList[i_DList].no_cells = D.no_cells;
							DList[i_DList].mass = D.mass;
							DList[i_DList].vol = D.vol;	
						}
						else 
						{
							DebugMessage("Error (calloc/realloc): No free memory "
										 "for DList available!\n");
							state = STATE_ERROR;
						}

						i_DList ++;
						free(D.cell_list);
					}
					else
					{
						setCellAsChecked(cell_checked_arr, cell_count, c);	
					}
				}			 
			}end_c_loop(c, ct)
		} 
	}
		

	Message("Droplets: %i\n", DList_length);
	printDList("mydroplets.txt", DList, DList_length);

	/*Free Memory*/
	if(DList_length > 0)
	{
		for(i=0; i<DList_length; ++i)
		{
			if(DList[i].cell_list != NULL)
			{
				free(DList[i].cell_list);
			}	
		}
		if(DList != NULL)
		{
			free(DList);
		} 
	}
	if(cell_checked_arr != NULL) 
	{
		free(cell_checked_arr);
	}
	#endif
 }
