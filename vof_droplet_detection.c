/*
UDF for continuos phase area detection (i.e. droplet detection) in Multiphase 
VoF (Volume of Fluid) ANSYS Fluent cases.

Settings:
    _PHASE_IDX : Phase index phase which should be detected
    _MIN_VOL_FRAC : Lower Limit for detection
    _DETECT_OVER_EDGES 1 : Detect connected areas over cell edges, 
                           if 0 only over cell faces 
    _FLUID_  1 : Id of Fluid Domain

WARNING: THIS IS AN EARLY VERSION THERE MAY BE INEXPECTED BUGS!

Current Limitations & Todos:
- Only double precision solver supported, due to file i/o

Ver: 0.5 (Christian Schubert)
- parallel working version (but only over face detection on parallel boundaries)
- several bug fixes

Copyright 2019-2020 Christian Schubert

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
#define _DETECT_OVER_EDGES 0

#define _FLUID_  1
/* ------------------------------------------------------------------------- */

#define _MYDEBUG 1
#define STATE_OK 1
#define STATE_ERROR 0
#define MIN_UDMI 1
#define UDMI_INT_TOL 1E-8
#define STRLENMAX 50

#define CELL_CHECKED 1
#define NO_MAX_CELL_EDGE_NEIGHBOR_CELLS 80
#define NO_MAX_CELL_FACE_NEIGHBOR_CELLS 30
/* ------------------------------------------------------------------------- */


int cfileexists(const char * filename){
/*
* Check if a file exist using fopen() function
* return 1 if the file exist otherwise return 0
*/
    FILE *file;
    if (file = fopen(filename, "r")){
        fclose(file);
        return 1;
    }
    return 0;
}

/* ------------------------------------------------------------------------- */

struct Cpa
{
    cell_t *cell_list;
    int    no_cells;
    real   com[ND_ND];                      /* center of mass (mass weighted average) */
    real   sumed_com_cell_weights[ND_ND];   /* intermediate values for com determination*/
    real   alpha_max;
    real   alpha_mean;
    real   mass;
    real   vol;
    int     no_parboundary_faces;
    face_t  *parboundary_faces_list;
    int    *boundary_id;
    int     no_boundaries;
    char    **parboundary_name_list;
    int     len_parboundary_name_list;
};

int initCpa(struct Cpa *d, cell_t c)	/* Initialize struct for single Cpa */	
{
    int state = STATE_OK;
    int i;

    (*d).no_cells = 1;
    (*d).mass = 0.0;
    (*d).vol  = 0.0;
    (*d).alpha_max = 0.0;
    (*d).alpha_mean  = 0.0;
    (*d).cell_list = NULL;
    (*d).cell_list = (cell_t*) calloc((*d).no_cells, sizeof(cell_t));
    (*d).cell_list[0] = c;

    for(i = 0; i < ND_ND; ++i)	/* ND_ND equals the number of dimensions: 2D -> ND_ND=2; 3D -> ND_ND=3 */	
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

    (*d).no_parboundary_faces = 0;
    (*d).parboundary_faces_list = NULL;
    (*d).parboundary_name_list = NULL;
    (*d).len_parboundary_name_list = 0;

    (*d).boundary_id  = NULL;
    (*d).no_boundaries = 0;

    return state;
}

int cpaCellsAppend(struct Cpa *d, cell_t val)
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


void updateCpaWeights(struct Cpa *d, real c_mass, real x_c[ND_ND])	/* Update center of mass */
{
    int i = 0;

    for(i = 0; i<ND_ND; ++i)
    {
        (*d).sumed_com_cell_weights[i] += x_c[i] * c_mass;
    }
}

void updateCpaBoundaryFaceID_List(struct Cpa *d, Thread *tf)
{
    int i = 0;
    int boundary_id_used = 0;
    int boundary_face_id = THREAD_ID(tf);


    if((*d).no_boundaries == 0)
    {
        (*d).boundary_id =  (int*) calloc(1, sizeof(int));
        (*d).boundary_id[0] = boundary_face_id;
        (*d).no_boundaries ++;
    }
    else
    {
        for(i = 0; i<(*d).no_boundaries; i++)
        {
            if(boundary_face_id == (*d).boundary_id[i])
            {
                boundary_id_used = 1;
                break;
            }
        }

        if(!boundary_id_used)
        {   
            (*d).no_boundaries ++;
            (*d).boundary_id =  (int*) realloc((*d).boundary_id, 
                                    (*d).no_boundaries*sizeof(int));
            (*d).boundary_id[(*d).no_boundaries-1] = boundary_face_id;
        }
    }
}


void updateCpaParBoundaryFaceID_List(struct Cpa *d,face_t newfaceid)
{
    int i = 0;
    int par_boundary_faceid_used = 0;


    if((*d).no_parboundary_faces == 0)
    {
        (*d).parboundary_faces_list =  (int*) calloc(1, sizeof(int));
        (*d).parboundary_faces_list[0] = newfaceid;
        (*d).no_parboundary_faces ++;
    }
    else
    {
        for(i = 0; i<(*d).no_parboundary_faces; i++)
        {
            if(newfaceid == (*d).parboundary_faces_list[i])
            {
                par_boundary_faceid_used = 1;
                break;
            }
        }

        if(!par_boundary_faceid_used)
        {   
            (*d).no_parboundary_faces ++;
            (*d).parboundary_faces_list =  (int*) realloc((*d).parboundary_faces_list, 
                                    (*d).no_parboundary_faces*sizeof(int));
            (*d).parboundary_faces_list[(*d).no_parboundary_faces-1] = newfaceid;
        }
    }
}

void updateParBoundaryNameList(struct Cpa *d, char pbname[])
{
    int i = 0;
    int cur_len;
    int par_boundary_name_used = 0;

    if((*d).len_parboundary_name_list == 0)
    {
        (*d).parboundary_name_list =  (char**) malloc(1 * sizeof(char*));

        (*d).parboundary_name_list[0] = (char*) malloc(STRLENMAX * sizeof(char));
        strcpy((*d).parboundary_name_list[0], pbname);
        (*d).len_parboundary_name_list ++;
    }
    else
    {
        for(i = 0; i<(*d).len_parboundary_name_list; ++i)
        {
            if(strcmp(pbname, (*d).parboundary_name_list[i]) == 0)
            {
                par_boundary_name_used = 1;
                break;
            } 
        }

        if(!par_boundary_name_used)
        {
            cur_len = (*d).len_parboundary_name_list;

            (*d).parboundary_name_list =  (char**) realloc((*d).parboundary_name_list, 
                                    (cur_len+1)*sizeof(char*));

            (*d).parboundary_name_list[cur_len] = (char*) malloc(STRLENMAX * sizeof(char));
            strcpy((*d).parboundary_name_list[cur_len], pbname);
            (*d).len_parboundary_name_list ++; 
        }
    }
}


void setFinalCpaCOM(struct Cpa *d)
{
    int i = 0;

    for(i = 0; i<ND_ND; ++i)
    {
        (*d).com[i] = (*d).sumed_com_cell_weights[i] / (*d).mass;
    }
}


int printCpaCells(char filesuffix[], struct Cpa *DList, int sizeDList, Thread *ct)
{
    FILE *fd = NULL;
    int state = STATE_OK;
    int i,j;
    char filename[100];

    sprintf(filename, "%i_cells_%s", myid, filesuffix);

    fd = fopen(filename, "a");

    if(fd == NULL)
    {
        Message("Error (printCpas()): Unable to open file %s "
                "for writing!\n", filename);
        state = STATE_ERROR;
    }
    else
    {
        fprintf(fd, "{\n");
        fprintf(fd, "%lf\n", CURRENT_TIME);
        fprintf(fd,"cell_count[cells]\n");
        for(i = 0; i < sizeDList; ++i)
        {
            fprintf(fd, "%i[", DList[i].no_cells);
            for(j = 0; j <DList[i].no_cells; ++j)
            {
                fprintf(fd, "%li", (long) C_ID(DList[i].cell_list[j],ct));
                if(j < (DList[i].no_cells -1))
                {
                    fprintf(fd, ", ");
                }
            }
            fprintf(fd, "]]\n");
        }
        fprintf(fd, "}\n");
    }
    
    return state;
}

int printCpas(char filesuffix[], struct Cpa *DList, int sizeDList)
{
    FILE *fd = NULL;
    int state = STATE_OK;
    int i,j;
    char filename[100];

    sprintf(filename, "%i_%s", myid, filesuffix);

    fd = fopen(filename, "a");

    if(fd == NULL)
    {
        Message("Error (printCpas()): Unable to open file %s "
                "for writing!\n", filename);
        state = STATE_ERROR;
    }
    else
    {
        fprintf(fd, "{\n");
        fprintf(fd, "%lf\n", CURRENT_TIME);
        fprintf(fd,"cell_count[com.(x)_in_m,com.(y),com.(z),mass_in_kg,"
                    "volume_in_m³,cellAveragedAlpha,maxAlpha,"
                    "[list_of_boundary_names],[list_of_processor_cells]]\n");

        for(i = 0; i < sizeDList; ++i)
        {
            fprintf(fd, "%i[", DList[i].no_cells);
            for(j = 0; j <ND_ND; ++j)
            {
                fprintf(fd, " %lf,", DList[i].com[j]);
            }
            fprintf(fd, "%lE,", DList[i].mass);
            fprintf(fd, "%lE,", DList[i].vol);
            fprintf(fd, "%lE,", DList[i].alpha_mean);
            fprintf(fd, "%lE,", DList[i].alpha_max);

            fprintf(fd, "[");
            if(DList[i].no_boundaries > 0 || DList[i].len_parboundary_name_list >0)
            {
                for(j = 0; j < DList[i].no_boundaries; ++j)
                {
                    fprintf(fd, "%i", DList[i].boundary_id[j]);

                    if(j < (DList[i].no_boundaries -1) || DList[i].len_parboundary_name_list>0)
                    {
                        fprintf(fd, ",");
                    }
                }

                for(j = 0; j < DList[i].len_parboundary_name_list; j++)
                {
                    fprintf(fd, "%s", DList[i].parboundary_name_list[j]);

                    if(j < (DList[i].len_parboundary_name_list-1))
                    {
                        fprintf(fd, ",");
                    }
                }
                fprintf(fd, "],");
            }
            else
            {
                fprintf(fd, "],");
            }

            fprintf(fd, "%i[", DList[i].no_parboundary_faces);

            if(DList[i].no_parboundary_faces > 0)
            {
                for(j = 0; j < (DList[i].no_parboundary_faces-1); ++j)
                {
                    fprintf(fd, "%i,", DList[i].parboundary_faces_list[j]);
                }
                fprintf(fd, "%i]]\n",DList[i].parboundary_faces_list[DList[i].no_parboundary_faces-1]);
            }
            else
            {
                fprintf(fd, "]]\n");
            }

        }
        fprintf(fd, "}\n");
    }

    if(fd != NULL)
    {
        fclose(fd);
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
        Message("Warning setCellAsChecked(): Index %i for checked cells range from 0 to %i out of bounds\n", i, cell_count);
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


int copyIntArray(int **toarr, int *fromarr, int arrsize)
{
    int i;
    int state = STATE_OK;

    *toarr = (int *) calloc(arrsize, sizeof(int));

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


int copyCharArrays(char ***toarr, char **fromarr, int arrsize, int chararrsize)
{
    int i;
    int state = STATE_OK;

    *toarr = (char **) malloc(arrsize * sizeof(char*));

    if(*toarr != NULL) 
    {
        for(i = 0; i < arrsize; ++i)
        {
            (*toarr)[i] = (char*) malloc(chararrsize * sizeof(char));
            strcpy((*toarr)[i], fromarr[i]);
        }
    }else 
    {
        DebugMessage("\nError: Kein freier Speicher für Copy Calloc vorhanden.\n");
        state = STATE_ERROR;
    }

    return state;
}

/* ------------------------------------------------------------------------- */

void resetIntArray(int *arr, int size, int val)
{
    int i = 0;

    for(i=0; i<size; ++i)
    {
        arr[i] = val;
    }
}


void getCellsFaceNeighborCells(
                                cell_t c0, 
                                Thread *c0t,
                                cell_t nocells_c0t, 
                                cell_t *cfncarr, 
                                int *cfncarr_size
                                )
{
    int n;
    int i=0;
    face_t f;
    Thread *tf;
    cell_t c0x, c1x;
    *cfncarr_size = 0;

    c_face_loop(c0, c0t, n)
    {
        f = C_FACE(c0,c0t,n);            /*return thread global face index*/
        tf = C_FACE_THREAD(c0,c0t,n);

        if(tf != NULL && !BOUNDARY_FACE_THREAD_P(tf))
        {
            c0x = F_C0(f, tf);
            c1x = F_C1(f, tf);
            
            if(i < NO_MAX_CELL_FACE_NEIGHBOR_CELLS)
            {
                if (c0x == c0)
                {   
                    if(c1x > 0 && c1x < nocells_c0t)
                    {
                        cfncarr[i] = c1x;
                        i += 1;
                    }
                }
                else if(c1x == c0)
                {
                    if(c0x > 0 && c0x < nocells_c0t)
                    {
                        cfncarr[i] = c0x;
                        i += 1;
                    }
                }	
            }
            else
            {
                Message("Error getCellsFaceNeighborCells(): Face index going to high!\n");
            } 
        }
    }

    *cfncarr_size = i;
}


void array1ElementFrequencyInArray2(int *arr1, int size1, int *arr2, int size2, int *freq_arr)
{
    int i;

    for(i=0; i<size1; ++i)
    {   
        freq_arr[i] = countValinIntArray(&(*arr2), size2, arr1[i]);
    }
}

int countValinIntArray(int *arr, int arr_size, int val)
{
    int i;
    int count = 0;

    for(i=0; i<arr_size; ++i)
    {
        if(arr[i] == val)
        {
            count ++;
        }
    }

    return count;
}


int occursInCellArray(cell_t *arr, int arr_size, cell_t val)
{
    int i;

    for(i=0; i<arr_size; ++i)
    {
        if(arr[i] == val)
        {
            return 1;
        }
    }

    return 0;
}


void getCellsEdgeNeighborCells(
                                cell_t c0,
                                Thread *c0t,
                                cell_t nocells_c0t,  
                                cell_t *cencarr,
                                int *cfnarr_size,
                                int *cencarr_size
                                )
{
    int i,j;
    cell_t ci;

    *cencarr_size = 0;

    cell_t c0_cfncarr[NO_MAX_CELL_FACE_NEIGHBOR_CELLS];
    cell_t ci_cfncarr[NO_MAX_CELL_FACE_NEIGHBOR_CELLS];
    int freq_c0_in_ci_cfncarr[NO_MAX_CELL_FACE_NEIGHBOR_CELLS];
    int c0_cfncarr_size = 0;
    int ci_cfncarr_size = 0;

    getCellsFaceNeighborCells(c0, c0t, nocells_c0t, c0_cfncarr, &c0_cfncarr_size);

    for(i=0; i<c0_cfncarr_size; ++i)
    {   
        if(*cencarr_size < NO_MAX_CELL_EDGE_NEIGHBOR_CELLS)
        {
            cencarr[*cencarr_size] = c0_cfncarr[i];
            *cencarr_size += 1;
        }
        else
        {
            Message("Error getCellsEdgeNeighborCells(): Edge index going to high 1!\n");
        }
    }

    *cfnarr_size = c0_cfncarr_size;

    for(i=0; i<c0_cfncarr_size; ++i)
    {   
        ci = c0_cfncarr[i];
        getCellsFaceNeighborCells(ci, c0t, nocells_c0t, ci_cfncarr, &ci_cfncarr_size);

        array1ElementFrequencyInArray2(ci_cfncarr, ci_cfncarr_size, c0_cfncarr, c0_cfncarr_size, freq_c0_in_ci_cfncarr);

        for(j=0; j<ci_cfncarr_size; ++j)
        {
            if(freq_c0_in_ci_cfncarr[j] >= 2)
            {
                if(*cencarr_size < NO_MAX_CELL_EDGE_NEIGHBOR_CELLS)
                {
                    cencarr[*cencarr_size] = ci_cfncarr[j];
                    *cencarr_size += 1;
                }
                else
                {
                    Message("Error getCellsEdgeNeighborCells(): Edge index going to high 2!\n");
                }
            }
        }
    }
}


face_t getGlobalFaceBetweenCells(cell_t c0, cell_t c1, Thread *ct)
{
    cell_t c0f, c1f;
    face_t f;
    Thread *tf;
    int n;

    c_face_loop(c0, ct, n)
    {
        f = C_FACE(c0,ct,n);          
        tf = C_FACE_THREAD(c0,ct,n);

        if(tf != NULL)
        {
            c0f = F_C0(f, tf);
            c1f= F_C1(f, tf);

            if ((c0f == c0 && c1f == c1) || (c0f == c1 && c1f == c0))
            {
                    return F_ID(f,tf);
            }
        }
    }

    return -1;
}
/*----------------------------------------------------------------------------*/

void cpa_detection()
{
/*
    
*/
    #if !RP_HOST
    cell_t c, cx, ci;                 /*Cell index*/
    Domain *domain =Get_Domain(1);          /*Fluid domain*/

    int state = STATE_OK;
    Thread *ct;                             /*Cell Thread pointer*/
    Thread *tf;                             /*Face Thread pointer*/
    Thread **pt;
    face_t f, fid;
    int fluid_IDs[] = {_FLUID_};                
    real x_c[ND_ND];
    real c_mass = 0;                               /*Local face index number*/
    int i,j,n;
    int cci;
    struct Cpa D;                       /*Single Cpa*/
    int dci;                                /*Cpa cell index*/	
    int nocells_ct = 0;
    int cell_count = 0;
    int *cell_checked_arr = NULL;           /*0=unchecked, 1=checked*/
    struct Cpa *DList = NULL;           /*Cpa List (dynamic allocation)*/
    int i_DList = 0;                        /* Current Cpa index in Cpa list */
    int DList_length = 0;                   /*Number Cpas in Cpa list */
    int c0fnc_array_size;
    char pbname[STRLENMAX];

    #if _DETECT_OVER_EDGES
    cell_t c0nc_array[NO_MAX_CELL_EDGE_NEIGHBOR_CELLS];  /*neighbor cells of c0*/
    #else
    cell_t c0nc_array[NO_MAX_CELL_FACE_NEIGHBOR_CELLS];
    #endif    

    int c0nc_array_size = 0;

    int no_fluid_IDs = sizeof(fluid_IDs)/sizeof(fluid_IDs[0]);
    

    /* Count domain cells and initialize UDMI */
    cell_count = 0;
    for (i = 0; i < no_fluid_IDs; i++)
    {
        ct = NULL;
        ct = Lookup_Thread(domain, fluid_IDs[i]);

        if(ct != NULL)
        {
            cell_count += THREAD_N_ELEMENTS_INT(ct);
            cell_count += THREAD_N_ELEMENTS_EXT(ct);
        }
        else
        {
            Message("Error CpaDetermination(): Cannot find fluid cell thread with id %i", fluid_IDs[i]);
        }
        
    }

    cell_checked_arr = (int*) calloc(cell_count,  sizeof(int));
    
    if(cell_checked_arr != NULL || cell_count == 0) 
    {
        Message("Found %i cells to check in myid: %i\n", cell_count, myid);
    }
    else
    {
        DebugMessage("Error (calloc): No free memory for cell_checked_arr available!\n");
        state = STATE_ERROR;
    }
    
    if(state == STATE_OK &&  N_UDM >= MIN_UDMI)
    {
        for (i = 0; i < no_fluid_IDs; i++)
        {
            Message("Looking for Droplets in domain with id %i, myid %i\n", fluid_IDs[i], myid);

            ct = NULL;
            ct = Lookup_Thread(domain, fluid_IDs[i]);

            if(ct != NULL)
            {   
                nocells_ct = THREAD_N_ELEMENTS_INT(ct);
                nocells_ct += THREAD_N_ELEMENTS_EXT(ct);

                begin_c_loop_int(c, ct)
                {
                    pt = THREAD_SUB_THREADS(ct);

                    if(pt == NULL)
                    {
                        Message("Error pt == NULL!\n");
                    }

                    if(c >= cell_count)
                    {
                        Message("Error c > cell_count!\n");
                    }

                    /* find initial droplet cell */
                    if (
                        pt != NULL
                        && (cell_checked_arr[c] != CELL_CHECKED) 
                        && (C_VOF(c, pt[_PHASE_IDX]) > _MIN_VOL_FRAC)
                        && (state == STATE_OK)
                    )
                    {
                        state = initCpa(&D, c); /*Initital cpa (droplet)*/
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

                                D.alpha_mean = D.alpha_mean + C_VOF(cx, pt[_PHASE_IDX]);

                                if(D.alpha_max < C_VOF(cx, pt[_PHASE_IDX]))
                                {
                                    D.alpha_max = C_VOF(cx, pt[_PHASE_IDX]);
                                }

                                updateCpaWeights(&D, c_mass, x_c);
                            
                                c_face_loop(cx, ct, n)
                                {
                                    f = C_FACE(cx,ct,n); /*return thread global face index*/
                                    tf = C_FACE_THREAD(cx,ct,n); 
                                    
                                    if(tf != NULL && BOUNDARY_FACE_THREAD_P(tf))
                                    {
                                        updateCpaBoundaryFaceID_List(&D, tf);
                                    }
                                }

                                if (C_UDMI(cx,ct,0) > UDMI_INT_TOL) 
                                {
                                    /*Partition Boundary interior cell*/
                                    /*add cell boundary id to boundary list ...
                                    currently only working for cell neighbors accross
                                    compute node boundaries ...*/
                                    c0nc_array_size = 0;
                                    getCellsFaceNeighborCells(cx, ct, nocells_ct, c0nc_array, &c0nc_array_size);

                                    for(cci=0; cci<c0nc_array_size; ++cci)
                                    {
                                        ci = c0nc_array[cci];

                                        if (C_PART(ci, ct) != myid && C_VOF(ci, pt[_PHASE_IDX]) > _MIN_VOL_FRAC)
                                        {
                                            fid = getGlobalFaceBetweenCells(cx, ci, ct);

                                            if (fid != -1)
                                            {
                                                updateCpaParBoundaryFaceID_List(&D, fid);

                                                sprintf(pbname, "procBoundary%ito%li", myid, (long) C_PART(ci, ct)); 
                                                updateParBoundaryNameList(&D, pbname);
                                            }
                                        }
                                    }
                                }

                                c0nc_array_size = 0;

                                #if _DETECT_OVER_EDGES
                                getCellsEdgeNeighborCells(cx, ct, nocells_ct, c0nc_array, &c0nc_array_size, &c0fnc_array_size);
                                #else
                                getCellsFaceNeighborCells(cx, ct, nocells_ct, c0nc_array, &c0nc_array_size);
                                #endif    

                                for(cci=0; cci<c0nc_array_size; ++cci)
                                {
                                    ci = c0nc_array[cci];
                                    pt = THREAD_SUB_THREADS(ct);

                                    if( (cell_checked_arr[ci] != CELL_CHECKED)
                                        && (C_PART(ci, ct) == myid)
                                        && (C_VOF(ci, pt[_PHASE_IDX]) > _MIN_VOL_FRAC
                                        && (ci != cx)
                                        )
                                    )
                                    {
                                        state = cpaCellsAppend(&D, ci);
                                        setCellAsChecked(cell_checked_arr, cell_count, ci);
                                    }
                                    else
                                    {
                                        setCellAsChecked(cell_checked_arr, cell_count, ci);
                                    }
                                }

                                setCellAsChecked(cell_checked_arr, cell_count, cx);
                            }

                            D.alpha_mean = D.alpha_mean/D.no_cells;

                            setFinalCpaCOM(&D);

                            if(DList_length == 0)
                            {
                                DList = (struct Cpa*) calloc(DList_length, 
                                                        sizeof(struct Cpa));
                            }
                            else
                            {
                                DList = (struct Cpa*) realloc(DList, 
                                        (DList_length) * sizeof(struct Cpa));  
                            }
                                
                            if(DList != NULL) 
                            {
                                state = copyCell_tArray(&(DList[i_DList].cell_list),
                                            D.cell_list, D.no_cells);
                                state = copyCell_tArray(&(DList[i_DList].parboundary_faces_list),
                                            D.parboundary_faces_list, D.no_parboundary_faces);
                                state = copyRealND_ND_Array(&(DList[i_DList].com), 
                                            D.com, ND_ND);
                                state = copyIntArray(&(DList[i_DList].boundary_id),
                                            D.boundary_id, D.no_boundaries);
                                state = copyCharArrays(&(DList[i_DList].parboundary_name_list),
                                            D.parboundary_name_list, D.len_parboundary_name_list, STRLENMAX);

                                DList[i_DList].no_boundaries = D.no_boundaries;       
                                DList[i_DList].no_cells = D.no_cells;
                                DList[i_DList].no_parboundary_faces = D.no_parboundary_faces;
                                DList[i_DList].mass = D.mass;
                                DList[i_DList].vol = D.vol;
                                DList[i_DList].alpha_mean = D.alpha_mean;
                                DList[i_DList].alpha_max = D.alpha_max;
                                DList[i_DList].len_parboundary_name_list = D.len_parboundary_name_list;
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
                }end_c_loop_int(c, ct)    
            }
            else
            {
                Message("Error DropletDetermination(): Cannot find fluid cell thread with id %i", fluid_IDs[i]);
            }
        } 
        Message("Found %i cpas in myid %i.\n", DList_length, myid);
        printCpas("cpa.txt", DList, DList_length);
        /*printCpaCells("cpa.txt", DList, DList_length, ct); */ /*For debug only*/
    }
    
    /*Free Memory*/
    if(DList_length > 0)
    {
        for(i=0; i<DList_length; ++i)
        {
            if(DList[i].cell_list != NULL)
            {
                free(DList[i].cell_list);
            }

            if(DList[i].parboundary_faces_list != NULL)
            {
                free(DList[i].parboundary_faces_list);
            }

            for(j=0;j<DList[i].len_parboundary_name_list;++j)
            {
                if(DList[i].parboundary_name_list[j] != NULL)
                {
                    free(DList[i].parboundary_name_list[j]);
                }
            }
            
            if(DList[i].parboundary_name_list != NULL)
            {
                free(DList[i].parboundary_name_list);
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

/*----------------------------------------------------------------------------*/

DEFINE_ON_DEMAND(CPAD_oD)
{
#if !RP_HOST
    cpa_detection();
#endif
}


DEFINE_EXECUTE_AT_END(CPAD_aE)
{
#if !RP_HOST
    cpa_detection();
#endif
}


DEFINE_ON_DEMAND(MARKPAR_oD)
{
#if !RP_HOST
    markParBoundCells();
#endif
}


void markParBoundCells()
{

    #if !RP_HOST
    int n, nocells_ct;
    cell_t c, c1, c0;                           /*Cell index*/
    Domain *domain =Get_Domain(1);          /*Fluid domain*/
    Thread *ct;                             /*Cell Thread pointer*/
    Thread *tf;                             /*Face Thread pointer*/
    face_t f;                               /*Face index*/        

    if(N_UDM > 0)
    {
        thread_loop_c(ct, domain) /*loops over all cell threads in domain*/
        {                        
            if(ct != NULL)
            {   
                nocells_ct = THREAD_N_ELEMENTS_EEXT(ct);
                nocells_ct += THREAD_N_ELEMENTS_INT(ct);

                begin_c_loop_ext(c, ct)
                {
                    C_UDMI(c,ct,0) = 0.0;
                }end_c_loop_ext(c, ct)

                begin_c_loop_ext(c, ct)
                {
                    c_face_loop(c, ct, n)
                    {
                        f = C_FACE(c,ct,n);            /*return thread global face index*/
                        tf = C_FACE_THREAD(c,ct,n);

                        if(tf != NULL)
                        {
                            c0 = F_C0(f, tf);
                            c1= F_C1(f, tf);
                            
                            if (c0 == c)
                            {   
                                if(c1> 0 && c1< nocells_ct)
                                {
                                    if(C_PART(c1,ct) == myid)
                                    {   
                                        C_UDMI(c1,ct,0) = myid+1;
                                    }
                                }
                            }
                            else if(c1 == c)
                            {
                                if(c0 > 0 && c0 < nocells_ct)
                                {
                                    if(C_PART(c0,ct) == myid)
                                    {   
                                        C_UDMI(c0,ct,0) = myid+1;
                                    }
                                }
                            }	
                        }
                    }
                }end_c_loop_ext(c, ct)
            }
        }
    }
    else
    {
        Message("Not enough UDMI Memory!\n");
    }
    #endif
}
