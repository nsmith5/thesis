/*
 * This file impliments I/0 using HDF5
 *
 * Just a warning: this is a hairy one. Read the HDF5 docs
 * for help on the details.
 */


#include <mpi.h>        // MPI runtime
#include <hdf5.h>       // HDF5 header
#include <stdlib.h>     // malloc etc
#include <unistd.h>     // access() to check if file exists
#include <stdbool.h>    // Booleans
#include <assert.h>     // Assertions
#include <string.h>     // String manipulation sprintf etc...

#include "binary.h"     // local libbinary header

#define FAIL -1
#define LEN 8

const int io_verbose = true;

void mpi_print (const char *str)
{
    int rank;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    if (rank == 0) printf("%s\n", str);
    MPI_Barrier (MPI_COMM_WORLD);
}

hid_t io_init_from_file (const char *filename)
{
    /* Open HDF5 file and pass back file ID */

    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    hid_t plist_id;
    hid_t file_id;
    herr_t status;

    plist_id = H5Pcreate (H5P_FILE_ACCESS);

    H5Pset_fapl_mpio (plist_id, comm, info);

    file_id = H5Fopen (filename, H5F_ACC_RDWR, plist_id);

    status = H5Pclose (plist_id);
    assert (status != FAIL);

    return file_id;
}


hid_t io_init_new_file (const char *filename)
{
    /*
     *  Create a file to save data to for this session
     */
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    hid_t plist_id;     /* Property list id */
    hid_t file_id;      /* File id */
    herr_t err;         /* Error status */

    plist_id = H5Pcreate (H5P_FILE_ACCESS);

    H5Pset_fapl_mpio (plist_id, comm, info);

    file_id = H5Fcreate (filename,
                         H5F_ACC_TRUNC,
                         H5P_DEFAULT,
                         plist_id );

    err = H5Pclose (plist_id);
    assert (err != FAIL);
    return file_id;
}

herr_t io_finalize (hid_t file_id)
{
    /* Close HDF5 file */
    return H5Fclose(file_id);
}

herr_t write_array_dataset (const char *name,
                            hid_t       group_id,
                            double     *arr,
                            state      *s)
{
     /* 
      * Write an N x N double array *arr* to
      * group with id *group_id*. Name it *name*.
      */
     hid_t dset_id, dataspace;
     hid_t memspace, plist_id;
     herr_t status;
     hsize_t size[2], count[2], offset[2];

     /* Describe the shape of the array */
     size[0] = s->N;
     size[1] = s->N;
     dataspace = H5Screate_simple(2, size, NULL);


      /* Create Dataset */
     dset_id = H5Dcreate(group_id,
                         name,
                         H5T_NATIVE_DOUBLE,
                         dataspace,
                         H5P_DEFAULT,
                         H5P_DEFAULT,
                         H5P_DEFAULT);

     /* Describe memory shape, property list and data shape */
     count[0] = 1;
     count[1] = s->N;
     offset[1] = 0;
     memspace = H5Screate_simple (2, count, NULL);

     /* Set up some of the MPI things */
     plist_id = H5Pcreate (H5P_DATASET_XFER);
     H5Pset_dxpl_mpio (plist_id, H5FD_MPIO_INDEPENDENT);

     /* Write data row by row in slabs */
     for (int row = 0; row < s->local_n0; row++)
     {
         offset[0] = s->local_0_start + row;

         status = H5Sselect_hyperslab (dataspace,
                                       H5S_SELECT_SET,
                                       offset,
                                       NULL,
                                       count,
                                       NULL);

         status = H5Dwrite (dset_id,
                            H5T_NATIVE_DOUBLE,
                            memspace,
                            dataspace,
                            plist_id,
                            arr + row*2*((s->N>>1) + 1));
     }

     /* Close everything you opened */
     status = H5Pclose (plist_id);
     status = H5Sclose (memspace);
     status = H5Dclose (dset_id);
     status = H5Sclose (dataspace);

     return status;
}

herr_t read_array_dataset (const char *name,
                           hid_t       group_id,
                           double     *arr,
                           state      *s)
{
    /* Opposite of *write_array_dataset* */
    hid_t dset_id, dataspace;
    hid_t memspace, plist_id;
    herr_t status;
    hsize_t count[2], offset[2];

    /* Open Dataset */
    dset_id = H5Dopen1 (group_id, name);

    /* Describe memory shape, property list and data shape */
    count[0] = 1;
    count[1] = s->N;
    offset[1] = 0;
    memspace = H5Screate_simple (2, count, NULL);

    /* Set up some of the MPI things */
    dataspace = H5Dget_space (dset_id);
    plist_id = H5Pcreate (H5P_DATASET_XFER);
    H5Pset_dxpl_mpio (plist_id, H5FD_MPIO_INDEPENDENT);

    /* Write data row by row in slabs */
    for (int row = 0; row < s->local_n0; row++)
    {
        offset[0] = s->local_0_start + row;

        status = H5Sselect_hyperslab (dataspace,
                                      H5S_SELECT_SET,
                                      offset,
                                      NULL,
                                      count,
                                      NULL);

        status = H5Dread (dset_id,
                          H5T_NATIVE_DOUBLE,
                          memspace,
                          dataspace,
                          plist_id,
                          arr + row*2*(s->N/2 + 1));
    }

    /* Close everything you opened */
    status = H5Pclose (plist_id);
    status = H5Sclose (memspace);
    status = H5Dclose (dset_id);
    status = H5Sclose (dataspace);

    return status;
}

herr_t write_double_attribute (const char *name,
                               hid_t       group_id,
                               double     *value)
{
     /* 
      * Write a double attribute to the group *group_id*
      * under the name *name*
      */
     hsize_t size = 1;
     herr_t status;
     hid_t attr_id, dataspace;

     dataspace = H5Screate_simple(1, &size, NULL);
     attr_id = H5Acreate2 (group_id,
                           name,
                           H5T_NATIVE_DOUBLE,
                           dataspace,
                           H5P_DEFAULT,
                           H5P_DEFAULT);

     status = H5Awrite (attr_id, H5T_NATIVE_DOUBLE, value);
     status = H5Aclose (attr_id);

     return status;
}

herr_t read_double_attribute (const char *name,
                              hid_t       group_id,
                              double     *value)
{
    /* Read double attribute from dataset 'group_id' */
    hid_t attr_id;
    herr_t status;

    attr_id = H5Aopen (group_id, name, H5P_DEFAULT);
    status = H5Aread (attr_id, H5T_NATIVE_DOUBLE, value);

    status = H5Aclose (attr_id);
    return status;
}

herr_t write_int_attribute (const char *name,
                            hid_t       group_id,
                            int        *value)
{
     /* ..... you get the idea... */
     hsize_t size = 1;
     herr_t status;
     hid_t attr_id, dataspace;

     dataspace = H5Screate_simple (1, &size, NULL);
     attr_id = H5Acreate2 (group_id,
                           name,
                           H5T_NATIVE_INT,
                           dataspace,
                           H5P_DEFAULT,
                           H5P_DEFAULT);

     status = H5Awrite (attr_id, H5T_NATIVE_INT, value);
     status = H5Aclose (attr_id);

     return status;
}

herr_t read_int_attribute (const char *name,
                           hid_t       group_id,
                           int        *value)
{
    /* Read integer attribute from dataset 'group_id' */
    hid_t attr_id;
    herr_t status;

    attr_id = H5Aopen (group_id, name, H5P_DEFAULT);
    status = H5Aread (attr_id, H5T_NATIVE_INT, value);

    status = H5Aclose (attr_id);

    return status;
}

herr_t save_state (state *s,
                   hid_t  file_id)
{
     /* Save entire state at a time step */
     hid_t group_id;
     herr_t status;
     int N_int;

     /* Make Group from simulation time `t` */
     char groupname[50];
     char step_str[10];
     sprintf(step_str, "%d", s->step);
     sprintf(groupname, "%0*d%s", LEN-(int)strlen(step_str), 0, step_str);

     group_id = H5Gcreate (file_id,
                           groupname,
                           H5P_DEFAULT,
                           H5P_DEFAULT,
                           H5P_DEFAULT);

     /* Write state attributes */

     status = write_double_attribute ("eta",    group_id, &s->eta);
     status = write_double_attribute ("chi",    group_id, &s->chi);
     status = write_double_attribute ("epsilon_0", group_id, &s->epsilon0);
     status = write_double_attribute ("sigma0", group_id, &s->sigma0);
     status = write_double_attribute ("sigma",  group_id, &s->sigma);
     status = write_double_attribute ("omega",  group_id, &s->omega);
     status = write_double_attribute ("Wc",     group_id, &s->Wc);

     status = write_double_attribute ("Mn",     group_id, &s->Mn);
     status = write_double_attribute ("Mc",     group_id, &s->Mc);

     status = write_double_attribute ("k0",     group_id, &s->k0);
     status = write_double_attribute ("alpha",     group_id, &s->alpha);
     status = write_double_attribute ("beta",     group_id, &s->beta);
     status = write_double_attribute ("rho",     group_id, &s->rho);
     status = write_double_attribute ("alphac",     group_id, &s->alphac);

     status = write_double_attribute ("dx",     group_id, &s->dx);
     status = write_double_attribute ("dt",     group_id, &s->dt);
     status = write_double_attribute ("Time",   group_id, &s->t);
     status = write_int_attribute ("Time Step", group_id, &s->step);

     N_int = (int)s->N;
     status = write_int_attribute ("N", group_id, &N_int);

     status = write_array_dataset ("Concentration", group_id, s->c->real, s);
     status = write_array_dataset ("Density",   group_id, s->n->real, s);

     status = H5Gclose (group_id);

     return status;
}

state* load_state (hid_t       file_id,
                   const char *datafile)
{
    state* s;
    ptrdiff_t N;
    int N_int;
    double dx, dt;
    hid_t group_id;

    group_id = H5Gopen2 (file_id, datafile, H5P_DEFAULT);

    read_int_attribute ("N", group_id, &N_int);
    read_double_attribute ("dx", group_id, &dx);
    read_double_attribute ("dt", group_id, &dt);

    N = (ptrdiff_t)N_int;

    s = create_state (N, dx, dt);

    /* Read state attributes */

    read_double_attribute ("eta",      group_id, &s->eta);
    read_double_attribute ("chi",      group_id, &s->chi);
    read_double_attribute ("epsilon_0",group_id, &s->epsilon0);
    read_double_attribute ("sigma0",   group_id, &s->sigma0);
    read_double_attribute ("sigma",    group_id, &s->sigma);
    read_double_attribute ("omega",    group_id, &s->omega);
    read_double_attribute ("Wc",       group_id, &s->Wc);

    read_double_attribute ("Mn", group_id, &s->Mn);
    read_double_attribute ("Mc", group_id, &s->Mc);

    read_double_attribute ("k0",     group_id, &s->k0);
    read_double_attribute ("alpha",  group_id, &s->alpha);
    read_double_attribute ("beta",   group_id, &s->beta);
    read_double_attribute ("rho",    group_id, &s->rho);
    read_double_attribute ("alphac", group_id, &s->alphac);

    read_double_attribute ("Time",   group_id, &s->t);
    read_int_attribute ("Time Step", group_id, &s->step);

    read_array_dataset ("Concentration", group_id, s->c->real, s);
    read_array_dataset ("Density",   group_id, s->n->real, s);

    set_C (s);
    set_propagators (s);

    H5Gclose (group_id);

    return s;
}

state* new_state_from_file (const char *filename)
{
    /* Read state out of a YAML file *filename* */
    FILE *ifp;
    char *mode = "r";
    state* s;
    char name[20];
    double value;
    int N;
    double dt, dx;
    double n0, c0;

    /* Setinal Values */
    dt = dx = c0 = n0 = -1.0;
    N = -1;

    // Open up the file
    ifp = fopen (filename, mode);
    assert (ifp != NULL);

    while (fscanf(ifp, "%s : %lf", name, &value) != EOF)
    {
        // Scan for the numerical details
        if (strcmp (name, "N") == 0)        N = (int)value;
        else if (strcmp(name, "dx") == 0)   dx = value;
        else if (strcmp(name, "dt") == 0)   dt = value;
    }

    if (N == -1 || dx == -1.0 || dt == -1.0)
    {
        my_error("Failed to load numerical values from YAML file");
    }

    // Allocate state with *ZERO* assertion that we
    // have legitimate values for N, dx, or dt.
    s = create_state (N, dx, dt);
    assert (s != NULL);
    
    // Rewind to the start of the file
    rewind (ifp);

    while (fscanf(ifp, "%s : %lf", name, &value) != EOF)
    {
        // Scan for all other variables from the YAML file
        if (strcmp (name, "eta") == 0)          s->eta = value;
        else if (strcmp (name, "Mc") == 0)      s->Mc = value;
        else if (strcmp (name, "Mn") == 0)      s->Mn = value;
        else if (strcmp (name, "chi") == 0)     s->chi = value;
        else if (strcmp (name, "epsilon0") == 0) s->epsilon0 = value;
        else if (strcmp (name, "sigma0") == 0)  s->sigma0 = value;
        else if (strcmp (name, "sigma") == 0)   s->sigma = value;
        else if (strcmp (name, "omega") == 0)   s->omega = value;
        else if (strcmp (name, "Wc") == 0)      s->Wc = value;
        else if (strcmp (name, "k0") == 0)      s->k0 = value;
        else if (strcmp (name, "alpha") == 0)   s->alpha = value;
        else if (strcmp (name, "beta") == 0)    s->beta = value;
        else if (strcmp (name, "rho") == 0)     s->rho = value;
        else if (strcmp (name, "alphac") == 0)  s->alphac = value;
        else if (strcmp (name, "n0") == 0)      n0 = value;
        else if (strcmp (name, "c0") == 0)      c0 = value;
    }

    if (n0 == -1.0 || c0 == -1.0)
    {
        my_error ("Failed to load initial density n0 or initial concentration c0 from YAML file");
    }

    // Initialize the fields as uniform from n0 and c0 parameters
    set_uniform (s->n->real, n0, s->local_n0, N);
    set_uniform (s->c->real, c0, s->local_n0, N);

    s->step = 0;
    s->t = 0.0;
    set_C (s);
    set_propagators (s);

    // Ultra dangerously return the state without any checking that
    // this process actually worked
    return s;
}
