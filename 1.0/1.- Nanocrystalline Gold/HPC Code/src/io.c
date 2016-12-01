/*
 * This file is part of Diffusion Equation MPI
 * Nathan Smith (c) 2016
 *
 * This file impliments I/0 using HDF5
 */


#include <mpi.h>
#include <hdf5.h>
#include <stdlib.h>
#include <unistd.h>     // access() to check if file exists
#include <stdbool.h>
#include <assert.h>
#include <string.h>
#include <yaml.h>

#include "error.h"
#include "state.h"
#include "io.h"

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
    return H5Fclose(file_id);
}

herr_t write_array_dataset (const char *name,
                            hid_t       group_id,
                            double     *arr,
                            state      *s)
{
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
     dataspace = H5Dget_space (dset_id);
     plist_id = H5Pcreate (H5P_DATASET_XFER);
     H5Pset_dxpl_mpio (plist_id, H5FD_MPIO_COLLECTIVE);

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
                            arr + row*2*(s->N/2 + 1));
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
    H5Pset_dxpl_mpio (plist_id, H5FD_MPIO_COLLECTIVE);

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

	MPI_Barrier(MPI_COMM_WORLD);
	printf("Loaded an array!\n");

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
    /* Read integer attribute from dataset 'group_id' */
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
     hid_t group_id;
     herr_t status;

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
     status = write_double_attribute ("kbT",    group_id, &s->kbT);
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

     status = write_array_dataset ("Concentration", group_id, s->c, s);
     status = write_array_dataset ("Density",   group_id, s->n, s);

     set_C (s);

     status = H5Gclose (group_id);

     return status;
}

herr_t load_state (state      *s,
                   hid_t       file_id,
                   const char *datafile)
{
    hid_t group_id;
    herr_t status;

    group_id = H5Gopen2 (file_id, datafile, H5P_DEFAULT);

    /* Read state attributes */

    status = read_double_attribute ("eta",      group_id, &s->eta);
    status = read_double_attribute ("chi",      group_id, &s->chi);
    status = read_double_attribute ("epsilon_0",group_id, &s->epsilon0);
    status = read_double_attribute ("sigma0",   group_id, &s->sigma0);
    status = read_double_attribute ("sigma",    group_id, &s->sigma);
    status = read_double_attribute ("omega",    group_id, &s->omega);
    status = read_double_attribute ("kbT",      group_id, &s->kbT);
    status = read_double_attribute ("Wc",       group_id, &s->Wc);

    status = read_double_attribute ("Mn", group_id, &s->Mn);
    status = read_double_attribute ("Mc", group_id, &s->Mc);

    status = read_double_attribute ("k0",     group_id, &s->k0);
    status = read_double_attribute ("alpha",  group_id, &s->alpha);
    status = read_double_attribute ("beta",   group_id, &s->beta);
    status = read_double_attribute ("rho",    group_id, &s->rho);
    status = read_double_attribute ("alphac", group_id, &s->alphac);

    status = read_double_attribute ("dx",     group_id, &s->dx);
    status = read_double_attribute ("dt",     group_id, &s->dt);
    status = read_double_attribute ("Time",   group_id, &s->t);
    status = read_int_attribute ("Time Step", group_id, &s->step);

    status = read_array_dataset ("Concentration", group_id, s->c, s);
    status = read_array_dataset ("Density",   group_id, s->n, s);

    set_C(s);

    status = H5Gclose (group_id);

    return status;

}
