void mpi_print (const char *str);

hid_t io_init_new_file (const char *filename);

hid_t io_init_from_file (const char *filename);

herr_t io_finalize (hid_t file_id);

herr_t save_state (state* s, hid_t file_id);

herr_t load_state (state *s, hid_t file_id, const char *datafile);
