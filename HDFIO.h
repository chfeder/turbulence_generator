#ifndef HDFIO_H
#define HDFIO_H

#include <hdf5.h>
#include <assert.h>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

// we use H5_HAVE_PARALLEL defined in hdf5.h to signal whether we have MPI support or not
#ifdef H5_HAVE_PARALLEL
#include <mpi.h>
#else
#ifndef MPI_Comm
#define MPI_Comm int
#endif
#ifndef MPI_COMM_NULL
#define MPI_COMM_NULL 0
#endif
#endif

/**
 * HDFIO class
 *
 * This class represents an HDF5 object/file and provides basic
 * HDF operations like opening a file, reading data form a dataset
 * of that file, creating new HDF files and datasets and writing
 * datasets into a created file. There are also some functions
 * for getting and setting different attributes, like the
 * dimensionaltity of the HDF dataset to be written or the
 * datasetnames inside an HDF5 file.
 *
 * @author Christoph Federrath
 * @version 2007-2020
 *
 */

class HDFIO
{
    private:
        std::string ClassSignature, Filename;
        int Rank; //dimensionality of the field
        // HDF5 stuff
        hid_t   File_id, Dataset_id, Dataspace_id;
        hsize_t HDFSize, HDFDims[4];
        herr_t  HDF5_status, HDF5_error;
        bool Debug;

        public:

        /**
         * Default constructor.
         *
         */
        HDFIO()
        : ClassSignature("HDFIO: "),
          Filename(), //create Nullstring
          Rank(1), //dimensionality 1 (x only)
          File_id(0), //set all to 0
          Dataset_id(0),
          Dataspace_id(0),
          HDFSize(0), //set Buffersize to 0
          HDF5_status(0), //set HDF5 status to 0
          HDF5_error(-1), //set to -1
          Debug(false)
        {
            for(int i = 0; i < 4; i++)
                HDFDims[i] = 0; //set Dimensions to (0,0,0,0)
        }

        // get function signature for printing
        std::string FuncSig(const std::string funcname)
        {
            return ClassSignature+funcname+": ";
        };

        /**
         * open an HDF5 file in serial mode
         * @param Filename HDF5 filename
         * @param read_write_char 'r': read only flag, 'w': write flag
         *
         */
        void open(const std::string Filename, const char read_write_char)
        {
            this->open(Filename, read_write_char, MPI_COMM_NULL); // open in serial mode
        }

        /**
         * open an HDF5 file (overloaded)
         * @param Filename HDF5 filename
         * @param read_write_char 'r': read only flag, 'w': write flag
         * @param comm: MPI communicator for parallel file I/O
         *
         */
        void open(const std::string Filename, const char read_write_char, MPI_Comm comm)
        {
            this->Filename = Filename;

            hid_t plist_id = H5P_DEFAULT;
#ifdef H5_HAVE_PARALLEL
            if (comm != MPI_COMM_NULL) {
                plist_id = H5Pcreate(H5P_FILE_ACCESS);
                H5Pset_fapl_mpio(plist_id, comm, MPI_INFO_NULL);
                /*
                hsize_t threshold = 0;
                hsize_t alignment = 0;
                HDF5_status = H5Pget_alignment(plist_id, &threshold, &alignment);
                std::cout<<"current threshold, alignment = "<<threshold<<" "<<alignment<<std::endl;
                HDF5_status = H5Pset_alignment(plist_id, (hsize_t)0, (hsize_t)(4 * 1024 * 1024));
                HDF5_status = H5Pget_alignment(plist_id, &threshold, &alignment);
                std::cout<<"new threshold, alignment = "<<threshold<<" "<<alignment<<std::endl;
                */
                assert( HDF5_status != HDF5_error );
            }
#endif
            switch (read_write_char)
            {
                case 'r':
                {
                    // open HDF5 file in read only mode
                    File_id = H5Fopen(Filename.c_str(), H5F_ACC_RDONLY, plist_id);
                    assert( File_id != HDF5_error );
                    break;
                }
                case 'w':
                {
                    // open HDF5 file in write mode
                    File_id = H5Fopen(Filename.c_str(), H5F_ACC_RDWR, plist_id);
                    assert( File_id != HDF5_error );
                    break;
                }
                default:
                {
                    // open HDF5 file in read only mode
                    File_id = H5Fopen(Filename.c_str(), H5F_ACC_RDONLY, plist_id);
                    assert( File_id != HDF5_error );
                    break;
                }
            }
#ifdef H5_HAVE_PARALLEL
            if (comm != MPI_COMM_NULL) H5Pclose(plist_id);
#endif
        }

        /**
         * close HDF5 file
         *
         */
        void close(void)
        {
            // close HDF5 file
            HDF5_status = H5Fclose(File_id);
            assert( HDF5_status != HDF5_error );
        }


        /**
         * read data from a dataset
         * @param *DataBuffer pointer to double/float/int array to which data is to be written
         * @param Datasetname datasetname
         * @param DataType (i.e. H5T_STD_I32LE)
         *
         */
        void read(void* const DataBuffer, const std::string Datasetname, const hid_t DataType)
        {
            this->read(DataBuffer, Datasetname, DataType, MPI_COMM_NULL); // serial mode
        }

        /**
         * read data from a dataset
         * @param *DataBuffer pointer to double/float/int array to which data is to be written
         * @param Datasetname datasetname
         * @param DataType (i.e. H5T_STD_I32LE)
         * @param comm: MPI communicator for parallel file I/O
         *
         */
        void read(void* const DataBuffer, const std::string Datasetname, const hid_t DataType, MPI_Comm comm)
        {
            // get dimensional information from dataspace and update HDFSize
            this->getDims(Datasetname);

            // open dataset
            Dataset_id = H5Dopen(File_id, Datasetname.c_str(), H5P_DEFAULT);
            assert( Dataset_id != HDF5_error );

            // open dataspace (to get dimensions)
            Dataspace_id = H5Dget_space(Dataset_id);
            assert( Dataspace_id != HDF5_error );

            /// create property list for collective dataset i/o
            hid_t plist_id = H5P_DEFAULT;
#ifdef H5_HAVE_PARALLEL
            if (comm != MPI_COMM_NULL) {
                plist_id = H5Pcreate(H5P_DATASET_XFER);
                H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
            }
#endif
            // read buffer /memspaceid //filespaceid
            HDF5_status = H5Dread(Dataset_id, DataType, H5S_ALL, H5S_ALL, plist_id, DataBuffer);
            assert( HDF5_status != HDF5_error );
#ifdef H5_HAVE_PARALLEL
            if (comm != MPI_COMM_NULL) {
                HDF5_status = H5Pclose(plist_id);
                assert( HDF5_status != HDF5_error );
            }
#endif
            HDF5_status = H5Sclose(Dataspace_id);
            assert( HDF5_status != HDF5_error );

            HDF5_status = H5Dclose(Dataset_id);
            assert( HDF5_status != HDF5_error );

        } //read

        /**
         * read attribute from a dataset
         * @param *DataBuffer pointer to double/float/int array to which data is to be written
         * @param Datasetname name of dataset
         * @param Attributename name of attribute to be read
         * @param DataType (i.e. H5T_STD_I32LE)
         *
         */
        void read_attribute(void* const DataBuffer, const std::string Datasetname, const std::string Attributename, const hid_t DataType)
        {
            this->read_attribute(DataBuffer, Datasetname, Attributename, DataType, MPI_COMM_NULL); // serial mode
        }

        /**
         * read attribute from a dataset
         * @param *DataBuffer pointer to double/float/int array to which data is to be written
         * @param Datasetname name of dataset
         * @param Attributename name of attribute to be read
         * @param DataType (i.e. H5T_STD_I32LE)
         * @param comm: MPI communicator for parallel file I/O
         *
         */
        void read_attribute(void* const DataBuffer, const std::string Datasetname, const std::string Attributename, const hid_t DataType, MPI_Comm comm)
        {
            // open dataset
            Dataset_id = H5Dopen(File_id, Datasetname.c_str(), H5P_DEFAULT);
            assert( Dataset_id != HDF5_error );

            /// create property list for collective dataset i/o
            hid_t plist_id = H5P_DEFAULT;
#ifdef H5_HAVE_PARALLEL
            if (comm != MPI_COMM_NULL) {
                plist_id = H5Pcreate(H5P_DATASET_XFER);
                H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
            }
#endif
            /// open attribute
            hid_t attr_id = H5Aopen(Dataset_id, Attributename.c_str(), plist_id);

            /// read attribute
            HDF5_status = H5Aread(attr_id, DataType, DataBuffer);
            assert( HDF5_status != HDF5_error );

#ifdef H5_HAVE_PARALLEL
            if (comm != MPI_COMM_NULL) {
                HDF5_status = H5Pclose(plist_id);
                assert( HDF5_status != HDF5_error );
            }
#endif
            HDF5_status = H5Aclose(attr_id);
            assert( HDF5_status != HDF5_error );

            HDF5_status = H5Dclose(Dataset_id);
            assert( HDF5_status != HDF5_error );

        } //read

        /**
         * read a slab of data from a dataset (overloaded)
         * @param *DataBuffer pointer to double/float/int array to which data is to be written
         * @param Datasetname datasetname
         * @param DataType (i.e. H5T_STD_I32LE)
         * @param offset (offset array)
         * @param count (count array, i.e. the number of cells[dim] to be selected)
         * @param out_rank (rank of the output array selected)
         * @param out_offset (output offset array)
         * @param out_count (output count array, i.e. the number of cells[dim] in output)
         *
         */
        void read_slab(void* const DataBuffer, const std::string Datasetname, const hid_t DataType,
                       const hsize_t offset[], const hsize_t count[],
                       const hsize_t out_rank, const hsize_t out_offset[], const hsize_t out_count[])
        {
            /// simply set total_out_count = out_count, in which case we selected a hyperslab that fits completely into the
            /// total size of the output array. Use the overloaded function below, if output goes into an offset hyperslab
            this->read_slab(DataBuffer, Datasetname, DataType, offset, count, out_rank, out_offset, out_count, out_count);
        };

        /**
         * read a slab of data from a dataset (overloaded)
         * @param *DataBuffer pointer to double/float/int array to which data is to be written
         * @param Datasetname datasetname
         * @param DataType (i.e. H5T_STD_I32LE)
         * @param offset (offset array)
         * @param count (count array, i.e. the number of cells[dim] to be selected)
         * @param out_rank (rank of the output array selected)
         * @param out_offset (output offset array)
         * @param out_count (output count array, i.e. the number of cells[dim] in output)
         * @param comm: MPI communicator for parallel file I/O
         *
         */
        void read_slab(void* const DataBuffer, const std::string Datasetname, const hid_t DataType,
                       const hsize_t offset[], const hsize_t count[],
                       const hsize_t out_rank, const hsize_t out_offset[], const hsize_t out_count[], MPI_Comm comm)
        {
            /// simply set total_out_count = out_count, in which case we selected a hyperslab that fits completely into the
            /// total size of the output array. Use the overloaded function below, if output goes into an offset hyperslab
            this->read_slab(DataBuffer, Datasetname, DataType, offset, count, out_rank, out_offset, out_count, out_count, comm);
        };

        /**
         * read a slab of data from a dataset
         * @param *DataBuffer pointer to double/float/int array to which data is to be written
         * @param Datasetname datasetname
         * @param DataType (i.e. H5T_STD_I32LE)
         * @param offset (offset array)
         * @param count (count array, i.e. the number of cells[dim] to be selected)
         * @param out_rank (rank of the output array selected)
         * @param out_offset (output offset array)
         * @param out_count (output count array, i.e. the number of cells[dim] in output)
         * @param total_out_count (total output count array, i.e. the total number of cells[dim] in output)
         *
         */
        void read_slab(void* const DataBuffer, const std::string Datasetname, const hid_t DataType, 
                       const hsize_t offset[], const hsize_t count[], 
                       const hsize_t out_rank, const hsize_t out_offset[], const hsize_t out_count[],
                       const hsize_t total_out_count[])
        {
            this->read_slab(DataBuffer, Datasetname, DataType, offset, count, out_rank, out_offset, out_count, total_out_count, MPI_COMM_NULL); // serial mode
        }

        /**
         * read a slab of data from a dataset
         * @param *DataBuffer pointer to double/float/int array to which data is to be written
         * @param Datasetname datasetname
         * @param DataType (i.e. H5T_STD_I32LE)
         * @param offset (offset array)
         * @param count (count array, i.e. the number of cells[dim] to be selected)
         * @param out_rank (rank of the output array selected)
         * @param out_offset (output offset array)
         * @param out_count (output count array, i.e. the number of cells[dim] in output)
         * @param total_out_count (total output count array, i.e. the total number of cells[dim] in output)
         * @param comm: MPI communicator for parallel file I/O
         *
         */
        void read_slab(void* const DataBuffer, const std::string Datasetname, const hid_t DataType,
                       const hsize_t offset[], const hsize_t count[],
                       const hsize_t out_rank, const hsize_t out_offset[], const hsize_t out_count[],
                       const hsize_t total_out_count[], MPI_Comm comm)
        {
            // get dimensional information from dataspace and update HDFSize
            this->getDims(Datasetname);

            // open dataset
            Dataset_id = H5Dopen(File_id, Datasetname.c_str(), H5P_DEFAULT);
            assert( Dataset_id != HDF5_error );

            // open dataspace (to get dimensions)
            Dataspace_id = H5Dget_space(Dataset_id);
            assert( Dataspace_id != HDF5_error );

            // select hyperslab
            HDF5_status = H5Sselect_hyperslab(Dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
            assert( HDF5_status != HDF5_error );

            // create memspace
            hid_t Memspace_id = H5Screate_simple(out_rank, total_out_count, NULL);
            HDF5_status = H5Sselect_hyperslab(Memspace_id, H5S_SELECT_SET, out_offset, NULL, out_count, NULL);
            assert( HDF5_status != HDF5_error );

            /// create property list for collective dataset i/o
            hid_t plist_id = H5P_DEFAULT;
#ifdef H5_HAVE_PARALLEL
            if (comm != MPI_COMM_NULL) {
                plist_id = H5Pcreate(H5P_DATASET_XFER);
                H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
            }
#endif
            // read buffer
            HDF5_status = H5Dread(Dataset_id, DataType, Memspace_id, Dataspace_id, plist_id, DataBuffer );
            assert( HDF5_status != HDF5_error );
#ifdef H5_HAVE_PARALLEL
            if (comm != MPI_COMM_NULL) {
                HDF5_status = H5Pclose(plist_id);
                assert( HDF5_status != HDF5_error );
            }
#endif
            HDF5_status = H5Sclose(Memspace_id);
            assert( HDF5_status != HDF5_error );

            HDF5_status = H5Sclose(Dataspace_id);
            assert( HDF5_status != HDF5_error );

            HDF5_status = H5Dclose(Dataset_id);
            assert( HDF5_status != HDF5_error );

        }//read_slab

        /**
         * overwrite a slab of data from a dataset
         * @param *DataBuffer pointer to double/float/int array to which data is to be written
         * @param Datasetname datasetname
         * @param DataType (i.e. H5T_STD_I32LE)
         * @param offset (offset array)
         * @param count (count array, i.e. the number of cells[dim] to be selected)
         * @param out_rank (rank of the output array selected)
         * @param out_offset (output offset array)
         * @param out_count (output count array, i.e. the number of cells[dim] in output)
         *
         */
        void overwrite_slab(void* const DataBuffer, const std::string Datasetname, const hid_t DataType, 
                            const hsize_t offset[], const hsize_t count[], 
                            const hsize_t out_rank, const hsize_t out_offset[], const hsize_t out_count[])
        {
            this->overwrite_slab(DataBuffer, Datasetname, DataType, 
                                 offset, count, out_rank, out_offset, out_count, MPI_COMM_NULL); // serial mode
        }

        /**
         * overwrite a slab of data from a dataset
         * @param *DataBuffer pointer to double/float/int array to which data is to be written
         * @param Datasetname datasetname
         * @param DataType (i.e. H5T_STD_I32LE)
         * @param offset (offset array)
         * @param count (count array, i.e. the number of cells[dim] to be selected)
         * @param out_rank (rank of the output array selected)
         * @param out_offset (output offset array)
         * @param out_count (output count array, i.e. the number of cells[dim] in output)
         * @param comm: MPI communicator for parallel file I/O
         *
         */
        void overwrite_slab(void* const DataBuffer, const std::string Datasetname, const hid_t DataType,
                            const hsize_t offset[], const hsize_t count[],
                            const hsize_t out_rank, const hsize_t out_offset[], const hsize_t out_count[],
                            MPI_Comm comm)
        {
            // get dimensional information from dataspace and update HDFSize
            this->getDims(Datasetname);

            // open dataset
            Dataset_id = H5Dopen(File_id, Datasetname.c_str(), H5P_DEFAULT);
            assert( Dataset_id != HDF5_error );

            // open dataspace (to get dimensions)
            Dataspace_id = H5Dget_space(Dataset_id);
            assert( Dataspace_id != HDF5_error );

            // select hyperslab
            HDF5_status = H5Sselect_hyperslab(Dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
            assert( HDF5_status != HDF5_error );

            // create memspace
            hid_t Memspace_id = H5Screate_simple(out_rank, out_count, NULL);
            HDF5_status = H5Sselect_hyperslab(Memspace_id, H5S_SELECT_SET, out_offset, NULL, out_count, NULL);
            assert( HDF5_status != HDF5_error );

            /// create property list for collective dataset i/o
            hid_t plist_id = H5P_DEFAULT;
#ifdef H5_HAVE_PARALLEL
            if (comm != MPI_COMM_NULL) {
                plist_id = H5Pcreate(H5P_DATASET_XFER);
                H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
            }
#endif
            // overwrite dataset
            HDF5_status = H5Dwrite(Dataset_id, DataType, Memspace_id, Dataspace_id, plist_id, DataBuffer);
            assert( HDF5_status != HDF5_error );
#ifdef H5_HAVE_PARALLEL
            if (comm != MPI_COMM_NULL) {
                HDF5_status = H5Pclose(plist_id);
                assert( HDF5_status != HDF5_error );
            }
#endif
            HDF5_status = H5Sclose(Memspace_id);
            assert( HDF5_status != HDF5_error );

            HDF5_status = H5Sclose(Dataspace_id);
            assert( HDF5_status != HDF5_error );

            HDF5_status = H5Dclose(Dataset_id);
            assert( HDF5_status != HDF5_error );

        }//overwrite_slab

        /**
         * overwrite data of an existing dataset (serial mode)
         * @param *DataBuffer pointer to double/float/int array containing data to be written
         * @param Datasetname datasetname
         * @param DataType (i.e. H5T_STD_I32LE)
         *
         */
        void overwrite(const void* const DataBuffer, const std::string Datasetname, const hid_t DataType)
        {
            this->overwrite(DataBuffer, Datasetname, DataType, MPI_COMM_NULL);
        }

        /**
         * overwrite data of an existing dataset (overloaded)
         * @param *DataBuffer pointer to double/float/int array containing data to be written
         * @param Datasetname datasetname
         * @param DataType (i.e. H5T_STD_I32LE)
         * @param comm: MPI communicator for parallel file I/O
         *
         */
        void overwrite(const void* const DataBuffer, const std::string Datasetname, const hid_t DataType, MPI_Comm comm)
        {
            // get dimensional information from dataspace and update HDFSize
            this->getDims(Datasetname);

            // open dataset
            Dataset_id = H5Dopen(File_id, Datasetname.c_str(), H5P_DEFAULT);
            assert( Dataset_id != HDF5_error );

            // open dataspace (to get dimensions)
            Dataspace_id = H5Dget_space(Dataset_id);
            assert( Dataspace_id != HDF5_error );

            /// create property list for collective dataset i/o
            hid_t plist_id = H5P_DEFAULT;
#ifdef H5_HAVE_PARALLEL
            if (comm != MPI_COMM_NULL) {
                plist_id = H5Pcreate(H5P_DATASET_XFER);
                H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
            }
#endif
            // overwrite dataset
            HDF5_status = H5Dwrite(Dataset_id, DataType, H5S_ALL, H5S_ALL, plist_id, DataBuffer);
            assert( HDF5_status != HDF5_error );
#ifdef H5_HAVE_PARALLEL
            if (comm != MPI_COMM_NULL) {
                HDF5_status = H5Pclose(plist_id);
                assert( HDF5_status != HDF5_error );
            }
#endif
            HDF5_status = H5Sclose(Dataspace_id);
            assert( HDF5_status != HDF5_error );

            HDF5_status = H5Dclose(Dataset_id);
            assert( HDF5_status != HDF5_error );

        }//overwrite

        /**
         * create new HDF5 file (serial mode)
         * @param Filename HDF5 filename
         *
         */
        void create(const std::string Filename)
        {
            this->create(Filename, MPI_COMM_NULL);
        }

        /**
         * create new HDF5 file (overloaded)
         * @param Filename HDF5 filename
         * @param comm: MPI communicator for parallel file I/O
         *
         */
        void create(const std::string Filename, MPI_Comm comm)
        {
            this->Filename = Filename;

            hid_t plist_id = H5P_DEFAULT;
#ifdef H5_HAVE_PARALLEL
            if (comm != MPI_COMM_NULL) {
                plist_id = H5Pcreate(H5P_FILE_ACCESS);
                H5Pset_fapl_mpio(plist_id, comm, MPI_INFO_NULL);
            }
#endif
            // create HDF5 file (overwrite)
            File_id = H5Fcreate(Filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
            assert( File_id != HDF5_error );

#ifdef H5_HAVE_PARALLEL
            if (comm != MPI_COMM_NULL) H5Pclose(plist_id);
#endif
        }


        /**
         * write HDF5 dataset (serial mode)
         * @param DataBuffer int/float/double array containing the data
         * @param Datasetname datasetname
         * @param Dimensions dataset dimensions
         * @param DataType (i.e. H5T_STD_I32LE)
         *
         */
        void write(const void* const DataBuffer, const std::string Datasetname,
                const std::vector<int> Dimensions, const hid_t DataType)
        {
            this->write(DataBuffer, Datasetname, Dimensions, DataType, MPI_COMM_NULL);
        }

        /**
         * write HDF5 dataset (overloaded)
         * @param DataBuffer int/float/double array containing the data
         * @param Datasetname datasetname
         * @param Dimensions dataset dimensions
         * @param DataType (i.e. H5T_STD_I32LE)
         * @param comm: MPI communicator for parallel file I/O
         *
         */
        void write(const void* const DataBuffer, const std::string Datasetname,
                const std::vector<int> Dimensions, const hid_t DataType, MPI_Comm comm)
        {
            this->setDims(Dimensions); // set dimensions
            this->write(DataBuffer, Datasetname, DataType, comm); // call write
        }

        /**
         * write HDF5 dataset (serial mode)
         * @param *DataBuffer pointer to int/float/double array containing the data
         * @param Datasetname datasetname
         * @param DataType (i.e. H5T_IEEE_F32BE, H5T_STD_I32LE, ...)
         *
         */
        void write(const void* const DataBuffer, const std::string Datasetname, const hid_t DataType)
        {
            this->write(DataBuffer, Datasetname, DataType, MPI_COMM_NULL);
        }

        /**
         * write HDF5 dataset (overloaded)
         * @param *DataBuffer pointer to int/float/double array containing the data
         * @param Datasetname datasetname
         * @param DataType (i.e. H5T_IEEE_F32BE, H5T_STD_I32LE, ...)
         * @param comm: MPI communicator for parallel file I/O
         *
         */
        void write(const void* const DataBuffer, const std::string Datasetname, const hid_t DataType, MPI_Comm comm)
        {
            // -------------- create dataspace
            Dataspace_id = H5Screate_simple(Rank, HDFDims, NULL);
            assert( Dataspace_id != HDF5_error );

            // -------------- create dataset
            Dataset_id = H5Dcreate(File_id, Datasetname.c_str(), DataType, Dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            assert( Dataset_id != HDF5_error );

            /// create property list for collective dataset i/o
            hid_t plist_id = H5P_DEFAULT;
#ifdef H5_HAVE_PARALLEL
            if (comm != MPI_COMM_NULL) {
                plist_id = H5Pcreate(H5P_DATASET_XFER);
                H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
            }
#endif
            // -------------- write dataset
            HDF5_status = H5Dwrite(Dataset_id, DataType, H5S_ALL, H5S_ALL, plist_id, DataBuffer);
            assert( HDF5_status != HDF5_error );
#ifdef H5_HAVE_PARALLEL
            if (comm != MPI_COMM_NULL) {
                HDF5_status = H5Pclose(plist_id);
                assert( HDF5_status != HDF5_error );
            }
#endif
            // -------------- close dataset
            HDF5_status = H5Dclose(Dataset_id);
            assert( HDF5_status != HDF5_error );

            // -------------- close dataspace
            HDF5_status = H5Sclose(Dataspace_id);
            assert( HDF5_status != HDF5_error );
        }


        /**
         * create empty HDF5 dataset
         * @param Datasetname datasetname
         * @param Dimensions dataset dimensions
         * @param DataType (i.e. H5T_IEEE_F32BE, H5T_STD_I32LE, ...)
         * @param comm: MPI communicator for parallel file I/O
         *
         */
        void create_dataset(const std::string Datasetname,
                            const std::vector<int> Dimensions, const hid_t DataType, MPI_Comm comm)
        {
            this->setDims(Dimensions); // set dimensions

            // -------------- create dataspace
            Dataspace_id = H5Screate_simple(Rank, HDFDims, NULL);
            assert( Dataspace_id != HDF5_error );

            // -------------- create dataset
            Dataset_id = H5Dcreate(File_id, Datasetname.c_str(), DataType, Dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            assert( Dataset_id != HDF5_error );

            // -------------- close dataset
            HDF5_status = H5Dclose(Dataset_id);
            assert( HDF5_status != HDF5_error );

            // -------------- close dataspace
            HDF5_status = H5Sclose(Dataspace_id);
            assert( HDF5_status != HDF5_error );
        }

        /**
         * delete HDF5 dataset
         * @param Datasetname datasetname
         *
         */
        void delete_dataset(const std::string Datasetname)
        {
            std::vector<std::string> dsets_in_file = getDatasetnames();
            bool dset_in_file = false; // see if the Datasetname exists in the file
            for (int i = 0; i < dsets_in_file.size(); i++) {
                if (dsets_in_file[i] == Datasetname) { dset_in_file = true; break; }
            }
            if (dset_in_file) { // if the Datasetname is in the file, delete it
                HDF5_status = H5Ldelete(File_id, Datasetname.c_str(), H5P_DEFAULT); // delete dataset
                assert( HDF5_status != HDF5_error );
            }
        }

        /**
         * delete all HDF5 datasets
         *
        */
        void delete_datasets(void)
        {
            std::vector<std::string> dsets_in_file = getDatasetnames();
            for (int i = 0; i < dsets_in_file.size(); i++) {
                this->delete_dataset(dsets_in_file[i]);
            }
        }

        /**
         * set the rank of a dataset
         * @param Rank the rank, dimensionality (0,1,2)
         *
         */
        void setRank(const int Rank)
        {
            this->Rank = Rank;
        }

        /**
         * set array dimension in different directions
         * @param dim dimension of the array
         *
         */
        void setDimX(const int dim_x)
        {
            HDFDims[0] = static_cast<hsize_t>(dim_x);
        }

        void setDimY(const int dim_y)
        {
            HDFDims[1] = static_cast<hsize_t>(dim_y);
        }

        void setDimZ(const int dim_z)
        {
            HDFDims[2] = static_cast<hsize_t>(dim_z);
        }

        void setDims(const std::vector<int> Dimensions)
        {
            Rank = Dimensions.size();
            for(int i = 0; i < Rank; i++)
                HDFDims[i] = static_cast<hsize_t>(Dimensions[i]);
        }

        /**
         * get the rank of the current dataset
         * @return int dimensionality
         *
         */
        int getRank(void) const
        {
            return Rank;
        }

        /**
         * get the rank of a dataset with name datasetname
         * @param Datasetname datasetname
         * @return int dimensionality
         *
         */
        int getRank(const std::string Datasetname)
        {
            // get dimensional information from dataspace and update HDFSize
            this->getDims(Datasetname);
            return Rank;
        }

        /**
         * get array dimension in different directions
         * @return dimension of the array
         *
         */
        int getDimX(void) const
        {
            return static_cast<int>(HDFDims[0]);
        }

        int getDimY(void) const
        {
            return static_cast<int>(HDFDims[1]);
        }

        int getDimZ(void) const
        {
            return static_cast<int>(HDFDims[2]);
        }

        /**
         * get dataset size of dataset with datasetname
         * @param Datasetname datasetname
         * @return size of the dataset
         *
         */
        int getSize(const std::string Datasetname)
        {
            // open dataset
            Dataset_id = H5Dopen(File_id, Datasetname.c_str(), H5P_DEFAULT);
            if (Dataset_id == HDF5_error)
            {
                std::cout << "HDFIO:  getSize():  CAUTION: Datasetname '" << Datasetname << "' does not exists in file '"<<this->getFilename()<<"'. Continuing..." << std::endl;
                return 0;
            }
            assert( Dataset_id != HDF5_error );

            // open dataspace
            Dataspace_id = H5Dget_space(Dataset_id);
            assert( Dataspace_id != HDF5_error );

            // get dimensional information from dataspace
            hsize_t HDFxdims[4], HDFmaxdims[4];
            Rank = H5Sget_simple_extent_dims(Dataspace_id, HDFxdims, HDFmaxdims);

            // from the dimensional info, calculate the size of the buffer.
            HDFSize = 1;
            for (int i = 0; i < Rank; i++) {
                HDFDims[i] = HDFxdims[i];
                HDFSize *= HDFDims[i];
            }

            HDF5_status = H5Sclose(Dataspace_id);
            assert( HDF5_status != HDF5_error );

            HDF5_status = H5Dclose(Dataset_id);
            assert( HDF5_status != HDF5_error );

            return static_cast<int>(HDFSize);
        }

        /**
         * get array dimension in different directions and update HDFSize
         * @param Datasetname datasetname for which dimensional info is read
         * @return dimension of the array
         *
         */
        std::vector<int> getDims(const std::string Datasetname)
        {
            // open dataset
            Dataset_id = H5Dopen(File_id, Datasetname.c_str(), H5P_DEFAULT);
            assert( Dataset_id != HDF5_error );

            // open dataspace
            Dataspace_id = H5Dget_space(Dataset_id);
            assert( Dataspace_id != HDF5_error );

            // get dimensional information from dataspace
            hsize_t HDFxdims[4], HDFmaxdims[4];

            Rank = H5Sget_simple_extent_dims(Dataspace_id, HDFxdims, HDFmaxdims);

            // from the dimensional info, calculate the size of the buffer.
            HDFSize = 1;
            for (int i = 0; i < Rank; i++) {
                HDFDims[i] = HDFxdims[i];
                HDFSize *= HDFDims[i];
            }

            HDF5_status = H5Sclose(Dataspace_id);
            assert( HDF5_status != HDF5_error );

            HDF5_status = H5Dclose(Dataset_id);
            assert( HDF5_status != HDF5_error );

            std::vector<int> ReturnDims(Rank);
            for(int i = 0; i < Rank; i++)
                ReturnDims[i] = static_cast<int>(HDFDims[i]);

            return ReturnDims;
        }

        /**
         * get HDF5 file ID
         * @return File_id
         *
         */
        hid_t getFileID(void)
        {
            return File_id;
        }

        /**
         * get HDF5 filename
         * @return filename
         *
         */
        std::string getFilename(void)
        {
            return Filename;
        }

        /**
         * get number of datasets in HDF5 file
         * @return the number of datasets
         *
         */
        int getNumberOfDatasets(void)
        {
            hsize_t *NumberofObjects = new hsize_t[1];
            H5Gget_num_objs(File_id, NumberofObjects);
            int returnNumber = static_cast<int>(NumberofObjects[0]);
            delete [] NumberofObjects;
            return returnNumber;
        }

        /**
         * get HDF5 datasetname
         * @param datasetnumber integer number identifying the dataset
         * @return datasetname
         *
         */
        std::string getDatasetname(const int datasetnumber)
        {
            char *Name = new char[256];
            H5Gget_objname_by_idx(File_id, datasetnumber, Name, 256);
            std::string returnName = static_cast<std::string>(Name);
            delete [] Name;
            return returnName;
        }

       /**
         * getDatasetnames
         * @return vector<string> datasetnames
         *
         */
       std::vector<std::string> getDatasetnames(void)
       {
           int nsets = this->getNumberOfDatasets();
           std::vector<std::string> ret(nsets);
           for (int i=0; i<nsets; i++)
           {
             ret[i] = this->getDatasetname(i);
             if (Debug) std::cout << "dataset in file: " << ret[i] << std::endl;
           }
           return ret;
       }

    // clear whitespace and/or terminators at the end of a string
    private: std::string Trim(const std::string input)
    {
        std::string ret = input;
        return ret.erase(input.find_last_not_of(" \n\r\t")+1);
    }

    /// functions that return FLASH integer, real, string, logical scalars as maps
    public: std::map<std::string, int> ReadFlashIntegerScalars()
    { return GetFLASHPropsInt("integer scalars"); };
    public: std::map<std::string, double> ReadFlashRealScalars()
    { return GetFLASHPropsReal("real scalars"); };
    public: std::map<std::string, std::string> ReadFlashStringScalars()
    { return GetFLASHPropsString("string scalars"); };
    public: std::map<std::string, bool> ReadFlashLogicalScalars()
    {   std::map<std::string, int> int_props = GetFLASHPropsInt("logical scalars");
        std::map<std::string, bool> ret; // define returned map
        for (std::map<std::string, int>::iterator it = int_props.begin(); it != int_props.end(); it++)
            ret[it->first] = (bool)it->second; // fill boolean return map
        return ret;
    };
    /// functions that return FLASH integer, real, string, logical runtime parameters as maps
    public: std::map<std::string, int> ReadFlashIntegerParameters()
    { return GetFLASHPropsInt("integer runtime parameters"); };
    public: std::map<std::string, double> ReadFlashRealParameters()
    { return GetFLASHPropsReal("real runtime parameters"); };
    public: std::map<std::string, std::string> ReadFlashStringParameters()
    { return GetFLASHPropsString("string runtime parameters"); };
    public: std::map<std::string, bool> ReadFlashLogicalParameters()
    {   std::map<std::string, int> int_props = GetFLASHPropsInt("logical runtime parameters");
        std::map<std::string, bool> ret; // define returned map
        for (std::map<std::string, int>::iterator it = int_props.begin(); it != int_props.end(); it++)
            ret[it->first] = (bool)it->second; // fill boolean return map
        return ret;
    };
    // returns map of FLASH integer (or logical) runtime parameters or scalars
    private: std::map<std::string, int> GetFLASHPropsInt(std::string dsetname)
    {
        int nProps = this->getDims(dsetname)[0];
        const unsigned int string_size = 79;
        hid_t h5string = H5Tcopy(H5T_C_S1);
        H5Tset_size(h5string, string_size);
        struct PropsStruct{ char name[string_size]; int value; };
        hid_t datatype = H5Tcreate(H5T_COMPOUND, sizeof(PropsStruct));
        H5Tinsert(datatype, "name", HOFFSET(PropsStruct, name), h5string);
        H5Tinsert(datatype, "value", HOFFSET(PropsStruct, value), H5T_NATIVE_INT);
        PropsStruct * pdat = new PropsStruct[nProps];
        this->read(pdat, dsetname, datatype);
        std::map<std::string, int> props;
        for (int i = 0; i < nProps; i++) {
            std::string pname = Trim((std::string)pdat[i].name);
            props[pname] = (int)pdat[i].value;
            if (Debug) std::cout<<FuncSig(__func__)<<"name, value = '"<<pname<<"', '"<<props[pname]<<"'"<<std::endl;
        }
        delete [] pdat;
        return props;
    };
    // returns map of FLASH real runtime parameters or scalars
    private: std::map<std::string, double> GetFLASHPropsReal(std::string dsetname)
    {
        int nProps = this->getDims(dsetname)[0];
        const unsigned int string_size = 79;
        hid_t h5string = H5Tcopy(H5T_C_S1);
        H5Tset_size(h5string, string_size);
        struct PropsStruct{ char name[string_size]; double value; };
        hid_t datatype = H5Tcreate(H5T_COMPOUND, sizeof(PropsStruct));
        H5Tinsert(datatype, "name", HOFFSET(PropsStruct, name), h5string);
        H5Tinsert(datatype, "value", HOFFSET(PropsStruct, value), H5T_NATIVE_DOUBLE);
        PropsStruct * pdat = new PropsStruct[nProps];
        this->read(pdat, dsetname, datatype);
        std::map<std::string, double> props;
        for (int i = 0; i < nProps; i++) {
            std::string pname = Trim((std::string)pdat[i].name);
            props[pname] = (double)pdat[i].value;
            if (Debug) std::cout<<FuncSig(__func__)<<"name, value = '"<<pname<<"', '"<<props[pname]<<"'"<<std::endl;
        }
        delete [] pdat;
        return props;
    };
    // returns map of FLASH string runtime parameters or scalars
    private: std::map<std::string, std::string> GetFLASHPropsString(std::string dsetname)
    {
        int nProps = this->getDims(dsetname)[0];
        const unsigned int string_size = 79;
        hid_t h5string = H5Tcopy(H5T_C_S1);
        H5Tset_size(h5string, string_size);
        struct PropsStruct{ char name[string_size]; char value[string_size]; };
        hid_t datatype = H5Tcreate(H5T_COMPOUND, sizeof(PropsStruct));
        H5Tinsert(datatype, "name", HOFFSET(PropsStruct, name), h5string);
        H5Tinsert(datatype, "value", HOFFSET(PropsStruct, value), H5Tcopy(h5string));
        PropsStruct * pdat = new PropsStruct[nProps];
        this->read(pdat, dsetname, datatype);
        std::map<std::string, std::string> props;
        for (int i = 0; i < nProps; i++) {
            std::string pname = Trim((std::string)pdat[i].name);
            props[pname] = Trim((std::string)pdat[i].value);
            if (Debug) std::cout<<FuncSig(__func__)<<": name, value = '"<<pname<<"', '"<<props[pname]<<"'"<<std::endl;
        }
        delete [] pdat;
        return props;
    };

}; // end: HDFIO
#endif
