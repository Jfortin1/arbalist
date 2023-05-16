#include "hdf5_utils.h"
#include "H5Cpp.h"
#include <algorithm>
#include <vector>
#include <string>

H5::DataSet create_1d_compressed_hdf5_dataset(H5::Group& location, const H5::DataType& dtype, const std::string& name, hsize_t length, int deflate_level, hsize_t chunk) {
    H5::DataSpace dspace(1, &length);
 	H5::DSetCreatPropList plist;

    if (deflate_level >= 0 && length) {
        plist.setDeflate(deflate_level);
        if (chunk > length) {
            plist.setChunk(1, &length);
        } else {
            plist.setChunk(1, &chunk);
        }
    }

    return location.createDataSet(name, dtype, dspace, plist);
}

void dump_sparse_matrix(const std::string& output_file, const std::string& output_group, std::vector<std::vector<int> >& collected, int nrow, int deflate_level, int chunk_dim) {
    // Dumping it all to HDF5.
    H5::FileAccPropList fapl(H5::FileAccPropList::DEFAULT.getId());
    fapl.setCache(0, 511, chunk_dim * 10 * sizeof(int), 1);
    H5::H5File fhandle(output_file, H5F_ACC_RDWR, H5::FileCreatPropList::DEFAULT, fapl);
    H5::Group ghandle = fhandle.openGroup(output_group);

    std::vector<hsize_t> gathered(collected.size() + 1);
    int max_count = 0;
    for (size_t i = 0, end = collected.size(); i < end; ++i) {
        auto& x = collected[i];
        std::sort(x.begin(), x.end());

        int nunique = 0;
        int last = -1;
        int count = 0;
        for (auto y : x) {
            if (y > last) {
                ++nunique;
                if (count > max_count) {
                    max_count = count;
                }
                count = 1;
                last = y;
            } else {
                ++count;
            }
        }

        gathered[i + 1] = gathered[i] + nunique;
        if (count > max_count) {
            max_count = count;
        }
    }

    {
        H5::DataSet phandle = create_1d_compressed_hdf5_dataset(ghandle, H5::PredType::NATIVE_UINT64, "indptr", gathered.size(), deflate_level, chunk_dim);
        hsize_t len = gathered.size();
        H5::DataSpace dataspace(1, &len), memspace(1, &len);
        dataspace.selectAll();
        memspace.selectAll();
        phandle.write(gathered.data(), H5::PredType::NATIVE_HSIZE, memspace, dataspace);
    }

    {
        H5::DataSet shandle = create_1d_compressed_hdf5_dataset(ghandle, H5::PredType::NATIVE_UINT64, "shape", 2, deflate_level, chunk_dim);
        hsize_t len = 2;
        H5::DataSpace dataspace(1, &len), memspace(1, &len);
        dataspace.selectAll();
        memspace.selectAll();
        hsize_t shape[2];
        shape[0] = nrow;
        shape[1] = collected.size();
        shandle.write(shape, H5::PredType::NATIVE_HSIZE, memspace, dataspace);
    }

    hsize_t total = gathered.back();
    H5::DataSet ihandle = create_1d_compressed_hdf5_dataset(ghandle, H5::PredType::NATIVE_UINT32, "indices", total, deflate_level, chunk_dim);

    const H5::PredType* dtype;
    if (max_count <= 255) {
        dtype = &(H5::PredType::NATIVE_UINT8);
    } else if (max_count <= 65535) {
        dtype = &(H5::PredType::NATIVE_UINT16);
    } else {
        dtype = &(H5::PredType::NATIVE_UINT32);
    }
    H5::DataSet dhandle = create_1d_compressed_hdf5_dataset(ghandle, *dtype, "data", total, deflate_level, chunk_dim);

    std::vector<int> ibuffer;
    std::vector<int> dbuffer;
    H5::DataSpace dataspace(1, &total), memspace;
    hsize_t sofar = 0;

    for (auto& x : collected) {
        ibuffer.clear();
        dbuffer.clear();

        // Saving the damn thing in reverse order, like an idiot.
        for (auto y : x) {
            if (ibuffer.empty() || y != ibuffer.back()) {
                ibuffer.push_back(y);
                dbuffer.push_back(1);
            } else {
                ++(dbuffer.back());
            }
        }

        hsize_t count = ibuffer.size();
        dataspace.selectHyperslab(H5S_SELECT_SET, &count, &sofar);
        memspace.setExtentSimple(1, &count);
        memspace.selectAll();
        dhandle.write(dbuffer.data(), H5::PredType::NATIVE_INT, memspace, dataspace);
        ihandle.write(ibuffer.data(), H5::PredType::NATIVE_INT, memspace, dataspace);

        sofar += count;
    }
}
