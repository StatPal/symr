/*
 * 
 * Adapted from https://answers.opencv.org/question/179162/how-can-i-load-4d-medical-image-of-type-nifti/
 * 
 * Precompile header using
g++ Read_files_2.hpp -I /usr/include/eigen3 -O3
 * - Don't do it now. taking huge gch file. 

 *
 * To do:
 	* Other format
 	* Include scl_slope, scl_inter
 *
 */


/* READ_DATA */
#ifndef READ_DATA
#define READ_DATA


#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <ostream>



/***************************************************
 * 
 * To read the data part
 * 
****************************************************/
template<typename T>
Eigen::MatrixXd read_data(char* const data_file, short dim[8], float& vox_offset, char will_write = 0){

	std::ifstream infile(data_file, std::ios::binary);
	if (infile.fail()) {
		std::cout << "Could not read file *.nii" << std::endl;
		exit(EXIT_FAILURE);
	}
	infile.seekg (vox_offset, std::ios::cur);

	int n = dim[3]*dim[2]*dim[1];
	long unsigned int num_voxels = dim[1] * dim[2] * dim[3] * dim[4] * dim[5] * dim[6] * dim[7];
	std::cout << "Reading " << num_voxels << " voxels" << std::endl;


	Eigen::MatrixXd img_mat;
	
	for (int i = 0; i < dim[7]; i++)
	{
		for (int j = 0; j < dim[6]; j++)
		{
			for (int k = 0; k < dim[5]; k++)
			{
				img_mat = Eigen::MatrixXd::Zero(n, dim[4]);
				for (int l = 0; l < dim[4]; l++)
				{
					std::vector<T> data(n, 0);
					infile.read(reinterpret_cast<char*>(&data[0]), sizeof(T)*n);
					
					for(int m = 0; m < n; ++m){
						img_mat(m, l) = data[ m ];
					}
					
					std::cout << l+1 << "-th set have size:" << img_mat.col(l).size() << 
					";  mean: " << img_mat.col(l).mean() << std::endl;
					
					if( will_write == 1){
						std::string file_name = "Output/output_new_";
						std::ostringstream oss;
						oss << i << "_" << j << "_" << k << "_" << l ;
						file_name += oss.str();
						file_name += ".txt";
						
						std::cout << "Writing " << dim[1] << "x" << dim[2] << "x" << dim[3] << " file " << file_name.c_str() << std::endl;
						std::cout << "mean: " << img_mat.col(l).mean() << std::endl;
						std::ofstream myfile;
						myfile.open(file_name.c_str());
						myfile << img_mat.col(l) ;
						myfile.close();
					}
				}
				std::cout << "-------------------------------------------\n-------------------------------------------\n";
				std::cout << std::flush;
			}
		}
	}
	
	return img_mat;
}



/***************************************************
 * 
 * To read the header and the data part 
 * help taken from https://brainder.org/2012/09/23/the-nifti-file-format/
 * and https://answers.opencv.org/question/179162/how-can-i-load-4d-medical-image-of-type-nifti/
 * 
****************************************************/
Eigen::MatrixXd Read_nift1(char* const data_file, short our_dim[8], char will_write = 0)
{
	//ifstream infile("../example.nii", std::ios::binary);
	std::ifstream infile(data_file, std::ios::binary);

	if (infile.fail())
	{
		std::cout << "Could not read file *.nii" << std::endl;
		exit(EXIT_FAILURE);
	}

	size_t bytes_read = 0;

	std::cout << "----------------\nReading header:\n----------------\n" ;

	int sizeof_header;
	infile.read(reinterpret_cast<char*>(&sizeof_header), sizeof(sizeof_header));
	bytes_read += infile.gcount();

	if (sizeof_header != 348)
	{
		std::cout << "Invalid header size: should be 348 bytes" << std::endl;
		exit(EXIT_FAILURE);
	}

	char data_type[10];
	infile.read(reinterpret_cast<char*>(&data_type), sizeof(data_type));
	bytes_read += infile.gcount();

	char db_name[18];
	infile.read(reinterpret_cast<char*>(&db_name), sizeof(db_name));
	bytes_read += infile.gcount();

	int extents;
	infile.read(reinterpret_cast<char*>(&extents), sizeof(extents));
	bytes_read += infile.gcount();

	short session_error;
	infile.read(reinterpret_cast<char*>(&session_error), sizeof(session_error));
	bytes_read += infile.gcount();

	char regular;
	infile.read(reinterpret_cast<char*>(&regular), sizeof(regular));
	bytes_read += infile.gcount();

	char dim_info;
	infile.read(reinterpret_cast<char*>(&dim_info), sizeof(dim_info));
	bytes_read += infile.gcount();

	short dim[8];
	infile.read(reinterpret_cast<char*>(&dim), sizeof(dim));
	bytes_read += infile.gcount();

	std::cout << dim[0] << " dimensions" << std::endl;
	std::cout << "Dim 1(X): " << dim[1] << std::endl;
	std::cout << "Dim 2(Y): " << dim[2] << std::endl;
	std::cout << "Dim 3(Z): " << dim[3] << std::endl;
	std::cout << "Dim 4(T): " << dim[4] << std::endl;
	std::cout << "Dim 5   : " << dim[5] << std::endl;
	std::cout << "Dim 6   : " << dim[6] << std::endl;
	std::cout << "Dim 7   : " << dim[7] << std::endl;

	float intent_p1;
	infile.read(reinterpret_cast<char*>(&intent_p1), sizeof(intent_p1));
	bytes_read += infile.gcount();

	float intent_p2;
	infile.read(reinterpret_cast<char*>(&intent_p2), sizeof(intent_p2));
	bytes_read += infile.gcount();

	float intent_p3;
	infile.read(reinterpret_cast<char*>(&intent_p3), sizeof(intent_p3));
	bytes_read += infile.gcount();

	short intent_code;
	infile.read(reinterpret_cast<char*>(&intent_code), sizeof(intent_code));
	bytes_read += infile.gcount();

	short datatype;
	infile.read(reinterpret_cast<char*>(&datatype), sizeof(datatype));
	bytes_read += infile.gcount();
	std::cout << "Datatype code: " << datatype << ";\t";

	short bitpix;
	infile.read(reinterpret_cast<char*>(&bitpix), sizeof(bitpix));
	bytes_read += infile.gcount();
	std::cout << "Bits per pixel: " << bitpix << std::endl;

	short slice_start;
	infile.read(reinterpret_cast<char*>(&slice_start), sizeof(slice_start));
	bytes_read += infile.gcount();

	float pixdim[8];
	infile.read(reinterpret_cast<char*>(&pixdim), sizeof(pixdim));
	bytes_read += infile.gcount();

	float vox_offset;
	infile.read(reinterpret_cast<char*>(&vox_offset), sizeof(vox_offset));
	bytes_read += infile.gcount();
	
	std::cout << "Bytes offset to data file: " << vox_offset << std::endl;
	
	float scl_slope;
	infile.read(reinterpret_cast<char*>(&scl_slope), sizeof(scl_slope));
	bytes_read += infile.gcount();

	float scl_inter;
	infile.read(reinterpret_cast<char*>(&scl_inter), sizeof(scl_inter));
	bytes_read += infile.gcount();

//Subrata - include these two into consideration.

	short slice_end;
	infile.read(reinterpret_cast<char*>(&slice_end), sizeof(slice_end));
	bytes_read += infile.gcount();

	char slice_code;
	infile.read(reinterpret_cast<char*>(&slice_code), sizeof(slice_code));
	bytes_read += infile.gcount();

	char xyzt_units;
	infile.read(reinterpret_cast<char*>(&xyzt_units), sizeof(xyzt_units));
	bytes_read += infile.gcount();

	float cal_max;
	infile.read(reinterpret_cast<char*>(&cal_max), sizeof(cal_max));
	bytes_read += infile.gcount();

	float cal_min;
	infile.read(reinterpret_cast<char*>(&cal_min), sizeof(cal_min));
	bytes_read += infile.gcount();

	float slice_duration;
	infile.read(reinterpret_cast<char*>(&slice_duration), sizeof(slice_duration));
	bytes_read += infile.gcount();

	float toffset;
	infile.read(reinterpret_cast<char*>(&toffset), sizeof(toffset));
	bytes_read += infile.gcount();

	int glmax;
	infile.read(reinterpret_cast<char*>(&glmax), sizeof(glmax));
	bytes_read += infile.gcount();

	int glmin;
	infile.read(reinterpret_cast<char*>(&glmin), sizeof(glmin));
	bytes_read += infile.gcount();

	char descrip[80];
	infile.read(reinterpret_cast<char*>(&descrip), sizeof(descrip));
	bytes_read += infile.gcount();

	char aux_file[24];
	infile.read(reinterpret_cast<char*>(&aux_file), sizeof(aux_file));
	bytes_read += infile.gcount();

	short qform_code;
	infile.read(reinterpret_cast<char*>(&qform_code), sizeof(qform_code));
	bytes_read += infile.gcount();

	short sform_code;
	infile.read(reinterpret_cast<char*>(&sform_code), sizeof(sform_code));
	bytes_read += infile.gcount();

	float quatern_b;
	infile.read(reinterpret_cast<char*>(&quatern_b), sizeof(quatern_b));
	bytes_read += infile.gcount();

	float quatern_c;
	infile.read(reinterpret_cast<char*>(&quatern_c), sizeof(quatern_c));
	bytes_read += infile.gcount();

	float quatern_d;
	infile.read(reinterpret_cast<char*>(&quatern_d), sizeof(quatern_d));
	bytes_read += infile.gcount();

	float qoffset_x;
	infile.read(reinterpret_cast<char*>(&qoffset_x), sizeof(qoffset_x));
	bytes_read += infile.gcount();

	float qoffset_y;
	infile.read(reinterpret_cast<char*>(&qoffset_y), sizeof(qoffset_y));
	bytes_read += infile.gcount();

	float qoffset_z;
	infile.read(reinterpret_cast<char*>(&qoffset_z), sizeof(qoffset_z));
	bytes_read += infile.gcount();

	float srow_x[4];
	infile.read(reinterpret_cast<char*>(&srow_x), sizeof(srow_x));
	bytes_read += infile.gcount();

	float srow_y[4];
	infile.read(reinterpret_cast<char*>(&srow_y), sizeof(srow_y));
	bytes_read += infile.gcount();

	float srow_z[4];
	infile.read(reinterpret_cast<char*>(&srow_z), sizeof(srow_z));
	bytes_read += infile.gcount();

	char intent_name[16];
	infile.read(reinterpret_cast<char*>(&intent_name), sizeof(intent_name));
	bytes_read += infile.gcount();

	char magic[4];
	infile.read(reinterpret_cast<char*>(&magic), sizeof(magic));
	bytes_read += infile.gcount();

	std::cout << "Read " << bytes_read << " bytes.\n-----------------------" << std::endl;
	if (bytes_read != 348) {
		std::cout << "Error reading header" << std::endl;
		exit(EXIT_FAILURE);		//Has read this many bytes
	}

	for(int i = 0; i < 8; ++i){
		our_dim[i] = dim[i];
	}

	switch(datatype){
		//case 1: return read_data<bool>(data_file, dim, vox_offset, will_write);
		//	break;
		case 2: return read_data<unsigned char>(data_file, dim, vox_offset, will_write);
			break;
		case 4: return read_data<signed short>(data_file, dim, vox_offset, will_write);
			break;
		case 8: return read_data<signed int>(data_file, dim, vox_offset, will_write);
			break;
		case 16: return read_data<float>(data_file, dim, vox_offset, will_write);
			break;
		//case 32: return read_data<complex>(data_file, dim, vox_offset, will_write);
		//	break;
		case 64: return read_data<double>(data_file, dim, vox_offset, will_write);
			break; 
		//case 128: return read_data<RGB>(data_file, dim, vox_offset, will_write);
		//	break; 
		case 256: return read_data<signed char>(data_file, dim, vox_offset, will_write);
			break; 
		case 512: return read_data<unsigned short>(data_file, dim, vox_offset, will_write);
			break; 
		case 768: return read_data<unsigned int>(data_file, dim, vox_offset, will_write);
			break; 
		case 1024: return read_data<signed long long>(data_file, dim, vox_offset, will_write);
			break; 
		case 1280: return read_data<unsigned long long>(data_file, dim, vox_offset, will_write);
			break; 
		case 1536: return read_data<long double>(data_file, dim, vox_offset, will_write);
			break; 
		//case 1792: return read_data<double pair>(data_file, dim, vox_offset, will_write);
		//	break;
		//case 2048: return read_data<long double pair>(data_file, dim, vox_offset, will_write);
		//	break;
		//case 2304: return read_data<RGBA>(data_file, dim, vox_offset, will_write);
		//	break; 
		default: std::cout << "Not supported yet" << std::endl;
			exit(EXIT_FAILURE);
			//https://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/datatype.html
	}
}

#endif /* !READ_DATA */
