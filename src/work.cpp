# include <stdio.h>
# include "mpi.h"
# include <gsl/gsl_vector.h>
# include <gsl/gsl_sf.h>
# include <gsl/gsl_math.h>
# include <gsl/gsl_randist.h>
# include <gsl/gsl_complex.h>
# include <gsl/gsl_complex_math.h>
# include <H5Cpp.h>
# include <fftw3-mpi.h>
//# include <fftw3.h>

//using namespace H5;
using namespace std;

const int VEC_RANK = 1;
int i;
const int ARRAY_SIZE = 251;

//int writedata(gsl_vector * data_vector, H5::H5File& file, H5std_string& name, H5::DataSpace& dataspace );
int writedata(gsl_vector * data_vector, H5::H5File& file, H5std_string& name, H5::DataSpace& dataspace );

int main(int argc, char **argv)
{
	// Set up things
	gsl_rng_env_setup();
	double tval;
	double xval;
	double yval;
	double rval;

	double in[ARRAY_SIZE];
	fftw_complex out[ARRAY_SIZE];

	H5std_string t_name("t");
	H5std_string x_name("x");
	H5std_string y_name("y");
	H5std_string r_name("r");
	H5std_string out_final_name("out_final");
	gsl_vector * t         = gsl_vector_alloc(ARRAY_SIZE);
	gsl_vector * x         = gsl_vector_alloc(ARRAY_SIZE);
	gsl_vector * y         = gsl_vector_alloc(ARRAY_SIZE);
	gsl_vector * r         = gsl_vector_alloc(ARRAY_SIZE);
	gsl_vector * out_final = gsl_vector_calloc(ARRAY_SIZE);

	gsl_rng * rgen    = gsl_rng_alloc(gsl_rng_ranlxd2);

	printf("seed = %lu\n", gsl_rng_default_seed);

	for (i = 0; i<ARRAY_SIZE; i++)
	{
		tval = i*0.001;
		xval = gsl_sf_sin(2.0*M_PI*50.0*tval) + gsl_sf_sin(2.0*M_PI*120.0*tval);
		rval = gsl_ran_gaussian(rgen,2.0);
		yval  = xval + (float)rval;
		gsl_vector_set(t, i, tval);
		gsl_vector_set(x, i, xval);
		gsl_vector_set(y, i, yval);
		gsl_vector_set(r, i, rval);
		in[i] = yval;
	}
	for (i = 0; i<11; i++)
	{
		printf("t: %0.3f\tx: %0.3f\tr: %0.3f\ty: %0.3f\n", 
				gsl_vector_get(t, i),
				gsl_vector_get(x, i), 
				gsl_vector_get(r, i), 
				gsl_vector_get(y, i));
	}


	const H5std_string FILE_NAME("output.h5");

	hsize_t dimsf[1];
	dimsf[0] = ARRAY_SIZE;
	H5::DataSpace dataspace( VEC_RANK, dimsf );
	H5::FloatType datatype( H5::PredType::NATIVE_FLOAT );
	H5::H5File file(FILE_NAME,H5F_ACC_TRUNC);

	fftw_plan p = fftw_plan_dft_r2c_1d(ARRAY_SIZE, in, out, FFTW_ESTIMATE);
	fftw_execute(p);

	gsl_vector_complex * out_gsl = gsl_vector_complex_alloc(ARRAY_SIZE);

	gsl_complex temp;
	for (i = 0; i<ARRAY_SIZE; i++)
	{
		GSL_SET_COMPLEX(&temp, out[i][0], out[i][1]);
		gsl_vector_set(out_final, i, gsl_complex_abs2(temp));
	}

	writedata(t, file, t_name, dataspace);
	writedata(x, file, x_name, dataspace);
	writedata(y, file, y_name, dataspace);
	writedata(r, file, r_name, dataspace);
	writedata(out_final, file, out_final_name, dataspace);
//        H5::DataSet dataset = file.createDataSet( t_name, datatype, dataspace );
//        dataset.write(data, H5::PredType::NATIVE_FLOAT );

//        H5::DataSet dataset = file.createDataSet( y_name, datatype, dataspace );
//        dataset.write(data, H5::PredType::NATIVE_FLOAT );
	return 0;
}

int writedata(gsl_vector * data_vector, H5::H5File& file, H5std_string& name, H5::DataSpace& dataspace )
{
	int i;
	float data[ARRAY_SIZE];
	for (i = 0; i<ARRAY_SIZE; i++)
		data[i] = gsl_vector_get(data_vector,i);

	H5::FloatType datatype( H5::PredType::NATIVE_FLOAT );
	H5::DataSet dataset = file.createDataSet( name, datatype, dataspace );
	dataset.write(data, H5::PredType::NATIVE_FLOAT );

	return 0;
}
