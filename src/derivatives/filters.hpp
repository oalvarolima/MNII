#ifndef FILTERS_HPP
#define FILTERS_HPP

#include <stb/stb_image.h>
#include <stb/stb_image_write.h>
#include <stb/stb_image_resize.h>
#include <Eigen/Dense>
#include <string>
#include <functional>
#include "../utils/Fn.hpp"

typedef Eigen::MatrixXd Matrix;

Matrix applyGaussianFilter(const Matrix& m);
Matrix grayScaleImgMatrix(std::string img_filename);
void writeImgOfMatrix(const Matrix& m, std::string imgName);
Matrix edgeDetector(const Matrix& m);
Matrix laplacianEdgeDetector(const Matrix& img);

#endif 
