#pragma once

#include <functional>
#include <stb/stb_image.h>
#include <stb/stb_image_resize.h>
#include <stb/stb_image_write.h>
#include <string>

#include "../utils/Fn.hpp"
#include "../eigen/mtxUtils.hpp"

Matrix applyGaussianFilter(const Matrix& m);
Matrix grayScaleImgMatrix(std::string img_filename);
void writeImgOfMatrix(const Matrix& m, std::string imgName);
Matrix edgeDetector(const Matrix& m);
Matrix laplacianEdgeDetector(const Matrix& img);