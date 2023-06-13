#include "filters.hpp"

Matrix appendZeros(const Matrix &m, uint32_t zerosAmount);

double applyKernelToPix(uint32_t row, uint32_t col, const Matrix &img, const Matrix &filterKernel);

void forEach(Matrix &m, const std::function<double(double)> &f);

Matrix grayScaleImgMatrix(std::string img_filename) {
    int width, height, channels;
    unsigned char *img = stbi_load(img_filename.c_str(), &width, &height, &channels, 0);
    if(!img) {
        LOG("Não foi possível carregar a imagem");
        return {};
    }

    unsigned char *p = img;
    Matrix grayScaleColorMtx(height, width);
    for(uint32_t row = 0; row < height; row++) {
        for(uint32_t col = 0; col < width; col++, p += channels) {
            grayScaleColorMtx(row, col) = (uint8_t)((*p + *(p + 1) + *(p + 2))/3.);
        }
    }

    stbi_image_free(img);

    return grayScaleColorMtx;
}

void writeImgOfMatrix(const Matrix &m, std::string imgName) {
    unsigned char *img = (unsigned char *)malloc(m.cols()*m.rows());
    unsigned char *p = img;
    for(uint32_t row = 0; row < m.rows(); row++) {
        for(uint32_t col = 0; col < m.cols(); col++, p++) {
            *p = (uint8_t)m(row, col);
        }
    }

    stbi_write_jpg(imgName.c_str(), m.cols(), m.rows(), 1, img, 100);

    free(img);
}

Matrix forEachPix(const Matrix &img, double(*filterFunc)(uint32_t, uint32_t, const Matrix &, const Matrix &),
                  const Matrix &filterKernel) {
    uint32_t halfKernelSize = filterKernel.rows() >> 1;
    Matrix img_withZerosInBoards = appendZeros(img, halfKernelSize);

    Matrix filteredImg(img.rows(), img.cols());
    for(uint32_t row = halfKernelSize; row < img_withZerosInBoards.rows() - halfKernelSize; row++) {
        for(uint32_t col = halfKernelSize; col < img_withZerosInBoards.cols() - halfKernelSize; col++) {
            filteredImg(row - halfKernelSize, col - halfKernelSize) = filterFunc(row, col, img_withZerosInBoards,
                                                                                 filterKernel);
        }
    }

    return filteredImg;
}

Matrix appendZeros(const Matrix &m, uint32_t zerosAmount) {
    Matrix horizontalZerosLine(zerosAmount, m.cols());
    horizontalZerosLine.fill(0);
    Matrix m_withZerosUpAndDown(m.rows() + 2*horizontalZerosLine.rows(), m.cols());
    m_withZerosUpAndDown << horizontalZerosLine, m, horizontalZerosLine;

    Matrix verticalZerosLine(m_withZerosUpAndDown.rows(), zerosAmount);
    verticalZerosLine.fill(0);
    Matrix m_withZerosOnBoards(m_withZerosUpAndDown.rows(), m_withZerosUpAndDown.cols() + 2*verticalZerosLine.cols());
    m_withZerosOnBoards << verticalZerosLine, m_withZerosUpAndDown, verticalZerosLine;

    return m_withZerosOnBoards;
}

Matrix applyGaussianFilter(const Matrix &img) {
    Matrix gaussianKernel(5, 5);
    gaussianKernel << 1, 4, 7, 4, 1,
            4, 16, 26, 16, 4,
            7, 26, 41, 26, 7,
            4, 16, 26, 16, 4,
            1, 4, 7, 4, 1;
    gaussianKernel /= gaussianKernel.sum();

    return forEachPix(img, applyKernelToPix, gaussianKernel);
}

double applyKernelToPix(uint32_t row, uint32_t col, const Matrix &img, const Matrix &filterKernel) {
    Matrix k = img.block(row - (filterKernel.rows() >> 1), col - (filterKernel.cols() >> 1), filterKernel.rows(),
                         filterKernel.cols());
    return (k.cwiseProduct(filterKernel)).sum();
}

Matrix edgeDetector(const Matrix &img) {
    Matrix sobelKernel(5, 5);
    sobelKernel << 2, 1, 0, -1, -2,
            2, 1, 0, -1, -2,
            4, 2, 0, -2, -4,
            2, 1, 0, -1, -2,
            2, 1, 0, -1, -2;
    sobelKernel /= sobelKernel.rows()*sobelKernel.cols();

    Matrix filteredOnY = forEachPix(img, applyKernelToPix, sobelKernel);
    Matrix filteredOnX = forEachPix(img, applyKernelToPix, sobelKernel.transpose());

    forEach(filteredOnX, [](double x) { return x*x; });
    forEach(filteredOnY, [](double x) { return x*x; });

    Matrix sum = filteredOnY + filteredOnX;
    forEach(sum, [](double x) { return sqrt(x); });
    //forEach(sum, [](double x){ return x > 22 ? 255. : 0.;});
    return sum;
}

void forEach(Matrix &m, const std::function<double(double)> &f) {
    for(uint32_t row = 0; row < m.rows(); row++) {
        for(uint32_t col = 0; col < m.cols(); col++) {
            m(row, col) = f(m(row, col));
        }
    }
}

bool equals(double n1, double n2, double EPS) {
    return fabs(n1 - n2) <= EPS;
}

Matrix laplacianEdgeDetector(const Matrix &img) {
    Matrix laplacianKernel(9, 9);
    laplacianKernel << 0, 1, 1, 2, 2, 2, 1, 1, 0,
            1, 2, 4, 5, 5, 5, 4, 2, 1,
            1, 4, 5, 3, 0, 3, 5, 4, 1,
            2, 5, 3, -12, -24, -12, 3, 5, 2,
            2, 5, 0, -24, -40, -24, 0, 5, 2,
            2, 5, 3, -12, -24, -12, 3, 5, 2,
            1, 4, 5, 3, 0, 3, 5, 4, 1,
            1, 2, 4, 5, 5, 5, 4, 2, 1,
            0, 1, 1, 2, 2, 2, 1, 1, 0;

    uint32_t sum = 0;
    forEach(laplacianKernel, [&sum](double x) {
        sum += fabs(x);
        return x;
    });
    laplacianKernel /= sum;

    Matrix filtered = forEachPix(img, applyKernelToPix, laplacianKernel);
    //forEach(filtered, [=](double x){ return x > 150 ? 255 : 0;});

    return filtered;
}
