#include "SpecialIntegrator.hpp"

const std::vector<std::string> weightFunctions = {"e^(-(x^2))", "e^x", "1/sqrt(1 - (x^2))"};
const std::vector<std::string> functionsToAdd = {"e^(x^2)", "e^(-x)", "sqrt(1 - (x^2))"};
const std::vector<std::vector<SpecialIntegrator::TableValues>> SpecialTableValues =
{
    {//Hermite
        {{0.886226925452758013649, 0.886226925452758013649}, {-0.7071067811865475244008, 0.7071067811865475244008}},
        {{0.295408975150919337883, 1.181635900603677351532, 0.295408975150919337883}, {-1.224744871391589049099, 0, 1.224744871391589049099}},
        {{0.081312835447245177143, 0.8049140900055128365061, 0.8049140900055128365061, 0.08131283544724517714303}, {-1.650680123885784555883, -0.5246476232752903178841, 0.5246476232752903178841,1.650680123885784555883}}
    },
    {//Laguerre
        {{0.8535533905932737622004, 0.1464466094067262377996}, {0.5857864376269049511983, 3.414213562373095048802}},
        {{0.71109300992917301545, 0.2785177335692408488014, 0.01038925650158613574897}, {0.4157745567834790833115, 2.29428036027904171982, 6.289945082937479196866}},
        {{0.603154104341633601636, 0.3574186924377996866415, 0.03888790851500538427244, 5.392947055613274501038E-4}, {0.3225476896193923118004, 1.745761101158346575687, 4.536620296921127983279, 9.395070912301133129234}}
    },
    {//Chebyshev
        {{1.570796326794896619231, 1.570796326794896619231}, {-0.7071067811865475244008, 0.7071067811865475244008}},
        {{1.047197551196597746154, 1.047197551196597746154, 1.047197551196597746154}, {-0.8660254037844386467637, 0, 0.8660254037844386467637}},
        {{0.7853981633974483096157, 0.7853981633974483096157, 0.785398163397448309616, 0.7853981633974483096157}, {-0.9238795325112867561282, -0.3826834323650897717285, 0.3826834323650897717285, 0.9238795325112867561282}}
    }
};

SpecialIntegrator::SpecialIntegrator(SpecialIntegrator::Type type, SpecialIntegrator::Degree degree)
: tableValues(SpecialTableValues[type][degree]), weightFunction(weightFunctions[type]), functionToAdd(functionsToAdd[type]) {}

double SpecialIntegrator::integrate(const Fn &f) {
    Fn weightFunction(this->weightFunction);
    Fn functionWithWeight(functionToAdd + "* (" + f.expression + ")");

    double sum = 0;
    for (uint32_t i = 0; i < tableValues.weights.size(); i++)
        sum += weightFunction.eval(tableValues.alfa[i]) * tableValues.weights[i] * functionWithWeight.eval(tableValues.alfa[i]);

    return sum;
}