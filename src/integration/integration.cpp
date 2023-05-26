#include "integration.hpp"

double integratorFunction(const Fn& f, double upper, double lower, uint32_t intervalsNum, 
                          double(intervalIntegrator)(const Fn&, double, double))
{
    double interval = (upper - lower) / intervalsNum;
    double result = 0;
    for( uint32_t n = 0; n < intervalsNum; n++, lower += interval )
    {
        result += intervalIntegrator(f, lower + interval, lower);
    }

    return result;
}

std::vector<itrInfo> integrate(const Fn& f, double upper, double lower, double EPS,
                               double(intervalIntegrator)(const Fn& f, double, double))
{
    std::vector<itrInfo> results;
    const uint32_t MaxItr = 1000000;
    double Error = EPS+EPS;
    double result, oldResult = 0;
    for( uint32_t intervalsNum = 1; Error >= EPS && intervalsNum <= MaxItr; intervalsNum++ )
    {
        result = integratorFunction(f, upper, lower, intervalsNum, intervalIntegrator);
        Error = fabs( (result - oldResult) / result );
        
        results.push_back( {intervalsNum, result, Error} );

        oldResult = result;
    }

    return results;
}

double NC::closed::two(const Fn& f, double upper, double lower)
{
    return .5*(upper - lower)*(f.eval(upper) + f.eval(lower));
}

double NC::closed::three(const Fn& f, double upper, double lower)
{
    double h = (upper - lower) / 2;
    return (1./3)*h*(f.eval(lower) + 4*f.eval(lower + h) + f.eval(upper));
}

double NC::closed::four(const Fn& f, double upper, double lower)
{
    double h = (upper - lower) / 3;
    return (3./8)*h*(f.eval(lower) + 3*f.eval(lower + h) + 3*f.eval(lower + 2*h) + f.eval(upper));
}

double NC::closed::five(const Fn& f, double upper, double lower)
{
    double h = (upper - lower) / 4;
    return (2./45)*h*(7*f.eval(lower) + 32*f.eval(lower + h) + 12*f.eval(lower + 2*h) + 32*f.eval(lower + 3*h) + 7*f.eval(upper));
}

double NC::opened::two(const Fn &f, double upper, double lower)
{
    double interval = (upper - lower) / 3;
    return (3./2)*interval*(f.eval(lower + interval) + f.eval(lower + 2*interval));
}

double NC::opened::three(const Fn &f, double upper, double lower)
{
    double interval = (upper - lower) / 4;
    return (4./3)*interval*(2*f.eval(lower + interval) - f.eval(lower + 2*interval) + 2*(f.eval(lower + 3*interval)));
}

double NC::opened::four(const Fn &f, double upper, double lower)
{
    double interval = (upper - lower) / 5;
    return (5./24)*interval*(11*f.eval(lower + interval) + f.eval(lower + 2*interval) + (f.eval(lower + 3*interval)) + 11*(f.eval(lower + 4*interval)));
}


double NC::opened::five(const Fn &f, double upper, double lower)
{
    double interval = (upper - lower) / 6;
    return (3./20)*interval*(22*f.eval(lower + interval) - 28*f.eval(lower + 2*interval) + 52*(f.eval(lower + 3*interval)) - 28*(f.eval(lower + 4*interval)) + 22*(f.eval(lower + 5*interval)));
}

double Gauss::Legendre::two(const Fn &f, double upper, double lower){
  auto x = [upper, lower](double alfa){ 
    return (upper+lower)/2 + (upper-lower)/2 * alfa;
  };
  std::vector<double> w = {1, 1};
  std::vector<double> alfa = {-0.5773502691896257, 0.5773502691896257};
  double somatorio = 0;
  for (int i=0; i<2; i++){
    somatorio += w[i] * f.eval(x(alfa[i]));
  }
  return (upper-lower)/2 * somatorio;
}

double Gauss::Legendre::three(const Fn &f, double upper, double lower){
  auto x = [upper, lower](double alfa){ 
    return (upper+lower)/2 + (upper-lower)/2 * alfa;
  };
  std::vector<double> w = {0.5555555555555556, 0.8888888888888889, 0.5555555555555556};
  std::vector<double> alfa = {-0.7745966692414834, 0, 0.7745966692414834};
  double somatorio = 0;
  for (int i=0; i<3; i++){
    somatorio += w[i] * f.eval(x(alfa[i]));
  }
  return (upper-lower)/2 * somatorio;
}

double Gauss::Legendre::four(const Fn &f, double upper, double lower){
  auto x = [upper, lower](double alfa){ 
    return (upper+lower)/2 + (upper-lower)/2 * alfa;
  };
  std::vector<double> w = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538};
  std::vector<double> alfa = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
  double somatorio = 0;
  for (int i=0; i<4; i++) {
    somatorio += w[i] * f.eval(x(alfa[i]));
  }
  return (upper-lower)/2 * somatorio;
}

double Gauss::Hermite::two(std::string f_inicial){
  const std::string nova_string = "e^(x^2) * (" + f_inicial + ")";
  Fn f = Fn(nova_string);
  Fn f2 = Fn("e^(-(x^2))");

  std::vector<double> w = {0.886226925452758013649, 0.886226925452758013649};
  std::vector<double> alfa = {-0.7071067811865475244008, 0.7071067811865475244008};
  double somatorio = 0;
  for (int i=0; i<2; i++){
    somatorio += f2.eval(alfa[i]) * w[i] * f.eval(alfa[i]);
  }
  return somatorio;
}

double Gauss::Hermite::three(std::string f_inicial){
  const std::string nova_string = "e^(x^2) * (" + f_inicial + ")";
  Fn f = Fn(nova_string);
  Fn f2 = Fn("e^(-(x^2))");

  std::vector<double> w = {0.295408975150919337883, 1.181635900603677351532, 0.295408975150919337883};
  std::vector<double> alfa = {-1.224744871391589049099, 0, 1.224744871391589049099};
  double somatorio = 0;
  for (int i=0; i<3; i++) {
    somatorio += f2.eval(alfa[i]) * w[i] * f.eval(alfa[i]);
  }
  return somatorio;
}

double Gauss::Hermite::four(std::string f_inicial){
  const std::string nova_string = "e^(x^2) * (" + f_inicial + ")";
  Fn f = Fn(nova_string);
  Fn f2 = Fn("e^(-(x^2))");

  std::vector<double> w = {0.081312835447245177143,
0.8049140900055128365061, 0.8049140900055128365061,
0.08131283544724517714303};
  std::vector<double> alfa = {-1.650680123885784555883, -0.5246476232752903178841, 0.5246476232752903178841,
1.650680123885784555883};
  double somatorio = 0;
  for (int i=0; i<4; i++){
    somatorio += f2.eval(alfa[i]) * w[i] * f.eval(alfa[i]);
  }
  return somatorio;
}

double Gauss::Laguerre::two(std::string f_inicial){
  const std::string nova_string = "e^x * (" + f_inicial + ")";
  Fn f = Fn(nova_string);
  Fn f2 = Fn("e^(-x)");

  std::vector<double> w = {0.8535533905932737622004, 0.1464466094067262377996};
  std::vector<double> alfa = {0.5857864376269049511983, 3.414213562373095048802};
  double somatorio = 0;
  for (int i=0; i<2; i++){
    somatorio += f2.eval(alfa[i]) * w[i] * f.eval(alfa[i]);
  }
  return somatorio;
}

double Gauss::Laguerre::three(std::string f_inicial){
  const std::string nova_string = "e^x * (" + f_inicial + ")";
  Fn f = Fn(nova_string);
  Fn f2 = Fn("e^(-x)");

  std::vector<double> w = {0.71109300992917301545,
0.2785177335692408488014, 0.01038925650158613574897};
  std::vector<double> alfa = {0.4157745567834790833115,
2.29428036027904171982, 6.289945082937479196866};
  double somatorio = 0;
  for (int i=0; i<3; i++){
    somatorio += f2.eval(alfa[i]) * w[i] * f.eval(alfa[i]);
  }
  return somatorio;
}

double Gauss::Laguerre::four(std::string f_inicial){
  const std::string nova_string = "e^x * (" + f_inicial + ")";
  Fn f = Fn(nova_string);
  Fn f2 = Fn("e^(-x)");

  std::vector<double> w = {0.603154104341633601636,
0.3574186924377996866415, 0.03888790851500538427244, 5.392947055613274501038E-4};
  std::vector<double> alfa = {0.3225476896193923118004,
1.745761101158346575687, 4.536620296921127983279, 9.395070912301133129234};
  double somatorio = 0;
  for (int i=0; i<4; i++){
    somatorio += f2.eval(alfa[i]) * w[i] * f.eval(alfa[i]);
  }
  return somatorio;
}

double Gauss::Chebyshev::two(std::string f_inicial){
  const std::string nova_string = "sqrt(1-x^2) * (" + f_inicial + ")";
  Fn f = Fn(nova_string);
  Fn f2 = Fn("1/sqrt(1-x^2)");

  std::vector<double> w = {1.570796326794896619231, 1.570796326794896619231};
  std::vector<double> alfa = {-0.7071067811865475244008, 0.7071067811865475244008};
  double somatorio = 0;
  for (int i=0; i<2; i++){
    somatorio += f2.eval(alfa[i]) * w[i] * f.eval(alfa[i]);
  }
  return somatorio;
}

double Gauss::Chebyshev::three(std::string f_inicial){
  const std::string nova_string = "sqrt(1-x^2) * (" + f_inicial + ")";
  Fn f = Fn(nova_string);
  Fn f2 = Fn("1/sqrt(1-x^2)");

  std::vector<double> w = {1.047197551196597746154,
1.047197551196597746154, 1.047197551196597746154};
  std::vector<double> alfa = {-0.8660254037844386467637,
0, 0.8660254037844386467637};
  double somatorio = 0;
  for (int i=0; i<3; i++){
    somatorio += f2.eval(alfa[i]) * w[i] * f.eval(alfa[i]);
  }
  return somatorio;
}

double Gauss::Chebyshev::four(std::string f_inicial){
  const std::string nova_string = "sqrt(1-x^2) * (" + f_inicial + ")";
  Fn f = Fn(nova_string);
  Fn f2 = Fn("1/sqrt(1-x^2)");

  std::vector<double> w = {0.7853981633974483096157, 0.7853981633974483096157, 0.785398163397448309616, 0.7853981633974483096157};
  std::vector<double> alfa = {-0.9238795325112867561282, -0.3826834323650897717285, 0.3826834323650897717285, 0.9238795325112867561282};
  double somatorio = 0;
  for (int i=0; i<4; i++) {
    somatorio += f2.eval(alfa[i]) * w[i] * f.eval(alfa[i]);
  }
  return somatorio;
}

void printResults(const std::vector<itrInfo> &results, const std::string &methodName)
{
    printf("\n                            %s\n\n", methodName.c_str());
    printf("N de intervalos               Integral                    EPS\n");
    for(auto [itrNum, Ivalue, EPS] : results)
    {
        printf("%s%u                     %s%.8f%s%.8f\n",
               std::string(8 - std::to_string(itrNum).length(), ' ').c_str(),
               itrNum,
               Ivalue < 0 ? "" : " ",
               Ivalue,
               std::string(26 - std::to_string(Ivalue).length(), ' ').c_str(),
               EPS);
    }
}
