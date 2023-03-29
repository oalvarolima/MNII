#include "Icalculator.hpp"

struct integralData
{
    double upper, lower;
    std::string expr;
    double EPS;
};

integralData getUserInput();

void runIcalculator()
{
    std::string userAnswer;
    LOG("\n\n=-=-=-=-=-=-=-=-=-= CALCULADOR NUMÉRICO DE INTEGRAIS =-=-=-=-=-=-=-=-=-=\n");
    while( userAnswer != "0" )
    {
        LOG("\nPara sair entre [0]");
        LOG("Para calcular o valor de uma integral entre [1]");
        getline(std::cin, userAnswer);
        if( userAnswer == "1" )
        {
            LOG("\nPara usar a Regra do trapézio entre [1]");
            LOG("Para usar 1/3 de Simpson entre [2]");
            LOG("Para usar 3/8 de Simpson entre [3]");

            getline(std::cin, userAnswer);
            integralData input = getUserInput();

            if( input.expr == "0" )
                continue;

            switch ( std::stoi(userAnswer) ) 
            {
                case 1:
                    printResults(integrate(Fn(input.expr), input.upper, input.lower, input.EPS, trapezoidalRule, 80, 1), "Regra do Trapézio");
                    break;
                case 2:
                    printResults(integrate(Fn(input.expr), input.upper, input.lower, input.EPS, simpson13, 80, 2), "1/3 de Simpson");
                    break;
                case 3:
                    printResults(integrate(Fn(input.expr), input.upper, input.lower, input.EPS, simpson38, 80, 3), "3/8 de Simpson");
                    break;
            }
        }
    }
}

integralData getUserInput()
{
    LOG("\nPara voltar [0]\nEntre a expressão limite_superior limite_inferior EPS");
    std::string rawInput, token;
    getline(std::cin, rawInput);
    std::stringstream input(rawInput);
    std::vector<std::string> tokens;

    while( std::getline(input, token, ' ') )
        tokens.push_back(token);

    double upper, lower, EPS;
    try 
    {
        if( tokens.size() < 4 )
            throw std::invalid_argument("");

        upper = std::stod(tokens[1]);
        lower = std::stod(tokens[2]);
        EPS = std::stod(tokens[3]);

        return {upper, lower, tokens[0], EPS};

    } catch (...) {
        if( tokens[0] == "0" )
            return {0, 0, "0", 0};

        LOG("Entrada errada, tente novamente.");
        return getUserInput();
    }
}
