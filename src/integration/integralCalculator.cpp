#include "integralCalculator.hpp"

struct integralData
{
    double upper, lower;
    std::string expr;
    double EPS;
};

integralData getUserInput();
void openedFormulas();
void closedFormulas();

void runIntegralCalc()
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
            LOG("\nPara sair entre [0]");
            LOG("Para calcular usando formulas fechadas entre [1]");
            LOG("Para calcular usando formulas abertas entre [2]");
            getline(std::cin, userAnswer);

            if( userAnswer == "0" )
                break;

            if( userAnswer == "1" )
                closedFormulas();
            else
                openedFormulas();
        }
    }
}

void closedFormulas()
{
    LOG("\nPara voltar entre [0]");
    LOG("Para usar a Regra do trapézio entre [1]");
    LOG("Para usar 1/3 de Simpson entre [2]");
    LOG("Para usar 3/8 de Simpson entre [3]");
    LOG("Para usar Regra de Boole entre [4]");


    std::string userAnswer;
    getline(std::cin, userAnswer);
    integralData input = getUserInput();

    if( input.expr == "0" )
        return;

    switch ( std::stoi(userAnswer) ) 
    {
        case 1:
            printResults(integrate(Fn(input.expr), input.upper, input.lower, input.EPS, trapezoidal), "Regra do Trapézio");
            break;
        case 2:
            printResults(integrate(Fn(input.expr), input.upper, input.lower, input.EPS, simpson13), "1/3 de Simpson");
            break;
        case 3:
            printResults(integrate(Fn(input.expr), input.upper, input.lower, input.EPS, simpson38), "3/8 de Simpson");
            break;
        case 4:
            printResults(integrate(Fn(input.expr), input.upper, input.lower, input.EPS, boole), "Regra de Boole");
            break;
    }
}

void openedFormulas()
{
    LOG("\nPara voltar entre [0]");
    LOG("Para usar Aberto de grau 1 entre [1]");
    LOG("Para usar Aberto de grau 2 entre [2]");
    LOG("Para usar Aberto de grau 3 entre [3]");
    LOG("Para usar Aberto de grau 4 entre [4]");


    std::string userAnswer;
    getline(std::cin, userAnswer);
    integralData input = getUserInput();

    if( input.expr == "0" )
        return;

    switch ( std::stoi(userAnswer) ) 
    {
        case 1:
            printResults(integrate(Fn(input.expr), input.upper, input.lower, input.EPS, opened2Points), "Aberto grau 1");
            break;
        case 2:
            printResults(integrate(Fn(input.expr), input.upper, input.lower, input.EPS, milne), "Aberto grau 2");
            break;
        case 3:
            printResults(integrate(Fn(input.expr), input.upper, input.lower, input.EPS, opened4points), "Aberto grau 3");
            break;
        case 4:
            printResults(integrate(Fn(input.expr), input.upper, input.lower, input.EPS, opened5points), "Aberto grau 4");
            break;
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
