/*It is a program that finds weights that satisfy n-dimensional AND, OR, and XOR gates through learning using perceptron.*/

#include <iostream>
#include <cmath>
#include <random>
#include <ctime>
#include <fstream>

#define SCOPE 5 //random scope
#define C 2 // learning rate
int errornum;
double tot_error;
bool is_XOR; 

struct Perc
{
    int layer;
    int num;
};
//Perceptron struct
class Output
{
public:
    double sigmoid(double input) //sigmoid function
    {
        double result;
        input = 0 - input;
        result = (double)1 / (1 + exp(input));
        return result;
    }
    void error(double error,int layer)//Function that exports the amount of change in error to a file.
    {
        using namespace std;
        errornum += 1;
        tot_error += error;
        int num = pow(2, layer);
        if (is_XOR) {
            if (errornum == num*3)
            {
                errornum = 0;
                tot_error = tot_error / layer;
                ofstream file("ERROR.txt", ios_base::out | ios_base::app);
                file << tot_error << "\n";
                tot_error = 0;
            }
        }
        else {
            if (errornum == num)
            {
                errornum = 0;
                tot_error = tot_error / layer;
                ofstream file("ERROR.txt", ios_base::out | ios_base::app);
                file << tot_error << "\n";
                tot_error = 0;
            }
        }
        
    }
    double calc(double input[], double weight[], int layer) //Find net
    {
        double net = 0;
        double result;
        for (int i = 0; i < layer + 1; i++)
        {
            net += input[i] * weight[i];

        }
        result = sigmoid(net);
        return result;
    }
};
class Learning //learning class
{
public:
    double* backdrop(double weight[], int layer, double target, double input[]) //return weight
    {
        Output output;

        double net;
        double num = pow(2, layer);

        net = output.calc(input, weight, layer); 
        for (int k = 0; k < layer + 1; k++)
        {
            weight[k] = weight[k] + ((target - net) * (1 - net) * net * C * input[k]); 
        }
        output.error(target - net,layer);
        return weight;
    }

};
class Gate //AND, OR, XOR
{
public:
    void AND(double weight[], int layer, int num)
    {
        Output output;
        Learning learning;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> dis(-SCOPE, SCOPE);

        double* input = new double[layer + 1]; 
        double& theta = weight[layer];//make last value of weight theta value
        input[layer] = -1; //input of theta is -1
        double* target = new double[num];
        using namespace std;

        ofstream file("AND.txt");

        for (int i = 0; i < layer + 1; i++) //Initialize the weight value randomly.
        {
            weight[i] = dis(gen);
        }
        while (true)
        {

            for (int i = 0; i < layer + 1; i++)
            {
                file << "\t" << weight[i];
            }
            file << "\n";

            double& theta = weight[layer];
            bool success = true;
            double net = 0;

            for (int i = 0; i < layer; i++)
            {
                input[i] = 1;
            }
            net = output.calc(input, weight, layer);
            if (net < 0.5)
            {
                success = false;
            }

            //The input value is the number of all cases where layers of inputs are 0,1, which is a duplicate permutation of "2π(layer)", Therefore, "2^(layer)=num" inputs exist.
            for (int i = 0; i < num - 1; i++)
            {
                int temp = i;
                for (int j = 0; j < layer; j++)
                {
                    input[j] = temp % 2;
                    temp = temp / 2;
                }

                net = 0;
                net = output.calc(input, weight, layer);
                if (net >= 0.5) 
                {
                    success = false;
                }
            }

            if (success) 
            {
                std::cout << "last weight: ";
                for (int i = 0; i < layer; i++)
                {
                    std::cout << weight[i] << " " << std::endl;
                }
                std::cout << "\n";
                std::cout << "last theta: " << weight[layer] << "\n";
                delete[] input;
                delete[] target;
                file.close();
                break;
            }

            for (int i = 0; i < num - 1; i++)
            {
                target[i] = 0;
            }
            target[num - 1] = 1;
            for (int i = 0; i < num; i++)
            {
                int temp = i;
                for (int j = 0; j < layer; j++)
                {
                    input[j] = temp % 2;
                    temp = temp / 2;
                }
                input[layer] = -1;
                double* w;
                w = learning.backdrop(weight, layer, target[i], input);
                for (int i = 0; i < layer + 1; i++)
                {
                    weight[i] = w[i];
                }
            }



        }


    }
    void OR(double weight[], int layer, int num)
    {
        Output output;
        Learning learning;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> dis(-SCOPE, SCOPE);

        double* input = new double[layer + 1];
        double& theta = weight[layer];
        input[layer] = -1;
        double* target = new double[num];
        using namespace std;

        ofstream file("OR.txt");
        for (int i = 0; i < layer + 1; i++) //Initialize the weight value randomly.
        {
            weight[i] = dis(gen);
        }
        while (true)
        {

            for (int i = 0; i < layer + 1; i++)
            {
                file << "\t" << weight[i];
            }
            file << "\n";
            double& theta = weight[layer];
            bool success = true;
            double net = 0;

            for (int i = 0; i < layer; i++) 
            {
                input[i] = 0;
            }
            net = output.calc(input, weight, layer);
            if (net >= 0.5)
            {
                success = false;
            }

            //The input value is the number of all cases where layers of inputs are 0,1, which is a duplicate permutation of "2π(layer)".Therefore, "2^(layer)=num" inputs exist.
            for (int i = 1; i < num; i++)
            {
                int temp = i;
                for (int j = 0; j < layer; j++)
                {
                    input[j] = temp % 2;
                    temp = temp / 2;
                }

                net = 0;
                net = output.calc(input, weight, layer); 
                if (net < 0.5) 
                {
                    success = false;
                }
            }

            if (success) 
            {
                std::cout << "last weight: ";
                for (int i = 0; i < layer; i++)
                {
                    std::cout << weight[i] << " " << std::endl;
                }
                std::cout << "\n";
                std::cout << "last theta: " << weight[layer] << "\n";
                delete[] input; 
                delete[] target;
                file.close();
                break;
            }

            for (int i = 1; i < num; i++)
            {
                target[i] = 1;
            }
            target[0] = 0;

            for (int i = 0; i < num; i++)
            {
                int temp = i;
                for (int j = 0; j < layer; j++)
                {
                    input[j] = temp % 2;
                    temp = temp / 2;
                }
                input[layer] = -1; 
                double* w;
                w = learning.backdrop(weight, layer, target[i], input);
                for (int i = 0; i < layer + 1; i++)
                {
                    weight[i] = w[i];
                }
            }

        }


    }
    void XOR(double weight[], int layer, int num)
    {
        Output output;
        Learning learning;

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> dis(-SCOPE, SCOPE);

        double* input = new double[layer + 1];
        double& theta = weight[layer];
        input[layer] = -1; 
        double* target = new double[num];
        double** h_weight = new double* [layer]; 
        using namespace std;
        double* h_net = new double[layer + 1];
        h_net[layer] = -1;

        ofstream file("XOR.txt");

        for (int i = 0; i < layer; i++) 
        {
            h_weight[i] = new double[layer + 1];
        }


        for (int i = 0; i < layer + 1; i++)
        {
            weight[i] = dis(gen);
        }
        for (int i = 0; i < layer; i++)
        {
            for (int j = 0; j < layer + 1; j++)
            {
                h_weight[i][j] = dis(gen);
            }
        }
        while (true)
        {
            for (int i = 0; i < layer + 1; i++)
            {
                file << "\t" << weight[i];
            }
            for (int i = 0; i < layer; i++)
            {
                for (int j = 0; j < layer + 1; j++)
                {
                    file << "\t" << h_weight[i][j];
                }
            }
            file << "\n";



            double& theta = weight[layer];
            bool success = true;
            double net = 0; 

            for (int i = 0; i < layer; i++) 
            {
                input[i] = 0;
            }
            for (int i = 0; i < layer; i++)
            {
                h_net[i] = output.calc(input, h_weight[i], layer);
            }
            net = output.calc(h_net, weight, layer); 
            if (net >= 0.5)
            {
                success = false;
            }
            for (int i = 0; i < layer; i++) 
            {
                input[i] = 1;
            }
            for (int i = 0; i < layer; i++)
            {
                h_net[i] = output.calc(input, h_weight[i], layer);
            }
            net = output.calc(h_net, weight, layer); // net를 계산합니다.
            if (net >= 0.5)
            {
                success = false;
            }

            //The input value is the number of all cases where layers of inputs are 0,1, which is a duplicate permutation of "2π(layer)".Therefore, there are "2 ^ (layer) = num" input.
            for (int i = 1; i < num - 1; i++)
            {
                int temp = i;
                for (int j = 0; j < layer; j++)
                {
                    input[j] = temp % 2;
                    temp = temp / 2;
                }

                net = 0;
                for (int i = 0; i < layer; i++)
                {
                    h_net[i] = output.calc(input, h_weight[i], layer);
                }
                net = output.calc(h_net, weight, layer); 
                if (net < 0.5)
                {
                    success = false;
                }
            }

            if (success) 
            {
                std::cout << "last weight: ";
                for (int i = 0; i < layer; i++)
                {
                    std::cout << weight[i] << " " << std::endl;
                }
                std::cout << "\n";
                std::cout << "last theta: " << weight[layer] << "\n";

                delete[] input; 
                delete[] target;
                delete[] h_net;
                file.close();
                for (int i = 0; i < layer; i++)
                {
                    delete[] h_weight[i];
                }
                delete[] h_weight;
                break;
            }
            for (int i = 0; i < num - 1; i++)
            {
                target[i] = 1;
            }
            target[num - 1] = 0;

            for (int i = 0; i < num; i++)
            {
                int temp = i;
                for (int j = 0; j < layer; j++)
                {
                    input[j] = temp % 2;
                    temp = temp / 2;
                }
                input[layer] = -1;
                double* w;
                w = learning.backdrop(h_weight[0], layer, target[i], input);
                for (int j = 0; j < layer + 1; j++)
                {
                    h_weight[0][j] = w[j];
                }
            }

            for (int i = 1; i < num - 1; i++)
            {
                target[i] = 0;
            }
            target[0] = 1;
            for (int i = 0; i < num; i++)
            {
                int temp = i;
                for (int j = 0; j < layer; j++)
                {
                    input[j] = temp % 2;
                    temp = temp / 2;
                }
                input[layer] = -1;
                double* w;
                for (int j = 1; j < layer; j++)
                {
                    w = learning.backdrop(h_weight[j], layer, target[i], input);
                    for (int k = 1; k < layer + 1; k++)
                    {
                        h_weight[j][k] = w[k];
                    }
                }

            }
            for (int i = 1; i < num - 1; i++)
            {
                target[i] = 1;
            }
            target[0] = 0;
            target[num - 1] = 0;
            for (int i = 0; i < num; i++)
            {
                switch (i)
                {
                case 0:
                    for (int k = 0; k < layer; k++)
                    {
                        input[k] = 1;
                    }
                    break;
                default:
                    for (int k = 1; k < layer; k++)
                    {
                        input[k] = 0;
                    }
                    input[0] = 1;
                    break;
                }
                if (i == (num - 1))
                {
                    for (int k = 0; k < layer; k++)
                    {
                        input[k] = 0;
                    }
                }
                input[layer] = -1;
                double* w;
                w = learning.backdrop(weight, layer, target[i], input); 
                for (int i = 0; i < layer + 1; i++)
                {
                    weight[i] = w[i];
                }
            }



        }


    }
};




int main()
{
    Perc data;
    Gate gate;
    int choice;
    errornum = 0;
    tot_error = 0;
    std::cout << "How many dimensions of layers do you want to configure?:";
    std::cin >> data.layer;

    std::cout << "Which gate are you looking for? (AND=1, OR=2, XOR=3)";
    std::cin >> choice;

    data.num = pow(2, data.layer);
    double* weight = new double[data.layer + 1]; 
    clock_t start, finish;
    double time;
    remove("ERROR.txt");
    switch (choice) 
    case 1:
        start = clock();
        gate.AND(weight, data.layer, data.num);
        finish = clock();
        break;
    case 2:
        start = clock();
        gate.OR(weight, data.layer, data.num);
        finish = clock();
        break;
    case 3:
        start = clock();
        is_XOR = true;
        gate.XOR(weight, data.layer, data.num);
        finish = clock();
        break;
    case 0:
        break;
    }

    time = double(finish - start);
    std::cout << "Running time is " << (time / CLOCKS_PER_SEC) << "seconds.\n" << std::endl;
    delete[] weight;
    system("pause");
}

