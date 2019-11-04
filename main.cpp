#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <random>
#include <iomanip>

const int NB_CHROMOSOME = 32;
const int TAILLE_CHROMOSOME = 20;
const int PROBA_MUTATION = 10;
const double M_PI = 3.141592653589793;
const int NB_ESSAI = 62500;

struct Particule
{
    int index;
    double fitness;
};

double f(const std::vector<double> &d)
{
	int j;
	double top = 0;

	for (j = 0; j<d.size(); j++)
	{
		top = top + (pow(d[j], (double)2) - 10 * cos(2 * M_PI*d[j]) + 10);
	}
	return top;
}

int f2(const std::vector<bool> &d)
{
    int sum = 0;
    int mult = 1;

    for(size_t i = 0; i < d.size(); i++)
    {
        if(d[i] == 0)
            sum += i+1;
        else
            mult *= (i+1);
    }

    return fabs(sum-36)+fabs(mult-360);
}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void genererChromosome(std::vector<double> &chromosome)
{
    chromosome.clear();
    for(int i = 0; i < TAILLE_CHROMOSOME; i++)
    {
        chromosome.push_back(fRand(-5.12, 5.12));
    }
}

std::vector<Particule> genererParticules(const std::vector<std::vector<double>> &population)
{
    std::vector<Particule> particules;
    for(int i = 0; i < population.size(); i++)
    {
        Particule p;
        p.index = i;
        p.fitness = f(population[i]);

        particules.push_back(p);
    }

    return particules;
}

std::vector<std::vector<double>> genererPopulation()
{
    std::vector<std::vector<double>> population;
    population.clear();
    for(int i = 0; i < NB_CHROMOSOME; i++)
    {
        population.push_back(std::vector<double>{});
        genererChromosome(population[i]);
    }

    return population;
}

void afficherPopulation(const std::vector<std::vector<double>> &pop)
{
    for(unsigned int i = 0; i < pop.size(); i++)
    {
        std::cout << i << ": ";
        for(unsigned int j = 0; j < pop[i].size(); j++)
        {
            std::cout << pop[i][j];
        }
        std::cout << std::endl;
    }
}


void crossover(const std::vector<double> &gagnant, std::vector<double> &perdant)
{
    int positionCroisement = rand() % TAILLE_CHROMOSOME;
    for(size_t i = positionCroisement; i < TAILLE_CHROMOSOME; i++)
    {
        perdant[i] = gagnant[i];
    }
}

void mutate(std::vector<double> &chromosome)
{
    double random = (rand() % 101); // Value between 0 - 100
    if(random <= PROBA_MUTATION)
    {
        int positionMutation = rand() % TAILLE_CHROMOSOME;
        chromosome[positionMutation] = fRand(-5.12, 5.12);
    }
}

bool team1Gagnante(const Particule &team1, const Particule &team2)
{
    return team1.fitness <= team2.fitness ? true : false;
}

void JouerTour(std::vector<std::vector<double>> &population, const std::vector<Particule> &particules)
{
    if(particules.size() > 1)
    {
        std::vector<Particule> Gagnant;
        for(int i = 0; i < particules.size() - 1; i+=2)
        {
            if(team1Gagnante(particules[i], particules[i+1]))
            {
                Gagnant.push_back(particules[i]);
                crossover(population[particules[i].index], population[particules[i+1].index]);
                mutate(population[particules[i+1].index]);
            }
            else
            {
                Gagnant.push_back(particules[i+1]);
                crossover(population[particules[i+1].index], population[particules[i].index]);
                mutate(population[particules[i].index]);
            }
        }
        JouerTour(population, Gagnant);
    }
}

double LCA ()
{
    std::vector<std::vector<double>> pop = genererPopulation();
    std::vector<Particule> particules;

    double m;
    int index;
    int essai = 0;

    do
    {
        essai++;
        particules = genererParticules(pop);
        JouerTour(pop, particules);

        m = particules[0].fitness;
        index = 0;
        for(int i = 1; i < particules.size(); i++)
        {
            if(m>particules[i].fitness)
            {
                index=i;
                m=particules[i].fitness;
            }
        }
    }while(essai != NB_ESSAI);

    for(int i = 0; i < pop[index].size(); i++)
    {
        std::cout << pop[index][i] << " ";
    }
    std::cout << std::endl << "Fitness: " << std::setprecision(10) << m << std::endl;
    return m;
}

int main()
{
    srand(time(0));
    double somme = 0.0;
    std::vector<double> rslt;
    for(int i = 0; i < 30; i++)
    {
        double x = LCA();
        rslt.push_back(x);
        somme += x;
    }

    double moyenne = somme / 30.0;
    std::cout << std::setprecision(10) << moyenne << std::endl;

    double somme2 = 0.0;
    for(int i = 0; i < 30; i++)
    {
        somme2 += (rslt[i] - moyenne)*(rslt[i] - moyenne);
    }
    somme2 /= 30;
    std::cout << std::setprecision(10) << sqrt(somme2) << std::endl;

    return 0;
}
