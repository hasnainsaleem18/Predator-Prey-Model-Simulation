# PredatorPrey.py
# Program to run simulation of the predator-prey relationship

def PredatorPrey(DT = 0.001, simLength = 12):
    numIterations = int(simLength/DT) + 1
    t = 0

    predator_population = 15
    predator_birth_fraction = 0.01
    predator_death_proportionality_constant = 1.06
    prey_population = 100
    prey_birth_fraction = 2
    prey_death_proportionality_constant = 0.02

    predator_births = (predator_birth_fraction * prey_population) * predator_population
    predator_deaths = predator_death_proportionality_constant * predator_population

    prey_births = prey_birth_fraction * prey_population
    prey_deaths = (prey_death_proportionality_constant * predator_population) * prey_population

    tLst = [t]
    predatorLst = [predator_population]
    preyLst = [prey_population]
    for i in range(1, numIterations):
        t = i * DT
        prey_population = prey_population + (prey_births - prey_deaths) * DT    #new population 
        predator_population = predator_population + (predator_births - predator_deaths) * DT #new population

        prey_births = prey_birth_fraction * prey_population   #slope
        prey_deaths = (prey_death_proportionality_constant * predator_population) * prey_population

        predator_births = (predator_birth_fraction * prey_population) * predator_population
        predator_deaths = predator_death_proportionality_constant * predator_population

        tLst.append(t)
        predatorLst.append(predator_population)
        preyLst.append(prey_population)

    return tLst, predatorLst, preyLst

def RK2(DT = 0.001, simLength = 12):
    numIterations = int(simLength/DT) + 1
    t = 0

    predator_population = 15
    predator_birth_fraction = 0.01
    predator_death_proportionality_constant = 1.06
    prey_population = 100
    prey_birth_fraction = 2
    prey_death_proportionality_constant = 0.02

    tLst = [t]
    predatorLst = [predator_population]
    preyLst = [prey_population]
    for i in range(1, numIterations):
        t = i * DT
       # slopes_1
        predator_births_1 = (predator_birth_fraction * prey_population) * predator_population  
        predator_deaths_1 = predator_death_proportionality_constant * predator_population

        prey_births_1 = prey_birth_fraction * prey_population
        prey_deaths_1 = (prey_death_proportionality_constant * predator_population) * prey_population

        # new population_1
        predator_population_1 = predator_population + (predator_births_1 - predator_deaths_1) * DT
        prey_population_1 = prey_population + (prey_births_1 - prey_deaths_1) * DT

        # slopes_2
        preadtor_births_2 = (predator_birth_fraction * prey_population_1) * predator_population_1
        predator_deaths_2 = predator_death_proportionality_constant * predator_population_1
        prey_briths_2 = prey_birth_fraction * prey_population_1
        prey_deaths_2 = (prey_death_proportionality_constant * predator_population_1) * prey_population_1

        # final new population
        predator_population = (predator_population) + ((predator_births_1 + preadtor_births_2 - predator_deaths_1 - predator_deaths_2) * DT) / 2
        prey_population = (prey_population) + ((prey_births_1 + prey_briths_2 - prey_deaths_1 - prey_deaths_2) * DT) / 2


        tLst.append(t)
        predatorLst.append(predator_population)
        preyLst.append(prey_population)

    return tLst, predatorLst, preyLst


def RK4(DT=0.001, simLength=12):
    numIterations = int(simLength / DT) + 1
    t = 0

    predator_population = 15
    predator_birth_fraction = 0.01
    predator_death_proportionality_constant = 1.06
    prey_population = 100
    prey_birth_fraction = 2
    prey_death_proportionality_constant = 0.02

    tLst = [t]
    predatorLst = [predator_population]
    preyLst = [prey_population]

    for i in range(1, numIterations):
        t = i * DT

        # Slopes_1
        predator_1 = DT * (predator_birth_fraction * prey_population * predator_population - predator_death_proportionality_constant * predator_population)
        prey_1 = DT * (prey_birth_fraction * prey_population - prey_death_proportionality_constant * predator_population * prey_population)

        # population_1 at midpoint_1
        predator_mid1 = predator_population + predator_1 / 2
        prey_mid1 = prey_population + prey_1 / 2

        #slopes_2
        predator_2 = DT * (predator_birth_fraction * prey_population * predator_mid1 - predator_death_proportionality_constant * predator_mid1)
        prey_2 = DT * (prey_birth_fraction * prey_mid1 - prey_death_proportionality_constant * predator_mid1 * prey_mid1)

        # population_2 at midpoint2
        predator_mid2 = predator_population + predator_2 / 2
        prey_mid2 = prey_population + prey_2 / 2

        #slopes_3
        predator_3 = DT * (predator_birth_fraction * prey_population * predator_mid2 - predator_death_proportionality_constant * predator_mid2)
        prey_3 = DT * (prey_birth_fraction * prey_mid2 - prey_death_proportionality_constant * predator_mid2 * prey_mid2)

        # New population_3
        predator_end = predator_population + predator_3
        prey_end = prey_population + prey_3

        # Slopes_4
        predator_4 = DT * (predator_birth_fraction * prey_population * predator_end - predator_death_proportionality_constant * predator_end)
        prey_4 = DT * (prey_birth_fraction * prey_end - prey_death_proportionality_constant * predator_end * prey_end)

        # Final new population
        predator_population += (predator_1 + 2 * predator_2 + 2 * predator_3 + predator_4) / 6
        prey_population += (prey_1 + 2 * prey_2 + 2 * prey_3 + prey_4) / 6

        tLst.append(t)
        predatorLst.append(predator_population)
        preyLst.append(prey_population)

    return tLst, predatorLst, preyLst
