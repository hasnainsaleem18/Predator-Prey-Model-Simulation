# sharkCompetition.py
# Program to run simulation of species competition between
# white tip sharks (WTS) and black tip sharks (BTS)

def sharkCompetition(DT = 0.001, simLength = 5):
    numIterations = int(simLength/DT) + 1
    t = 0

    WTS_population = 20
    BTS_population = 15

    WTS_birth_fraction = 1
    WTS_death_proportionality_constant = 0.27
    WTS_births = WTS_population * WTS_birth_fraction
    WTS_deaths = (WTS_death_proportionality_constant * BTS_population) * WTS_population

    BTS_birth_fraction = 1
    BTS_death_proportionality_constant = 0.2
    BTS_births = BTS_birth_fraction * BTS_population
    BTS_deaths = (BTS_death_proportionality_constant * WTS_population)*BTS_population

    tLst = [t]
    WTSLst = [WTS_population]
    BTSLst = [BTS_population]
    for i in range(1, numIterations):
        t = i * DT
        BTS_population = BTS_population + (BTS_births - BTS_deaths) * DT
        WTS_population = WTS_population + (WTS_births - WTS_deaths) * DT

        BTS_births = BTS_birth_fraction * BTS_population
        BTS_deaths = (BTS_death_proportionality_constant * WTS_population) * BTS_population

        WTS_births = WTS_population * WTS_birth_fraction
        WTS_deaths = (WTS_death_proportionality_constant * BTS_population) * WTS_population
        tLst.append(t)
        WTSLst.append(WTS_population)
        BTSLst.append(BTS_population)

    return tLst, WTSLst, BTSLst


def RK2(DT=0.001, simLength=5):
    numIterations = int(simLength/DT) + 1
    t = 0

    WTS_population = 20
    BTS_population = 15

    WTS_birth_fraction = 1
    WTS_death_proportionality_constant = 0.27

    BTS_birth_fraction = 1
    BTS_death_proportionality_constant = 0.2

    tLst = [t]
    WTSLst = [WTS_population]
    BTSLst = [BTS_population]

    for i in range(1, numIterations):
        t = i * DT

        #slopes_1
        WTS_births_1 = WTS_population * WTS_birth_fraction
        WTS_deaths_1 = (WTS_death_proportionality_constant * BTS_population) * WTS_population
        BTS_births_1 = BTS_birth_fraction * BTS_population
        BTS_deaths_1 = (BTS_death_proportionality_constant * WTS_population) * BTS_population

        #new population_1
        WTS_population_1 = WTS_population + (WTS_births_1 - WTS_deaths_1) * DT
        BTS_population_1 = BTS_population + (BTS_births_1 - BTS_deaths_1) * DT

        #slopes_2
        WTS_births_2 = WTS_population_1 * WTS_birth_fraction
        WTS_deaths_2 = (WTS_death_proportionality_constant * BTS_population_1) * WTS_population_1
        BTS_births_2 = BTS_birth_fraction * BTS_population_1
        BTS_deaths_2 = (BTS_death_proportionality_constant * WTS_population_1) * BTS_population_1

        #final new population
        WTS_population += (((WTS_births_1 - WTS_deaths_1) + (WTS_births_2 - WTS_deaths_2)) * DT) / 2
        BTS_population += (((BTS_births_1 - BTS_deaths_1) + (BTS_births_2 - BTS_deaths_2)) * DT) / 2

        tLst.append(t)
        WTSLst.append(WTS_population)
        BTSLst.append(BTS_population)

    return tLst, WTSLst, BTSLst


def RK4(DT=0.001, simLength=5):
    numIterations = int(simLength / DT) + 1
    t = 0

    WTS_population = 20
    BTS_population = 15

    WTS_birth_fraction = 1
    WTS_death_proportionality_constant = 0.27

    BTS_birth_fraction = 1
    BTS_death_proportionality_constant = 0.2

    tLst = [t]
    WTSLst = [WTS_population]
    BTSLst = [BTS_population]

    for i in range(1, numIterations):
        t = i * DT
        dt_2 = DT / 2.0

        #populaton_1
        WTS_1 = WTS_birth_fraction * WTS_population - WTS_death_proportionality_constant * BTS_population * WTS_population
        BTS_1 = BTS_birth_fraction * BTS_population - BTS_death_proportionality_constant * WTS_population * BTS_population

        #populaton_2
        WTS_2 = WTS_birth_fraction * (WTS_population + WTS_1 * dt_2) - WTS_death_proportionality_constant * (BTS_population + BTS_1 * dt_2) * (WTS_population + WTS_1 * dt_2)
        BTS_2 = BTS_birth_fraction * (BTS_population + BTS_1 * dt_2) - BTS_death_proportionality_constant * (WTS_population + WTS_1 * dt_2) * (BTS_population + BTS_1 * dt_2)

        #populaton_3
        WTS_3 = WTS_birth_fraction * (WTS_population + WTS_2 * dt_2) - WTS_death_proportionality_constant * (BTS_population + BTS_2 * dt_2) * (WTS_population + WTS_2 * dt_2)
        BTS_3 = BTS_birth_fraction * (BTS_population + BTS_2 * dt_2) - BTS_death_proportionality_constant * (WTS_population + WTS_2 * dt_2) * (BTS_population + BTS_2 * dt_2)

        #populaton_4
        WTS_4 = WTS_birth_fraction * (WTS_population + WTS_3 * DT) - WTS_death_proportionality_constant * (BTS_population + BTS_3 * DT) * (WTS_population + WTS_3 * DT)
        BTS_4 = BTS_birth_fraction * (BTS_population + BTS_3 * DT) - BTS_death_proportionality_constant * (WTS_population + WTS_3 * DT) * (BTS_population + BTS_3 * DT)

        #final population
        WTS_population += (WTS_1 + 2 * WTS_2 + 2 * WTS_3 + WTS_4) * (DT / 6)
        BTS_population += (BTS_1 + 2 * BTS_2 + 2 * BTS_3 + BTS_4) * (DT / 6)

        tLst.append(t)
        WTSLst.append(WTS_population)
        BTSLst.append(BTS_population)

    return tLst, WTSLst, BTSLst









