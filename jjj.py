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

        # Slopes_1
        WTS_births_1 = WTS_population * WTS_birth_fraction
        WTS_deaths_1 
        WTS_deaths_1 = (WTS_death_proportionality_constant * BTS_population) * WTS_population
        BTS_births_1 = BTS_birth_fraction * BTS_population
        BTS_deaths_1 = (BTS_death_proportionality_constant * WTS_population) * BTS_population

        # Slopes at the midpoint of the interval
        WTS_mid = WTS_population + 0.5 * DT * WTS_births_1
        BTS_mid = BTS_population + 0.5 * DT * BTS_births_1
        WTS_births_mid = WTS_mid * WTS_birth_fraction
        WTS_deaths_mid = (WTS_death_proportionality_constant * BTS_mid) * WTS_mid
        BTS_births_mid = BTS_birth_fraction * BTS_mid
        BTS_deaths_mid = (BTS_death_proportionality_constant * WTS_mid) * BTS_mid

        # Slopes at the end of the interval
        WTS_end = WTS_population + DT * WTS_births_mid
        BTS_end = BTS_population + DT * BTS_births_mid
        WTS_births_end = WTS_end * WTS_birth_fraction
        WTS_deaths_end = (WTS_death_proportionality_constant * BTS_end) * WTS_end
        BTS_births_end = BTS_birth_fraction * BTS_end
        BTS_deaths_end = (BTS_death_proportionality_constant * WTS_end) * BTS_end

        # Update populations using the RK4 method
        WTS_population += DT / 6 * (WTS_births_1 + 2 * WTS_births_mid + 2 * WTS_births_end)
        BTS_population += DT / 6 * (BTS_births_1 + 2 * BTS_births_mid + 2 * BTS_births_end)

        tLst.append(t)
        WTSLst.append(WTS_population)
        BTSLst.append(BTS_population)

    return tLst, WTSLst, BTSLst




