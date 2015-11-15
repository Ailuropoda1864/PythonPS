# Problem Set 3: Simulating the Spread of Disease and Virus Population Dynamics 

import numpy
import random
import pylab
import copy



''' 
Begin helper code
'''

class NoChildException(Exception):
    """
    NoChildException is raised by the reproduce() method in the SimpleVirus
    and ResistantVirus classes to indicate that a virus particle does not
    reproduce. You can use NoChildException as is, you do not need to
    modify/add any code.
    """

'''
End helper code
'''

#
# PROBLEM 2
#
# Enter your definitions for the SimpleVirus and Patient classes in this box.
class SimpleVirus(object):

    """
    Representation of a simple virus (does not model drug effects/resistance).
    """
    def __init__(self, maxBirthProb, clearProb):
        """
        Initialize a SimpleVirus instance, saves all parameters as attributes
        of the instance.        
        maxBirthProb: Maximum reproduction probability (a float between 0-1)        
        clearProb: Maximum clearance probability (a float between 0-1).
        """
        self.maxBirthProb = float(maxBirthProb)
        self.clearProb = float(clearProb)


    def doesClear(self):
        """ Stochastically determines whether this virus particle is cleared from the
        patient's body at a time step. 
        returns: True with probability self.getClearProb and otherwise returns
        False.
        """

        return random.random() < self.clearProb

    
    def reproduce(self, popDensity):
        """
        Stochastically determines whether this virus particle reproduces at a
        time step. Called by the update() method in the Patient and
        TreatedPatient classes. The virus particle reproduces with probability
        self.maxBirthProb * (1 - popDensity).
        
        If this virus particle reproduces, then reproduce() creates and returns
        the instance of the offspring SimpleVirus (which has the same
        maxBirthProb and clearProb values as its parent).         

        popDensity: the population density (a float), defined as the current
        virus population divided by the maximum population.         
        
        returns: a new instance of the SimpleVirus class representing the
        offspring of this virus particle. The child should have the same
        maxBirthProb and clearProb values as this virus. Raises a
        NoChildException if this virus particle does not reproduce.               
        """

        if random.random() < self.maxBirthProb * (1 - popDensity):
            return SimpleVirus(self.maxBirthProb, self.clearProb)
        else:
            raise NoChildException



class Patient(object):
    """
    Representation of a simplified patient. The patient does not take any drugs
    and his/her virus populations have no drug resistance.
    """    

    def __init__(self, viruses, maxPop):
        """
        Initialization function, saves the viruses and maxPop parameters as
        attributes.

        viruses: the list representing the virus population (a list of
        SimpleVirus instances)

        maxPop: the maximum virus population for this patient (an integer)
        """

        self.viruses = viruses
        self.maxPop = maxPop


    def getTotalPop(self):
        """
        Gets the size of the current total virus population. 
        returns: The total virus population (an integer)
        """

        return len(self.viruses)


    def update(self):
        """
        Update the state of the virus population in this patient for a single
        time step. update() should execute the following steps in this order:
        
        - Determine whether each virus particle survives and updates the list
        of virus particles accordingly.   
        
        - The current population density is calculated. This population density
          value is used until the next call to update() 
        
        - Based on this value of population density, determine whether each 
          virus particle should reproduce and add offspring virus particles to 
          the list of viruses in this patient.                    

        returns: The total virus population at the end of the update (an
        integer)
        """
        #step 1
        viruses_copy = self.viruses[:]
        for virus in viruses_copy:
            if virus.doesClear():
                self.viruses.remove(virus)
                
        #step 2
        popDensity = self.getTotalPop()/float(self.maxPop)
        
        #step 3
        viruses_copy2 = self.viruses[:]
        for virus in viruses_copy2:
            try:
                offspring = virus.reproduce(popDensity)
                self.viruses.append(offspring)
            except NoChildException:
                continue
                
        return self.getTotalPop()


#
# PROBLEM 3
#
def simulationWithoutDrug(numViruses, maxPop, maxBirthProb, clearProb,
                          numTrials):
    """
    Run the simulation and plot the graph for problem 3 (no drugs are used,
    viruses do not have any drug resistance).    
    For each of numTrials trial, instantiates a patient, runs a simulation
    for 300 timesteps, and plots the average virus population size as a
    function of time.

    numViruses: number of SimpleVirus to create for patient (an integer)
    maxPop: maximum virus population for patient (an integer)
    maxBirthProb: Maximum reproduction probability (a float between 0-1)        
    clearProb: Maximum clearance probability (a float between 0-1)
    numTrials: number of simulation runs to execute (an integer)
    """
    Patients = [Patient([SimpleVirus(maxBirthProb, clearProb) \
                for v in range(numViruses)], maxPop) \
                for t in range(numTrials)]
    num_timesteps = 300
    num_virus = [0 for n in range(num_timesteps)]
    for patient in Patients:
        for timestep in range(num_timesteps):
            num_virus[timestep] += patient.update()
    for timestep in range(num_timesteps):
        num_virus[timestep] /= float(numTrials)
    pylab.plot(num_virus, label='virus')
    pylab.xlabel('Number of time steps')
    pylab.ylabel('Average size of the virus population')
    pylab.title('Simulation (no drug treatment)')
    pylab.legend(loc = 'best')
    pylab.show()
    
## Uncomment this line to run the simulation          
#simulationWithoutDrug(100, 1000, 0.1, 0.05,100)


#
# PROBLEM 4
#
class ResistantVirus(SimpleVirus):
    """
    Representation of a virus which can have drug resistance.
    """   

    def __init__(self, maxBirthProb, clearProb, resistances, mutProb):
        """
        Initialize a ResistantVirus instance, saves all parameters as attributes
        of the instance.

        maxBirthProb: Maximum reproduction probability (a float between 0-1)       

        clearProb: Maximum clearance probability (a float between 0-1).

        resistances: A dictionary of drug names (strings) mapping to the state
        of this virus particle's resistance (either True or False) to each drug.
        e.g. {'guttagonol':False, 'srinol':False}, means that this virus
        particle is resistant to neither guttagonol nor srinol.

        mutProb: Mutation probability for this virus particle (a float). This is
        the probability of the offspring acquiring or losing resistance to a drug.
        """

        SimpleVirus.__init__(self, maxBirthProb, clearProb)
        self.resistances = resistances
        self.mutProb = float(mutProb)


    def isResistantTo(self, drug):
        """
        Get the state of this virus particle's resistance to a drug. This method
        is called by getResistPop() in TreatedPatient to determine how many virus
        particles have resistance to a drug.       

        drug: The drug (a string)

        returns: True if this virus instance is resistant to the drug, False
        otherwise.
        """
        
        try:
            return self.resistances[drug]
        except:
            return False # assume no resistance if drug is not in resistances dictionary


    def reproduce(self, popDensity, activeDrugs):
        """
        Stochastically determines whether this virus particle reproduces at a
        time step. Called by the update() method in the TreatedPatient class.

        A virus particle will only reproduce if it is resistant to ALL the drugs
        in the activeDrugs list. For example, if there are 2 drugs in the
        activeDrugs list, and the virus particle is resistant to 1 or no drugs,
        then it will NOT reproduce.

        Hence, if the virus is resistant to all drugs
        in activeDrugs, then the virus reproduces with probability:      

        self.maxBirthProb * (1 - popDensity).                       

        If this virus particle reproduces, then reproduce() creates and returns
        the instance of the offspring ResistantVirus (which has the same
        maxBirthProb and clearProb values as its parent). The offspring virus
        will have the same maxBirthProb, clearProb, and mutProb as the parent.

        For each drug resistance trait of the virus (i.e. each key of
        self.resistances), the offspring has probability 1-mutProb of
        inheriting that resistance trait from the parent, and probability
        mutProb of switching that resistance trait in the offspring.       

        For example, if a virus particle is resistant to guttagonol but not
        srinol, and self.mutProb is 0.1, then there is a 10% chance that
        that the offspring will lose resistance to guttagonol and a 90%
        chance that the offspring will be resistant to guttagonol.
        There is also a 10% chance that the offspring will gain resistance to
        srinol and a 90% chance that the offspring will not be resistant to
        srinol.

        popDensity: the population density (a float), defined as the current
        virus population divided by the maximum population       

        activeDrugs: a list of the drug names acting on this virus particle
        (a list of strings).

        returns: a new instance of the ResistantVirus class representing the
        offspring of this virus particle. The child should have the same
        maxBirthProb and clearProb values as this virus. Raises a
        NoChildException if this virus particle does not reproduce.
        """

        for drug in activeDrugs:
            if not self.isResistantTo(drug):
                raise NoChildException
        if random.random() < self.maxBirthProb * (1 - popDensity):
            offspring_resistances = copy.deepcopy(self.resistances)
            for drug in self.resistances.keys():
                if random.random() < self.mutProb: # mutation occurs (i.e., not resistant -> resistant, resistant -> not resistant)
                    offspring_resistances[drug] = not self.resistances[drug]
            return ResistantVirus(self.maxBirthProb, self.clearProb, offspring_resistances, self.mutProb)
        else:
            raise NoChildException         


            

class TreatedPatient(Patient):
    """
    Representation of a patient. The patient is able to take drugs and his/her
    virus population can acquire resistance to the drugs he/she takes.
    """

    def __init__(self, viruses, maxPop):
        """
        Initialization function, saves the viruses and maxPop parameters as
        attributes. Also initializes the list of drugs being administered
        (which should initially include no drugs).              

        viruses: The list representing the virus population (a list of
        virus instances)

        maxPop: The  maximum virus population for this patient (an integer)
        """

        Patient.__init__(self, viruses, maxPop)
        self.Rx = []


    def addPrescription(self, newDrug):
        """
        Administer a drug to this patient. After a prescription is added, the
        drug acts on the virus population for all subsequent time steps. If the
        newDrug is already prescribed to this patient, the method has no effect.

        newDrug: The name of the drug to administer to the patient (a string).

        postcondition: The list of drugs being administered to a patient is updated
        """

        if newDrug not in self.Rx:
            self.Rx.append(newDrug)



    def getResistPop(self, drugResist):
        """
        Get the population of virus particles resistant to the drugs listed in
        drugResist.       

        drugResist: Which drug resistances to include in the population (a list
        of strings - e.g. ['guttagonol'] or ['guttagonol', 'srinol'])

        returns: The population of viruses (an integer) with resistances to all
        drugs in the drugResist list.
        """
        num = 0
        for virus in self.viruses:
            for drug in drugResist:
                if not virus.isResistantTo(drug):
                    break
            else:
                num += 1
        return num
                


    def update(self):
        """
        Update the state of the virus population in this patient for a single
        time step. update() should execute these actions in order:

        - Determine whether each virus particle survives and update the list of
          virus particles accordingly

        - The current population density is calculated. This population density
          value is used until the next call to update().

        - Based on this value of population density, determine whether each 
          virus particle should reproduce and add offspring virus particles to 
          the list of viruses in this patient.
          The list of drugs being administered should be accounted for in the
          determination of whether each virus particle reproduces.

        returns: The total virus population at the end of the update (an
        integer)
        """

        #step 1
        viruses_copy = self.viruses[:]
        for virus in viruses_copy:
            if virus.doesClear():
                self.viruses.remove(virus)
                
        #step 2
        popDensity = self.getTotalPop() / float(self.maxPop)
        
        #step 3
        viruses_copy2 = self.viruses[:]
        for virus in viruses_copy2:
            try:
                offspring = virus.reproduce(popDensity, self.Rx)
                self.viruses.append(offspring)
            except NoChildException:
                continue
                
        return self.getTotalPop()


##test 
#virus1 = ResistantVirus(1.0, 0.0, {"drug1": True}, 0.0)
#virus2 = ResistantVirus(1.0, 0.0, {"drug1": False, "drug2": True}, 0.0)
#virus3 = ResistantVirus(1.0, 0.0, {"drug1": True, "drug2": True}, 0.0)
#patient = TreatedPatient([virus1, virus2, virus3], 100)
#print patient.getResistPop(['drug1']), "2"
#print patient.getResistPop(['drug2']), "2"
#print patient.getResistPop(['drug1','drug2']), "1"
#print patient.getResistPop(['drug3']), "0"
#print patient.getResistPop(['drug1', 'drug3']), "0"
#print patient.getResistPop(['drug1','drug2', 'drug3']), "0"

#
# PROBLEM 5
#
def simulationWithDrug(numViruses, maxPop, maxBirthProb, clearProb, resistances,
                       mutProb, numTrials):
    """
    Runs simulations and plots graphs for problem 5.

    For each of numTrials trials, instantiates a patient, runs a simulation for
    150 timesteps, adds guttagonol, and runs the simulation for an additional
    150 timesteps.  At the end plots the average virus population size
    (for both the total virus population and the guttagonol-resistant virus
    population) as a function of time.

    numViruses: number of ResistantVirus to create for patient (an integer)
    maxPop: maximum virus population for patient (an integer)
    maxBirthProb: Maximum reproduction probability (a float between 0-1)        
    clearProb: maximum clearance probability (a float between 0-1)
    resistances: a dictionary of drugs that each ResistantVirus is resistant to
                 (e.g., {'guttagonol': False})
    mutProb: mutation probability for each ResistantVirus particle
             (a float between 0-1). 
    numTrials: number of simulation runs to execute (an integer)
    
    """

    Patients = [TreatedPatient([ResistantVirus(maxBirthProb, clearProb, resistances, mutProb) \
                for v in range(numViruses)], maxPop) \
                for t in range(numTrials)]
    timesteps_before_drug = 150
    timesteps_after_drug = 150
    total_timesteps = timesteps_before_drug + timesteps_after_drug
    tot_virus = [0 for n in range(total_timesteps)]
    guttagonol_res_virus = tot_virus[:]
    
    for patient in Patients:
        #before drug administration
        for timestep in range(timesteps_before_drug):
            tot_virus[timestep] += patient.update()
            guttagonol_res_virus[timestep] += patient.getResistPop(['guttagonol'])
        #administrate drug    
        patient.addPrescription('guttagonol')
        #after drug administration
        for timestep in range(timesteps_after_drug):
            tot_virus[timesteps_before_drug + timestep] += patient.update()
            guttagonol_res_virus[timesteps_before_drug + timestep] += patient.getResistPop(['guttagonol'])
        
    for timestep in range(total_timesteps):
        tot_virus[timestep] /= float(numTrials)
        guttagonol_res_virus[timestep] /= float(numTrials)
        
    pylab.plot(tot_virus, label='total virus')
    pylab.plot(guttagonol_res_virus, label='guttagonol-resistant virus')
    pylab.xlabel('Number of time steps')
    pylab.ylabel('Average size of the virus population')
    pylab.title('Simulation (treatment with guttagonol)')
    pylab.legend(loc = 'best')
    pylab.show()

## Uncomment this line to run the simulation 
#simulationWithDrug(100, 1000, 0.1, 0.05, {'guttagonol': False}, 0.005, 50)


# After completion of problem set
# Just for fun: adding more parameters to simulationWithDrug
def simulationWithDrug_v2(numViruses, maxPop, maxBirthProb, clearProb, resistances,
                       mutProb, numTrials, stepBeforeDrug = 150, stepAfterDrug = 150,
                       Rx = ['guttagonol']):
    """
    Runs simulations and plots graphs for problem 5.

    For each of numTrials trials, instantiates a patient, runs a simulation for
    stepBeforeDrug timesteps, adds Rx, and runs the simulation for an additional
    stepAfterDrug timesteps.  At the end plots the average virus population size
    (for both the total virus population and the drug(s)-resistant virus
    population) as a function of time.

    numViruses: number of ResistantVirus to create for patient (an integer)
    maxPop: maximum virus population for patient (an integer)
    maxBirthProb: Maximum reproduction probability (a float between 0-1)        
    clearProb: maximum clearance probability (a float between 0-1)
    resistances: a dictionary of drugs that each ResistantVirus is resistant to
                 (e.g., {'guttagonol': False})
    mutProb: mutation probability for each ResistantVirus particle
             (a float between 0-1). 
    numTrials: number of simulation runs to execute (an integer)
    stepBeforeDrug: number of timesteps of simulation to run 
                    before drug administration (an integer)
    stepAfterDrug: number of timesteps of simulation to run 
                    after drug administration (an integer)
    Rx: a list of drugs to be prescribed
    
    """

    Patients = [TreatedPatient([ResistantVirus(maxBirthProb, clearProb, resistances, mutProb) \
                for v in range(numViruses)], maxPop) \
                for t in range(numTrials)]
    total_timesteps = stepBeforeDrug + stepAfterDrug
    virus_count = [0 for n in range(total_timesteps)]
    virus = {'total virus': virus_count[:]}
    for drug in Rx:
        virus[drug] = virus_count[:]
    
    for patient in Patients:
        #before drug administration
        for timestep in range(stepBeforeDrug):
            virus['total virus'][timestep] += patient.update()
            for drug in Rx:
                virus[drug][timestep] += patient.getResistPop(Rx)
        #administrate drug    
        for drug in Rx:
            patient.addPrescription(drug)
        #after drug administration
        for timestep in range(stepAfterDrug):
            virus['total virus'][stepBeforeDrug + timestep] += patient.update()
            for drug in Rx:
                virus[drug][stepBeforeDrug + timestep] += patient.getResistPop(Rx)
        
    for timestep in range(total_timesteps):
        for value in virus.values():
            value[timestep] /= float(numTrials)
    
    for key in virus.keys():        
        pylab.plot(virus[key], label = key)
    pylab.xlabel('Number of time steps')
    pylab.ylabel('Average size of the virus population')
    pylab.title('Simulation with treatment')
    pylab.legend(loc = 'best')
    pylab.show()

## Uncomment the lines below to run the simulation 
#simulationWithDrug_v2(100, 1000, 0.1, 0.05, {'guttagonol': False}, 0.005, 50)
#simulationWithDrug_v2(100, 1000, 0.1, 0.05, {'drug1': False, 'drug2': True, 'drug3': True}, 
#                        0.005, 30, stepAfterDrug = 500, Rx = ['drug2', 'drug3'])