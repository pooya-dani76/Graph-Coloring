from numpy import random
import numpy as np
import os
import sys
import msvcrt
import networkx as nx
import matplotlib.pyplot as plt

##################################### Default Values ###################################
InitialPopulationCoef = 2  # Value of c in N=cn formula for number of initial population
P_Crossover = 0.8  # Probebility Of Crossover
P_Mutation = 0.1  # Probebility Of Mutation
NumberOfGenerations = 200
Population = []  # Population List
InterMediatePopulation = []
Colors = ['red' , 'blue' , 'green' , 'yellow' , 'purple' , 'brown' , 'grey' , 'orange' , 'black' ]


class Chromosome:   # Storage a Chromosome with its Objective Function
    """
    Implements a Chromosome with its Objective Function
    """
    ObjectiveFunction = 0
    def __init__(self, chromosome):
        self.chromosome = chromosome

def ProbebilityAutheticator(Number):  # Detect An Entered Value is Between 0 and 1 Or not
    Authenticate = True
    if Number > 1 or Number < 0:
        Authenticate = False           

    return Authenticate            

def IntNumberAutheticator(Number):  # Detect An Entered Value is integer or not
    Authenticate = True
    for i in Number:
        if ord(i) > 57 or ord(i) < 48:
            Authenticate = False
            break
    return Authenticate

def MatrixElementsAuthenticator(Element):  # Check The Entered Matrix Value is Valid Or Not
    """
    This Function Check The Entered Matrix Value is Valid Or Not (0 , 1 is Valid)
    """
    if IntNumberAutheticator(Element) == True:
        if Element == '1' or Element == '0':
            return True

        else:
            return False
    else:
        print('Invalid Format!\n')

while True:    # Get Number of Vertex
    os.system('cls')
    NumberOfVetexTemp = input('Enter Number of Vertex of Graph : \n')
    if IntNumberAutheticator(NumberOfVetexTemp) == True:
        if int(NumberOfVetexTemp) > 0:
            NumberOfVetex = int(NumberOfVetexTemp)
            break
        else:
            print('Error : A Graph Cannot has {} Vertex!\n'.format(NumberOfVetexTemp))
            msvcrt.getch()
    else:
        print('Invalid Format For Number of Vertex!\n Number of Vertex Should be an Integer Value...\n')
        msvcrt.getch()

while True:    # Get Number of Colors
    os.system('cls')
    NumberOfColorsTemp = input('Enter Number Of Exist Color : \n')
    if IntNumberAutheticator(NumberOfColorsTemp) == True:
        if int(NumberOfColorsTemp) > 0:
            NumberOfColors = int(NumberOfColorsTemp)
            break
        else:
            print('Error : Color Number Cannot be {}!\n'.format(NumberOfColorsTemp))
            msvcrt.getch()
    else:
        print('Invalid Format For Number of Colors!\n Number of Colors Should be an Integer Value...\n')
        msvcrt.getch()

while True:    # Get Initial Population Coef
    os.system('cls')
    print('Do You Want Enter Desire Initial Population Coef or Use Default Value for It?\n')
    print('1.Enter Initial Population Coef Manually\n')
    print('Press Other key To Continue With Default Value(Default Value = {})\n'.format(InitialPopulationCoef))
    Selector = msvcrt.getch()
    if Selector.decode("utf-8") == '1':
        while True:   
            InitialPopulationCoefTemp = input('Enter Initial Population Coef : \n')
            if IntNumberAutheticator(InitialPopulationCoefTemp) == True:
                if int(InitialPopulationCoefTemp) > 0:
                    InitialPopulationCoef = int(InitialPopulationCoefTemp)
                    break
                else:
                    print('Error : Initial Population Coef Cannot be {}!\n'.format(InitialPopulationCoefTemp))
                    msvcrt.getch()
            else:
                print('Invalid Initial Population Coef!\n Initial Population Coef Should be an Integer Value...\n')
                msvcrt.getch()
        break   

    else: 
        break    

while True:    # Get Number Of Generations
    os.system('cls')
    print('Do You Want Enter Number Of Generations or Use Default Value for It?\n')
    print('1.Enter Number Of Generations Manually\n')
    print('Press Other Key To Continue With Default Value(Default Value = {})\n'.format(NumberOfGenerations))
    Selector = msvcrt.getch()
    if Selector.decode("utf-8") == '1':
        while True:   
            NumberOfGenerationsTemp = input('Enter Number Of Generations : \n')
            if IntNumberAutheticator(NumberOfGenerationsTemp) == True:
                if int(NumberOfGenerationsTemp) > 0:
                    NumberOfGenerations = int(NumberOfGenerationsTemp)
                    break
                else:
                    print('Error : Number Of Generations Cannot be {}!\n'.format(NumberOfGenerationsTemp))
                    msvcrt.getch()
            else:
                print('Invalid Number Of Generations!\n Number Of Generations Should be an Integer Value...\n')
                msvcrt.getch()
        break   

    else: 
        break 

while True:    # Get Probebility Of Crossover
    os.system('cls')
    print('Do You Want Enter Probebility Of Crossover or Use Default Value for It?\n')
    print('1.Enter Probebility Of Crossover Manually\n')
    print('Press Other Key To Continue With Default Value(Default Value = {})\n'.format(P_Crossover))
    Selector = msvcrt.getch()
    if Selector.decode("utf-8") == '1':
        while True: 
            P_CrossoverTemp = input('Enter Probebility Of Crossover  : \n')
            if ProbebilityAutheticator(float(P_CrossoverTemp)) == True:
                P_Crossover = float(P_CrossoverTemp)
                break
            else:
                print('Error : Probebility Of Crossover  Cannot be {}!\n'.format(P_CrossoverTemp))
                msvcrt.getch()
        break 
    else:
        break 

while True:    # Get Probebility Of Mutation
    os.system('cls')
    print('Do You Want Enter Probebility Of Mutation or Use Default Value for It?\n')
    print('1.Enter Probebility Of Mutation Manually\n')
    print('Press Other Key To Continue With Default Value(Default Value = {})\n'.format(P_Mutation))
    Selector = msvcrt.getch()
    if Selector.decode("utf-8") == '1':
        while True: 
            P_MutationTemp = input('Enter Probebility Of Mutation  : \n')
            if ProbebilityAutheticator(float(P_MutationTemp)) == True:
                P_Mutation = float(P_MutationTemp)
                break
            else:
                print('Error : Probebility Of Mutation  Cannot be {}!\n'.format(P_MutationTemp))
                msvcrt.getch()
        break
    else:
        break  

# Storage best Chromosome
Best = Chromosome(random.randint(NumberOfColors , size=NumberOfVetex))
Best.ObjectiveFunction = ((NumberOfVetex * (NumberOfVetex-1)) / 2) + NumberOfColors

# Create Adjacency Matrix for Graph
AdjacencyMatrix = np.zeros((NumberOfVetex, NumberOfVetex))
os.system('cls')

if NumberOfColors >= NumberOfVetex:
    print('You Can Assign A Color To Each Vertex...\n')
    print('No Need To A Genetic Algorithm ;)\n')
    msvcrt.getch()
    sys.exit()

##################################### Functions ###################################
def NumberOfEdges():   # Calculate Number of Edges
    count = 0
    for i in range(0, NumberOfVetex):
        for j in range(i+1, NumberOfVetex):
            if AdjacencyMatrix[i][j] == 1:
                count += 1
    return count

def VertexAuthenticator(Vertex):  # Check The Entered Vertex Number is Exist Or Not
    """
    This Function Check The Entered Vertex Number is Exist Or Not
    """
    if Vertex < NumberOfVetex and Vertex >= 0:
        return True

    else:
        return False

def Get_A_Connection():  # Get , Autheticate & Add A Connection To Neighborhood Matrix
    """
    This Function Gets , Autheticate & Add A Connection To Neighborhood Matrix
    """
    print('Enter 2 Vertex Number has been Connected:\n')
    while True:
        OriginTemp = input('Origin Vertex Number : \n')
        if IntNumberAutheticator(OriginTemp) == True:
            if VertexAuthenticator(int(OriginTemp)) == True:
                Origin = int(OriginTemp)
                break
            else:
                print('Your Entered Vertex Number is Invalid!\n')
                print('Valid Vertex Numbers is:\n')
                for a in range(NumberOfVetex):
                    print('{} \t'.format(a))
                print('\n\n')
                msvcrt.getch()
                os.system('cls')
        else:
            print('Invalid Format!\n')
            msvcrt.getch()
            os.system('cls')
    while True:
        DestinationTemp = input('Destination Vertex Number : \n')
        if IntNumberAutheticator(DestinationTemp) == True:
            if VertexAuthenticator(int(DestinationTemp)) == True:
                Destination = int(DestinationTemp)
                break
            else:
                print('Your Entered Vertex Number is Invalid!\n')
                print('Valid Vertex Numbers is:\n')
                for a in range(NumberOfVetex):
                    print('{} \t'.format(a))
                print('\n\n')
                msvcrt.getch()
                os.system('cls')
        else:
            print('Invalid Format!\n')
            msvcrt.getch()
            os.system('cls')

    if VertexAuthenticator(Origin) == True and VertexAuthenticator(Destination) == True:
        AdjacencyMatrix[Origin][Destination] = 1
        AdjacencyMatrix[Destination][Origin] = 1
        print('Connection {}->{}  Added! \n'.format(Origin, Destination))
        print('Press Any Key To Continue...')
        msvcrt.getch()

def GetConnections():  # Provide A Menu To Add Connections To The Input Graph
    """
    This Function Provide A Menu To Add Connections To The Input Graph
    """
    while True:
        os.system('cls')
        print('Select One :\n')
        print('1.Add Connection\n')
        print('2.Continue Algorithm...\n')

        Selector = msvcrt.getch()
        os.system('cls')

        if Selector.decode("utf-8") == '1':
            Get_A_Connection()

        elif Selector.decode("utf-8") == '2':
            break
        else:
            print('Invalid Option!\n')
            msvcrt.getch()

def GetAdjacencyMatrix():   # Gets Upper Triangle of Neighborhood Matrix as A Graph
    """
    Gets Upper Triangle of Neighborhood Matrix as A Graph
    """
    print('Enter Neighborhood Matrix(Enter Just Upper Triangle of The Matrix):\n')
    for i in range(0, NumberOfVetex):
        for j in range(i+1, NumberOfVetex):
            while True:
                Temp = input('Enter Element ({} , {})\n'.format(i, j))
                if MatrixElementsAuthenticator(Temp) == True:
                    AdjacencyMatrix[i][j] = int(Temp)
                    break
                else:
                    print(
                        'Error : Neighborhood Matrix Elements Can be 0(Not Connected 2 Vertex) or 1(Connected 2 Vertex)\n')

def ColorNumbersOF(chromosome):  # Find Number Of Colors Used in a Solution
    Colors = []
    for i in chromosome:
        flag = 0
        for j in Colors:
            if j == i:
                flag = 1
        if flag == 0:
            temp = i
            Colors.append(temp) 
    return Colors.__len__()        

def ObjectiveFunction(Chromosome):  # Calculate Conflicts Two Neighbor Vertex
    Count = 0
    for i in range(0, NumberOfVetex):
        for j in range(i+1, NumberOfVetex):
            if AdjacencyMatrix[i][j] == 1:
                if Chromosome[i] == Chromosome[j]:
                    Count += 1
    return Count

def CreateInitialPopulation():
    for i in range(InitialPopulationCoef * NumberOfVetex):
        Chrom = Chromosome(random.randint(NumberOfColors, size=NumberOfVetex))
        Chrom.ObjectiveFunction = ObjectiveFunction(Chrom.chromosome) + ColorNumbersOF(Chrom.chromosome)
        Population.append(Chrom) 

def Selection():  # Implement FPS For Selection
    """
    Implements FSP Method 
    """ 
    InterMediatePopulation.clear()
    RouletteWheel = [0]
    ChanceBound = 0
    SumOFs = 0
    M = NumberOfEdges() + 1
    for i in Population:
        SumOFs = SumOFs + i.ObjectiveFunction
    Total = ((NumberOfVetex * InitialPopulationCoef) * M) - SumOFs
    for j in Population:
        SelectionChance = ChanceBound + ((M - j.ObjectiveFunction)/Total)
        ChanceBound = SelectionChance
        RouletteWheel.append(SelectionChance)
    for K in range(NumberOfVetex * InitialPopulationCoef):
        RollDice = random.rand()
        for m in range(NumberOfVetex * InitialPopulationCoef):
            if RollDice >= RouletteWheel[m] and RollDice < RouletteWheel[m+1]:
                TempChromosome = Chromosome(Population[m].chromosome)
                TempChromosome.ObjectiveFunction = Population[m].ObjectiveFunction
                InterMediatePopulation.append(TempChromosome)
                break

def Crossover():  # Implement 1.pt Crossover
    """
    Implement 1.pt Crossover
    """
    ChromosomeCount = 0
    while ChromosomeCount < InitialPopulationCoef*NumberOfVetex:
        IsCrossover = random.rand()
        if (InitialPopulationCoef*NumberOfVetex) % 2 == 1 and ChromosomeCount == (InitialPopulationCoef*NumberOfVetex) - 1:
            break
        if IsCrossover <= P_Crossover:
            CutPoint = random.randint(NumberOfVetex)
            Temp1 = InterMediatePopulation[ChromosomeCount].chromosome.copy()
            Temp2 = InterMediatePopulation[ChromosomeCount+1].chromosome.copy()
            InterMediatePopulation[ChromosomeCount].chromosome = np.concatenate((Temp1[:CutPoint+1] , Temp2[CutPoint+1:]))
            InterMediatePopulation[ChromosomeCount+1].chromosome = np.concatenate((Temp2[:CutPoint+1] , Temp1[CutPoint+1:]))
  
        ChromosomeCount += 2         

def Mutation():   # Implement Mutation with Change a Gene
    """
    Implement Mutation
    """
    for i in range(NumberOfVetex * InitialPopulationCoef):
        IsMutation = random.rand()
        if IsMutation < P_Mutation:
            MutationPoint = random.randint(NumberOfVetex)
            while True:
                newGene = random.randint(NumberOfColors)
                if newGene != InterMediatePopulation[i].chromosome[MutationPoint]:
                    break
            InterMediatePopulation[i].chromosome[MutationPoint] = newGene

def Presevation(Generation):  # Implement Elitism by Presevation Method
    """
    Implement Elitism by Presevation Method
    """
    for i in Generation:
        if i.ObjectiveFunction <= Best.ObjectiveFunction:
            if i.ObjectiveFunction == Best.ObjectiveFunction and i.ObjectiveFunction - ColorNumbersOF(i.chromosome) < Best.ObjectiveFunction - ColorNumbersOF(Best.chromosome) :
                Best.chromosome = i.chromosome
                Best.ObjectiveFunction = i.ObjectiveFunction 
            if i.ObjectiveFunction < Best.ObjectiveFunction:   
                Best.chromosome = i.chromosome
                Best.ObjectiveFunction = i.ObjectiveFunction      

def Replacment():  # Implement Generational + Elitism
    """
    Implement Generational + Elitism
    """
    Population.clear()
    for i in InterMediatePopulation:
        Temp = Chromosome(i.chromosome)
        Temp.ObjectiveFunction = ObjectiveFunction(i.chromosome) + ColorNumbersOF(i.chromosome)
        Population.append(Temp)

def GraphDraw(chromosome): # Draw Result Colored Graph
    """
    Draw Result Colored Graph 
    """
    Graph = nx.Graph()
    GraphColors = []

    for k in range(NumberOfVetex):
        Graph.add_node(k)

    for i in range(0, NumberOfVetex):
        for j in range(i+1, NumberOfVetex):
            if AdjacencyMatrix[i][j] == 1:
                Graph.add_edge(i , j)          

    for m in chromosome:
        GraphColors.append(Colors[m]) 

    nx.draw(Graph , node_color = GraphColors , with_labels = True , node_size = 700)  
    plt.show()      
            
def Menu():  # Provide Menu Of The Program
    """
    Provide Menu Of The Program
    """
    os.system('cls')
    print('{}×{} Neighborhood Matrix Created!\n\n'.format(
        NumberOfVetex, NumberOfVetex))
    print('Enter Graph : \n\n')
    print('1.Enter Graph by Neighborhood Matrix\n')
    print('2.Enter Graph by Define Connections')

    Selector = msvcrt.getch()
    os.system('cls')

    if Selector.decode("utf-8") == '1':
        GetAdjacencyMatrix()

    elif Selector.decode("utf-8") == '2':
        GetConnections()

    else:
        print('Invalid Option!\n')

##################################### Main Program ################################### 
Menu()
CreateInitialPopulation()

for GenerationNumber in range(NumberOfGenerations):

    os.system('cls')
    print('Generation {}'.format(GenerationNumber)) 
    for i in Population:
        print(i.chromosome)  
            
    Presevation(Population)

    Selection()

    Crossover()

    Mutation()

    Replacment()

print('Answer Is :\n')
for i in range(NumberOfVetex):
    print('Vertex {} = {}'.format(i , Colors[Best.chromosome[i]]))

if Best.ObjectiveFunction - ColorNumbersOF(Best.chromosome) > 0:
    print('Sorry I Cannot Find a Better Solution or There is No Any Better Solution!\n')

GraphDraw(Best.chromosome)
