# Graph-Coloring
Task is coloring a graph with n vertices with minimum number of color possible

## 🔧 Requirements:
      pip install -r Requirements.txt

## Representation
  Arrays with length n that every element represents a vertex calls chromosome
  Inside of each chromosome element placed a number between 0 & maximum number of color is accessable (detemined by user)
  
## Inputs:
### Graph
  There is 2 ways for giving graph input to program:
   - Adjacency Matrix 
   - Getting Connections one by one

### Initial Population coefficient
  Number of population in every generation is multiple of number of vertices.
  
  In this input you determine this coefficient.
  
  Number of population = Population coefficient × Number of vertices

### Number of Generation
  An integer value -> default is 200

### Crossover Rate
  Probability of crossover occur -> default is 0.8
  
### Mutation Rate 
  Probability of mutation occur -> default is 0.1
   
  ⚠️ For prevent error occur , every input should be autheticate by program!
   
## Objective Function(OF):  
   Objevtive function for this problem has 2 parts:
   - Any 2 connected vertices should not have same color.
   - Number of colors used should be minimum
   So , for Implement this purpose we propose this OF :
   
        > OF = (Number Of conflicts) + (Number of colors used)
        
   And OF should be minimum   
   
   OF for every chromosome has been stored with itself in a class
   
## Initiat Population:
   Generate chromosomes randomly 
   
   As mentioned above number of population in every generation is multiple of number of vertices.
   
## Selection:
  We used roulette wheel method for implement selection 
  
  Selection chance of every chromosome is depend on its OF (Fittn)
  
   > Selection Chance = (Chromosome OF) / (All  Population OF)
  
  All selected chromosome append to intermediate population

## Crossover:
  We used 1-pt crossover in this algorithm
  
  Cut point has been determine randomly
  
## Mutation:
  Mutation point has been determine randomly too.
  
  Another valid value will be replace to the mutation point 
  
## Replacement:
  We use Generational + Elitism method.
  
  Store best chromosome of previous population (if better from previous best) and replace new population with old population.
  
  This method garantees find global minima
  
## Output:
  If Number of color be greater than or be equal with number of vertices , so you no need to run this algorithm & program knows this!
  
  In otherwise program determines color for each vertex in the best possible state.
  
  If there is impossible coloring with no conflict , program shows best state with minimum conflict.
  
  Output shows you best chromosome and colored corresponding graph.
  
## Example:
  Coloring graph with 6 vertices with default values (mentioned above)
    
  <p align = "center">
  <img src = "https://github.com/pooya-dani76/Graph-Coloring/blob/main/Examples/E1-C.PNG">
  </p>
   
   
  <p align = "center">
  <img src = "https://github.com/pooya-dani76/Graph-Coloring/blob/main/Examples/E1-G.png">
  </p>
  
## Run Example:
  - Run `EC.py` file
  - Enter number of vertices
  - Enter algorithm parameter(such as crossove rate & etc)
  - Enter graph
  - Get answer & enjoy!
  
Good Luck Guys! :trollface: 
    
