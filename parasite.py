import numpy as np

input = [
    {
    "room": 1,
    "grid": [
        [0, 3],
        [0, 1]
        ],
    "interestedIndividuals": ["0,0"]
    },
    {
    "room": 2,
    "grid": [
        [0, 3, 2],
        [0, 1, 1],
        [1, 0, 0]
    ],
    "interestedIndividuals": [
    "0,2", "2,0", "1,2"
    ]
    }
]

VACANT = 0
HEALTHY = 1
VACCINATED = 2
INFECTED = 3

def ortho_neighbours(tup):
    yield (tup[0]+1,tup[1])
    yield (tup[0]  ,tup[1]-1)
    yield (tup[0]  ,tup[1]+1)
    yield (tup[0]-1,tup[1])
    
def all_neighbours(tup):
    for n in ortho_neighbours(tup):
        yield n
    yield (tup[0]+1,tup[1]+1)
    yield (tup[0]+1,tup[1]-1)
    yield (tup[0]-1,tup[1]-1)
    yield (tup[0]-1,tup[1]+1)

def isFullyInfected(grid, infect_time):
    for idx, person in np.ndenumerate(grid):
        if(person==HEALTHY and infect_time[idx]==-1):
            return False
    return True

def simulate(grid, allow_diag=False):
    neighbours =  all_neighbours if allow_diag else ortho_neighbours
    infect_time = grid.copy()
    infect_time.fill(-1)
    infected_set = set()
    
    # set each infected square as 0 in the time array
    for idx,square in np.ndenumerate(grid):
        if (square == INFECTED):
            infect_time[idx] = 0
            infected_set.add(idx)
            
    MAX_TICKS = 100
    for tick in range(1,MAX_TICKS):
        to_be_infected = set()
        for idx in infected_set:
            for neighbour in neighbours(idx):
                try:
                    if(grid[neighbour]==HEALTHY and infect_time[neighbour]==-1 
                       and neighbour[0]>=0 and neighbour[1]>=0): # because negative indices are valid in python
                        infect_time[neighbour] = tick
                        to_be_infected.add(neighbour)
                except IndexError:
                    pass
        if(len(to_be_infected)==0): break
        else: infected_set = infected_set.union(to_be_infected)
    return (infect_time,infected_set)

for room in input:
    # initialize a grid with same dimensions as room grid
    grid = np.array(room["grid"])
    
    infect_time,infected_set = simulate(grid)
    # part 1
    p1 = {}
    for person in room["interestedIndividuals"]:
        p1[person] = infect_time[eval(person)]
    print(f'p1:{p1}')
    
    # part 2
    p2 = -1
    if (isFullyInfected(room,infect_time)): p2 = np.amax(infect_time)
    
    # Part 3
    infect_time, infected_set = simulate(grid, allow_diag=True)
    p3 = -1
    if (isFullyInfected(room,infect_time)): p3 = np.amax(infect_time)
    
    # Part 4
    energy=0
    infect_time, infected_set = simulate(grid)
    while(not isFullyInfected(grid,infect_time)):
        energy+=1
        # expand infection through vacant spaces
        for idx,infected_square in np.ndenumerate(infect_time):
            for neighbour in ortho_neighbours(idx):
                if(neighbour[0]<0 or neighbour[1]<0): continue
                try:
                    grid[idx]=INFECTED
                except IndexError:
                    pass

    output = {
        "room": room["room"],
        "p1": p1,
        "p2": p2,
        "p3": p3,
        "p4": energy
    }
    
    print(output)