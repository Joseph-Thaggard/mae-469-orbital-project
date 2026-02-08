# Package imports

#import matplotlib 
 
# HIGHLY simplistic to not account for input files yet
# Minimum descriptors of orbits using predefined characteristics
# Eventually use "add_object" or some other method to add an object 
#   class instead of manually defining them. 
# Use lists instead to give characteristics (or lists of lists)

## Simulation parameters

dt = 1
csv_path = "/Users/joe/Desktop/MAE-469"

# Define simulation as sun-centered cartesian with Sol at [0,0,0]
# List syntax: Name, mass, position, velocity
# Mass unit is kg
Sol = ["Sol",1.989e30,[0,0,0],[0,0,0]]
Mercury = ["Mercury",3.301e23,[],[]] 
Venus = ["Venus",4.86e24,[],[]]
Earth = ["Earth",5.972e24,[],[]]
Mars = ["Mars",6.42e23,[],[]]
Jupiter = ["Jupiter",1.898e27,[],[]] 
Saturn = ["Saturn",5.683e26,[],[]]
Neptune = ["Neptune",1.024e26,[],[]]
Uranus = ["Uranus",8.68e25,[],[]]
Spaceship = ["Spaceship",1000,[],[]]

# Compile simulation objects

sim_objects = [Sol,Mercury,Venus,Earth,Mars,Jupiter,Saturn,Neptune,Uranus,Spaceship]

# Force calculator is distance-based just for gravity. 
# Considers gravity for all objects on all objects

def calculate_force(sim_objects):
    # Load properties for force calculation
    positions = []
    mass = []
    for object in sim_objects:
        positions.append(object[2])
        mass.append(object[1])
    print(positions)
    print(mass)
    # Calculate forces
    for object in sim_objects:
        print(object[0])
        #mass = object[1]
        #print(mass)
        initial_position = object[2]
        xi = initial_position[0]
        yi = initial_position[1]
        zi = initial_position[2]
        #print(initial_position)
        dr = [] # 3D distance by √(dx^2+dy^2+dz^2)
        for i in range(1,len(mass)):
            print(mass[i])
            print(positions[i])

    return sim_objects


def diagnostics(sim_objects):
    print(sim_objects)

calculate_force(sim_objects)
