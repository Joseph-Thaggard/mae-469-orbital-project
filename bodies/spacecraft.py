class Spacecraft: 
    def __init__(self,name,mass,position,velocity, burn_schedule=[]):
        self.name=name
        self.mass=mass
        self.position=position
        self.velocity=velocity
        self.mu=6.67430e-11*mass
        self.burn_schedule=burn_schedule
        

