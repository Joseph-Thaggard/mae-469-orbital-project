class Body: 
    def __init__(self,name,mass,radius,position,velocity):
        self.name=name
        self.mass=mass
        self.radius=radius
        self.position=position
        self.velocity=velocity
        self.mu=6.67430e-11*mass
        

