#Justin Ryan-Final
#Import required libraries
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import wx 

#Create a class "charge" to hold the properties of our moving particle
class Charge:
    def __init__(self):
        '''Charged particle with velocity, charge, and mass.
        '''
        #The initial location of the charge will be the origin
        self.location = np.array([0,0,0])
        #Collect the magnitude of the velocity from our Speed dialogue box
        speed = float(initialVelocityMagnitude.GetValue())*10**6
        #Depending on the velocity orientation chosen, apply this speed value
        #To a component of velocity and assign that velocity to the charge
        if initialVelocityDirection.GetSelection() == 0: #Perpindicular
            self.velocity = np.array([speed,speed,0.])
        elif initialVelocityDirection.GetSelection() == 1: #Parallel
            self.velocity = np.array([0.,0.,speed])
        elif initialVelocityDirection.GetSelection() == 2: #45 Degrees
            self.velocity = np.array([speed,speed,speed])
        #Based on the Choice of charge, assign a mass and charge to particle
        if chargeType.GetSelection() == 0:
            self.charge = 1.602e-19
            self.mass = 9.11e-31
        elif chargeType.GetSelection() == 1:
            self.charge = -1.602e-19
            self.mass = 9.11e-31
        elif chargeType.GetSelection() == 2:
            self.charge = float(customCharge.GetValue())
            self.mass = float(customMass.GetValue())
            
#Define a force function to be used in acceleration calculation
def force(particle):
    '''Calculates the current magnetic force on the charge as a Vector.
    Args:
        param1 (:obj:'Charge'): Charge object
    '''
    #Useing equation F=q(vxB) where v and B are vectors. Computed using Numpy.Cross
    vCrossB = np.cross(particle.velocity,bField(particle.location))
    return particle.charge*vCrossB
    
#Define a function for acceleration (dv/dt) to be used in Time Step analysis
def acceleration(particle):
    '''Returns a value for acceleration of the particle.
    Args:
        param1 (:obj:'Charge'): Particle undergoing simulation'''
    f = force(particle) #pull force from function
    return f/particle.mass #from Newton's Second Law

#Define a function to compute the magnetic field at a given position        
def bField(position):
    '''Returns the b-field at the current position.
    Args:
        param1 (:obj:'numpy.ndarray'): Position of particle
    '''
    #Assuming constant Magnetic field to achieve cyclotron motion
    bStrength = magStrengthInput.GetValue()/10.0
    return np.array([0.,0.,bStrength]) #Return B as a vector

#Based on our initial conditions for position and velocity, step through trajectory    
def trajectory(particle):
    '''Plots calculates the trajectory of the object using time step estimations.'''
    #Due to high velocity of particles, need very small timesteps
    deltaT = 1e-15
    #Continue stepping through for many steps
    for i in range(20000):
        #r_next = r + dr*dt
        particle.location = particle.location + deltaT*particle.velocity
        #Append this location to the list of locations
        path.append(particle.location)
        dv = acceleration(particle) #Pull acceleration from function
        #v_next = v + dv*dt
        particle.velocity = particle.velocity + dv*deltaT
    
#Create a unit vector formula to be used in determining the angle between the 
#velocity and b-field
def unit_vector(vector):
    '''returns the unit vector of the inputted vector
    Arg:
        param1 (:obj:'numpy.ndarray'): 3D vector
    '''
    return vector/np.linalg.norm(vector)
    
def angle(v1,v2):
    '''returns the angle between two vectors.
    Args:
        param1 (:obj:'numpy.ndarray'): 3D Vector
        param2 (:obj:'numpy.ndarray'): 3D Vector
    '''
    unitv1 = unit_vector(v1)
    unitv2 = unit_vector(v2)
    #From numpy documentation on finding angle between two vectors
    #np.clip ensures that the rounding in dot product doesn't exceed arccos domain
    return np.linalg.norm(np.arccos(np.clip(np.dot(unitv1,unitv2),-1.0,1.0)))
    
def radius(particle):
    '''Calculates the radius of the cyclotron motion for the inputted charge.
    Args:
        param1 (:obj:'numpy.ndarray'): Charge object
    '''
    m = particle.mass
    b=np.linalg.norm(bField(particle.location))
    speed = np.linalg.norm(particle.velocity)
    vperp = speed*np.sin(angle(particle.velocity,b))
    q = particle.charge
    return abs(m*vperp/q/b)
        
#Define a function that will be called upon clicking the "Simulate" Button
def simulate(event):
    #We want the path array to remain global so that it can be cleared and 
    #filled for each simulation that takes place while the program is running
    global path,particle
    path = [np.array([0.,0.,0.])]
    #Initialize the particle
    particle = Charge()
    #Bfield string for plot title
    bString = str(magStrengthInput.GetValue()/10.0)
    #Calculate Trajectory
    trajectory(particle)
    #Plotting
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    #Plots the start point of the particle
    ax.scatter(0,0,0,c='r',label='Start')
    #Since the "path" list is composed of vectors, use zip to extract the components
    #into 3 separate lists
    x,y,z = zip(*path)
    ax.plot(x,y,z)
    plt.legend(loc='lower left')
    plt.title(chargeType.GetStringSelection()+' Moving '+initialVelocityDirection.GetStringSelection()+' to '+bString+' B Field')
    plt.show()
    #Print the radius of motion
    print('\nRadius of Motion (m): '+str(radius(particle)))

#GUI - initialize application and window
app = wx.App()
win = wx.Frame(None, title='Charge Trajectory',size=(500,400))
win.Show()

#Create individual Window Components
#RadioBox for selecting charge type
charges = ['Positron','Electron','Custom']
chargeType = wx.RadioBox(win,-1,'Charge', pos=(5,25),choices=charges,style=wx.RA_SPECIFY_ROWS)

#Create input boxes that will take in values for mass and charge if the user selects "custom" charge type
customMassText = wx.StaticText(win,-1,'Mass (if custom)', pos=(120,30))
customMass = wx.TextCtrl(win,pos=(120,50))
customChargeText = wx.StaticText(win,-1,'Charge (if custom)', pos=(120,100))
customCharge = wx.TextCtrl(win,pos=(120,120))

#I used a slider for magnetic field strength as a proof-of-concept. Range from small to very strong
magStrengthText = wx.StaticText(win,-1,'Magnetic Field Strength (0-10 T)',pos=(250,30))
magStrengthInput = wx.Slider(win,-1,pos=(250,50))
magStrengthInput.SetRange(0,100)
magStrengthInput.SetValue(100)

#Allow the user to input speed, affects radius of motion
velocityText = wx.StaticText(win,-1,'Speed of Particle (x10^6 m/s)',pos=(5,230))
initialVelocityMagnitude = wx.TextCtrl(win,pos=(5,260),size=(40,25))
initialVelocityMagnitude.SetValue('1')

#Allow the user to select velocity orientation relative to magnetic field. Demonstrates cross product nature
directions = ['Perpindicular','Parallel','45 Degrees']
initialVelocityDirection = wx.RadioBox(win,-1,'Velocity Relative to B', pos=(250,110),choices=directions,style=wx.RA_SPECIFY_ROWS)

#Create simulate button and bind it to the simulate function
simTrajectory = wx.Button(win,label='Simulate',pos=(300,270))
simTrajectory.Bind(wx.EVT_BUTTON, simulate)

#Force the program to continue running until the app is terminated.
app.MainLoop()