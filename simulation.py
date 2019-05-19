"""
2-D Physics Simulation Application
by Alexander Barnhart
Version 1.1

A script that utilizes the pygame and thorpy packages to create an interactive 2-Dimensional
physics simulation. The pygame and thorpy packages must be installed in order for the script
to run on a machine. The script was written for Python 3 and is best suited to run in this
environment. To execute from the command line type: "python3 [filepath]/simulation.py"

Functions:
	pythag: Solves for the hypotenuse of a triangle using the Pythagorean theorem.
	dot_product: Computes the dot product of two vectors.
	vec_subtract: Computes the difference of two vectors.
	point_distance: Finds the distance between two points.
	is_mouse_in_window: Checks if the mouse is inside the simulation window's bounds.
	does_collide: Determines if two balls collide with eachother at a moment. [Used only for object placement].
	do_they_collide: Another function that checks if two balls collide.
	does_hit_ramp: determines if a ball collides with a ramp.
	is_within_ramp: determines if the ball is within a ramp's vertical and horizontal bounds.
	bounce_wall: Calculates change in momentum after a ball collides with a wall.
	bounce_ball: Calculates change in momentum after two balls collide.
	bounce_ramp: calculate change in momentum after a ball collides with a ramp.
	handle_collisions: Checks if two balls are colliding in the moment and calls the appropriate function.
	set_demo_one:	Sets demo 1 for demonstration purposes.
	set_demo_two:	Sets demo 2 for demonstration purposes.
	set_demo_three:	Sets demo 3 for demonstration purposes.
	set_demo_four:	Sets demo 4 for demonstration purposes.
	set_demo_five:	Sets demo 5 for demonstration purposes.
	in_bounds: determines if the mouse is within a rectangle, represented by an array of four points.
	toggle_sim: Switches the sim_on boolean between True and False.
	clear_sim: Removes all objects from the simulation.
	reset_sim: Resets all objects to their original positions.
	placing_obj: Switches the user_placing_object boolean between True and False.
	undo_placement: Removes the most recently added object.
	placing_ball: Tells the program the user is placing a ball.
	placing_ramp: Tells the program the user is placing a ramp.
	draw_sim_buttons: Renders the buttons that controls the simulation and its objects.
	draw_text: Renders the text on the application.
	draw_settings_area: Renders the simulation and object settings windows.
	set_value: handles input checking for the user input fields. 
"""

import pygame, thorpy, sys, math

#-------class definitions-------
class Ball:
	"""
	Class that represents a ball object in the simulation.
	Ball Methods:
		__init__(x,y,mass,rad,vang,vmag,aang,amag,fang,fmag,color)
			- Constructor that sets initial attribute values.	
		set_force(Fx, Fy)
			- Sets the x and y force vectors for the object.
		set_mass(m)
			- Sets the mass value for the object.
		update_values()
			- Updates the values of the object with respect to time.
		set_gravity(grav)
			- Updates the gravity value of the object.
		reset()
			- Resets all the values to their initial states as defined by the constructor.

		Ball Properties:
			_color: The object's color (Red[0-255],Green[0-255],Blue[0-255]) 
			_initx: Initial x position
			_inity: Initial y position
			_x: Current x position
			_y: Current y position
			_vang: Velocity vector angle
			_vmag: Velocity vector magnitude
			_initvx: Initial velocity vector magnitude for x axis
			_initvy: Initlal velocity vector magnitude for y axis
			_vx: Current x axis velocity	
			_vy: Current y axis velocity
			_mass: Object's mass
			_gravity: Current gravity value
			_normal_force: The ball's normal force value
			_aang: Angle of acceleration vector
			_amag: Magnitude of acceleration vector
			_initax: Initial acceleration vector magnitude for x axis
			_initay: Initial acceleration vector magnitude for y axis
			_ax: Current x axis acceleration
			_ay: Current y axis acceleration
			_Fang: Force vector angle
			_Fmag: Force magnitude angle
			_initFx: Initial force vector magnitude for x axis
			_initFy: Initial force vector magnitude for y axis
			_Fx: Current x axis force
			_Fy: Current y axis force
			_dragFy: The y axis force of wind resistance
			_dragFx: The x axis force of wind resistance
			_mu_s: The static friction coefficient
			_mu_k: The kinetic friction coefficient
			_radius: Radius of the circle(ball)

	"""
	def __init__(self, x, y, mass, rad, vang, vmag, aang, amag, fang, fmag, color):
		"""
		__init__(x, y, mass, rad, vang, vmag, aang, amag, fang, fmag, color) constructor method
		Method Arguments:
			x (integer)	-- required
				The object's x position.
			y (integer)	-- required
				The object's y position.
			mass (float) -- required 
				The object's mass.
			rad (integer)	-- required
				The object's radius.
			vang (float) -- required
				The object's velocity vector angle.
			vmag (float) -- required
				The object's velocity vector magnitude.
			aang (float) -- required
				The object's acceleration vector angle.
			amag (float) -- required
				The object's acceleration vector magnitude.
			fang (float) -- required
				The object's force vector angle.
			fmag (float) -- required
				The object's force vector magnitude.
			color (tuple of integers between 0 and 255, ex:(133,0,255)) -- required
				The object's color.
		"""
		self._color = color					
		self._initx = x							
		self._inity = y							
		self._x = self._initx				
		self._y = self._inity	
		self._vang = vang						
		self._vmag = vmag/10						
		self._initvx = self._vmag*math.cos(self._vang*(math.pi/180.0)) 
		self._initvy = self._vmag*math.sin(self._vang*(math.pi/180.0))	
		self._vx = self._initvx			
		self._vy = self._initvy			
		self._mass = mass						
		self._gravity = -9.8				
		self._normal_force = (self._mass * self._gravity) / 100000.0
		self._aang = fang						
		if mass != 0:
			self._amag = (fmag/mass)/100000
		else:
			self._amag = 0 
		self._initax = self._amag*math.cos(self._aang*(3.141592654/180)) 
		self._initay = self._amag*math.sin(self._aang*(3.141592654/180)) 
		self._ax = self._initax			
		self._ay = self._initay			
		self._Fang = fang						
		self._Fmag = fmag						
		self._initFx = self._Fmag*math.cos(self._Fang*(3.141592654/180))	
		self._initFy = self._Fmag*math.sin(self._Fang*(3.141592654/180))	
		self._Fx = self._initFx			
		self._Fy = self._initFy
		self._dragFy = 0.0
		self._dragFx = 0.0			
		self._mu_s = 0.0	
		self._mu_k = 0.0		
		self._radius = rad					

	def set_force(self,Fx, Fy, dragFx, dragFy):		
		"""
		set_force(Fx, Fy) method
		Method Arguments:
			Fx (float) -- required
				Force vector x direction.
			Fy (float) -- required
				Force vector y direction.
		"""
		self._Fx = Fx
		self._Fy = Fy
		self._dragFy = dragFy
		self._dragFx = dragFx
		self._ay = self._ay + (self._Fy/self._mass)
		self._ax = self._ax + (self._Fx/self._mass)

	def set_mass(self,m):			
		"""
		set_mass(m) method
		Method Arguments:
			m (float) -- required
				The object's mass.
		"""
		self._mass = m

	def update_values(self, air_density):	
		"""
		update_values() method
		Method Arguments:
			None
		"""
		self._ax = self._Fx / self._mass
		self._ay = self._Fy / self._mass
		if self._Fx != 0.0:
			self._Fx = 0.0
		if self._Fy != 0.0:
			self._Fy = 0.0
		self._dragFx = (0.5*air_density*(self._vx*self._vx)*0.47*500)/100000 
		self._dragFy = (0.5*air_density*(self._vy*self._vy)*0.47*500)/100000
		self._vx = self._vx + self._ax
		if self._vx > 0:
			self._vx = self._vx - self._dragFx/self._mass
		else:
			self._vx = self._vx + self._dragFx/self._mass
		self._vy = self._vy + (self._ay + self._gravity)
		if self._vy > 0:
			self._vy = self._vy - self._dragFy/self._mass
		else:
			self._vy = self._vy + self._dragFy/self._mass
		if abs(self._vx) < 0.08:
			self._vx = 0.0
		if abs(self._vy) < 0.08:
			self._vy = 0.0
		self._x = self._x + self._vx
		self._y = self._y - self._vy

	def set_gravity(self, grav):	
		"""
		set_gravity(grav) method
		Method Arguments:
			grav (float) -- required
				The value of gravity that affects the object in the simulation.
		"""
		self._gravity = grav/100.0
		self._normal_force = float(self._mass * self._gravity)

	def set_friction(self, mu_s, mu_k):
		self._mu_s = mu_s
		self._mu_k = mu_k

	def reset(self):		
		"""
		reset() method
		Method Arguments:
			None
		"""
		self._x = self._initx
		self._y = self._inity
		self._vx = self._initvx
		self._vy = self._initvy
		self._ax = self._initax
		self._ay = self._initay
		self._Fx = self._initFx
		self._Fy = self._initFy
		self._dragFx = 0.0
		self._dragFy = 0.0

class Ramp:
	"""
	Class that represents a ramp object in the simulation.
	Ball Methods:
		__init__(x,y,height,width,color)
			- Constructor that sets initial attribute values.

		Ramp Properties:
			_color: The object's color (Red[0-255],Green[0-255],Blue[0-255]) 
			_x: Current x position of the right angle
			_y: Current y position of the right angle
			_height: Ramp's vertical height
			_width: Ramp's horizontal width
			_hypot = the length of the hypotenuse
			_angle = the value of the angle of the slope. [arctan(height/width)]
			_slope = the slope value of the hypotenuse
			_pointlist = the three points put into a list. Used for drawing the object.

	"""
	def __init__(self, x, y, height, width, color):
		"""
		__init__(x, y, height, width, color) constructor method
		Method Arguments:
			x (integer)	-- required
				The object's x position.
			y (integer)	-- required
				The object's y position.
			height (integer) -- required
				The ramp's vertical height.
			width (integer) -- required
				The ramp's horizontal width.
			 (3-tuple ocolorf integers between 0 and 255, ex:(133,0,255)) -- required
				The object's color.
		"""
		self._color = color
		self._x = x
		self._y = y
		self._height = height
		self._width = width
		self._hypot = math.hypot(width,height)
		self._angle = math.atan2(height,width)
		self._slope = ((y+height)-y)/(x-(x+width))
		self._pointlist = [(x, y), (x+self._width, y), (x, y-self._height)]

#-------global variables-------
sim_on = False	#if the simulation is running
demos_on = False
user_placing_object = False	#if the user is placing the object
object_held = 0		#current object being placed(0 = ball, 1 = ramp)
clock = pygame.time.Clock()		#pygame clock object used for timing
meter = 100				#what i'm using to measure as a "meter" of pixels
air_density = 1.5
elastic_constant = 0.2
#--------simulation window bounds---------
upper_bound = 50			#floor
lower_bound = 550			#ceiling
left_bound = 50				#left wall
right_bound = 750			#right wall
#-----------------------------------------
balls = []						#list of objects in simulation
ramps = []						#list of ramps in the simulation
ramp_edge = ''				#used for ramp collision
add_buffer = []
gravity = -0.0098			#gravity value, can be adjusted
is_elastic = False		#if the collisions are elastic
radian_const = (math.pi/180.0)	#constant used for radian/degree conversion
click_points = [0,0]				#(x,y) coordinates of the last mouse click
selected_color = (0,0,0)		#color to apply to object about to be placed
#-------functions-------
def pythag(a,b):	
	"""
	Computes and returns the hypotenuse of a triangle using the Pythagorean Theorem.
	Parameters:
		a (float) -- required
			A side of the triangle.
		b (float) -- required
			A side of the triangle.

	Returns:
		The square root of a^2 + b^2. A single float value.
	
	"""
	return math.sqrt((a*a) + (b*b))

def dot_product(vec1, vec2):	
	"""
	Computes the dot product of two vectors.
	Parameters:
		vec1 (tuple(float,float)) -- required
			A vector to multiply.
		vec2 (tuple(float,float)) -- required
			A vector to multiply.

	Returns:
		The dot product of vec1 and vec2. A single float value.
	"""
	return (vec1[0]*vec2[0])+(vec1[1]*vec2[1])

def vec_subtract(vec1,vec2):	
	"""
	Computes the difference of two vectors.
	Parameters:
		vec1 (tuple(float,float)) -- required
			A vector to subtract.
		vec2 (tuple(float,float)) -- required
			A vector to subtract.

	Returns:
		The difference of the two vectors. A tuple(vector) with two elements.
	"""
	return (vec1[0]-vec2[0],vec1[1]-vec2[1])
 
def point_distance(x1, x2, y1, y2):	
	"""
	Computes the distance between two points.
	Parameters:
		x1 (integer) -- required
			The x coordinate of the first point.
		x2 (integer) -- required
			The x coordinate of the second point.
		y1 (integer) -- required
			The y coordinate of the first point.
		y2 (integer) -- required
			The y coordinate of the second point.

	Returns:
		The distance between the two points in space. A single float value.
	"""
	tempx = x2-x1
	tempx = tempx * tempx
	tempy = y2-y1
	tempy = tempy * tempy
	return math.sqrt(tempx + tempy)


def is_mouse_in_window(x, y):
	"""
	Determines if the mouse is within the bounds of the simulation window.
	Parameters:
		x (integer) -- required
			x position of the mouse.
		y (integer) -- required
			y position of the mouse.

	Returns:
		A boolean. True if the mouse is within the window, False if the mouse is not.
	"""
	if (x <= right_bound and x >= left_bound) and (y <= lower_bound and y >= upper_bound):
		return True
	else:
		return False

def does_collide(x, y):
	"""
	Determines if an object collides with any other object. Only used for initial object
  placement.
	Parameters:
		x (integer) -- required
			x position of an object.
		y (integer) -- required
			y position of an object.

	Returns:
		A boolean. True if the object collides with another, False if the object does not.
	"""
	global balls
	if len(balls) == 0:
		return False
	for ball in balls:
		if point_distance(x, ball._x, y, ball._y) < (20 + ball._radius):
			return True
		else:
			continue
	return False

def do_they_collide(ball1, ball2):
	"""
	Determines if one ball collides with another. Used for collision handling.
	Parameters:
		ball1 (Ball) -- required
			A Ball object.
		ball2 (Ball) -- required
			A Ball object.

	Returns:
		A boolean. True if the two objects collide, False if the objects do not.
	"""
	if point_distance(ball1._x, ball2._x, ball1._y, ball2._y) < (ball1._radius + ball2._radius):
		return True
	else:
		return False

def does_hit_ramp(ball, ramp):
	"""
	Determines if a ball collides with a ramp.
	Parameters:
		ball (Ball) -- required
			A Ball object.
		ramp (Ramp) -- required
			A Ramp object.

	Returns:
		True if it collides, False it doesn't.
	"""	
	global ramp_edge
	y = ball._y
	x = ball._x
	m = ramp._slope
	x1 = ramp._x+ramp._width
	x2 = ramp._x
	y1 = ramp._y
	y2 = ramp._y-ramp._height
	A = y1-y2
	B = x2-x1
	C = ((x1-x2)*y1) + (x1*(y2-y1))
	#first, check the three end points
	pt1 = (x1,y1)	#ah
	pt2	= (x2,y2)	#oh
	pt3 = (ramp._x,ramp._y) #oa
	dist = abs((A*x) + (B*y) + C)/math.sqrt(A*A + B*B)
	if is_within_ramp(ball,ramp) == True:
		if point_distance(ball._x,pt1[0],ball._y,pt1[1]) <= ball._radius:
			ramp_edge = 'cah'
			return True
		elif point_distance(ball._x,pt2[0],ball._y,pt2[1]) <= ball._radius:
			ramp_edge = 'coh'
			return True
		elif point_distance(ball._x,pt3[0],ball._y,pt3[1]) <= ball._radius:
			ramp_edge = 'coa'
			return True
		if ball._radius >= dist:
			ramp_edge = 'h'
			return True
		if ball._x - ball._radius < ramp._x and ball._x > ramp._x or ball._x + ball._radius > ramp._x and ball._x < ramp._x:
				ramp_edge = 'o'
				return True
		if ball._y - ball._radius < ramp._y and ball._y > ramp._y or ball._y + ball._radius > ramp._y and ball._y < ramp._y:
				ramp_edge = 'a'
				return True
		else:
			return False
	else:
		return False

def is_within_ramp(ball, ramp):
	"""
	Determines if the ball is within the horizontal and vertical bounds of the ramp. Used to
  calculate collisions.
	Parameters:
		ball (Ball) -- required
			A Ball object.
		ramp (Ramp) -- required
			A Ramp object.

	Returns:
		Nothing.
	"""
	if ramp._height > 0:
		if (ball._x+ball._radius) > min(ramp._x, ramp._x+ramp._width) and (ball._x-ball._radius) < max(ramp._x, ramp._x+ramp._width) and (ball._y - ball._radius) < ramp._y and (ball._y + ball._radius) > (ramp._y - ramp._height):
			return True
		else:
			return False
	elif ramp._height < 0:
		if (ball._x+ball._radius) > min(ramp._x, ramp._x+ramp._width) and (ball._x-ball._radius) < max(ramp._x, ramp._x+ramp._width) and (ball._y + ball._radius) > ramp._y and (ball._y - ball._radius) < (ramp._y-ramp._height):
			return True
		else:
			return False
	
def bounce_wall(ball, wall):	
	"""
	Calculates the resulting momentum for a ball that has collided with a wall. Used for
	collision handling. Updates the velocity of a ball object in accordance with momentum
	conservation.
	Parameters:
		ball (Ball) -- required
			A Ball object.
		wall (integer) -- required
			Numeric representation of one of the four walls 
			(left wall = 1, top wall = 2, right wall = 3, bottom wall = 4).

	Returns:
		Nothing.
	"""
	global elastic_constant
	mov = 0.0
	rad = ball._radius
	FF = ball._normal_force * float(ball._mu_k)
	mass = ball._mass
	#-----------------------------------------------------------------
	if wall == 1:			#hits the leftmost wall
		ball._vx = ball._vx*-1
		mov = left_bound - (ball._x-rad)
		ball._x = ball._x + mov
		if ball._vx >= 0:
			ball._vx = ball._vx - abs(ball._vx * elastic_constant)
		else:
			ball._vx = ball._vx + abs(ball._vx * elastic_constant)
	#-----------------------------------------------------------------
	elif wall == 3:		#hits the rightmost wall
		ball._vx = ball._vx*-1
		mov = (ball._x+rad) - right_bound
		ball._x = ball._x - mov
		if ball._vx >= 0:
			ball._vx = ball._vx - abs(ball._vx * elastic_constant)
		else:
			ball._vx = ball._vx + abs(ball._vx * elastic_constant)
	#-----------------------------------------------------------------
	elif wall == 2:		#hits the top wall (ceiling)
		FFa = FF/mass
		if ball._vx > 0:
			ball._vx = ball._vx - FFa
		else:
			ball._vx = ball._vx + FFa
		ball._vy = ball._vy*-1
		mov = upper_bound - (ball._y-rad)
		ball._y = ball._y + mov
		if ball._vy >= 0:
			ball._vy = ball._vy - abs(ball._vy * elastic_constant)
		else:
			ball._vy = ball._vy + abs(ball._vy * elastic_constant)
	#-----------------------------------------------------------------
	elif wall == 4:		#hits the bottom wall (floor)
		FFa = FF/mass
		if ball._vx > 0:
			ball._vx = ball._vx + FFa
		else:
			ball._vx = ball._vx - FFa
		ball._vy = ball._vy*-1
		mov = (ball._y+rad) - lower_bound
		ball._y = ball._y - mov
		if ball._vy >= 0:
			ball._vy = ball._vy - abs(ball._vy * elastic_constant)
		else:
			ball._vy = ball._vy + abs(ball._vy * elastic_constant)
	if abs(ball._vx) < 0.05:
		ball._vx = 0.0
	if abs(ball._vy) < 0.05:
		ball._vy = 0.0 
	
def bounce_ball(ball1, ball2):	
	"""
	Calculates the resulting momentums of two ball objects after a collision. Updates the
	velocity of both ball objects in accordance with momentum conservation.
	Parameters:
		ball1 (Ball) -- required
			A Ball object.
		ball2 (Ball) -- required
			A Ball object.

	Returns:
		Nothing.

	Additional Notes:
		The following was used as a reference for the equations:
			https://imada.sdu.dk/~rolf/Edu/DM815/E10/2dcollisions.pdf
	"""
	global radian_const,elastic_constant

	m1 = ball1._mass			#ball1's mass
	v1x = ball1._vx				#ball1's x-velocity
	v1y	= ball1._vy				#ball1's y-velocity
	v1 = pythag(v1x,v1y)	#ball1's velocity vector
	x1 = ball1._x					#ball1's x position
	y1 = ball1._y					#ball1's y position
	#-----------------
	m2 = ball2._mass			#ball2's mass
	v2x = ball2._vx				#ball2's x-velocity
	v2y	= ball2._vy				#ball2's y-velocity
	v2 = pythag(v2x,v2y)	#ball2's velocity vector
	x2 = ball2._x					#ball2's x position
	y2 = ball2._y					#ball2's y position
	phi = 0.0
	dist = ball1._radius + ball2._radius
	phi = math.acos((x2-x1)/pythag((x2-x1),(y2-y1)))
	N = ((x2-x1),(y1-y2))
	temp = math.sqrt((N[0]*N[0])+(N[1]*N[1]))
	UN = (N[0]/temp,N[1]/temp)
	UT = (-UN[1],UN[0])
	V1 = (v1x,v1y)
	V2 = (v2x,v2y)
	V1N = dot_product(UN,V1)
	V1T = dot_product(UT,V1)
	V2N = dot_product(UN,V2)
	V2T = dot_product(UT,V2)
	new_V1N = ((V1N*(m1-m2))+(2.0*m2*V2N))/(m1+m2)
	new_V2N = ((V2N*(m2-m1))+(2.0*m1*V1N))/(m1+m2)
	new_V1N = (new_V1N*UN[0],new_V1N*UN[1])
	new_V1T = (V1T*UT[0],V1T*UT[1])
	new_V2N = (new_V2N*UN[0],new_V2N*UN[1])
	new_V2T = (V2T*UT[0],V2T*UT[1])
	new_V1 = (new_V1N[0]+new_V1T[0],new_V1N[1]+new_V1T[1])
	new_V2 = (new_V2N[0]+new_V2T[0],new_V2N[1]+new_V2T[1])	
	
	ball1._vx = new_V1[0]
	ball1._vy = new_V1[1]
	ball2._vx = new_V2[0]
	ball2._vy = new_V2[1]

	if ball1._vx >= 0:
		ball1._vx = ball1._vx - abs(ball1._vx*elastic_constant)
	else:
		ball1._vx = ball1._vx + abs(ball1._vx * elastic_constant)
	if ball1._vy >= 0:
		ball1._vy = ball1._vy - abs(ball1._vy * elastic_constant)
	else:
		ball1._vy = ball1._vy + abs(ball1._vy * elastic_constant)

	if ball2._vx >= 0:
		ball2._vx = ball2._vx - abs(ball2._vx*elastic_constant)
	else:
		ball2._vx = ball2._vx + abs(ball2._vx * elastic_constant)
	if ball2._vy >= 0:
		ball2._vy = ball2._vy - abs(ball2._vy * elastic_constant)
	else:
		ball2._vy = ball2._vy + abs(ball2._vy * elastic_constant)

	if abs(ball1._vx) < 0.05:
		ball1._vx = 0.0
	if abs(ball1._vy) < 0.05:
		ball1._vy = 0.0
	if abs(ball2._vx) < 0.05:
		ball2._vx = 0.0
	if abs(ball2._vy) < 0.05:
		ball2._vy = 0.0

	mov = (ball1._radius+ball2._radius) - pythag((x2-x1),(y2-y1))
	movy = math.sin(phi)*mov
	movx = math.cos(phi)*mov
	if(ball2._y > ball1._y):
		ball2._y = ball2._y + movy+1
	else:
		ball1._y = ball1._y + movy+1

	if(ball2._x > ball1._x):
		ball2._x = ball2._x + movx+1
	else:
		ball1._x = ball1._x + movx+1
	
	FF1 = ball1._normal_force * float(ball1._mu_s)	#static friction force of ball 1
	F1x = (ball1._mass * abs(v1x-new_V1[0]))		#horizontal force of ball 1
	FF2 = ball2._normal_force * float(ball2._mu_s)	#static friction force of ball 2
	F2x = (ball2._mass * abs(v2x-new_V2[0]))	#horizontal force of ball 2
	if v1 == 0.0:		#ball 1 is at rest
		if abs(FF1) > abs(F2x):
			ball1._vx = 0.0
			ball2._vx = -(v2x*elastic_constant)
	elif v2 == 0.0:		#ball 2 is at rest
		if abs(FF2) > abs(F1x):
			ball2._vx = 0.0
			ball1._vx = -(v1x*elastic_constant)
					
def bounce_ramp(ball, ramp):
	"""
	Determines resulting momentum for a ball when it collides with a ramp.
	Parameters:
		ball (Ball) -- required
			A Ball object.
		ramp (Ramp) -- required
			A Ramp object.

	Returns:
		Nothing.
	"""
	global ramp_edge, radian_const
	FF = ball._normal_force * float(ball._mu_k)
	v = (ball._vx, ball._vy)
	vmag = pythag(v[0], v[1])
	vang = math.atan2(v[1], v[0])
	nmag = pythag(ramp._width,ramp._height)
	
	if ramp_edge == 'cah':
		n = (ramp._width/nmag, -ramp._height/nmag)
		v_n = (v[0]*n[0]) + (v[1]*n[1])
		f = (n[0]*v_n, n[1]*v_n)
		f = (2*f[0],2*f[1])
		new_v = (v[0]-f[0], v[1]-f[1])
		ball._vx = new_v[0]	
		ball._vy = new_v[1]
	if ramp_edge == 'coh':
		n = (-ramp._width/nmag, ramp._height/nmag)
		v_n = (v[0]*n[0]) + (v[1]*n[1])
		f = (n[0]*v_n, n[1]*v_n)
		f = (2*f[0],2*f[1])
		new_v = (v[0]-f[0], v[1]-f[1])
		ball._vx = new_v[0]	
		ball._vy = new_v[1]
	if ramp_edge == 'coa':
		n = (-ramp._width/nmag, -ramp._height/nmag)
		v_n = (v[0]*n[0]) + (v[1]*n[1])
		f = (n[0]*v_n, n[1]*v_n)
		f = (2*f[0],2*f[1])
		new_v = (v[0]-f[0], v[1]-f[1])
		ball._vx = new_v[0]	
		ball._vy = new_v[1]
	if ramp_edge == 'h':
		nmag = pythag(ramp._width,ramp._height)
		n = (ramp._height/nmag, ramp._width/nmag)
		v_n = (v[0]*n[0]) + (v[1]*n[1])
		f = (n[0]*v_n, n[1]*v_n)
		f = (2*f[0],2*f[1])
		new_v = (v[0]-f[0], v[1]-f[1])
		ball._vx = new_v[0]	
		ball._vy = new_v[1]
	if ramp_edge == 'o':
		ball._vx = ball._vx*-1
	if ramp_edge == 'a':
		ball._vy = ball._vy*-1
	if ramp._x < ball._x:
		ball._x = ball._x + 1
	else:
		ball._x = ball._x - 1
	if ramp._y < ball._y:
		ball._y = ball._y + 1
	else:
		ball._y = ball._y - 1
 
	if ball._vx >= 0:
		ball._vx = ball._vx - abs(ball._vx * elastic_constant)
	else:
		ball._vx = ball._vx + abs(ball._vx * elastic_constant)
	if ball._vy >= 0:
		ball._vy = ball._vy - abs(ball._vy * elastic_constant)
	else:
		ball._vy = ball._vy + abs(ball._vy * elastic_constant)

	FFa = FF/ball._mass
	if ball._vx > 0:
		ball._vx = ball._vx + FFa
	else:
		ball._vx = ball._vx - FFa
	if ball._vy > 0:
		ball._vy = ball._vy + FFa
	else:
		ball._vy = ball._vy - FFa

	if abs(ball._vx) < 0.05:
		ball._vx = 0.0
	if abs(ball._vy) < 0.05:
		ball._vy = 0.0

def handle_collisions():		
	"""
	Detects collisions of any kind within the simulation. Calls the appropriate functions for
	calculating resulting momentum.
	Parameters:
		None.

	Returns:
		Nothing.
	"""
	global balls, ramps, lower_bound, upper_bound, left_bound, right_bound
	for i in range(0,len(balls)):
		ball = balls[i]
		#first, check if any object is against the walls.
		if ball._y > lower_bound - ball._radius:
			bounce_wall(ball, 4)
		if ball._y < upper_bound + ball._radius:
			bounce_wall(ball, 2)
		if ball._x < left_bound + ball._radius:
			bounce_wall(ball, 1) 
		if ball._x > right_bound - ball._radius:
			bounce_wall(ball, 3)
		#then, check if any balls collide.
		for j in range(i,len(balls)):
			otherball = balls[j]
			if do_they_collide(ball, otherball) == True and i != j:
				bounce_ball(ball,otherball)
		#then, check if the ball collide with ramp.
	for g in range(0,len(balls)):
		ball = balls[g]
		for k in range(0,len(ramps)):
			ramp = ramps[k]
			if does_hit_ramp(ball,ramp) == True:
				bounce_ramp(ball,ramp)

def set_demo_one():
	"""
	Sets up demo 1 in the simulation. Used for demonstration.
	Parameters:
		None.

	Returns:
		Nothing.
	"""
	global balls,ramps,gravity,air_density,sim_on
	balls.clear()
	ramps.clear()
	ball1 = Ball(100, 300, 30, 30.0, 0, 40, 0, 0, 0, 0, Black)
	ball2 = Ball(700, 337, 30, 30.0, 180, 40, 0, 0, 0, 0, Red)
	balls.append(ball1)
	balls.append(ball2)
	gravity_entry.set_value('0.0')
	drag_entry.set_value('0.0')
	elastic_radio.set_value(True)	
	sim_on = False
	

def set_demo_two():
	"""
	Sets up demo 2 in the simulation. Used for demonstration.
	Parameters:
		None.

	Returns:
		Nothing.
	"""
	global balls,ramps,gravity,air_density,sim_on
	balls.clear()
	ramps.clear()
	ball = Ball(100, 300, 20, 20.0, 0, 0, 0, 0, 0, 0, Black)
	ball._mu_k = 0.6
	ramp = Ramp(50,550,100,500,(0,0,0))
	balls.append(ball)
	ramps.append(ramp)
	gravity_entry.set_value('-9.8')
	drag_entry.set_value('0.0')
	elastic_radio.set_value(False)
	sim_on = False

def set_demo_three():
	"""
	Sets up demo 3 in the simulation. Used for demonstration.
	Parameters:
		None.

	Returns:
		Nothing.
	"""
	global balls,ramps,gravity,air_density,sim_on
	balls.clear()
	ramps.clear()
	ball = Ball(100, 100, 20, 20.0, 0, 90, 0, 0, 0, 0, Black)
	ball._mu_k = 0.6
	balls.append(ball)
	gravity_entry.set_value('-9.8')
	drag_entry.set_value('0.0')
	elastic_radio.set_value(False)
	sim_on = False

def set_demo_four():
	"""
	Sets up demo 4 in the simulation. Used for demonstration.
	Parameters:
		None.

	Returns:
		Nothing.
	"""
	global balls,ramps,gravity,air_density,sim_on
	balls.clear()
	ramps.clear()
	gravity_entry.set_value('0.0')
	drag_entry.set_value('0.0')
	elastic_radio.set_value(True)
	ball = Ball(100, 100, 20, 20, 0, 90, 0, 0, 0, 0, Black)
	balls.append(ball)
	ball = Ball(150, 275, 50, 40, 300, 40, 0, 0, 0, 0, Blue)
	balls.append(ball)
	ball = Ball(500, 275, 10, 15, 45, 120, 0, 0, 0, 0, Red)
	balls.append(ball)
	ball = Ball(600, 400, 25, 30, 90, 10, 0, 0, 0, 0, Green)
	balls.append(ball)
	ball = Ball(300, 300, 90, 50, 180, 3, 0, 0, 0, 0, Cyan)
	balls.append(ball)
	sim_on = False

def set_demo_five():
	"""
	Sets up demo 3 in the simulation. Used for demonstration.
	Parameters:
		None.

	Returns:
		Nothing.
	"""
	global balls,ramps,gravity,air_density,sim_on
	balls.clear()
	ramps.clear()
	ball = Ball(375, 275, 20, 20.0, 0, 0, 0, 0, 0, 0, Red)
	balls.append(ball)
	ball = Ball(350, 245, 20, 20.0, 0, 0, 0, 0, 0, 0, Yellow)
	balls.append(ball)
	ball = Ball(395, 245, 20, 20.0, 0, 0, 0, 0, 0, 0, Green)
	balls.append(ball)
	ball = Ball(375, 210, 20, 20.0, 0, 0, 0, 0, 0, 0, Blue)
	balls.append(ball)
	ball = Ball(330, 210, 20, 20.0, 0, 0, 0, 0, 0, 0, Black)
	balls.append(ball)
	ball = Ball(420, 210, 20, 20.0, 0, 0, 0, 0, 0, 0, Cyan)
	balls.append(ball)
	ball = Ball(375, 500, 20, 20.0, 90, 150, 0, 0, 0, 0, Black)
	balls.append(ball)
	gravity_entry.set_value('0.0')
	drag_entry.set_value('15.0')
	elastic_radio.set_value(True)
	sim_on = False

def in_bounds(bounds_array):
	"""
	Returns True if the mouse is in the bounds of a given array that represents four sides of a rectangle.
	Parameters:
		bounds_array (list) -- required
			- must be in the following format: [left bound, right bound, top bound, bottom bound]

	Returns:
		boolean
	"""
	global click_points
	if click_points[0] > bounds_array[0] and click_points[0] < bounds_array[1] and click_points[1] > bounds_array[2] and click_points[1] < bounds_array[3]:
		return True
	else:
		return False

#-------color RGB constants-------
White = (230,230,230)
Black = (0,0,0)
Red = (255,0,0)
Green = (0,200,0)
Blue = (0,0,255)
Yellow = (255,255,0)
Cyan = (0,200,200)
Grey = (130,130,130)

back_rectangle = pygame.Rect(970,170,42,150)
black_square = pygame.Rect(980,180,20,20)
red_square = pygame.Rect(980,202,20,20)
green_square = pygame.Rect(980,224,20,20)
blue_square = pygame.Rect(980,246,20,20)
yellow_square = pygame.Rect(980,268,20,20)
cyan_square = pygame.Rect(980,290,20,20)

#bounds_format = [left,right,top,bottom]
black_bounds = [980, 1000, 180, 200]
red_bounds = [980, 1000, 202, 222]
green_bounds = [980, 1000, 224, 244]
blue_bounds = [980, 1000, 246, 266]
yellow_bounds = [980, 1000, 268, 288]
cyan_bounds = [980, 1000, 290, 310]

pygame.init()
pygame.font.init()

#-------initialize screen-------
pygame.key.set_repeat(300,30)
screen_size = (1100,600)
screen = pygame.display.set_mode(screen_size, pygame.DOUBLEBUF)
screen.fill(Grey)		
UPDATE = pygame.USEREVENT+1				#initialize UPDATE as a self-defined event
pygame.time.set_timer(UPDATE, 1)		#update every 1 millisecond

font1 = pygame.font.SysFont("meera", 20)
font2 = pygame.font.SysFont("dejavuserif", 10)
sim_area_title = font1.render("Simulation Settings", True, Black)
sim_window_title = font1.render("Simulation Window", True, Black)
obj_settings_title = font1.render("Object Settings", True, Black)
F_letter = font2.render("F", True, Black)

#-------button/input init-------
def toggle_sim():
	"""
	Function called when the user clicks the 'Start/Stop Simulation' buttons. Switches the sim_on
 	boolean variable to True if False and to False if True.
	Parameters:
		None.

	Returns:
		Nothing.
	"""
	global sim_on
	if sim_on == False:
		sim_on = True
	else:
		sim_on = False

def clear_sim():
	"""
	Function called when the user clicks the 'Clear' button. Clears the balls array.
	Parameters:
		None.

	Returns:
		Nothing.
	"""
	global balls,sim_on,ramps,add_buffer
	balls.clear()
	ramps.clear()
	add_buffer.clear()
	sim_on = False

def reset_sim():
	"""
	Function called when the user clicks the 'Reset' button. Resets all objects to their
	original positions.
	Parameters:
		None.

	Returns:
		Nothing.
	"""
	global balls
	global sim_on
	for ball in balls:
		ball.reset()
	sim_on = False

def placing_obj():
	"""
	Function called when the user clicks the 'Place Object' button. Switches the
	user_placing_object button to True if False and to False if True.
	Parameters:
		None.

	Returns:
		Nothing.
	"""
	global user_placing_object
	if user_placing_object == False:
		user_placing_object = True	
	else:
		user_placing_object = False	
	
def undo_placement():
	"""
	Function called when the user clicks the 'Undo' button. Removes the last element pushed into 
	the balls array.
	Parameters:
		None.

	Returns:
		Nothing.
	"""
	global balls, ramps, add_buffer
	if len(add_buffer) > 0:
		temp = add_buffer.pop()
		if temp == 'ball':
			balls.pop()
		elif temp == 'ramp':
			ramps.pop()

def placing_ball():
	"""
	Changes object_held to 0, which means the user is 'holding a ball'.
	Parameters:
		None.

	Returns:
		Nothing.
	"""
	global object_held
	object_held = 0

def placing_ramp():
	"""
	Changes object_held to 1, which means the user is 'holding a ramp'.
	Parameters:
		None.

	Returns:
		Nothing.
	"""
	global object_held
	object_held = 1

#initalize the start simulation button
start_button = thorpy.make_button("Start Simulation", func=toggle_sim)
start_box = thorpy.Box.make(elements=[start_button])
start_menu = thorpy.Menu(start_box)
for element in start_menu.get_population():
	element.surface = screen
start_box.set_topleft((648,550))
start_box.update()

#initialize the stop simulation button
stop_button = thorpy.make_button("Stop Simulation", func=toggle_sim)
stop_box = thorpy.Box.make(elements=[stop_button])
stop_menu = thorpy.Menu(stop_box)
for element in stop_menu.get_population():
	element.surface = screen
stop_box.set_topleft((648,550))
stop_box.update()

#intialize the clear button
clear_button = thorpy.make_button("Clear", func=clear_sim)
clear_box = thorpy.Box.make(elements=[clear_button])
clear_menu = thorpy.Menu(clear_box)
for element in clear_menu.get_population():
	element.surface = screen
clear_box.set_topleft((48,550))
clear_box.update()

#intialize the reset button
reset_button = thorpy.make_button("Reset", func=reset_sim)
reset_box = thorpy.Box.make(elements=[reset_button])
reset_menu = thorpy.Menu(reset_box)
for element in reset_menu.get_population():
	element.surface = screen
reset_box.set_topleft((98,550))
reset_box.update()

#initialize place object, redo, in object control box
place_obj_button = thorpy.make_button("Place Object", func=placing_obj)
place_obj_box = thorpy.Box.make(elements=[place_obj_button])
place_obj_menu = thorpy.Menu(place_obj_box)
for element in place_obj_menu.get_population():
	element.surface = screen
place_obj_box.set_topleft((600,14))
place_obj_box.update()

#object to indicate that you are placing a ball
ball_button = thorpy.make_button("Ball", func=placing_ball)
ball_box = thorpy.Box.make(elements=[ball_button])
ball_menu = thorpy.Menu(ball_box)
for element in ball_menu.get_population():
	element.surface = screen
ball_box.set_topleft((300,14))
ball_box.update()

#object to indicate that you are placing a ball
ramp_button = thorpy.make_button("Ramp", func=placing_ramp)
ramp_box = thorpy.Box.make(elements=[ramp_button])
ramp_menu = thorpy.Menu(ramp_box)
for element in ball_menu.get_population():
	element.surface = screen
ramp_box.set_topleft((400,14))
ball_box.update()

#cancel placement button
cancel_obj_button = thorpy.make_button("Cancel", func=placing_obj)
cancel_obj_box = thorpy.Box.make(elements=[cancel_obj_button])
cancel_obj_menu = thorpy.Menu(cancel_obj_box)
for element in cancel_obj_menu.get_population():
	element.surface = screen
cancel_obj_box.set_topleft((600,14))
cancel_obj_box.update()

#initialize undo button
undo_button = thorpy.make_button("Undo", func=undo_placement)
undo_box = thorpy.Box.make(elements=[undo_button])
undo_menu = thorpy.Menu(undo_box)
for element in undo_menu.get_population():
	element.surface = screen
undo_box.set_topleft((700,14))
undo_box.update()

#BEGIN MAIN SETTINGS AREA
 
#simulation settings area input fields
gravity_entry = thorpy.Inserter.make(name="      Gravity", value="-9.8")
drag_entry = thorpy.Inserter.make(name=   "Air Density", value="1.5")
elastic_radio = thorpy.Checker.make("Elastic Collsions", type_="radio")
sim_settings_box = thorpy.Box.make(elements=[gravity_entry, drag_entry,elastic_radio])
sim_settings_menu = thorpy.Menu(sim_settings_box)
for element in sim_settings_menu.get_population():
	element.surface = screen
sim_settings_box.fit_children(margins=(10, 10))
sim_settings_box.set_topleft((800, 50))
sim_settings_box.update()

#Begin intialization of object settings area
#--------
mass_entry = thorpy.Inserter.make(name="         Mass", value="10.0")
radius_entry = thorpy.Inserter.make(name="      Radius", value="20.0")
obj_settings_text = thorpy.make_text("Attributes             ", 18, Black)
obj_settings_box = thorpy.Box.make(elements=[obj_settings_text ,mass_entry, radius_entry])
obj_settings_menu = thorpy.Menu(obj_settings_box)
for element in obj_settings_menu.get_population():
	element.surface = screen
obj_settings_box.fit_children(margins=(10, 10))
obj_settings_box.set_topleft((800, 170))
obj_settings_box.update()

#force magnitude and angle entry fields
force_mag_entry = thorpy.Inserter.make(name="Magnitude", value="0.0")
force_angle_entry = thorpy.Inserter.make(name="         Angle", value="0")
force_title = thorpy.make_text("Force                  ", 18, Black)
force_box = thorpy.Box.make(elements=[force_title, force_mag_entry, force_angle_entry])
force_menu = thorpy.Menu(force_box)
for element in force_menu.get_population():
	element.surface = screen
force_box.fit_children(margins=(10, 10))
force_box.set_topleft((800, 270))
force_box.update()

#velocity magnitude and angle entry fields
velocity_mag_entry = thorpy.Inserter.make(name="Magnitude", value="0.0")
velocity_angle_entry = thorpy.Inserter.make(name="         Angle", value="0")
velocity_title = thorpy.make_text("Velocity              ", 18, Black)
velocity_box = thorpy.Box.make(elements=[velocity_title, velocity_mag_entry, velocity_angle_entry])
velocity_menu = thorpy.Menu(velocity_box)
for element in velocity_menu.get_population():
	element.surface = screen
velocity_box.fit_children(margins=(10, 10))
velocity_box.set_topleft((800, 370))
velocity_box.update()

#friction coefficient entry fields
mu_s_entry = thorpy.Inserter.make(name="         Static", value="0.0")
mu_k_entry = thorpy.Inserter.make(name="        Kinetic", value="0.0")
friction_title = thorpy.make_text("Friction              ", 18, Black)
friction_box = thorpy.Box.make(elements=[friction_title, mu_s_entry, mu_k_entry])
friction_menu = thorpy.Menu(friction_box)
for element in friction_menu.get_population():
	element.surface = screen
friction_box.fit_children(margins=(10, 10))
friction_box.set_topleft((800, 470))
friction_box.update()

#ramp entry fields
ramp_height_entry = thorpy.Inserter.make(name="Height", value="100.0")
ramp_width_entry = thorpy.Inserter.make(name="Width", value="100.0")
ramp_title = thorpy.make_text("Ramp Settings     ", 18, Black)
ramp_settings_box = thorpy.Box.make(elements=[ramp_title, ramp_height_entry, ramp_width_entry])
ramp_settings_menu = thorpy.Menu(ramp_settings_box)
for element in ramp_settings_menu.get_population():
	element.surface = screen
ramp_settings_box.fit_children(margins=(10, 10))
ramp_settings_box.set_topleft((800, 170))
ramp_settings_box.update()

def draw_sim_buttons():
	"""
	Renders the simulation and object control buttons. This is called every time the screen
	updates.
	Parameters:
		None.

	Returns:
		Nothing.
	"""
	global sim_on
	clear_box.blit()
	reset_box.blit()
	if user_placing_object == False:
		place_obj_box.blit()
	else:
		cancel_obj_box.blit()
		ball_box.blit()
		ramp_box.blit()
	undo_box.blit()
	if sim_on == False:
		start_box.blit()
	else:
		stop_box.blit()

def draw_text():
	"""
	Renders the text in the application. This is called every time the screen updates.
	Parameters:
		None.

	Returns:
		Nothing.
	"""
	screen.blit(sim_area_title, (800,25))
	screen.blit(sim_window_title, (50,25))
	if user_placing_object == True:
		screen.blit(obj_settings_title, (800, 150))

def draw_settings_area():
	"""
	Renders the object and simulation settings boxes. This is called every time the screen
	updates.
	Parameters:
		None.

	Returns:
		Nothing.
	"""
	global click_points,black_bounds,red_bounds,green_bounds,blue_bounds,yellow_bounds,cyan_bounds,selected_color
	sim_settings_box.blit()
	if user_placing_object == True:
		select_rectangle = pygame.Rect(980,180,20,20)
		selected_color = (0,0,0)
		#adding color choices
		pygame.draw.rect(screen, (180,180,180), back_rectangle)
		pygame.draw.rect(screen, Black, black_square)
		pygame.draw.rect(screen, Red, red_square)
		pygame.draw.rect(screen, Green, green_square)
		pygame.draw.rect(screen, Blue, blue_square)
		pygame.draw.rect(screen, Yellow, yellow_square)
		pygame.draw.rect(screen, Cyan, cyan_square)
		if in_bounds(black_bounds) == True:
			pygame.draw.rect(screen, (255,255,255), black_square, 3)
			selected_color = (0,0,0)
		elif in_bounds(red_bounds) == True:
			pygame.draw.rect(screen, (255,255,255), red_square, 3)
			selected_color = (255,0,0)
		elif in_bounds(green_bounds) == True:
			pygame.draw.rect(screen, (255,255,255), green_square, 3)
			selected_color = (0,255,0)
		elif in_bounds(blue_bounds) == True:
			pygame.draw.rect(screen, (255,255,255), blue_square, 3)
			selected_color = (0,0,255)
		elif in_bounds(yellow_bounds) == True:
			pygame.draw.rect(screen, (255,255,255), yellow_square, 3)
			selected_color = (255,255,0)
		elif in_bounds(cyan_bounds) == True:
			pygame.draw.rect(screen, (255,255,255), cyan_square, 3)
			selected_color = (0,255,255)
		else:
			pygame.draw.rect(screen, (255,255,255), select_rectangle, 3)
			
		if object_held == 0:
			obj_settings_box.blit()
			force_box.blit()
			velocity_box.blit()
			friction_box.blit()
		if object_held == 1:
			ramp_settings_box.blit()
	else:
		obj_settings_box.update()
		force_box.update()

def set_value(entry):
	"""
	Aids in setting user parameters by checking for errors before setting the values in the simulation.
	Parameters:
		entry (String).

	Returns:
		float.
	"""
	if entry != "":
		return float(entry)
	else:
		return 0.0

#-------begin main display loop-------
	
while True:
	clock.tick(60)	#60 fps limiter, prevents flicker
	screen.fill(Grey)		
	pygame.draw.rect(screen, White, (50, 50, 700, 500), 0)	#simulation window
	draw_sim_buttons()
	draw_text()
	draw_settings_area()

	#begin event loop
	for event in pygame.event.get():
		if event.type == pygame.QUIT:
			sys.exit()
	
		#when the mouse is clicked
		if event.type == pygame.MOUSEBUTTONDOWN:
			if event.button == 1:		#left click
				temp = pygame.mouse.get_pos()		#mouse's current position
				click_points = [temp[0], temp[1]]
				#only draw the ball if it has enough space to do so and it is in the sim window
				if is_mouse_in_window(temp[0], temp[1]) == True and does_collide(temp[0],temp[1]) == False and user_placing_object == True: 
					if object_held == 0:
						#get all the values from the text entry fields, create the ball object
						#and add it to balls[]
						color = selected_color
						mass = set_value(mass_entry.get_value())
						radius = set_value(radius_entry.get_value())
						Fang = set_value(force_angle_entry.get_value())
						Fmag = set_value(force_mag_entry.get_value())
						vang = set_value(velocity_angle_entry.get_value())
						vmag = set_value(velocity_mag_entry.get_value())
						new_ball = Ball(temp[0], temp[1], mass, radius, vang, vmag, 0, 0, Fang, Fmag, color)
						#---------IMPORTANT----------						
						new_ball.set_friction(mu_s_entry.get_value(), mu_k_entry.get_value())
						#----------------------------
						balls.append(new_ball)
						add_buffer.append('ball')
						user_placing_object = False
					if object_held == 1:
						color = selected_color
						height = set_value(ramp_height_entry.get_value())
						width = set_value(ramp_width_entry.get_value())
						new_ramp = Ramp(temp[0], temp[1], height, width, color)
						ramps.append(new_ramp)
						add_buffer.append('ramp')
			elif event.button == 3:		#right click
				user_placing_object = False

		#when any key is pressed
		if event.type == pygame.KEYDOWN:
			if event.key == pygame.K_q or event.key == pygame.K_ESCAPE:		#Q or ESC key
				sys.exit()
			if event.key == pygame.K_d and demos_on == True:		#Q or ESC key
				demos_on = False
				sim_on = False
			if event.key == pygame.K_d and demos_on == False:		#Q or ESC key
				demos_on = True
				balls.clear()
				ramps.clear()
				sim_on = False
			if demos_on == True:
				if event.key == pygame.K_1:		#1 key
					if sim_on == False:
						set_demo_one()
				if event.key == pygame.K_2:		#2 key
					if sim_on == False:
						set_demo_two()
				if event.key == pygame.K_3:		#3 key
					if sim_on == False:
						set_demo_three()
				if event.key == pygame.K_4:		#4 key
					if sim_on == False:
						set_demo_four()
				if event.key == pygame.K_5:		#5 key
					if sim_on == False:
						set_demo_five()

		if event.type == UPDATE:		#update the screen's visuals
			pos = pygame.mouse.get_pos()
			if sim_on == True:
				is_elastic = elastic_radio.get_value()
				air_density = set_value(drag_entry.get_value())
				if is_elastic == False:
					elastic_constant = 0.2
				else:
					elastic_constant = 0.0
			#if the mouse is in the sim window and the user is placing the object, draw the outline
			#of the object
			if is_mouse_in_window(pos[0],pos[1]) == True and user_placing_object == True:
				if object_held == 0:
					pygame.draw.circle(screen, Black, pos, int(set_value(radius_entry.get_value())), 1)
					pygame.draw.line(screen, Black, pos, ((80+set_value(radius_entry.get_value()))*math.cos(set_value(velocity_angle_entry.get_value())*(3.14/180))+pos[0], -(80+set_value(radius_entry.get_value()))*math.sin(set_value(velocity_angle_entry.get_value())*(3.14/180))+pos[1]), 1)
				elif object_held == 1:
					pygame.draw.polygon(screen, Black, [pos, (pos[0]+int(set_value(ramp_width_entry.get_value())), pos[1]), (pos[0], pos[1]-int(set_value(ramp_height_entry.get_value())))], 1)
			handle_collisions()

			for ball in balls:
				if sim_on == True:
					gravity = set_value(gravity_entry.get_value())
				ball.set_gravity(gravity)
				ptx = int(ball._x)		#ball's x position
				pty = int(ball._y)		#ball's y position
				rad = int(ball._radius)		#radius of the ball being drawn
				f_ang = int(ball._Fang)	#angle of the force vector
				F_letter = font2.render(("F = "+ str(int(ball._Fmag)) +" N"), True, Black)
				v_ang = int(ball._vang)	#angle of the force vector
				v_letter = font2.render(("v = "+ str(int(ball._vmag*10)) +" m/s"), True, Black)
				#draw ball on screen
				if sim_on == False:
					if ball._Fmag != 0:
						#draw force vector line
						pygame.draw.line(screen, Black, (ptx, pty), ((80+rad)*math.cos(f_ang*(3.14/180))+ptx, -(80+rad)*math.sin(f_ang*(3.14/180))+pty), 1)	
						#draw force vector label
						screen.blit(F_letter, ((80+rad)*math.cos(f_ang*radian_const)+ptx - 5, -(80+rad)*math.sin(f_ang*(3.14/180))+pty + 5))

					if ball._vmag != 0:
						#draw velocity vector line
						pygame.draw.line(screen, Black, (ptx, pty), ((80+rad)*math.cos(v_ang*(3.14/180))+ptx, -(80+rad)*math.sin(v_ang*(3.14/180))+pty), 1)	
						#draw velocity vector label
						screen.blit(v_letter, ((80+rad)*math.cos(v_ang*(3.14/180))+ptx - 5, -(80+rad)*math.sin(v_ang*(3.14/180))+pty + 5))
				pygame.draw.circle(screen, ball._color, (ptx, pty), rad)
				if sim_on == True:
					ball.update_values(air_density)
			for ramp in ramps:
				pygame.draw.polygon(screen, ramp._color, ramp._pointlist)
		
		#GUI button/text entry event reactions
		if sim_on == False:			
			start_menu.react(event)
		else:
			stop_menu.react(event)		
		clear_menu.react(event)
		reset_menu.react(event)
		sim_settings_menu.react(event)
		if user_placing_object == True:
			cancel_obj_menu.react(event)
			ramp_menu.react(event)
			ball_menu.react(event)
			if object_held == 0:
				obj_settings_menu.react(event)
				force_menu.react(event)
				velocity_menu.react(event)
				friction_menu.react(event)
			elif object_held == 1:
				ramp_settings_menu.react(event)
		else:
			place_obj_menu.react(event)
		undo_menu.react(event)
	
		pygame.display.update()	#update the screen
#-------end main display loop-------
pygame.quit()
#-------end program-------
