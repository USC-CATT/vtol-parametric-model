#!/usr/bin/env python
# coding: utf-8

# # SymPy Background
# 
# We have been stuck on the linearization of our equations of motion for a long time. What we would like to do is to use the _SymPy_ module from python to do some of the heavy lifting for us, and then we can verify its correctness manually. To do this, first we want to perform the parametric linearization on a simpler problem, and confirm that we know how to go from the nonlinear equations to the linearized ones at various trim points.

# ## Simpler problem
# The simpler problem we want to solve is
# $$
# \ddot{x} + 2\dot{x}^3 + 2x = 0
# $$
# 
# Which we can rearrange into the state space equation
# \begin{align}
# \dot{x}_1 &= (x_1 + x_2)^2 \\
# \dot{x}_2 &= (x_1 + 1)(x_2 + 2) + u
# \end{align}

# This system has equilibrium points at 
# $$
# x_2 = 0
# $$
# $$
# u = 2x_1
# $$

# In[3]:


from sympy import *
x1, x2, t, u = symbols("x1 x2 t u")


# In[4]:


f1 =Function("f1")
f2 = Function("f2")


# f1 = x2
# f2 = -2*x1 - 2*x2**3 + u

# f2 - f1

# In[5]:


f2


# In[6]:


df1dx1 = f1.diff(x1)
df1dx1


# In[7]:


df1dx2 = f1.diff(x2)
df1dx2


# In[8]:


df1du = f1.diff(u)
df1du


# In[9]:


df2dx1 = f2.diff(x1)
df2dx1


# In[10]:


df2dx2 = f2.diff(x2)
df2dx2


# In[11]:


df2du = f2.diff(u)
df2du


# In[12]:


A = Matrix([[f.diff(x) for x in [x1, x2]] for f in [f1, f2]])
A


# In[13]:


A.subs(x2, 0)


# In[14]:


B = Matrix([[f.diff(u) for u in [u]] for f in [f1, f2]])
B


# # Applying Our system
# We now want to turn to the system of equations we have interest in, which are much more complex than the ones we were just considering. First we need to define our variables.

# ### Relevant Papers
# Alaez et al: good diagram of coordinate systems and roudh draft of equations (ours are slightly modified from theirs)
# 
# Cook et al 2021: good diagram of our system, this is also the official way the matlab simulator is calculating things so we should /eventually/ check that its equivalent
# 
# time derivative of rotation matrices (not super relevant, useful to know how certain coordinate transoform derivatives work)
# 
# acheson 2025: official definition of GUAM system, has propeller numbers we should probably use
# 
# 

# In[15]:


# Shortcut functions for sin and cosine
x, y, z = symbols('x y z')
c = Lambda(x, cos(x))
s = Lambda(x, sin(x))


# ## Variable Definitions

# Position Variable $\vec{p}_{B/N}^N$

# In[16]:


pn, pe, pd = symbols("p_n p_e p_d")
pos = Matrix([[pn], [pe], [pd]])
pos


# Orientation $\vec{\Phi}_{B/N}^N$

# In[17]:


phi, tht, psi = symbols("phi theta psi")
PHI = Matrix([[phi], [tht], [psi]])
PHI


# Velocity $\vec{v}_{B/E}^B$

# In[18]:


u, v, w = symbols("u v w")
vel = Matrix([[u], [v], [w]])
vel


# Angular Velocity $\vec{\omega}_{B/E}^B$

# In[19]:


p, q, r = symbols("P Q R")
omg = Matrix([[p], [q], [r]])
omg


# ## Dynamic Derivative Definitions

# $\dot{\vec{p}}_{B/N}^N$

# In[20]:


C_bfn = Matrix([[c(tht)*c(psi), c(tht)*s(psi), -s(tht)],
                [-c(phi)*s(psi)+s(phi)*s(tht)*c(psi), c(phi)*c(psi)+s(phi)*s(tht)*s(psi), s(phi)*c(tht)],
                [s(phi)*s(psi)+c(phi)*s(tht)*c(psi), -s(phi)*c(psi)+c(phi)*s(tht)*s(psi), c(phi)*c(tht)]])
C_bfn


# In[21]:


dpos_brn_n = C_bfn.T*vel
dpos_brn_n


# $\dot{\vec{\Phi}}_{B/N}^N$

# In[22]:


H_phi = Matrix([[1, s(phi)*tan(tht), c(phi)*tan(tht)],
                [0, c(phi), -s(phi)],
                [0, s(phi)/c(tht), c(phi)/c(tht)]])
H_phi


# In[23]:


dPHI_brn_n = H_phi*omg
dPHI_brn_n


# $\dot{\vec{v}}_{B/E}^B$

# In[24]:


m, g = symbols("m g")
gvec = Matrix([[0],[0],[g]])
C_bfn*gvec


# In[25]:


omg.cross(vel)


# In[26]:


wn, we, wd = symbols("w_n w_e w_d")
v_wre_e = Matrix([[wn],[we],[wd]])
v_rel = vel - C_bfn*v_wre_e
v_rel


# In[27]:


V_T = v_rel.norm()


# In[28]:


alpha = atan2(v_rel[2], v_rel[0])


# In[29]:


beta = asin(v_rel[1]/V_T)


# In[30]:


alpha


# In[31]:


C_wfb = Matrix([[c(alpha)*c(beta), s(beta), s(alpha)*c(beta)],
                [-c(alpha)*s(beta), c(beta), -s(alpha)*s(beta)],
                [-s(alpha), 0, c(alpha)]])


# In[32]:


rho = symbols("\\rho")


# Define Wing features

# In[33]:


ail1, ail2, elev, rudd = symbols("\\delta{}a_1 \\delta{}a_2 \\delta{}e \\delta{}r")
i1, i2, i3, i4 = symbols("i_1 i_2 i_3 i_4")
xs1, ys1, zs1 = symbols("x_{s1} y_{s1} z_{s1}")
xvecs1 = Matrix([[xs1],[ys1],[zs1]])
xp1, yp1, zp1 = symbols("x_{p1} y_{p1} z_{p1}")
xvecp1 = Matrix([[xp1],[yp1],[zp1]])


# Front left wing

# In[34]:


b1, c1 = symbols("b_1 c_1")
S1 = b1*c1
AR1 = b1/c1
Cl01, dCldu1, cu1 = symbols("C_{L01} \\frac{\\partial{}C_{L1}}{\\partial\\delta{}u} c_{u1}")
Cla1 = 2*pi*AR1/(2+AR1)
CL1 = Lambda((x, y), Cl01 + Cla1*x + dCldu1*y*cu1/c1)
alpha_s1 = alpha + i1


# In[35]:


CL1(x, y)


# In[36]:


L1 = 1/2*rho*v_rel.dot(v_rel)*S1*CL1(alpha_s1, ail1)


# In[37]:


Cd01, Cda1, a01, e1, dCddu1 = symbols("C_{D01} C_{d\\alpha{}1} \\alpha_{01} e_1 \\frac{\\partial{}C_{D1}}{\\partial\\delta{}u_1}")
CD1 = Lambda((x, y), Cd01 + Cda1*(x-a01)**2 + CL1(x, y)**2/pi/e1/AR1 + dCddu1*y*cu1/c1)


# In[38]:


CD1(x, y)


# In[39]:


D1 = 1/2*rho*v_rel.dot(v_rel)*S1*CD1(alpha_s1, ail1)


# In[40]:


F_a1_w = Matrix([[-D1],[0],[-L1]])


# In[41]:


F_a1_b = C_wfb.T*F_a1_w
F_a1_b


# In[42]:


Thp1, tau1 = symbols("T_{p1} \\tau_1")
n1n, n1e, n1d = symbols("n_{1N} n_{1E} n_{1D}")
n_p1 = Matrix([[n1n],[n1e],[n1d]])
s_p1 = symbols("s_1")


# In[43]:


F_p1_b = Thp1*n_p1


# In[44]:


M_p1_b = tau1*s_p1*n_p1


# Aside:
# The motors can and should be described separately, using the following equations of motion
# \begin{align}
# \tau &= \frac{K_T}{R_m}(T - K_v\omega_m \\
# \dot{\omega}_m &= \frac{-K_TK_V}{R_mJ_m}\omega_m - \frac{C_D}{J_m}\omega_m^2 + \frac{K_T}{R_mJ_m}T \\
# T_p &= C_{Lm}\omega^2
# \end{align}
# With $\omega_m$ as the rotational speed of the motor, $\tau$ the torque generated by the motor, $T$ the requested thrust input to the motor system, and $T$ the actual thrust output from the motor system. For now, we will just model the relevant "outputs" of $\tau$ and $T_p$
# 
# s represents the sign of the direction the motor is rotating (+1 for CCW, -1 for CW)

# In[45]:


Cm01, Cma1, dCmdu1 = symbols("C_{M01} C_{M\\alpha{}1} \\frac{\\partial{}C_M1}{\\partial\\delta{}u_1}") 


# In[46]:


CM1 = Lambda((x, y), Cm01 + Cma1*x + dCmdu1*y*cu1/c1)


# In[47]:


M1 = 1/2*rho*S1*c1*v_rel.dot(v_rel)*CM1(alpha_s1,ail1)


# In[48]:


M_a1_b = Matrix([[0],[M1],[0]])


# In[49]:


M_f1_b = xvecs1.cross(F_a1_b)


# In[50]:


M_T1_b = xvecp1.cross(F_p1_b)


# In[51]:


F_ext_b = (F_a1_b) + (F_p1_b)


# In[52]:


M_ext_b = (M_a1_b + M_f1_b) + (M_T1_b + M_p1_b)


# $\dot{\vec{v}}_{B/E}^B$

# In[53]:


dvel_bre_b = 1/m*F_ext_b + C_bfn*gvec - omg.cross(vel)


# $\dot{\vec{\omega}}_{B/E}^B$

# In[54]:


Ixx, Iyy, Izz, Ixy, Iyz, Ixz = symbols("I_{xx} I_{yy} I_{zz} I_{xy} I_{yz} I_{xz}")
J = Matrix([[Ixx, Ixy, Ixz],
            [Ixy, Iyy, Iyz],
            [Ixz, Iyz, Izz]])


# In[55]:


domg_bre_b = J.inv()*(M_ext_b - omg.cross(J*omg))


# ## Linearization
# 
# We now have 12 state variables $(p_N, p_E, p_D, \phi, \theta, \psi, u, v, w, p, q, r)$ and a certain number of input variables, which for now since we only defined one wing and one propeller is $(a, e, r, T_p, \tau, w_N, w_E, w_D)$. We want to take each of our equations and linearize them against these variables, to get the linearized equations.

# In[56]:


dpos_brn_n.diff(pos)


# In[57]:


dpos_brn_n.subs({tht: 0, psi: 0}).diff(phi)


# In[58]:


dpos_brn_n.subs({phi: 0, psi: 0}).diff(tht)


# In[59]:


dpos_brn_n.subs({phi: pi/6, tht: pi/4, psi: 0}).diff(vel)


# In[60]:


dpos_brn_n.diff(PHI)


# In[61]:


dpos_brn_n.diff(vel)


# In[62]:


dpos_brn_n.diff(omg)


# In[63]:


surfaces = Matrix([[ail1], [ail2], [elev], [rudd]])


# In[64]:


props = Matrix([[Thp1], [tau1]])


# In[65]:


wind = Matrix([[wn],[we],[wd]])


# In[66]:


dpos_brn_n.diff(surfaces)


# In[67]:


dpos_brn_n.diff(props)


# In[68]:


dpos_brn_n.diff(wind)


# In[69]:


dPHI_brn_n.diff(pos)


# In[70]:


dPHI_brn_n.diff(PHI)


# In[71]:


dPHI_brn_n.diff(vel)


# In[72]:


dPHI_brn_n.diff(omg)


# In[73]:


dPHI_brn_n.diff(surfaces)


# In[74]:


dPHI_brn_n.diff(props)


# In[75]:


dPHI_brn_n.diff(wind)


# In[76]:


dvel_bre_b.diff(pos)


# In[77]:


dvel_bre_b.diff(PHI)


# In[78]:


dvel_bre_b.diff(vel)


# In[79]:


dvel_bre_b.diff(omg)


# In[80]:


dvel_bre_b.diff(surfaces)


# In[81]:


dvel_bre_b.diff(props)


# In[82]:


dvel_bre_b.diff(wind)


# In[83]:


domg_bre_b.diff(pos)


# In[84]:


domg_bre_b.diff(PHI)


# In[85]:


domg_bre_b.diff(vel)


# In[86]:


domg_bre_b.diff(omg)


# In[87]:


domg_bre_b.diff(surfaces)


# In[88]:


domg_bre_b.diff(props)


# In[89]:


domg_bre_b.diff(wind)

