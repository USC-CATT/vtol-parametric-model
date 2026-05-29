from sympy import *

# Shortcut functions for sin and cosine
x, y, z = symbols('x y z')
c = Lambda(x, cos(x))
s = Lambda(x, sin(x))

#Variable Definitions

# Position Variable
pn, pe, pd = symbols("p_n p_e p_d")
pos = Matrix([[pn], [pe], [pd]])
print("pos =")
print(pos)

# Orientation
phi, tht, psi = symbols("phi theta psi")
PHI = Matrix([[phi], [tht], [psi]])
print("PHI =")
print(PHI)

# Velocity
u, v, w = symbols("u v w")
vel = Matrix([[u], [v], [w]])
print("vel =")
print(vel)

# Angular Velocity
p, q, r = symbols("p q r")
omg = Matrix([[p], [q], [r]])
print("omg =")
print(omg)


#Kinematics

C_bfn = Matrix([
    [c(tht)*c(psi), c(tht)*s(psi), -s(tht)],
    [-c(phi)*s(psi)+s(phi)*s(tht)*c(psi), c(phi)*c(psi)+s(phi)*s(tht)*s(psi), s(phi)*c(tht)],
    [s(phi)*s(psi)+c(phi)*s(tht)*c(psi), -s(phi)*c(psi)+c(phi)*s(tht)*s(psi), c(phi)*c(tht)]
])

dpos_brn_n = C_bfn.T * vel

H_phi = Matrix([
    [1, s(phi)*tan(tht), c(phi)*tan(tht)],
    [0, c(phi), -s(phi)],
    [0, s(phi)/c(tht), c(phi)/c(tht)]
])

dPHI_brn_n = H_phi * omg


#Velocity Derivative Setup

m, g = symbols("m g")
gvec = Matrix([[0], [0], [g]])

wn, we, wd = symbols("w_n w_e w_d")
v_wre_e = Matrix([[wn], [we], [wd]])

v_rel = vel - C_bfn * v_wre_e
V_T = v_rel.norm()

alpha = atan2(v_rel[2], v_rel[0])
beta = asin(v_rel[1]/V_T)

C_wfb = Matrix([
    [c(alpha)*c(beta), s(beta), s(alpha)*c(beta)],
    [-c(alpha)*s(beta), c(beta), -s(alpha)*s(beta)],
    [-s(alpha), 0, c(alpha)]
])

rho = symbols("\\rho")


# wing 1 aerodynamics

ail1, ail2, elev, rudd = symbols("\\delta{}a_1 \\delta{}a_2 \\delta{}e \\delta{}r")

i1 = symbols("i_1")
xs1, ys1, zs1 = symbols("x_{s1} y_{s1} z_{s1}")
xvecs1 = Matrix([[xs1], [ys1], [zs1]])

b1, c1 = symbols("b_1 c_1")
S1 = b1*c1
AR1 = b1/c1

Cl01, dCldu1, cu1 = symbols("C_{L01} \\frac{\\partial{}C_{L1}}{\\partial\\delta{}u} c_{u1}")
Cla1 = 2*pi*AR1/(2 + AR1)

CL1 = Lambda((x, y), Cl01 + Cla1*x + dCldu1*y*cu1/c1)
alpha_s1 = alpha + i1

L1 = Rational(1, 2)*rho*v_rel.dot(v_rel)*S1*CL1(alpha_s1, ail1)

Cd01, Cda1, a01, e1, dCddu1 = symbols(
    "C_{D01} C_{d\\alpha{}1} \\alpha_{01} e_1 \\frac{\\partial{}C_{D1}}{\\partial\\delta{}u_1}"
)

CD1 = Lambda(
    (x, y),
    Cd01 + Cda1*(x - a01)**2 + CL1(x, y)**2/(pi*e1*AR1) + dCddu1*y*cu1/c1
)

D1 = Rational(1, 2)*rho*v_rel.dot(v_rel)*S1*CD1(alpha_s1, ail1)

F_a1_w = Matrix([[-D1], [0], [-L1]])
F_a1_b = C_wfb.T * F_a1_w

Cm01, Cma1, dCmdu1 = symbols(
    "C_{M01} C_{M\\alpha{}1} \\frac{\\partial{}C_M1}{\\partial\\delta{}u_1}"
)

CM1 = Lambda((x, y), Cm01 + Cma1*x + dCmdu1*y*cu1/c1)

M1 = Rational(1, 2)*rho*S1*c1*v_rel.dot(v_rel)*CM1(alpha_s1, ail1)

M_a1_b = Matrix([[0], [M1], [0]])
M_f1_b = xvecs1.cross(F_a1_b)


# Force Function

def prop_force_moment(Th, tau, spin, xvec, nvec):
    """
    Th   = thrust magnitude
    tau  = reaction torque magnitude
    spin = +1 or -1 spin direction
    xvec = prop position vector from CG, body frame
    nvec = thrust direction unit vector, body frame
    """
    F = Th * nvec
    M = xvec.cross(F) + tau * spin * nvec
    return F, M


# Define each propeller force

Thp = symbols('T_p1:9')          # T_p1 ... T_p8
taup = symbols('tau_p1:9')       # tau_p1 ... tau_p8
spinp = symbols('s_p1:9')        # s_p1 ... s_p8

xp = symbols('x_p1:9')
yp = symbols('y_p1:9')
zp = symbols('z_p1:9')

xvecp = [
    Matrix([[xp[i]], [yp[i]], [zp[i]]])
    for i in range(8)
]


n_lift = Matrix([[0], [0], [-1]])

F_lift_props_b = Matrix([[0], [0], [0]])
M_lift_props_b = Matrix([[0], [0], [0]])

for i in range(8):
    F_i, M_i = prop_force_moment(
        Thp[i],
        taup[i],
        spinp[i],
        xvecp[i],
        n_lift
    )

    F_lift_props_b += F_i
    M_lift_props_b += M_i


# Forward propeller 

Thc, tauc, spinc = symbols('T_c tau_c s_c')
xc, yc, zc = symbols('x_c y_c z_c')

xvecc = Matrix([[xc], [yc], [zc]])

# Assuming forward body x direction
n_cruise = Matrix([[1], [0], [0]])

F_cruise_b, M_cruise_b = prop_force_moment(
    Thc,
    tauc,
    spinc,
    xvecc,
    n_cruise
)


# Total forces and moments

F_ext_b = F_a1_b + F_lift_props_b + F_cruise_b

M_ext_b = M_a1_b + M_f1_b + M_lift_props_b + M_cruise_b


# Equations of motion

dvel_bre_b = (1/m)*F_ext_b + C_bfn*gvec - omg.cross(vel)

Ixx, Iyy, Izz, Ixy, Iyz, Ixz = symbols("I_{xx} I_{yy} I_{zz} I_{xy} I_{yz} I_{xz}")

J = Matrix([
    [Ixx, Ixy, Ixz],
    [Ixy, Iyy, Iyz],
    [Ixz, Iyz, Izz]
])

domg_bre_b = J.inv() * (M_ext_b - omg.cross(J*omg))


# State matrices

state = Matrix.vstack(pos, PHI, vel, omg)

surfaces = Matrix([[ail1], [ail2], [elev], [rudd]])

lift_thrusts = Matrix([[Thp[i]] for i in range(8)])
lift_torques = Matrix([[taup[i]] for i in range(8)])

cruise_inputs = Matrix([[Thc], [tauc]])

wind = Matrix([[wn], [we], [wd]])

inputs = Matrix.vstack(
    surfaces,
    lift_thrusts,
    lift_torques,
    cruise_inputs,
    wind
)


# State derivative

xdot = Matrix.vstack(
    dpos_brn_n,
    dPHI_brn_n,
    dvel_bre_b,
    domg_bre_b
)


# Linearization

A = xdot.jacobian(state)
B = xdot.jacobian(inputs)

print("State vector shape:", state.shape)
print("Input vector shape:", inputs.shape)
print("A matrix shape:", A.shape)
print("B matrix shape:", B.shape)