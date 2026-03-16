import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from scipy.integrate import solve_ivp

plt.rcParams['font.family'] = 'serif'           
plt.rcParams['font.size'] = 12                  
plt.rcParams['axes.titlesize'] = 14             
plt.rcParams['axes.labelsize'] = 12             
plt.rcParams['axes.linewidth'] = 1.5            
plt.rcParams['lines.linewidth'] = 2.5           
plt.rcParams['xtick.direction'] = 'in'          
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.size'] = 6            
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['xtick.major.width'] = 1.5         
plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['legend.frameon'] = True           
plt.rcParams['legend.edgecolor'] = 'black'
plt.rcParams['grid.alpha'] = 0.5                
plt.rcParams['grid.linestyle'] = '--'           


#ZMZ-409 ENGINE INITIAL DATA
D = 0.0955      # Cylinder bore, m
S = 0.094       # Piston stroke, m
CR = 9.0        # Compression ratio
L_k = 0.297     # Crank radius to connecting rod length ratio
R_crank = S / 2.0
V_h = (np.pi * D**2 / 4.0) * S  # Swept volume
V_c = V_h / (CR - 1.0)          # Clearance volume
gamma = 1.3                     # Specific heat ratio

#CYLINDER KINEMATICS
def kinematics(phi_deg):
    phi_rad = np.radians(phi_deg)
    #Instantaneous cylinder volume
    x = R_crank * (1 - np.cos(phi_rad) + (1/L_k) * (1 - np.sqrt(1 - (L_k * np.sin(phi_rad))**2)))
    V = V_c + (np.pi * D**2 / 4.0) * x
    #Volume derivative with respect to crank angle
    dV_dphi_rad = (np.pi * D**2 / 4.0) * R_crank * np.sin(phi_rad) * (1 + (L_k * np.cos(phi_rad)) / np.sqrt(1 - (L_k * np.sin(phi_rad))**2))
    dV_dphi = dV_dphi_rad * (np.pi / 180.0) 
    return V, dV_dphi

#WORKING CYCLE CALCULATION FUNCTION
def calculate_cycle(ethanol_percent, spark_advance):
    #Interpolating fuel properties based on ethanol fraction
    Hu = 44.0 - (44.0 - 38.84) * (ethanol_percent / 30.0)
    l0 = 14.7 - (14.7 - 12.99) * (ethanol_percent / 30.0)
    phi_z = 60.0 - (60.0 - 50.0) * (ethanol_percent / 30.0) # Faster combustion for ethanol
    
    phi_0 = spark_advance
    m_wiebe = 2.0
    
    def engine_diff_eq(phi_deg, P):
        V, dV_dphi = kinematics(phi_deg)
        
        #Wiebe heat release function
        if phi_0 <= phi_deg <= (phi_0 + phi_z):
            term = (phi_deg - phi_0) / phi_z
            dx_dphi = 6.908 * (m_wiebe + 1) * (term**m_wiebe) * np.exp(-6.908 * term**(m_wiebe + 1)) / phi_z
        else:
            dx_dphi = 0.0
            
        m_fuel = (V_h * 1.15) / l0
        Q_cycle = m_fuel * Hu * 1e6
        dQ_dphi = Q_cycle * dx_dphi
        
        P_val = P[0] if isinstance(P, np.ndarray) else P
        
        #First Law of Thermodynamics
        dP_dphi = ((gamma - 1) / V) * dQ_dphi - gamma * (P_val / V) * dV_dphi
        return [dP_dphi]

    phi_eval = np.linspace(-180, 180, 360)
    sol = solve_ivp(engine_diff_eq, (-180, 180), [100000], t_eval=phi_eval, method='RK45', max_step=1.0)
    
    P_bar = sol.y[0] / 100000 
    V_liters = [kinematics(ang)[0] * 1000 for ang in phi_eval]
    return phi_eval, P_bar, V_liters

#INTERFACE AND PLOTS SETUP
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
plt.subplots_adjust(bottom=0.25)

#Baseline Gasoline (strict gray dashed line)
phi_base, P_base, V_base = calculate_cycle(0, -10)
line_base_phi, = ax1.plot(phi_base, P_base, color='#555555', linestyle='--', label='Baseline (Gasoline, SA -10°)', zorder=1)
line_base_v, = ax2.plot(V_base, P_base, color='#555555', linestyle='--', label='Baseline (Gasoline)', zorder=1)

#Ethanol blend (deep burgundy academic color)
line_p_phi, = ax1.plot(phi_base, P_base, color='#900C3F', label='Ethanol-Gasoline Blend', zorder=2)
line_p_v, = ax2.plot(V_base, P_base, color='#900C3F', label='Ethanol-Gasoline Blend', zorder=2)

ax1.set_title('In-Cylinder Pressure Diagram (P-θ)')
ax1.set_xlabel('Crank Angle, CAD')
ax1.set_ylabel('In-Cylinder Pressure P, bar')
ax1.grid(True)
ax1.legend(loc='upper right')
ax1.set_ylim(0, 80)
ax1.set_xlim(-180, 180)

ax2.set_title('Indicator Diagram (P-V)')
ax2.set_xlabel('Cylinder Volume V, L')
ax2.set_ylabel('In-Cylinder Pressure P, bar')
ax2.grid(True)
ax2.set_ylim(0, 80)

#SLIDERS (INTERACTIVE)
ax_eth = plt.axes([0.15, 0.1, 0.65, 0.03])
ax_spark = plt.axes([0.15, 0.05, 0.65, 0.03])
slider_eth = Slider(ax_eth, 'Ethanol Fraction (%)', 0.0, 30.0, valinit=0.0, valstep=1.0)
slider_spark = Slider(ax_spark, 'Spark Advance (deg bTDC)', -30.0, 0.0, valinit=-10.0, valstep=1.0)

def update(val):
    eth = slider_eth.val
    spark = slider_spark.val
    _, P_new, V_new = calculate_cycle(eth, spark)
    line_p_phi.set_ydata(P_new)
    line_p_v.set_ydata(P_new)
    fig.canvas.draw_idle()

slider_eth.on_changed(update)
slider_spark.on_changed(update)

plt.show()
