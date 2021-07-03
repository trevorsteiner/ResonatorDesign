from simphony.library import siepic
from simphony.netlist import Subcircuit
from simphony.simulation import SweepSimulation
import matplotlib.pyplot as plt

# Declare the models used in the circuit
gc = siepic.ebeam_gc_te1550()                           # grating coupler
y = siepic.ebeam_y_1550()                               # y-branch
wg150 = siepic.ebeam_wg_integral_1550(length=150e-6)    # 150 micron waveguide
wg50 = siepic.ebeam_wg_integral_1550(length=50e-6)      # 50 micron waveguide

# Create the circuit, add all individual instances
circuit = Subcircuit('MZI')
e = circuit.add([
    (gc, 'input'),
    (gc, 'output'),
    (y, 'splitter'),
    (y, 'recombiner'),
    (wg150, 'wg_long'),
    (wg50, 'wg_short'),
])

# You can set pin names individually:
circuit.elements['input'].pins['n2'] = 'input'
circuit.elements['output'].pins['n2'] = 'output'

# Or you can rename all the pins simultaneously:
circuit.elements['splitter'].pins = ('in1', 'out1', 'out2')
circuit.elements['recombiner'].pins = ('out1', 'in2', 'in1')

# Circuits can be connected using the elements' string names:
circuit.connect_many([
    ('input', 'n1', 'splitter', 'in1'),
    ('splitter', 'out1', 'wg_long', 'n1'),
    ('splitter', 'out2', 'wg_short', 'n1'),
    ('recombiner', 'in1', 'wg_long', 'n2'),
    ('recombiner', 'in2', 'wg_short', 'n2'),
    ('output', 'n1', 'recombiner', 'out1'),
])

# Run a simulation on the netlist.
simulation = SweepSimulation(circuit, 1500e-9, 1600e-9)
result = simulation.simulate()

f, s = result.data('input', 'output')
plt.plot(f, s)
plt.title("MZI")
plt.tight_layout()
plt.show()