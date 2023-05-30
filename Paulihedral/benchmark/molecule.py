from qiskit_nature.drivers import UnitsType, Molecule
from qiskit_nature.drivers.second_quantization import (
    ElectronicStructureDriverType,
    ElectronicStructureMoleculeDriver,
)
from qiskit_nature.problems.second_quantization import ElectronicStructureProblem
from qiskit_nature.converters.second_quantization import QubitConverter
from qiskit_nature.mappers.second_quantization import JordanWignerMapper, ParityMapper,BravyiKitaevMapper
from .mypauli import pauliString


H2 = [['H', [0., 0, 0]],['H', [0, 0, -1.5]]]
He2 = [['He', [0., 0, 0]],['He', [0, 0, 2.6740000]]]
HF = [['F', [0., 0, 0]],['H', [0, 0, 0.9153800]]]
HCl = [['Cl', [0., 0., 0.]],['H', [0, 0, 1.2744500]]]
ClF = [['F', [0., 0., 0.]],['Cl', [0., 0., 1.6303300]]]
LiH = [['Li', [0,0,0]],['H', [0,0,0.76]]]
N2 = [['N', [0., 0, 0]],['N', [0, 0, -1.5]]]
NH2 = [['N', [0.,0.,0.]],['H', [0.,0.,1.0240500]],['H', [0.9971600,0.,-0.2331400]]]  # NH2 charge = 0, multiplicity=2
NH3 = [['N', [0.,0.,0.1128900]],['H', [0., 0.9380200,-0.2634100]],
    ['H', [0.8123500,-0.4690100,-0.2634100]],
    ['H', [-0.8123500,-0.4690100,-0.2634100]]]
H2O = [['O', [0.,0.,0.]],['H', [0.95,-0.55,0.]],['H', [-0.95,-0.55,0.]]]
H2S = [['S', [0., 0, 0]],['H', [0, 0, -1.5]],['H', [0, 0, 1.5]]]
CH2 = [['C', [0.,0.,0.]],['H', [0.,0.,1.1077100]],['H', [1.0837800,0.,-0.2289900]]]
CH4 = [['C', [0.,1.0,1.0]],['H', [0.051054399,0.051054399,0.051054399]],
    ['H', [1.948945601,1.948945601,0.051054399]],
    ['H', [0.051054399,1.948945601,1.948945601]],
    ['H', [1.948945601,0.051054399,1.948945601]]]
MgO = [['Mg', [0., 0, 0]],['O', [0, 0, -1.5]]]
CO2 = [['O', [1.4, 0., 0.]],['C', [0., 0., 0.]],['O', [-1.4, 0., 0.]]]
NaCl = [['Na', [0., -1.5, -1.5]],['Cl', [1.5, -1.5, -1.5]]]
KOH = [['K', [3.4030, 0.2500, 0.0]],['O', [2.5369, -0.2500, -1.5]],['H', [2.0,0.06,0.0]]]
FeO = [['Fe', [0., 0, 0]],['O', [0, 0, -1.5]]]

def get_qubit_op(geo):
    molecule = Molecule(geometry=geo, charge=0, multiplicity=2)
    driver = ElectronicStructureMoleculeDriver(
    molecule, basis="sto3g", driver_type=ElectronicStructureDriverType.PYSCF)
    es_problem = ElectronicStructureProblem(driver)
    second_q_op = es_problem.second_q_ops()
    # JordanWigner mapping
    qubit_converter = QubitConverter(mapper=BravyiKitaevMapper())
    qubit_op = qubit_converter.convert(second_q_op[0])
    # Parity mapping
    '''qubit_converter = QubitConverter(mapper=ParityMapper(), two_qubit_reduction=True)
    qubit_op = qubit_converter.convert(second_q_op[0], num_particles=es_problem.num_particles)'''

    return qubit_op

def gene_molecule_oplist(atom_geo):
    qubit_op = get_qubit_op(atom_geo)
    oplist = []
    for i in qubit_op.primitive:
        oplist.append([pauliString(str(i.paulis[0]), coeff=i.coeffs[0])])
    return oplist

