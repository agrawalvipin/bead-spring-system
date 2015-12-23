#!/usr/bin/env python
import numpy
from numpy.linalg import norm
from scipy import constants

class BeadSpringSystem:
    """ Stores bead positions and connectivities and writes to data file.
      Attributes:
        beads (list): each element stores [mol, bead_type, [x,y,z]].
        bonds (list): each element stores [bond_type, i, j].
        angles (list): each elemen tstores [angle_type, i, j, k].
        box_length (float): stores the edge length of the simulation box.
        bead_masses (list): each element stores the bead mass (float).
        angle_types (list): each element stores the angle type (int).
    """
    def __init__(self, box_length=1.0):
        """ Initialize an empty system. """
        self.beads      = []
        self.bonds      = []
        self.angles     = []
        self.box_length = box_length
        self.bead_masses = []
        self.angle_types = []

    def write_to_lammps(self, lmpdatafile):
        """ Writes an input file to LAMMPS data format, see:
            http://lammps.sandia.gov/doc/read_data.html
        """
        with open(lmpdatafile, 'w') as f:
            # Get number of unique bead and bond types.
            num_bead_types = len(set([b[1] for b in self.beads]))
            num_bond_types = len(set([b[0] for b in self.bonds]))
            num_angle_types= len(set([b[0] for b in self.angles]))
            f.write('LAMMPS 2005 data file for soft beads\n\n')
            f.write('% 7d atoms\n'%len(self.beads))
            f.write('% 7d bonds\n'%len(self.bonds))
            if (num_angle_types>0):
                f.write('% 7d angles\n\n'%len(self.angles))
            f.write('% 4d atom types\n' %num_bead_types)
            f.write('% 4d bond types\n' %num_bond_types)
            if (num_angle_types>0):
                f.write('% 4d angle types\n\n' %num_angle_types)
            f.write(' %15.9f %15.9f xlo xhi\n'%   (0.0, self.box_length))
            f.write(' %15.9f %15.9f ylo yhi\n'%   (0.0, self.box_length))
            f.write(' %15.9f %15.9f zlo zhi\n\n'% (0.0, self.box_length))
            f.write('Masses\n\n')
            for i,m in enumerate(self.bead_masses):
                f.write('%d %f\n'%(i+1,m))
            f.write('\nAtoms\n\n')
            for index,bead in enumerate(self.beads):
                bead = tuple([index+1] + bead[0:2] + list(bead[2]))
                f.write(' % 6d % 6d % 3d %.3f %.3f %.3f\n' %bead)
            f.write('\nBonds\n\n')
            for index,bond in enumerate(self.bonds):
                bond = tuple([index+1] + list(bond))
                f.write('% 6d % 3d % 6d % 6d\n' %bond)
            if (num_angle_types>0):
                f.write('\n')
                f.write('Angles\n\n')
                for index,angle in enumerate(self.angles):
                    angle = tuple([index+1] + list(angle))
                    f.write('% 6d % 3d % 6d % 6d %6d\n' %angle)

def create(param, box_length):
    """ Creates a bead spring system. 
    Arguments:
        param:   class of input parameters defined by command line arguments.
        box_length:  the size of the sides of the cubic system box.
    """
    system = BeadSpringSystem(box_length)
    # Define atom_types such that for each bead there is a unique number.
    unique_types = sorted(set(''.join(param.block)))
    system.atom_types = {t:i+1 for i,t in enumerate(unique_types)}

    # Gets bond types for each combination of pairs of beads.
    system.bond_types = {}
    for blk in param.block:
        for a,b in zip(blk, blk[1:]):
            pair = min(a,b) + max(a,b)
            if pair not in system.bond_types:
                system.bond_types[pair] = len(system.bond_types)+1

    angle_types = {}
    if param.angles is not None:
        for blk in range(len(param.block)):
            fullblock = param.nblk[blk]*param.block[blk]
            for i in range(len(fullblock)-2):
                angle = fullblock[i:i+3]
                if angle[0] > angle[-1]:
                    angle = angle[::-1]

                if angle not in angle_types:
                    angle_types[angle]=len(angle_types)+1

    # If there are more than one blocks then we need to check last, first.
    for blk in range(len(param.block)):
        if param.nblk[blk] > 1:
            pair = ''.join(sorted([param.block[blk][0], param.block[blk][-1]]))
            if pair not in system.bond_types:
                system.bond_types[pair] = len(system.bond_types)+1

    system.angle_types = angle_types if param.angles is not None else None
    print 'Block type is ', param.block
    print 'Atom types are:', system.atom_types
    print 'Bond types are:', system.bond_types
    if param.angles is not None: print 'Angle types are:', system.angle_types

    # Given the index of an atom in a chain, return the type, blocktype j
    def get_atom_type(i,j):
        atom = param.block[j][i%len(param.block[j])]
        return system.atom_types[atom]

    # Given the index, i of the first atom in a chain, blocktype j
    # return the bond type between i and i+1.
    def get_bond_type(i,j):
        itype = param.block[j][(i)  %len(param.block[j])]
        jtype = param.block[j][(i+1)%len(param.block[j])]

        pair = min(itype,jtype) + max(itype,jtype)
        return system.bond_types[pair], pair

    # return the angle type between i-1, i and i+1. block type j
    def get_angle_type(i,j):
        fullblock = param.nblk[j]*param.block[j]
        itype = fullblock[(i-1)%len(fullblock)]
        jtype = fullblock[  (i)%len(fullblock)]
        ktype = fullblock[(i+1)%len(fullblock)]
        triplet = min(itype,ktype) + jtype+ max(itype,ktype)
        return angle_types[triplet], triplet

    def random_unit():
        """ Returns a random vector with unit length. """
        r = numpy.random.rand(3)-0.5
        return r/numpy.linalg.norm(r)
        
    chain_size = [len(b)*n for b,n in zip(param.block, param.nblk)]
    nmol = 0
    for blki in range(len(param.block)):
        for molecule in range(param.nchain[blki]):
            first_bead = system.box_length*numpy.random.rand(3)
            system.beads.append([nmol, get_atom_type(0,blki), first_bead])
            for z in range(chain_size[blki] - 1):
                # What bond type am I and what is the length
                btype,bpair = get_bond_type(z,blki)
                ran_vect  = param.bond_length[bpair] * random_unit()
                next_bead = system.beads[-1][2] + ran_vect

                # Check if angle ijk is not near zero degrees.
                if z > 0:
                    while True:
                        btypeij,bpairij = get_bond_type(z-1,blki)
                        btypejk,bpairjk = get_bond_type(z,blki)
                        rij     = system.beads[-2][2] - system.beads[-1][2]
                        rjk     = next_bead - system.beads[-1][2]
                        cos_ijk = numpy.dot(rij,rjk)/norm(rij)/norm(rjk)

                        if param.angles is not None:
                            atype_ijk,triplet = get_angle_type(z,blki)
                            theta_ijk = param.angles[triplet]
                            max_cos = numpy.cos((theta_ijk+15.0)*numpy.pi/180.0)
                            min_cos = numpy.cos((theta_ijk-15.0)*numpy.pi/180.0)
                            if cos_ijk > max_cos and cos_ijk < min_cos: 
                                break
                        elif cos_ijk < 0.8: 
                            break
                        ran_vect  = param.bond_length[bpairjk] * random_unit()
                        next_bead = system.beads[-1][2] + ran_vect

                system.bonds.append([btype, len(system.beads), len(system.beads)+1])
                if param.angles is not None and z > 0:
                    atype_ijk,triplet = get_angle_type(z,blki)
                    n = len(system.beads)
                    system.angles.append([atype_ijk,n-1,n,n+1])
                system.beads.append([nmol, get_atom_type(z+1,blki), next_bead])
            nmol += 1
    return system

def make_system(path, param, density):
    """ Makes chains of atoms in a box.
    Arguments:
      path:    where to write the lammps input file.
      param:   class of input parameters defined by command line arguments.
      density: target system density (in g/cc) used to compute box size.
    """
    print 'Bead masses are: '
    for k,v in param.bead_mass.iteritems():
        print '  {}: {}'.format(k,v)

    mass = 0.0
    for blk,nb,nc in zip(param.block, param.nblk, param.nchain):
        mass += nb*nc*sum(param.bead_mass[b] for b in blk)
    numbead = sum(nb*nc*len(b) for nb,nc,b in 
                    zip(param.nblk, param.nchain, param.block))
    print '# of beads:    {:8d}'.format(numbead)
    print 'Density is:    {:8.4f} g/cc'.format(density)
    density *= 1e-24*constants.N_A 
    print 'Density is:    {:8.4f} g/mol/A^3'.format(density)
    volume = mass / density 

    box_length = volume**(1.0/3.0)
    print 'Volume is:     {:8.1f} A^3'.format(volume)
    print 'box length is: {:8.3f} A'.format(box_length)
    system = create(param, box_length)

    for ai in sorted(set(''.join(param.block))):
        system.bead_masses.append(param.bead_mass[ai])

    system.write_to_lammps(path)
    return system.bond_types,system.atom_types,system.angle_types

def main():
    """ Standalone mode (called from command line. """
    import optparse 
    parser = optparse.OptionParser()
    parser.add_option('', '--num_chains', dest='nchains', default=40,
                      help='sets number of chains')
    parser.add_option('', '--num-blocks', dest='nblocks', default=14,
                      help='sets number of blocks')
    parser.add_option('', '--block', dest='blockstr', default='S',
                      help='sets beads in block')
    parser.add_option('', '--filename',  dest='path', default='bead_system.lammps',
                      help='specifies output file')
    parser.add_option('','--bond_length', dest='r0', default= 5,
                      help='sets the bond length')
    parser.add_option('','--angle', dest='t0', default=None,
                      help='sets the angle')
    parser.add_option('','--density', type='float', dest='rho', default=1.0,
                      help='sets the system density')
    parser.add_option('','--bead_masses', dest='masses', default= 72.10,
                      help='sets the bead masses')
    parser.add_option('','--seed', dest='seed', default=-1, type=int,
                        help='Sets the seed for the random number generator.')

    opt,args = parser.parse_args()
    if opt.seed > 0: numpy.random.seed(opt.seed)

    bonds = opt.r0.strip().split(',')
    angles = opt.t0.strip().split(',') if opt.t0 is not None else None
    masses = opt.masses.strip().split(',')

    class SystemParameters: pass
    """ Holds input parameters passed in from the command line.
    Attributes:
      angles: (dict: string->float): equilibrium bond angles.
      bead_mass (dict: string->float): stores the mass of each bead type.
      block (list): repeating motif of beads within a chain (per chain type).
      bond_length: (dict: string->float): equilibrium bond lengths.
      nblk (list): number of repeated blocks in a chain (per chain type).
      nchain (list): number of chains in the system (per chain type).
    """
    p = SystemParameters()
    p.nchain = [int(s) for s in opt.nchains.strip().split(',')]
    p.nblk = [int(s) for s in opt.nblocks.strip().split(',')]
    p.block = opt.blockstr.strip().split(',')
    p.bead_mass = dict((n,float(v)) for n,v in (a.split('=') for a in masses))
    p.bond_length = dict((n,float(v)) for n,v in (a.split('=') for a in bonds))
    if angles:
        p.angles = dict((n,float(v)) for n,v in (a.split('=') for a in angles))
    else:
        p.angles = None
    make_system(opt.path, p, opt.rho)

# If called at top level.
if __name__ == '__main__':
    main()
