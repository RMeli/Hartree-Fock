"""
    Copyright (C) 2015 Rocco Meli

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from basis import *

import numpy as np

# H2
H2 = [Atom("H",(0,0,0),1,["1s"]),Atom("H",(0,0,1.4),1,["1s"])]
sto3g_H2 = STO3G(H2)

# HeH+
HeH = [Atom("He",(0,0,1.4632),2,["1s"]),Atom("H",(0,0,0),1,["1s"])]
sto3g_HeH = STO3G(HeH)

# He
He = [Atom("He",(0,0,0),2,["1s"])]
sto3g_He = STO3G(He)

# B
Be = [Atom("Be",(0,0,0),4,["1s","2s"])]
sto3g_Be = STO3G(Be)

# N2
N2 = [  Atom("N",(0,0,0),7,["1s","2s","2p"]),
        Atom("N",(0,0,2.074),7,["1s","2s","2p"])]
sto3g_N2 = STO3G(N2)

# HF
HF = [  Atom("H",(0,0,0),1,["1s"]),
        Atom("F",(0,0,1.807),9,["1s","2s","2p"])]
sto3g_HF = STO3G(HF)

# B
O = [Atom("O",(0,0,0),8,["1s","2s","2p"])]
sto3g_O = STO3G(O)

# H2O
# Szabo and Ostlund
H2O = [ Atom("H",(1.809*np.sin(104.52/180*np.pi/2),0,0),1,["1s"]),
        Atom("H",(-1.809*np.sin(104.52/180*np.pi/2),0,0),1,["1s"]),
        Atom("O",(0,1.809*np.cos(104.52/180*np.pi/2),0),8,["1s","2s","2p"])]
# Mathematica Journal
#H2O = [ Atom("H",(0,1.43233673,-0.96104039),1,["1s"]),
#        Atom("H",(0,-1.43233673,-0.96104039),1,["1s"]),
#        Atom("O",(0,0,0.24026010),8,["1s","2s","2p"])]
sto3g_H2O = STO3G(H2O)
