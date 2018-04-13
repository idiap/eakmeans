#Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
#Written by James Newling <jnewling@idiap.ch>

#eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.

from IPython.core.debugger import Tracer



rawheader = r"""Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.

"""

hashheader = r"""#Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
#Written by James Newling <jnewling@idiap.ch>

#eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.

"""


cppheader = r"""/*
%s
*/

"""%(rawheader, )

import os

import commands


hfiles = commands.getstatusoutput("find . -name \"*.h\"")[1].split("\n")
cppfiles = commands.getstatusoutput("find . -name \"*.cpp\"")[1].split("\n")
hppfiles = commands.getstatusoutput("find . -name \"*.hpp\"")[1].split("\n")
cppheaderable = hfiles + cppfiles + hppfiles

makefiles = commands.getstatusoutput("find . -name \"Makefile\"")[1].split("\n")
pyfiles =  commands.getstatusoutput("find . -name \"*.py\"")[1].split("\n")
pyxfiles = commands.getstatusoutput("find . -name \"*.pyx\"")[1].split("\n")
pyxbldfiles = commands.getstatusoutput("find . -name \"*.pyxbld\"")[1].split("\n")
hashheaderable = makefiles + pyfiles + pyxfiles + pyxbldfiles



#allfiles = commands.getstatusoutput("find . -type f")[1].split("\n")
#for f in allfiles:
    #if f not in cppheaderable and f not in hashheaderable and f not in rawheaderable:
        #print f
        #if "dat" not in f:
            #os.remove(f)
        ##if ".so" in f:
            ##print f
            ##os.remove(f)
        ##if ".o" in f:
            ##print f

            ##os.remove(f)


if True == True:
    for files, header in zip([cppheaderable, hashheaderable], [cppheader, hashheader]): #rawheader,  rawheaderable,
        for fn in files:

            if fn:
                print "headering ", fn, "..."

                filly = open(fn, "r")
                lines = filly.read()
                filly.close()

                filly = open(fn, "w")
                filly.write(header)
                filly.write(lines)
                filly.close()



