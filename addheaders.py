#Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
#Written by James Newling <jnewling@idiap.ch>

#eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.

from IPython.core.debugger import Tracer


old_rawheader = r"""Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.

"""

old_hashheader = r"""#Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
#Written by James Newling <jnewling@idiap.ch>

#eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.

"""

old_cppheader = r"""/*
%s
*/

"""%(old_rawheader, )

new_rawheader = r"""
Copyright (c) 2015-2018 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <>
All rights reserved.

eakmeans is a library for exact and approximate k-means written in C++ and
Python. This file is part of eakmeans. See file COPYING for more details.

This file is part of eakmeans.

eakmeans is free software: you can redistribute it and/or modify
it under the terms of the 3-Clause BSD Licence. See
https://opensource.org/licenses/BSD-3-Clause for more details.

eakmeans is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See file
COPYING for more details.
"""

new_hashheader = '\n'.join(['# ' + l for l in new_rawheader.split('\n')]) + '\n'
new_cppheader = "/*%s*/\n\n" % new_rawheader

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


for files, old_header, new_header in zip(
        [cppheaderable, hashheaderable],
        [old_cppheader, old_hashheader],
        [new_cppheader, new_hashheader]
        ):
    for fn in files:

        if fn:
            print "headering ", fn, "..."

            filly = open(fn, "r")
            lines = filly.read()
            filly.close()

            if lines.startswith(old_header):
                lines = lines[len(old_header):]

            filly = open(fn, "w")
            filly.write(new_header)
            filly.write(lines)
            filly.close()
