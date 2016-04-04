#Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
#Written by James Newling <jnewling@idiap.ch>

#eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.

import sys
import os
import commands
import shutil

names = commands.getstatusoutput('find .. -name "*.hpp" -type "f"')[1].split("\n")

bobs = []
for n in names:
	shutil.copy(n, 	n.split("/")[-1])
	
	#if "whileprototying" not in n and "test" not in n and "experiments" not in n and "junk" not in n:
		#if "util" in n and "main" in n:
			#pass
			##print "-------->  ", n  
		#else:
			#bobs.append(n)

#for b in bobs:
	#shutil.copy(b, 	b.split("/")[-1])
	#print b

		#print n

	#bobs.append(n.split("/")[-1])
#bobs.sort()

#for b in bobs:
	#print b
##for n in names:
	##if "arrutilv2l0" not in n and "kmeansstandalone" not in n:
		##shutil.copy(n, 	n.split("/")[-1])
