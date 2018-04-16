# 
# Copyright (c) 2015-2018 Idiap Research Institute, http://www.idiap.ch/
# Written by James Newling <james.newling@gmail.com>
# All rights reserved.
# 
# eakmeans is a library for exact and approximate k-means written in C++ and
# Python. This file is part of eakmeans. See file COPYING for more details.
# 
# This file is part of eakmeans.
# 
# eakmeans is free software: you can redistribute it and/or modify
# it under the terms of the 3-Clause BSD Licence. See
# https://opensource.org/licenses/BSD-3-Clause for more details.
# 
# eakmeans is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See file
# COPYING for more details.
# 
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
