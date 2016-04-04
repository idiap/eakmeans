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
