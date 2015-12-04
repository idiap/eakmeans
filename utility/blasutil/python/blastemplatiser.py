#EAKMeans is a fast Exact K-means library written in C++ with 
#command-line interface, shared library + header files and 
#Python bindings

#Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
#Written by James Newling <jnewling@idiap.ch>

#This file is part of EAKMeans.

#EAKMeans is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License version 3 as
#published by the Free Software Foundation.

#EAKMeans is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with EAKMeans. If not, see <http://www.gnu.org/licenses/>.

#python script to take cblas.h and make template functions, templated over floating type.


import sys
from IPython.core.debugger import Tracer;


print "This script should be run with caution, it automatically overwrites and regenerates C++ template code. It should not be necessary to run as an end-user. Comment out the following exit line if you wish to proceed."
sys.exit(1)


#change to location of cblas.h file:
blash_fn = '/idiap/user/jnewling/openblas/include/cblas.h'

#change to your preferred typename for floating point precision types:
template_float = "TFloatType"

#change to your preferred namespace:
templatenamespace = 'wblas'

#name of resulting header template file:
rootdir = '/idiap/home/jnewling/libraries/blasutil/'
target_header = rootdir + 'include/blastemplates.h'
target_source = rootdir + 'src/blastemplates.cpp'


blash_file = open(blash_fn)
cblash_lines = blash_file.readlines()
blash_file.close()


def get_signature(function_declaration):
	"""
	decompose a function declaration into
	rtype, ftype, fname, cvparamters : 
	(return type), {s, d, c, z}, "axpy", (const volatile pointer reference qualified function parameters)
	"""
	
	prebracket, postbracket = function_declaration.split('(')
	inbracket = postbracket.split(')')[0]
	rtype, fname = prebracket.split()
	cvparameters = inbracket.split(',')
	cvparamteres = [x.strip().replace("  ", " ") for x in cvparameters]
	fname = fname.split('cblas_')[1]
	if fname[0] in ['i', 'I']:
		ftype = fname[0:2]
		fname = fname[2::]
	else:
		ftype = fname[0]
		fname = fname[1::]
	
	return [rtype, ftype, fname, cvparameters]
	
def get_fragment(target = "saxpy"):
	"""
	return function declaration in cblas.h containing target
	"""
	target = "_" + target + "("
	function_lines = []
	
	nlines = len(cblash_lines)
	current_line_index = 0
	
	cblash_function_declaration = ""
	while current_line_index < nlines:
		current_line = cblash_lines[current_line_index].strip().replace(" (", "(")
		if current_line:
			if target in current_line:
				cblash_function_declaration = current_line
			
			else:
				if cblash_function_declaration:
					cblash_function_declaration += current_line
			
		if cblash_function_declaration:
			if ";" in cblash_function_declaration:
				function_lines.append(cblash_function_declaration)
				cblash_function_declaration = ""
					
		current_line_index += 1
	
	if len(function_lines) == 0:
		sys.exit("No function declaration containing %s found"%(target,))
	
	if len(function_lines) > 1:
		"More than 1 function declaration containing %s found:"%(target, )
		for l in function_lines:
			print l
			sys.stdout.flush()
		sys.exit("abort as %d targets found for %s"%(len(function_lines), target))
	
	
	#for f in function_lines:
		#f.replace("  ", " ")

	space_stripped = function_lines[0]
	while "  " in space_stripped:
		space_stripped = space_stripped.replace("  ", " ")
	return space_stripped

def get_generic(astring):
	return astring.replace('float', template_float).replace('double', template_float).replace("  ", " ")

def get_generic_signature(f):
	"""
	replace double and float with FloatType throughout
	"""
	if f[1][0] in ['i', 'I']:
		generic_signature = [get_generic(f[0]), 'ix', f[2], [get_generic(x) for x in f[3]]]
	else:
		generic_signature = [get_generic(f[0]), 'x', f[2], [get_generic(x) for x in f[3]]]
	
	return generic_signature		

def get_generic_function_declaration(g_sig):
	if g_sig[1][0] == 'i':
		generic_function_declaration = g_sig[0] + " i_" + g_sig[2] + "("

	else:
		generic_function_declaration = g_sig[0] + " " + g_sig[2] + "("
		
	for cvparmi, cvparm in enumerate(g_sig[3]):
		cvparm = cvparm.strip()
		while "  " in cvparm:
			cvparm.replace("  ", " ")
			
		generic_function_declaration += cvparm
		if cvparmi < len(g_sig[3]) - 1:
			generic_function_declaration += ", "
	generic_function_declaration += ");"
	
		
	return generic_function_declaration
		

def get_parameter_names(cvparameters):
	return [p.split()[-1].split("*")[-1].split("&")[-1] for p in cvparameters]
	
def get_function_call(function_line): 
	function_name = function_line.split()[1].split("(")[0]
	cvparms = [x.strip() for x in function_line.split("(")[1].split(")")[0].split(",")]
	p_names = get_parameter_names(cvparms)
	function_call = function_name + "("
	for cvparmi, cvparm in enumerate(p_names):
		function_call += cvparm
		if cvparmi < len(p_names) - 1:
			function_call += ", "
	function_call += ");"
	
	function_type = function_line.split()[0]
	if function_type != 'void':
		function_call = 'return ' + function_call
		#function_call += '\n\treturn value;'
	
	return function_call

	
		
	

def get_h_cpp_template_frags(generic_target = 'gemm', type_targets = 'sd'):
	

	if 'i_' not in generic_target:
		function_lines = [get_fragment('%s%s'%(x,generic_target)) for x in type_targets]
	else:
		function_lines = [get_fragment(generic_target.replace('_', x)) for x in type_targets]
	
	function_signatures = [get_signature(l) for l in function_lines]
	generic_signatures = [get_generic_signature(f) for f in function_signatures]
	
	generic_function_declarations = [get_generic_function_declaration(g_sig) for g_sig in generic_signatures]
	
		
	if len(generic_function_declarations) != len(type_targets):
		sys.exit("logic error")
	
	
	for i in range(len(type_targets)):
		if (generic_function_declarations[i] != generic_function_declarations[0]):
			print generic_function_declarations[i]
			print generic_function_declarations[0]
			Tracer()()
			sys.exit("function declarations do not match for different float types")
	
	header_string = "\ntemplate <typename %s>\n"%(template_float,) + generic_function_declarations[0]
		
	cpp_string = ""		
	for il, l in enumerate(function_lines):
		cpp_string += "\ntemplate <> \n"
		if 'i_' not in generic_target:
			cpp_string += l[0:-1].replace('cblas_%s%s'%(type_targets[il], generic_target), generic_target)
		else:
			cpp_string += l[0:-1].replace('cblas_%s'%(generic_target.replace('_', type_targets[il])), generic_target)
		cpp_string += "{\n\t"
		cpp_string += "%s\n"%(get_function_call(l))
		cpp_string += "}\n"
	
	return header_string, cpp_string
	

def get_all_float_targets():
	targets = []
	for l in cblash_lines:
		if ("cblas_s" in l):
			targets.append(l.split("cblas_s")[1].split("(")[0].strip())
		if ("cblas_is" in l): #is : index of single (used in for example ISAMAX)
			targets.append('i_' + l.split("cblas_is")[1].split("(")[0].strip())
			
	print len(targets)
	return targets



def make_h_cpp():
	h_file = "#ifndef BLASTEMPLATES_H \n#define BLASTEMPLATES_H \n#include <cblas.h> \n\nnamespace %s{"%(templatenamespace,)
	cpp_file = "#include \"blastemplates.h\"\n#include <cblas.h>\n\nnamespace %s{"%(templatenamespace,)
	
	all_targets = get_all_float_targets()
	
	#sdsdot not needed (use sdot)  
	#scasum is for complex, not single
	#scnrm2 not sure what this is
	
	
	for toremo in ['dsdot', 'casum', 'cnrm2']: 
		
		if toremo in all_targets:
			all_targets.remove(toremo)
	
	print all_targets
	for target in all_targets:
		#print "---------------------"
		#print target
		bla = get_h_cpp_template_frags(target, 'sd')
		#print "---------------------"
		h_file += bla[0]
		h_file += "\n"
		cpp_file += bla[1]
	

	h_file += "\n} \n\n#endif\n\n\n"
	cpp_file += "\n}\n\n\n"
	print h_file
	filly = open(target_header, 'w')
	filly.write(h_file)
	filly.close()
	filly = open(target_source, 'w')
	filly.write(cpp_file)
	filly.close()
	
	print cpp_file
	


