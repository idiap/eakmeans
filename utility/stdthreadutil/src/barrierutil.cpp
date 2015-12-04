/*
EAKMeans is a fast Exact K-means library written in C++ with 
command-line interface, shared library + header files and 
Python bindings

Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

This file is part of EAKMeans.

EAKMeans is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3 as
published by the Free Software Foundation.

EAKMeans is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with EAKMeans. If not, see <http://www.gnu.org/licenses/>.



*/

#include "barrierutil.h"
#include <iostream>

namespace stdthreadutil{


void btask(const size_t & ti, const size_t & nthreads, size_t & completions, std::mutex & workend_mutex,std::condition_variable & condvar, const std::function<void()> & task, const std::function<void()> & endtask) {
			
	/* do work. If mutex needed in task, use lambda function to wrap it in before passing to barrieredtask */
	task();

	/* lock a mutex and signal that work is complete */
	std::unique_lock<std::mutex> lk (workend_mutex);
	++ completions;	
	

	/* if this thread is the last thread to complete its work, do the final bit of work (combining results or whatever it may be) and then signal to the other threads that they may continue*/
	if (completions == nthreads){
		endtask();
		completions = 0;
		condvar.notify_all();
	}
	
	/* if this thread is not the last one to finish its work, wait for the final thread to finish and then continue */
	else{
		condvar.wait(lk, [&completions] {return (completions == 0);});
	}
	lk.unlock();
}

void btasks(
const size_t & ti, const size_t & nthreads, std::vector<size_t> & section_completions, std::vector<std::mutex> & sectionend_mutexes, std::vector<std::condition_variable> & section_condvars, const std::vector<std::function<void(size_t)>> & section_tasks, const std::vector<std::function<void()>> &  sectionend_tasks){
	for (size_t si = 0; si < section_tasks.size(); ++si){
		btask(ti, nthreads, section_completions[si], sectionend_mutexes[si], section_condvars[si], [&section_tasks, si, ti](){section_tasks[si](ti);}, sectionend_tasks[si]);
	}
}

void btask_rbtasks(size_t ti, size_t nthreads,std::vector<size_t> & section_completions, std::vector<std::mutex> & sectionend_mutexes, std::vector<std::condition_variable> & section_condvars, const std::function<void(size_t)> & initialisation_task,const std::function<void()> & initialisationend_task,  const std::vector<std::function<void(size_t)>> & section_tasks, const std::vector<std::function<void()>> &  sectionend_tasks, const std::function<bool()> & getiscomplete){
	size_t nsections = section_tasks.size();

	/* piggy back on last task's mutex & condvar to perform initialisations*/
	btask(ti, nthreads, section_completions[nsections - 1], sectionend_mutexes[nsections -1], section_condvars[nsections -1], [&initialisation_task, ti](){initialisation_task(ti);}, initialisationend_task);
	
	while (!getiscomplete()){
		btasks(ti, nthreads, section_completions, sectionend_mutexes, section_condvars, section_tasks, sectionend_tasks);
	}
}

void btasks_rbtasks(
size_t ti, size_t nthreads,std::vector<size_t> & x_completions, std::vector<std::mutex> & xend_mutexes, std::vector<std::condition_variable> & x_condvars, const std::vector< std::function<void(size_t)>> & initialisation_tasks,const std::vector< std::function<void()>> & initialisationend_tasks,  const std::vector< std::function<void(size_t)>> & section_tasks,const std::vector< std::function<void()>> &  sectionend_tasks, const std::function<bool()> & getiscomplete){

	/* if there is only one initialisation task, reduce problem to btask_rbtasks */
	if (initialisation_tasks.size() == 1){
		btask_rbtasks(ti, nthreads, x_completions, xend_mutexes, x_condvars, initialisation_tasks[0], initialisationend_tasks[0], section_tasks, sectionend_tasks, getiscomplete);
	}

	else{
		btasks(ti, nthreads, x_completions, xend_mutexes, x_condvars, initialisation_tasks, initialisationend_tasks);
		while (!getiscomplete()){
			btasks(ti, nthreads, x_completions, xend_mutexes, x_condvars, section_tasks, sectionend_tasks);
		}
	}
}




/* launch barriered tasks */
int launch_btasks(size_t nthreads, const std::vector<std::function<void(size_t)>> & section_tasks,const std::vector<std::function<void()>> & sectionend_tasks)
{

	size_t nsections = section_tasks.size();

	if (nsections != sectionend_tasks.size()){
		throw std::runtime_error("section tasks and sectionend tasks are not the same, they should be");	
	}
	
	if (nsections == 1){
		throw std::runtime_error("with just one section lapping may occur, very bad. We suggest adding another section with trivial work");
	}
	
	
	
	std::vector<std::thread> threads;	
	std::vector<size_t> section_completions(nsections, 0);
	std::vector<std::mutex> sectionend_mutexes (nsections); 
	std::vector<std::condition_variable> section_condvars (nsections);

	for (size_t ti = 0; ti < nthreads; ++ti){
		threads.emplace_back(btasks, ti, nthreads, std::ref(section_completions), std::ref(sectionend_mutexes), std::ref(section_condvars), std::ref(section_tasks), std::ref(sectionend_tasks));
	}
	
	for (auto & t : threads){
		t.join();
	}
	
	return 0;
}

int launch_btasks_rbtasks(size_t nthreads, const std::vector<std::function<void(size_t)>> & initialisation_tasks,const std::vector<std::function<void()>> & initialisationend_tasks,  const std::vector<std::function<void(size_t)>> & section_tasks,const std::vector<std::function<void()>> & sectionend_tasks, const std::function<bool()> & getiscomplete, const std::function<void()> & closing_task ){
	if (initialisation_tasks.size() != initialisationend_tasks.size()){
		throw std::runtime_error("initialisation tasks and initialisationend tasks are not the same size in btasks_rbtasks, they should be");	
	}
	
	if (section_tasks.size() != sectionend_tasks.size()){
		throw std::runtime_error("section tasks and section end tasks are are not the same size in btasks_rbtasks, they should be");	
	}
	
	if (section_tasks.size() == 1){
		throw std::runtime_error("with just one section lapping may occur, very bad. We suggest adding another section with trivial work (in btasks_rbtasks)");
	}
	
	size_t nsections = section_tasks.size();
	size_t ninittasks = initialisation_tasks.size();
	size_t nx = std::max(nsections, ninittasks);

	
	std::vector<std::thread> threads;	
	std::vector<size_t> section_completions (nx, 0);
	std::vector<std::mutex> sectionend_mutexes (nx); 
	std::vector<std::condition_variable> section_condvars (nx);

	for (size_t ti = 0; ti < nthreads; ++ti){
		threads.emplace_back(btasks_rbtasks, ti, nthreads, std::ref(section_completions), std::ref(sectionend_mutexes), std::ref(section_condvars), std::ref(initialisation_tasks), std::ref(initialisationend_tasks), std::ref(section_tasks), std::ref(sectionend_tasks), std::ref(getiscomplete));
	}	
	
	for (auto & t : threads){
		t.join();
	}
	
	closing_task();
	
	return 0;
}

int launch_btask_rbtasks(size_t nthreads, const std::function<void(size_t)> & initialisation_task,const std::function<void()> & initialisationend_task,  const std::vector<std::function<void(size_t)>> & section_tasks,const std::vector<std::function<void()>> & sectionend_tasks, const std::function<bool()> & getiscomplete,const std::function<void()> & closing_task ){
	return launch_btasks_rbtasks(nthreads, {initialisation_task}, {initialisationend_task}, section_tasks, sectionend_tasks, getiscomplete, closing_task);
}




}
