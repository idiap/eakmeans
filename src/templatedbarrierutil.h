/*
Copyright (c) 2015-2018 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <james.newling@gmail.com>
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
*/

#ifndef BARRIERUTIL_H
#define BARRIERUTIL_H

#include <thread>
#include <mutex> 
#include <vector>
#include <condition_variable>

#include <chrono>
#include <iostream>

namespace stdthreadutil{

/* btask for `barriered task'
 * used for situation when 
 * (1) several threads perform same task
 * (2) when all finished, one thread performs finishing task
 * (3) all threads are released */
 
template <typename TInt> 
inline void btask(
/* identity of thread, not used in function but left as parameter as potentially useful for debugging */
const TInt & ti, 
/* number of threads performing the task*/
const TInt & nthreads, 
/* reference to number of threads which have completed the task*/
TInt & completions, 
/* when work of this thread is complete, use the following tools (shared by all workers on this task) to notify the others */
std::mutex & workend_mutex,
std::condition_variable & condvar, 
/* task to perform, and task to perform at end if this thread finishes last*/
const std::function<void()> & task, 
const std::function<void()> & endtask)
{
			
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


/* btasks for `barrirered tasks'
 * used in situation when several threads of equal status perform
 for (# tasks) { do task | do end task if last to complete | }
 */
template <typename TInt> 
inline void btasks(
const TInt & ti,
const TInt & nthreads,
std::vector<TInt> & section_completions,
std::vector<std::mutex> & sectionend_mutexes,
std::vector<std::condition_variable> & section_condvars, 
const std::vector<std::function<void(TInt)>> & section_tasks,
const std::vector<std::function<void()>> &  sectionend_tasks
){
	for (TInt si = 0; si < section_tasks.size(); ++si){
		btask(ti, nthreads, section_completions[si], sectionend_mutexes[si], section_condvars[si], [&section_tasks, si, ti](){section_tasks[si](ti);}, sectionend_tasks[si]);
	}
}

/* btask_rbasks for `barriered task then repeat barriered tasks' 
 * used in situation when 
 * several threads of equal status perform
 * do initialisation task | do intitialisation end task if last to complete | 
 * while (condition is true) for (# tasks) do task | do end task if last to complete |
 */
template <typename TInt>
void btask_rbtasks(
TInt ti, 
TInt nthreads,
std::vector<TInt> & section_completions, 
std::vector<std::mutex> & sectionend_mutexes, 
std::vector<std::condition_variable> & section_condvars, 
const std::function<void(TInt)> & initialisation_task,
const std::function<void()> & initialisationend_task,  
const std::vector<std::function<void(TInt)>> & section_tasks,
const std::vector<std::function<void()>> &  sectionend_tasks, 
const std::function<bool()> & getiscomplete){

	TInt nsections = section_tasks.size();

	/* piggy back on last task's mutex & condvar to perform initialisations*/
	btask(ti, nthreads, section_completions[nsections - 1], sectionend_mutexes[nsections -1], section_condvars[nsections -1], [&initialisation_task, ti](){initialisation_task(ti);}, initialisationend_task);
	
	while (!getiscomplete()){
		btasks(ti, nthreads, section_completions, sectionend_mutexes, section_condvars, section_tasks, sectionend_tasks);
	}
}

/* barriered tasks (inititialisation) then while condition repeat barriered tasks  */
template <typename TInt>
void btasks_rbtasks(
TInt ti, 
TInt nthreads,
std::vector<TInt> & x_completions, 
std::vector<std::mutex> & xend_mutexes, 
std::vector<std::condition_variable> & x_condvars, 
const std::vector< std::function<void(TInt)>> & initialisation_tasks,
const std::vector< std::function<void()>> & initialisationend_tasks,  
const std::vector< std::function<void(TInt)>> & section_tasks,
const std::vector< std::function<void()>> &  sectionend_tasks, 
const std::function<bool()> & getiscomplete){


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
template <typename TInt>
int launch_btasks(
TInt nthreads, 
const std::vector<std::function<void(TInt)>> & section_tasks,
const std::vector<std::function<void()>> & sectionend_tasks
)
{

	TInt nsections = section_tasks.size();

	if (nsections != sectionend_tasks.size()){
		throw std::runtime_error("section tasks and sectionend tasks are not the same, they should be");	
	}
	
	if (nsections == 1){
		throw std::runtime_error("with just one section lapping may occur, very bad. We suggest adding another section with trivial work");
	}
	
	
	
	std::vector<std::thread> threads;	
	std::vector<TInt> section_completions(nsections, 0);
	std::vector<std::mutex> sectionend_mutexes (nsections); 
	std::vector<std::condition_variable> section_condvars (nsections);

	for (TInt ti = 0; ti < nthreads; ++ti){
		threads.emplace_back(btasks<TInt>, ti, nthreads, std::ref(section_completions), std::ref(sectionend_mutexes), std::ref(section_condvars), std::ref(section_tasks), std::ref(sectionend_tasks));
	}
	
	for (auto & t : threads){
		t.join();
	}
	
	return 0;
}

template <typename TInt>
int launch_btasks_rbtasks(
TInt nthreads, 
const std::vector<std::function<void(TInt)>> & initialisation_tasks,
const std::vector<std::function<void()>> & initialisationend_tasks,  
const std::vector<std::function<void(TInt)>> & section_tasks,
const std::vector<std::function<void()>> & sectionend_tasks, 
const std::function<bool()> & getiscomplete,
const std::function<void()> & closing_task 
)
{
	if (initialisation_tasks.size() != initialisationend_tasks.size()){
		throw std::runtime_error("initialisation tasks and initialisationend tasks are not the same size in btasks_rbtasks, they should be");	
	}
	
	if (section_tasks.size() != sectionend_tasks.size()){
		throw std::runtime_error("section tasks and section end tasks are are not the same size in btasks_rbtasks, they should be");	
	}
	
	if (section_tasks.size() == 1){
		throw std::runtime_error("with just one section lapping may occur, very bad. We suggest adding another section with trivial work (in btasks_rbtasks)");
	}
	
	TInt nsections = section_tasks.size();
	TInt ninittasks = initialisation_tasks.size();
	TInt nx = std::max(nsections, ninittasks);

	
	std::vector<std::thread> threads;	
	std::vector<TInt> section_completions (nx, 0);
	std::vector<std::mutex> sectionend_mutexes (nx); 
	std::vector<std::condition_variable> section_condvars (nx);

	for (TInt ti = 0; ti < nthreads; ++ti){
		threads.emplace_back(btasks_rbtasks<TInt>, ti, nthreads, std::ref(section_completions), std::ref(sectionend_mutexes), std::ref(section_condvars), std::ref(initialisation_tasks), std::ref(initialisationend_tasks), std::ref(section_tasks), std::ref(sectionend_tasks), std::ref(getiscomplete));
	}	
	
	for (auto & t : threads){
		t.join();
	}
	
	closing_task();
	
	return 0;
}

template <typename TInt>
int launch_btask_rbtasks(
TInt nthreads, 
const std::function<void(TInt)> & initialisation_task,
const std::function<void()> & initialisationend_task,  
const std::vector<std::function<void(TInt)>> & section_tasks,
const std::vector<std::function<void()>> & sectionend_tasks, 
const std::function<bool()> & getiscomplete,
const std::function<void()> & closing_task 
){
	return launch_btasks_rbtasks(nthreads, {initialisation_task}, {initialisationend_task}, section_tasks, sectionend_tasks, getiscomplete, closing_task);
}




}


#endif
