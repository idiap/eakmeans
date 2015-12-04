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

#ifndef BARRIERUTIL_H
#define BARRIERUTIL_H

#include <thread>
#include <mutex> 
#include <vector>
#include <condition_variable>
#include <chrono>

namespace stdthreadutil{
	

/* btask for `barriered task'
 * used for situation when 
 * (1) several threads perform same task
 * (2) when all finished, one thread performs finishing task
 * (3) all threads are released */
 
void btask(
/* identity of thread, not used in function but left as parameter as potentially useful for debugging */
const size_t & ti, 
/* number of threads performing the task*/
const size_t & nthreads, 
/* reference to number of threads which have completed the task*/
size_t & completions, 
/* when work of this thread is complete, use the following tools (shared by all workers on this task) to notify the others */
std::mutex & workend_mutex,
std::condition_variable & condvar, 
/* task to perform, and task to perform at end if this thread finishes last*/
const std::function<void()> & task, 
const std::function<void()> & endtask);


/* btasks for `barrirered tasks'
 * used in situation when several threads of equal status perform
 for (# tasks) { do task | do end task if last to complete | }
 */
inline void btasks(
const size_t & ti,
const size_t & nthreads,
std::vector<size_t> & section_completions,
std::vector<std::mutex> & sectionend_mutexes,
std::vector<std::condition_variable> & section_condvars, 
const std::vector<std::function<void(size_t)>> & section_tasks,
const std::vector<std::function<void()>> &  sectionend_tasks
);

/* btask_rbasks for `barriered task then repeat barriered tasks' 
 * used in situation when 
 * several threads of equal status perform
 * do initialisation task | do intitialisation end task if last to complete | 
 * while (condition is true) for (# tasks) do task | do end task if last to complete |
 */
void btask_rbtasks(
size_t ti, 
size_t nthreads,
std::vector<size_t> & section_completions, 
std::vector<std::mutex> & sectionend_mutexes, 
std::vector<std::condition_variable> & section_condvars, 
const std::function<void(size_t)> & initialisation_task,
const std::function<void()> & initialisationend_task,  
const std::vector<std::function<void(size_t)>> & section_tasks,
const std::vector<std::function<void()>> &  sectionend_tasks, 
const std::function<bool()> & getiscomplete);

/* barriered tasks (inititialisation) then while condition repeat barriered tasks  */
void btasks_rbtasks(
size_t ti, 
size_t nthreads,
std::vector<size_t> & x_completions, 
std::vector<std::mutex> & xend_mutexes, 
std::vector<std::condition_variable> & x_condvars, 
const std::vector< std::function<void(size_t)>> & initialisation_tasks,
const std::vector< std::function<void()>> & initialisationend_tasks,  
const std::vector< std::function<void(size_t)>> & section_tasks,
const std::vector< std::function<void()>> &  sectionend_tasks, 
const std::function<bool()> & getiscomplete);

/* launch barriered tasks */
int launch_btasks(
size_t nthreads, 
const std::vector<std::function<void(size_t)>> & section_tasks,
const std::vector<std::function<void()>> & sectionend_tasks
);


int launch_btasks_rbtasks(
size_t nthreads, 
const std::vector<std::function<void(size_t)>> & initialisation_tasks,
const std::vector<std::function<void()>> & initialisationend_tasks,  
const std::vector<std::function<void(size_t)>> & section_tasks,
const std::vector<std::function<void()>> & sectionend_tasks, 
const std::function<bool()> & getiscomplete,
const std::function<void()> & closing_task 
);



int launch_btask_rbtasks(
size_t nthreads, 
const std::function<void(size_t)> & initialisation_task,
const std::function<void()> & initialisationend_task,  
const std::vector<std::function<void(size_t)>> & section_tasks,
const std::vector<std::function<void()>> & sectionend_tasks, 
const std::function<bool()> & getiscomplete,
const std::function<void()> & closing_task 
);


}


#endif
