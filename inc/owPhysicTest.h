/*
 * owPhysicTest.h
 *
 *  Created on: Apr 30, 2014
 *      Author: serg
 */

#ifndef OWPHYSICTEST_H_
#define OWPHYSICTEST_H_

#include "owPhysicsFluidSimulator.h"

// Run the energy conservation test for the specified number of iterations.
// If iterations <= 0 the default value inside the implementation is used.
void test_energy_conservation(int argc, char **argv, int iterations = 0);


#endif /* OWPHYSICTEST_H_ */
