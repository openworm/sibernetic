/*******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) 2011, 2013 OpenWorm.
 * http://openworm.org
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the MIT License
 * which accompanies this distribution, and is available at
 * http://opensource.org/licenses/MIT
 *
 * Contributors:
 *     	OpenWorm - http://openworm.org/people.html
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 *******************************************************************************/

#ifndef OW_PHYSICS_CONSTANT_H
#define OW_PHYSICS_CONSTANT_H

#include <math.h>

#include "owOpenCLConstant.h"
/** Main physical constants contain here
 */

#ifndef M_PI
#define M_PI 3.1415927f
#endif


const float rho0 = 1000.0f;                         // Standard value of liquid density for water (kg/m^3)

const float mass = 0.54e-13f;                       // Mass for one particle (kg).
                                                    // Some facts about C. elegans:
                                                    // Adult worm mass = 3.25e-06 grams = 3.25e-09 kg
                                                    // worm density is around 1000 kg/m3
                                                    // Adult worm length =  1 mm =   1000 um =    1e-03 m
                                                    // Adult worm broad diameter = 60..80 um = 6..8e-05 m // we'll consider it to be equal 80 um (radius = 40 um)
                                                    // Adult worm volume = 0.0033 mm3
                                                    // 1000*40*40*Cw = 0.0033
                                                    // then Cw = 2.0625
                                                    // so, if we need a worm body model composed of, for example, 1e+05 particles
                                                    // each particle's mass should be 3.25e-09 / 1e+05 = 3.25e-14 kg
                                                    // and length of the worm will be (calculation follows):
                                                    // n - number of particles per 1 um
                                                    // (1000*n)*(40*n)*(40*n)*Cw = 1e+5 particles
                                                    // then n^3 = 0.303, n = 0.311
                                                    // then worm length = (1000*n) = 311 particles, radius = (40*n) = 12 particles
                                                    // So, in this case (1e+5 particles) we need r0 = 3.2 um = 3.2e-6 m
                                                    // and particle mass = 3.25e-14 kg
                                                    // NOTE: we use this value of mass because we are oriented on modeling of
                                                    // C. elegans's body model. But we you can use your own value of mass
                                                    // TODO: make it as an input parameter

const float timeStep = 5.0e-06f;                    // Time step of simulation (s)
                                                    // NOTE: "For numerical stability and convergence, several time step
                                                    // constraints must be satisfied. The Courant-Friedrich-Levy
                                                    // (CFL) condition for SPH (dt <= lambda_v*(h/v_max))
                                                    // states that the speed of numerical propagation must be
                                                    // higher than the speed of physical propagation, where v_max = max(||v_i(t)||)
                                                    // is the maximum magnitude of the velocity throughout the simulation. lambda_v is a constant factor, e. g.
                                                    // lambda_v = 0.4. In other words, a particle i must not move more than its smoothing length h in one time step. Fur-
                                                    // thermore, high accelerations might influence the simulation
                                                    // results negatively. Therefore, the time step must also satisfy
                                                    // dt <= lambda_f*sqrt(h/F_max)
                                                    // where F_max=max(||F_i(t)||) denotes the magnitude of
                                                    // the maximum force per unit mass for all particles throughout
                                                    // the simulation. lambda_f = 0.25."
                                                    // For more info [1, page 5].
                                                    // NOTE: actually it depends on mass too for bigger value of mass it possible
                                                    // to use bigger value of time step. Dependence on mass could be described by following
                                                    // mass influent on simulation scale and due to the fact that we're simulating incompressible
                                                    // liquid start configuration of particles should satisfy condition that density in
                                                    // every particle of configuration <= rh0. So if we decrease mass we should decrease
                                                    // distance among particles than it leads to that we should decrease time step.
                                                    // TODO: find dependence and make choice automatically
                                                    // [1] M. Ihmsen, N. Akinci, M. Gissler, M. Teschner, Boundary Handling and Adaptive Time-stepping for PCISPH Proc. VRIPHYS, Copenhagen, Denmark, pp. 79-88, Nov 11-12, 2010.
                                                    // ATTENTION! too large values can lead to 'explosion' of elastic matter objects
							/*0.00411**/
const float simulationScale = 0.0041f*pow(mass,1.f/3.f)/pow(0.00025f,1.f/3.f);//pow(mass,1.f/3.f)/pow(rho0,1.f/3.f); // Simulation scale coefficient. It means that N * simulationScale
                                                                   // converts from simulation scale to meters N / simulationScale convert from meters simulation scale
                                                                   // If you want to take real value of distance in meters you need multiple on simulation scale
                                                                   // NOTE: simulationScale depends from mass of particle. If we place one particle
                                                                   // into volume with some size of side we want that density in this value is equal to rho0

const float h = 3.34f;                              // Smoothed radius value. This is dimensionless invariant parameter.
                                                    // For taken real value in meter you need multiple this on simulationScale.
                                                    // h is a spatial distance, over which their properties are "smoothed" by a kernel function [1].
                                                    // [1] https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics

const float hashGridCellSize = 2.0f * h;            // All bounding box is divided on small spatial cells with size of side == h. Size of side for one spatial cell
                                                    // This require for spatial hashing and => searching a neighbors
const float r0 = 0.5f * h;                          // Standard distance between two boundary particle == equilibrium distance between 2 particles [1]
                                                    // [1] M. Ihmsen, N. Akinci, M. Gissler, M. Teschner, Boundary Handling and Adaptive Time-stepping for PCISPH Proc. VRIPHYS, Copenhagen, Denmark, pp. 79-88, Nov 11-12, 2010.

const float viscosity = 0.00005f;                   // liquid viscosity value //why this value? Dynamic viscosity of water at 25 C = 0.89e-3 Pa*s
const double beta = timeStep*timeStep*mass*mass*2/(rho0*rho0); // B. Solenthaler's dissertation, formula 3.6 (end of page 30)

const double Wpoly6Coefficient = 315.0 / ( 64.0 * M_PI * pow( (double)(h*simulationScale), 9.0 ) ); // Wpoly6Coefficient for kernel Wpoly6 [1]
                                                                                                    // [1] Solenthaler (Dissertation) page 17 eq. (2.20)

const double gradWspikyCoefficient= -45.0 / ( M_PI * pow( (double)(h*simulationScale), 6.0 ) );     // gradWspikyCoefficient for kernel gradWspiky [1]
                                                                                                    // [1] Solenthaler (Dissertation) page 18 eq. (2.21)

const double divgradWviscosityCoefficient = - gradWspikyCoefficient;                                // divgradWviscosityCoefficient for kernel Viscous [1]
                                                                                                    // [1] Solenthaler (Dissertation) page 18 eq. (2.22)
/* We' re using Cartesian coordinate system
				y|
				 |___x
				 /
			   z/
 */
const float gravity_x = 0.0f;                       // Value of vector Gravity component x
const float gravity_y = -9.8f;                      // Value of vector Gravity component y
const float gravity_z = 0.0f;                       // Value of vector Gravity component z
const int maxIteration = 3;                         // Number of iterations for Predictive-Corrective scheme

const float mass_mult_Wpoly6Coefficient = (float) ( (double)mass * Wpoly6Coefficient );                       // Conversion of double value to float. For work with only 1st precision arithmetic.
const float mass_mult_gradWspikyCoefficient = (float) ( (double)mass * gradWspikyCoefficient );               // It needs for work with devices don't support double precision.
const float mass_mult_divgradWviscosityCoefficient = (float) ( (double)mass * divgradWviscosityCoefficient ); // Also it helps to increase performance and memory consumption

/** Following parameters need for decreasing repeating calculation of this values
 */
const float hashGridCellSizeInv = 1.0f / hashGridCellSize; // Inverted value for hashGridCellSize
const float simulationScaleInv = 1.0f / simulationScale;   // Inverted value for simulationScale
const float _hScaled = h * simulationScale;         // scaled smoothing radius
const float _hScaled2 = _hScaled*_hScaled;          // squared scaled smoothing radius

const float surfTensCoeff = mass_mult_Wpoly6Coefficient * simulationScale;
//const float surfTensCoeff = -1.5e-09f * 0.3f* (float)(Wpoly6Coefficient * pow(h*simulationScale*h*simulationScale/2.0,3.0)) * simulationScale; // Surface coefficient. Actually it is -1.5e-09f * 0.3f
                                                                                                                                               // But for decreasing number of repeating calculation we suppose that
                                  /*5.00e-05*/                                                                                                            // surfTensCoeff = -1.5e-09f * 0.3f* (float)(Wpoly6Coefficient * pow(h*simulationScale*h*simulationScale/2.0,3.0)) * simulationScale
const float elasticityCoefficient = 3.00e-05f / mass; // Elasticity coefficient. Actually it isn't
                                                      // elasticity coefficient (elasticity coefficient = 1.95e-05f)
                                                      // But for decreasing number of repeating calculation we suppose that  elasticityCoefficient = 1.95e-05f / mass
#endif // #ifndef OW_PHYSICS_CONSTANT_H
