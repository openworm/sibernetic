/*******************************************************************************
 * The MIT License (MIT)
 *
 * Copyright (c) 2015 OpenWorm.
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

#ifndef OW_VTKEXPORT_H
#define OW_VTKEXPORT_H

#include "owConfigProperty.h"
#include "owOpenCLConstant.h"
#include "owPhysicsConstant.h"

/*
  Function exportState() serves to export the current state to the VTK file.
  Following is exported:
  - Location of the particles
  - Types of the particles (Liquid, elastic, boundary as defined in owOpenCLConstant.h)
  - Velocity of the particles
  - Connections (CellType == 1)
  - Membranes (CellType == 2)
  - Muscle numbers (Positive for muscles, -1 for elastic connections,
      -2 for membranes)
  - Muscle activation (Positive for muscles, -1 for elastic connections,
      -2 for membranes)

  State in the time step N is saved to the file buffers/state_N.vtp.


  Possible improvements:
  TODO: Write data in binary to save some space

 */


namespace owVtkExport {
	extern bool isActive;

	/** Export the state in one time step to a VTK file
	 *
	 *  @param iteration
	 *  Iteration number
	 *  @param config
	 *  Pointer to the owConfigProperty object it includes information about
	 *  @param position
	 *  Pointer to the position buffer
	 *  @param connections
	 *  Pointer to the connection buffer
	 *  @param velocity
	 *  Pointer to the velocity buffer
	 *  @param membranes
	 *  Pointer to the membranes buffer
	 *  @param muscleActivationSignal
	 *  Pointer to the muscle activation signal buffer
	 */
	void exportState(int iteration, owConfigProperty * config, float * position,
					 float * connections, float * velocity, int * membranes,
					 float * muscleActivationSignal);
}

#endif // #ifndef OW_VTKEXPORT_H
