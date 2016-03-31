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

#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <iostream>

#include "owVtkExport.h"

namespace owVtkExport {
	bool isActive = false;
	static int numConns = -1;
	static bool isBigEndian;

	void countExistingConnections(float * connections, owConfigProperty * config) {
		numConns = 0;
		for (unsigned int i = 0; i < MAX_NEIGHBOR_COUNT * config->numOfElasticP; ++i) {
			if (connections[4 * i + 0] != -1) {
				numConns += 1;
			}
		}
	}

	void setEndianness() {
		union {
			uint32_t i;
			char c[4];
		} bint = {0x01020304};

		isBigEndian = (bint.c[0] == 1);
	}

	void printParticles(std::ofstream & outFile, float * position,
						float * velocity, owConfigProperty * config) {

		outFile << "<Points>\n";
		outFile << "<DataArray NumberOfComponents=\"3\" type=\"Float32\""
				<< " format=\"ascii\">";
		for (int i = 0; i < config->getParticleCount(); i++) {
				outFile << position[i * 4 + 0] << " "
						<< position[i * 4 + 1] << " "
						<< position[i * 4 + 2] << " ";
		}
		outFile << "</DataArray>\n";
		outFile << "</Points>\n";

		outFile << "<PointData>\n";
		outFile << "<DataArray Name=\"ParticleType\" type=\"Int32\" format=\"ascii\">";
		for (int i = 0; i < config->getParticleCount(); i++) {
			outFile << (int) position[i * 4 + 3] << " ";
		}
		outFile << "</DataArray>\n";

		outFile << "<DataArray Name=\"Velocity\" type=\"Float32\""
				<< " NumberOfComponents=\"3\" format=\"ascii\">";
		for (int i = 0; i < config->getParticleCount(); i++) {
			outFile << velocity[i * 4 + 0] << " "
					<< velocity[i * 4 + 1] << " "
					<< velocity[i * 4 + 2] << " ";
		}
		outFile << "</DataArray>\n";
		outFile << "</PointData>\n";
	}

	void printMembranes(std::ofstream & outFile, int * membranes,
						owConfigProperty * config) {
		outFile << "<Polys>\n";
		outFile << "<DataArray type=\"Int32\" Name=\"connectivity\""
				<< " format=\"ascii\">";
		for (unsigned int i = 0; i < config->numOfMembranes; ++i) {
			outFile << membranes[i * 3 + 0] << " "
					<< membranes[i * 3 + 1] << " "
					<< membranes[i * 3 + 2] << " ";
		}
		outFile << "</DataArray>\n";
		outFile << "<DataArray type=\"Int32\" Name=\"offsets\""
				<< " format=\"ascii\">";
		for (unsigned int i = 0; i < config->numOfMembranes; ++i) {
			outFile << 3*(i + 1) << " ";
		}
		outFile << "</DataArray>\n";
		outFile << "</Polys>\n";
	}

	void printConnections(std::ofstream & outFile, float * connections,
						  owConfigProperty * config) {
		outFile << "<Lines>\n";
		outFile << "<DataArray type=\"Int32\" Name=\"connectivity\""
				<< " format=\"ascii\">";
		for (unsigned int i = 0; i < MAX_NEIGHBOR_COUNT * config->numOfElasticP; ++i) {
			if (connections[4 * i + 0] != -1) {
				outFile << i / MAX_NEIGHBOR_COUNT << " "
						<< (int) connections[4 * i + 0] << " ";
			}
		}
		outFile << "</DataArray>\n";
		outFile << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
		for (int i = 0; i < numConns; i++) {
			outFile << 2*(i + 1) << " ";
		}
		outFile << "</DataArray>\n";
		outFile << "</Lines>\n";
	}


	void printCellTypes(std::ofstream & outFile, owConfigProperty * config) {
		outFile << "<DataArray type=\"Float32\" Name=\"CellType\" format=\"ascii\">";
		for (int i = 0; i < numConns; i++) {
			outFile << "1.0 ";
		}
		for (unsigned int i = 0; i < config->numOfMembranes; ++i) {
			outFile << "2.0 ";
		}
		outFile << "</DataArray>";
	}

	void printMuscleNumbers(std::ofstream & outFile, float * connections,
							owConfigProperty * config) {
		outFile << "<DataArray type=\"Float32\" Name=\"MuscleNumber\""
				<< " format=\"ascii\">";
		for (unsigned int i = 0; i < MAX_NEIGHBOR_COUNT * config->numOfElasticP; ++i) {
			if (((int) connections[4 * i + 0]) >= 0) {
				int muscleId = floor(connections[4 * i + 2]) - 1;
				if (muscleId >= 0) {
					outFile << connections[4*i + 2] - floor(connections[4*i + 2])
							<< " ";
				}
				else {
					outFile << "-1.0 ";
				}
			}
		}
		for (unsigned int i = 0; i < config->numOfMembranes; ++i) {
			outFile << "-2.0 ";
		}
		outFile << "</DataArray>";
	}

	void printMuscleActivation(std::ofstream & outFile, float * connections,
							   float * muscleActivationSignal,
							   owConfigProperty * config) {
		outFile << "<DataArray type=\"Float32\" Name=\"MuscleActivation\""
				<< " format=\"ascii\">";
		for (unsigned int i = 0; i < MAX_NEIGHBOR_COUNT * config->numOfElasticP; ++i) {
			if (((int) connections[4 * i + 0]) >= 0) {
				int muscleId = floor(connections[4 * i + 2]) - 1;
				if (muscleId >= 0) {
					outFile << muscleActivationSignal[muscleId] << " ";
				}
				else {
					outFile << "-1.0 ";
				}
			}
		}
		for (unsigned int i = 0; i < config->numOfMembranes; ++i) {
			outFile << "-2.0 ";
		}
		outFile << "</DataArray>";
	}


	template <typename T> std::string to_string(T val, int width, char fill) {
		std::stringstream stream;
		stream << std::setw(width) << std::setfill(fill) << val;
		return stream.str();
	}

	void exportState(int iteration, owConfigProperty * config, float * position,
					 float * connections, float * velocity, int * membranes,
					 float * muscleActivationSignal) {

		std::string filename = config->getLoadPath() + std::string("state_")
			+ to_string(iteration, 8, '0') + std::string(".vtp");

		std::ofstream outFile(filename.c_str());
		if (!outFile) {
			throw std::runtime_error("Cannot create VTK file.");
		}

		if (numConns == -1) {
			/* Initialization of variables */
			countExistingConnections(connections, config);
			setEndianness();
		}

		outFile << "<VTKFile type=\"PolyData\" version=\"0.1\""
				<< " byte_order="
				<< (isBigEndian ? "\"BigEndian\"" : "\"LittleEndian\"")
				<< ">\n";
		outFile << "<PolyData>\n";
		outFile << "<Piece"
				<< " NumberOfPoints=\"" << config->getParticleCount() << "\""
				<< " NumberOfPolys=\""	<< config->numOfMembranes << "\""
				<< " NumberOfLines=\""	<< numConns << "\""
				<< ">\n";

		printParticles(outFile, position, velocity, config);
		printMembranes(outFile, membranes, config);
		printConnections(outFile, connections, config);

		outFile << "<CellData>\n";
		printCellTypes(outFile, config);
		printMuscleNumbers(outFile, connections, config);
		printMuscleActivation(outFile, connections, muscleActivationSignal, config);
		outFile << "</CellData>\n";

		outFile << "</Piece>\n";
		outFile << "</PolyData>\n";
		outFile << "</VTKFile>\n";

		outFile.close();
	}
}
