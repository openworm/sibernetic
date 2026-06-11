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

#include <algorithm>
#include <cstring>
#include <regex>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>

#if defined(__APPLE__) || defined(__MACOSX)
#include <mach/mach.h>
#include <mach/mach_time.h>
#endif

#include "owHelper.h"
#include "owPhysicsConstant.h"

// .npy v1.0 binary format constants for position_buffer.npy / position_buffer_meta.npy
// Header layout: 8-byte magic+version, 2-byte HEADER_LEN, 118-byte dict string = 128 bytes total.
// The row count in position_buffer.npy is a 10-char space-padded field at byte offset 61.
static const int NPY_HDR_LEN    = 118;
static const int NPY_HDR_TOTAL  = 128;
static const int NPY_NROW_OFFSET = 61;

/** owHelpre class constructor
 */
owHelper::owHelper(void) { refreshTime(); }
/** owHelpre class destructor
 */
owHelper::~owHelper(void) {}
/** Update time variable
 *
 *  Using precision system functions for getting exact system time.
 *  Initialization of start time value.
 *  For more info:
 *  Windows QueryPerformanceFrequency -
 * http://msdn.microsoft.com/ru-ru/library/windows/desktop/ms644905(v=vs.85).aspx
 *  Linux clock_gettime -
 * http://man7.org/linux/man-pages/man2/clock_gettime.2.html MacOS -
 * https://developer.apple.com/library/mac/documentation/Darwin/Conceptual/KernelProgramming/services/services.html
 */
void owHelper::refreshTime() {
#if defined(_WIN32) || defined(_WIN64)
  QueryPerformanceFrequency(&frequency); // get ticks per second
  QueryPerformanceCounter(&t1);          // start timer
  t0 = t1;
#elif defined(__linux__)
  clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
  t0 = t1;
#elif defined(__APPLE__)
  t1 = mach_absolute_time();
  t0 = t1;
#endif
}

/** Read physical parametr ftrom string
 *
 */
std::pair<std::string, float> readPhysParam(std::string s) {
  std::regex pattern("^\\s*(\\w+)\\s*:\\s*(\\d+(\\.\\d*([eE]?[+-]?\\d+)?)?)"
                     "\\s*(//.*)?$");
  std::smatch matches;
  if (std::regex_search(s, matches, pattern)) {
    return std::make_pair(matches[1].str(), std::stof(matches[2].str()));
  }
  return std::make_pair("", float());
}

int read_position = 0;
/** TODO make description
 *
 */
void findValf(std::string &str, char delimiter, size_t &start, float &val) {
  size_t end = str.find(delimiter, start + 1);
  if (end != std::string::npos) {
    val = std::atof(str.substr(start, end - start)
                        .c_str()); // TODO make check if this a number
    start = end;
  } else // it means usualy that we reached the end of string and there are no
         // \t
    val = std::atof(
        str.substr(start).c_str()); // TODO make check if this a number
}
/** TODO Documentation
 */
void findVali(std::string &str, char delimiter, size_t &start, int &val) {
  size_t end = str.find(delimiter, start + 1);
  if (end != std::string::npos) {
    val = std::atoi(str.substr(start, end - start).c_str());
    start = end;
  } else // it means usualy that we reached the end of string and there are no
         // \t
    val = std::atoi(str.substr(start).c_str());
}
/** Preparing initial data before load full configuration
 *
 *  Before load configuration data from file (initial position and velocity,
 *  connection data if is and membrane data if is) Sibernetic
 *  should allocate a memory in RAM. Method starts with reading position file
 *  first 6 lines in file correspond to dimensions of boundary box, than it
 * reads all file till the end an calculate numbers of elastic, liquid and
 * boundary particles and total number too. Also this read membranes file for
 * counting membranes numbers.
 *
 *  @param config
 *  pointer to owConfigProperty object it includes information about
 *  boundary box dimensions
 */
void owHelper::preLoadConfiguration(owConfigProperty *config) {
  int p_count = 0;
  std::string file_name = config->getConfigPath() + config->getConfigFileName();
  std::string inputStr;
  std::ifstream configFile(file_name.c_str(), std::ios_base::binary);
  float x, y, z, p_type;
  char delimiter = '\t';
  size_t pos;
  LOADMODE mode = NOMODE;
  if (configFile.is_open()) {
    while (configFile.good()) {
      std::getline(configFile, inputStr);
      inputStr.erase(std::remove(inputStr.begin(), inputStr.end(), '\r'),
                     inputStr.end());
      inputStr.erase(std::remove(inputStr.begin(), inputStr.end(), '\n'),
                     inputStr.end());
      if (inputStr.compare("[simulation box]") == 0) {
        configFile >> config->xmin;
        configFile >> config->xmax;
        configFile >> config->ymin;
        configFile >> config->ymax;
        configFile >> config->zmin;
        configFile >> config->zmax;
        read_position = configFile.tellg();
        mode = NOMODE;
        continue;
      }
      if (inputStr.compare("[physical parameters]") == 0) {
        mode = PHYS_CONST_READ_MODE;
        continue;
      }
      if (inputStr == "[position]") {
        mode = POSITION;
        continue;
      }
      if (inputStr == "[velocity]") {
        mode = VELOCITY;
        continue;
      }
      if (inputStr == "[membranes]") {
        mode = MEMBRANE;
        continue;
      }
      if (inputStr == "[particleMemIndex]") {
        break;
      }
      switch (mode) {
      case POSITION: {
        p_type = -1.1f; // reinitialize
        pos = 0;
        findValf(inputStr, delimiter, pos, x);
        findValf(inputStr, delimiter, pos, y);
        findValf(inputStr, delimiter, pos, z);
        findValf(inputStr, delimiter, pos, p_type);
        p_count++;
        switch ((int)p_type) {
        case LIQUID_PARTICLE:
          config->numOfLiquidP++;
          break;
        case ELASTIC_PARTICLE:
          config->numOfElasticP++;
          break;
        case BOUNDARY_PARTICLE:
          config->numOfBoundaryP++;
          break;
        }
        break;
      }
      case MEMBRANE: {
        config->numOfMembranes++;
        break;
      }
      case PHYS_CONST_READ_MODE: {
        auto p = readPhysParam(inputStr);
        if (p.first.compare("")) {
          config->setConst(p);
        }
        break;
      }
      default:
        continue;
      }
    }
  } else
    throw std::runtime_error("Could not open file configuration file: " +
                             file_name);
  configFile.close();
  config->setParticleCount(p_count);
}
/** Load full configuration
 *
 *  Load configuration data from file (initial position and velocity,
 *  connection data if is and membrane data if is). Method starts with
 *  reading position file than velocity elastic connection and membranes data
 * files
 *
 *  @param position_cpp
 *  pointer to position_cpp buffer
 *  @param velocity_cpp
 *  pointer to velocity_cpp buffer
 *  @param elasticConnections
 *  reference on pointer to elasticConnections buffer.
 *  In this function we allocate memory for elasticConnections.
 *  @param membraneData_cpp
 *  pointer to membraneData_cpp buffer
 *  @param particleMembranesList_cpp
 *  pointer to particleMembranesList_cpp buffer
 *  @param config
 *  pointer to owConfigProperty object it includes information about
 */
void owHelper::loadConfiguration(float *position_cpp, float *velocity_cpp,
                                 float *&elasticConnections,
                                 int *membraneData_cpp,
                                 int *&particleMembranesList_cpp,
                                 owConfigProperty *config) {
  std::string file_name = config->getConfigPath() + config->getConfigFileName();
  std::string inputStr;
  std::ifstream configFile(file_name.c_str(), std::ios_base::binary);
  char delimiter = '\t';
  size_t pos;
  LOADMODE mode = NOMODE;
  int i = 0;
  if (configFile.is_open()) {
    configFile.seekg(read_position);
    while (configFile.good()) {
      // configFile >> inputStr;
      std::getline(configFile, inputStr);
      inputStr.erase(std::remove(inputStr.begin(), inputStr.end(), '\r'),
                     inputStr.end());
      inputStr.erase(std::remove(inputStr.begin(), inputStr.end(), '\n'),
                     inputStr.end());
      if (inputStr == "[position]") {
        mode = POSITION;
        continue;
      }
      if (inputStr == "[velocity]") {
        i = 0;
        mode = VELOCITY;
        continue;
      }
      if (inputStr == "[connection]") {
        i = 0;
        mode = CONNECTION;
        continue;
      }
      if (inputStr == "[membranes]") {
        i = 0;
        mode = MEMBRANE;
        continue;
      }
      if (inputStr == "[particleMemIndex]") {
        i = 0;
        mode = PMEMINDEX;
        continue;
      }
      if (inputStr == "[end]")
        break;
      switch (mode) {
      case POSITION: {
        float x, y, z, p_type;
        p_type = -1.1f; // reinitialize
        pos = 0;
        findValf(inputStr, delimiter, pos, x);
        findValf(inputStr, delimiter, pos, y);
        findValf(inputStr, delimiter, pos, z);
        findValf(inputStr, delimiter, pos, p_type);
        position_cpp[4 * i + 0] = x;
        position_cpp[4 * i + 1] = y;
        position_cpp[4 * i + 2] = z;
        position_cpp[4 * i + 3] = p_type;
        i++;
        break;
      }
      case VELOCITY: {
        float x, y, z, p_type;
        p_type = -1.1f; // reinitialize
        pos = 0;
        findValf(inputStr, delimiter, pos, x);
        findValf(inputStr, delimiter, pos, y);
        findValf(inputStr, delimiter, pos, z);
        findValf(inputStr, delimiter, pos, p_type);
        velocity_cpp[4 * i + 0] = x;
        velocity_cpp[4 * i + 1] = y;
        velocity_cpp[4 * i + 2] = z;
        velocity_cpp[4 * i + 3] = p_type;
        i++;
        break;
      }
      case CONNECTION: {
        pos = 0;
        float jd, rij0, val1, val2;
        findValf(inputStr, delimiter, pos, jd);
        findValf(inputStr, delimiter, pos, rij0);
        findValf(inputStr, delimiter, pos, val1);
        findValf(inputStr, delimiter, pos, val2);
        elasticConnections[4 * i + 0] = jd;
        elasticConnections[4 * i + 1] =
            rij0 * config->getConst("simulationScale");
        elasticConnections[4 * i + 2] = val1;
        elasticConnections[4 * i + 3] = val2;
        i++;
        break;
      }
      case MEMBRANE: {
        pos = 0;
        int id, jd, kd;
        findVali(inputStr, delimiter, pos, id);
        findVali(inputStr, delimiter, pos, jd);
        findVali(inputStr, delimiter, pos, kd);
        membraneData_cpp[3 * i + 0] = id;
        membraneData_cpp[3 * i + 1] = jd;
        membraneData_cpp[3 * i + 2] = kd;
        i++;
        break;
      }
      case PMEMINDEX: {
        int id;
        pos = 0;
        findVali(inputStr, delimiter, pos, id);
        particleMembranesList_cpp[i] = id;
        i++;
        break;
      }
      default:
        continue;
      }
    }
  } else
    throw std::runtime_error("Could not open file configuration file: " + file_name);
  configFile.close();
  std::cout << "Configuration has been loaded" << std::endl;
}
/** Load configuration from simulation to files
 *
 *  This method is required for work with "load config to file" mode.
 *  In this mode information about simulation's evolution is storing into a file
 *  on every step (every time it reads data block with size = c).
 *  If Sibernetic runs in this mode it means that
 *  no calculation on OpenCL device runs.
 *
 *  @param position
 *  pointer to position buffer
 *  @param config
 *  pointer to owConfigProperty object it includes information about
 *  @param connections
 *  reference on pointer to elasticConnections buffer.
 *  @param membranes
 *  pointer to membranes buffer
 *  @param firstIteration
 *  if true it means that we first time record information
 *  to a file and on first iteration it put to
 *  the file info about dimensions of boundary box
 *  NOTE: next 2 parameters are an experimental
 *  @param filter_p
 *  pointer to filter particle buffer, if you need storing only
 *  a bunch of particles not all of them
 *  @param size
 *  size of filter_p array
 */
// Write a 128-byte .npy v1.0 header to f.
// shape_str is the Python shape literal, e.g. "(11,)" or "(%10ld, 4)"-formatted.
static void write_npy_header(FILE *f, const char *shape_str) {
  const unsigned char magic[8] = {0x93,'N','U','M','P','Y',0x01,0x00};
  fwrite(magic, 1, 8, f);
  unsigned short hlen = (unsigned short)NPY_HDR_LEN;
  fwrite(&hlen, 2, 1, f);
  char hdr[119] = {};
  int n = snprintf(hdr, sizeof(hdr),
      "{'descr': '<f4', 'fortran_order': False, 'shape': %s, }", shape_str);
  memset(hdr + n, ' ', 117 - n);
  hdr[117] = '\n';
  fwrite(hdr, 1, NPY_HDR_LEN, f);
}

void owHelper::loadConfigurationToFile(float *position,
                                       owConfigProperty *config,
                                       float *connections, int *membranes,
                                       bool firstIteration, int *filter_p,
                                       int size) {
  std::string basePath = config->getLoadPath();

  if (firstIteration) {
    // position_buffer_meta.npy — shape (11,) float32
    // [xmin xmax ymin ymax zmin zmax numElastic numLiquid numBoundary timeStep logStep]
    std::string metaFileName = basePath + "/position_buffer_meta.npy";
    FILE *metaFile = fopen(metaFileName.c_str(), "wb");
    if (!metaFile)
      throw std::runtime_error("Cannot create position_buffer_meta.npy");
    write_npy_header(metaFile, "(11,)");
    float meta[11] = {
        config->xmin,                   config->xmax,
        config->ymin,                   config->ymax,
        config->zmin,                   config->zmax,
        (float)config->numOfElasticP,   (float)config->numOfLiquidP,
        (float)config->numOfBoundaryP,  config->getTimeStep(),
        (float)config->getLogStep()
    };
    fwrite(meta, sizeof(float), 11, metaFile);
    fclose(metaFile);

    // position_buffer.npy — shape (N, 4) float32, row count updated in-place each frame.
    // Row count stored as space-padded 10-char field at byte offset NPY_NROW_OFFSET.
    std::string posFileName = basePath + "/position_buffer.npy";
    FILE *posFile = fopen(posFileName.c_str(), "wb");
    if (!posFile)
      throw std::runtime_error("Cannot create position_buffer.npy");
    char shape[32];
    snprintf(shape, sizeof(shape), "(%10ld, 4)", 0L);
    write_npy_header(posFile, shape);
    fclose(posFile);
  }

  // Open for read+write, seek to end, append rows, then update row count.
  std::string posFileName = basePath + "/position_buffer.npy";
  FILE *posFile = fopen(posFileName.c_str(), "r+b");
  if (!posFile)
    throw std::runtime_error("Cannot open position_buffer.npy for update");

  char rowBuf[11] = {0};
  fseek(posFile, NPY_NROW_OFFSET, SEEK_SET);
  fread(rowBuf, 1, 10, posFile);
  long currentRows = atol(rowBuf);

  long newRows = 0;
  fseek(posFile, 0, SEEK_END);
  if (size == 0) {
    for (int i = 0; i < config->getParticleCount(); i++) {
      if ((int)position[4 * i + 3] != BOUNDARY_PARTICLE || firstIteration) {
        fwrite(&position[i * 4], sizeof(float), 4, posFile);
        newRows++;
      }
    }
  } else {
    for (int index = 0; index < size; index++) {
      int i = filter_p[index];
      fwrite(&position[i * 4], sizeof(float), 4, posFile);
    }
    newRows = size;
  }

  char newRowBuf[11];
  snprintf(newRowBuf, sizeof(newRowBuf), "%10ld", currentRows + newRows);
  fseek(posFile, NPY_NROW_OFFSET, SEEK_SET);
  fwrite(newRowBuf, 1, 10, posFile);
  fclose(posFile);

  if (firstIteration) {
    std::string connectionFileName =
        config->getLoadPath() + std::string("/connection_buffer.txt");
    std::ofstream connectionFile(connectionFileName.c_str(),
                                 std::ofstream::trunc);
    if (!connectionFile)
      throw std::runtime_error("There was a problem with creation of "
                               "connection data file for logging. Check the "
                               "path.");
    int con_num = MAX_NEIGHBOR_COUNT * config->numOfElasticP;
    for (int i = 0; i < con_num; ++i)
      connectionFile << connections[4 * i + 0] << "\t" << connections[4 * i + 1]
                     << "\t" << connections[4 * i + 2] << "\t"
                     << connections[4 * i + 3] << "\n";
    connectionFile.close();
    std::string membraneFileName =
        config->getLoadPath() + std::string("/membranes_buffer.txt");
    std::ofstream membranesFile(membraneFileName.c_str(), std::ofstream::trunc);
    if (!membranesFile)
      throw std::runtime_error("There was a problem with creation of membrane "
                               "data file for logging. Check the path.");
    membranesFile << config->numOfMembranes << "\n";
    for (unsigned int i = 0; i < config->numOfMembranes; i++)
      membranesFile << membranes[3 * i + 0] << "\t" << membranes[3 * i + 1]
                    << "\t" << membranes[3 * i + 2] << "\n";
    membranesFile.close();
  }
}

void owHelper::loadPressureToFile(float *pressure_buffer,
                                  std::vector<size_t> & shell_particles,
                                  float *position_buffer,
                                  int iteration, owConfigProperty * config)
{
  std::ofstream pressureFile;
  std::string pressureFileName =
      config->getLoadPath() + std::string("/pressure_buffer.txt");
  if (!iteration) {
    pressureFile.open(pressureFileName.c_str(), std::ofstream::trunc);
    if (!pressureFile)
      throw std::runtime_error("There was a problem with creation of position "
                               "file for logging check the path.");
  } else {
    pressureFile.open(pressureFileName.c_str(), std::ofstream::app);
    if (!pressureFile)
      throw std::runtime_error("There was a problem with creation of position "
                               "file for logging Check the path.");
  }
  pressureFile << "[Iteration " << iteration << "]\n";
  for (int i = 0; (unsigned)i < shell_particles.size(); ++i) {
    int id = shell_particles[i];
    pressureFile << "Particle:\t" << id << "\n";
    pressureFile << "\tPosition:\t";
    pressureFile << position_buffer[id * 4 + 0] << "\t"
                 << position_buffer[id * 4 + 1] << "\t"
                 << position_buffer[id * 4 + 2] << "\t"
                 << position_buffer[id * 4 + 3] << "\n";
    pressureFile << "\tPressure:\t" << pressure_buffer[id] << "\n";
  }
  pressureFile.close();
}
/** Load configuration from simulation to files
 *
 *  Make configuration file
 *  @param position
 *  pointer to position buffer
 *  @param config
 *  pointer to owConfigProperty object it includes information about
 *  @param connections
 *  reference on pointer to elasticConnections buffer.
 *  @param membranes
 *  pointer to membranes buffer
 */
void owHelper::loadConfigurationToFile(float *position, float *velocity,
                                       float *connections, int *membranes,
                                       int *particleMemIndex,
                                       const char *filename,
                                       owConfigProperty *config) {
  std::ofstream configFile(filename, std::ofstream::trunc);
  if (!configFile)
    throw std::runtime_error("There was a problem with creation of file for "
                             "saving configuration. Check the path.");
  configFile << "[simulation box]"
             << "\n";
  configFile << config->xmin << "\n";
  configFile << config->xmax << "\n";
  configFile << config->ymin << "\n";
  configFile << config->ymax << "\n";
  configFile << config->zmin << "\n";
  configFile << config->zmax << "\n";
  configFile << "[position]\n";
  for (int i = 0; i < config->getParticleCount(); i++)
    configFile << position[i * 4 + 0] << "\t" << position[i * 4 + 1] << "\t"
               << position[i * 4 + 2] << "\t" << position[i * 4 + 3] << "\n";
  configFile << "[velocity]\n";
  for (int i = 0; i < config->getParticleCount(); i++)
    configFile << velocity[i * 4 + 0] << "\t" << velocity[i * 4 + 1] << "\t"
               << velocity[i * 4 + 2] << "\t" << velocity[i * 4 + 3] << "\n";
  configFile << "[connection]\n";
  int con_num = MAX_NEIGHBOR_COUNT * config->numOfElasticP;
  for (int i = 0; i < con_num; i++)
    configFile << connections[4 * i + 0] << "\t"
               << connections[4 * i + 1] / simulationScale << "\t"
               << connections[4 * i + 2] << "\t" << connections[4 * i + 3]
               << "\n";
  configFile << "[membranes]\n";
  for (unsigned int i = 0; i < config->numOfMembranes; ++i)
    configFile << membranes[3 * i + 0] << "\t" << membranes[3 * i + 1] << "\t"
               << membranes[3 * i + 2] << "\n";
  configFile << "[particleMemIndex]\n";
  int particleMemIndexCount =
      config->numOfElasticP * MAX_MEMBRANES_INCLUDING_SAME_PARTICLE;
  for (int i = 0; i < particleMemIndexCount; ++i)
    configFile << particleMemIndex[i] << "\n";
  configFile << "[end]";
  configFile.close();
}
// This function needed for visualiazation buffered data
long position_index = 0;
FILE *positionFile_bin = nullptr;

/** Load configuration from file to simulation
 *
 *  This method is required for work with "load config from file" mode.
 *  In this mode information about simulation's evolution is being taken from
 * file on every step (every time it reads data block with size =
 * PARTICLE_COUNT). If Sibernetic runs in this mode it means that no
 * calculation on OpenCL device runs.
 *
 *  @param position
 *  pointer to position buffer
 *  @param connections
 *  reference on pointer to elasticConnections buffer.
 *  @param membranes
 *  pointer to membranes buffer
 *  @param config
 *  pointer to owConfigProperty object it includes information about
 *  @param iteration
 *  if iteration==0 it means that we first time record information
 *  to a file and on first iteration it put to
 *  the file info about dimensions of boundary box
 */
bool owHelper::loadConfigurationFromFile(float *&position, float *&connections,
                                         int *&membranes,
                                         owConfigProperty *config,
                                         int iteration) {
  if (iteration == 0) {
    // Read metadata from position_buffer_meta.npy
    std::string metaFileName =
        config->getLoadPath() + std::string("/position_buffer_meta.npy");
    FILE *metaFile = fopen(metaFileName.c_str(), "rb");
    if (!metaFile)
      throw std::runtime_error("Cannot open position_buffer_meta.npy");
    fseek(metaFile, NPY_HDR_TOTAL, SEEK_SET);
    float meta[11];
    fread(meta, sizeof(float), 11, metaFile);
    fclose(metaFile);
    config->xmin           = meta[0];  config->xmax           = meta[1];
    config->ymin           = meta[2];  config->ymax           = meta[3];
    config->zmin           = meta[4];  config->zmax           = meta[5];
    config->numOfElasticP  = (int)meta[6];
    config->numOfLiquidP   = (int)meta[7];
    config->numOfBoundaryP = (int)meta[8];
    config->setTimeStep(meta[9]);
    config->setLogStep((int)meta[10]);
    config->setParticleCount(config->numOfElasticP + config->numOfLiquidP +
                             config->numOfBoundaryP);
    position = new float[4 * config->getParticleCount()];

    std::string posFileName =
        config->getLoadPath() + std::string("/position_buffer.npy");
    positionFile_bin = fopen(posFileName.c_str(), "rb");
    if (!positionFile_bin)
      throw std::runtime_error("Cannot open position_buffer.npy");
    fseek(positionFile_bin, NPY_HDR_TOTAL, SEEK_SET);
  }

  if (!positionFile_bin)
    return false;

  // Frame 0 includes boundary particles; subsequent frames do not.
  int N = (iteration == 0)
      ? config->getParticleCount()
      : (config->numOfElasticP + config->numOfLiquidP);

  size_t nread = fread(position, sizeof(float), 4 * N, positionFile_bin);
  if (nread < (size_t)(4 * N)) {
    fclose(positionFile_bin);
    positionFile_bin = nullptr;
    return false;
  }

  if (iteration == 0) {
    config->setParticleCount(config->numOfElasticP + config->numOfLiquidP);
    unsigned int i = 0;
    std::string connectionFileName =
        config->getLoadPath() + std::string("/connection_buffer.txt");
    std::ifstream connectionFile(connectionFileName.c_str());
    connections = new float[MAX_NEIGHBOR_COUNT * config->numOfElasticP * 4];
    if (connectionFile.is_open()) {
      float jd, rij0, val1, val2;
      while (connectionFile.good() &&
             i < MAX_NEIGHBOR_COUNT * (unsigned)config->numOfElasticP) {
        connectionFile >> jd >> rij0 >> val1 >> val2;
        connections[4 * i + 0] = jd;
        connections[4 * i + 1] = rij0;
        connections[4 * i + 2] = val1;
        connections[4 * i + 3] = val2;
        i++;
      }
    }
    connectionFile.close();
    std::string membraneFileName =
        config->getLoadPath() + std::string("/membranes_buffer.txt");
    std::ifstream membranesFile(membraneFileName.c_str());
    if (membranesFile.is_open()) {
      int m_count = 0;
      membranesFile >> m_count;
      int mi = 0;
      membranes = new int[3 * m_count];
      /* TODO: skip first two membranes? */
      while (membranesFile.good() && mi < m_count) {
        membranesFile >> membranes[3 * mi + 0] >> membranes[3 * mi + 1] >>
            membranes[3 * mi + 2];
        mi++;
      }
    }
    membranesFile.close();
  }
  return true;
}
/** Print value of elapsed time from last handling to watch_report method.
 *
 *  This function is required for logging time consumption info.
 *
 *  @param str
 *  represents output string format.
 */
void owHelper::watch_report(const char *str) {
#if defined(_WIN32) || defined(_WIN64)
  QueryPerformanceCounter(&t2);
  printf(str, (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart);
  t1 = t2;
  elapsedTime = (t2.QuadPart - t0.QuadPart) * 1000.0 / frequency.QuadPart;
#elif defined(__linux__)
  clock_gettime(CLOCK_MONOTONIC_RAW, &t2);
  time_t sec = t2.tv_sec - t1.tv_sec;
  long nsec;
  if (t2.tv_nsec >= t1.tv_nsec) {
    nsec = t2.tv_nsec - t1.tv_nsec;
  } else {
    nsec = 1000000000 - (t1.tv_nsec - t2.tv_nsec);
    sec -= 1;
  }
  printf(str, (float)sec * 1000.f + (float)nsec / 1000000.f);
  t1 = t2;
  elapsedTime = (float)(t2.tv_sec - t0.tv_sec) * 1000.f +
                (float)(t2.tv_nsec - t0.tv_nsec) / 1000000.f;
#elif defined(__APPLE__)
  uint64_t elapsedNano;
  static mach_timebase_info_data_t sTimebaseInfo;

  if (sTimebaseInfo.denom == 0) {
    (void)mach_timebase_info(&sTimebaseInfo);
  }

  t2 = mach_absolute_time();
  elapsedNano = (t2 - t1) * sTimebaseInfo.numer / sTimebaseInfo.denom;
  printf(str, (float)elapsedNano / 1000000.f);
  t1 = t2;
  elapsedNano = (t2 - t0) * sTimebaseInfo.numer / sTimebaseInfo.denom;
  elapsedTime = (float)elapsedNano / 1000000.f;
#endif
}
