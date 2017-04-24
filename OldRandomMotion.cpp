#include <iostream>
#include <fstream> // files, reading and writing
#include <sstream>
#include <iterator>
#include <vector>
#include <stdlib.h> // rand
#define _USE_MATH_DEFINES // for pi
#include <math.h>   // sqrt
#include <stdio.h>
//#include <iomanip> // set precison to print floats
#include <time.h>

//#include </home/alexandre/local/FlowVR/include/flowvr/moduleapi.h>
#include <flowvr/moduleapi.h>
//#include <flowvr/flowvr_vrpn_stamps.h>

// Include GLEW
//#include <GL/glew.h>

#include </people/kouyoumdjian/local/glfw/include/GLFW/glfw3.h>
//#include <GLFW/glfw3.h> // glfwGetTime
//#include <glfw3.h> // glfwGetTime
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

using namespace std;
using namespace flowvr;

/*
const int XDIM = 10;
const int YDIM = 10;
const int ZDIM = 10;
*/

const float XSTART = -2.4f;
const float YSTART = 0.0f;
const float ZSTART = 0.0f;
const float ZSTART3D = -1.35f;

const float XEND = 4.8f;
const float YEND = 2.7f;
const float ZEND = 1.35f;


const float MARGIN = 0.65f; // 0.50f is theoretically enough, but arithmetic precision issues may arise in practice
const float RADIUS = 0.10f; // Sphere radius in renderer; dirty as shit.

const float CAMERA_MARGIN = 2.5f; // Bit of padding between the camera and the simulation space, to limit occlusion to reasonable levels.
//const float SPACING = 0.6f;

/*
const int XNBT = int(XDIM/SPACING) + 1;
const int YNBT = int(YDIM/SPACING) + 1;
const int ZNBT = int(ZDIM/SPACING) + 1;
*/

//const int NB_TARGETS = XNBT * YNBT * ZNBT;
//const int CHOSEN_NB_TARGETS_2D = 450;
int CHOSEN_NB_TARGETS_2D;

//const int XNBT_2D = (int) sqrt(CHOSEN_NB_TARGETS_2D);
//const int YNBT_2D = XNBT_2D;
const float GRID_SIZE = 4.8f;
//const float SPACING = (GRID_SIZE - MARGIN)/XNBT_2D;

// Same thing as XNBT_2D^2 for a square grid. Not necessarily equal to CHOSEN_NB_TARGETS_2D due to rounding errors.
//const int NB_TARGETS_2D = XNBT_2D * YNBT_2D;

// Used by MoveTargets for the version with a 3D grid of targets
const float SPEED           = 0.000006f;
const float AMPLITUDE       = GRID_SIZE/3.0f;
const float SQ_AMPLITUDE    = AMPLITUDE*AMPLITUDE;

const float SPRING_FACTOR   = 0.0002f;
const float SHRINK_FACTOR   = 0.89f; // 12 iterations
const float MAX_EXP_GAP     = 60.0f;

/*
const glm::mat4 MVPC = glm::mat4(
	0.130323,   0,          0,          -6.51613,   // first column
	0,          0.130323,   0,          -6.51613,   // second column
	0,          0,          -0.1002,    9.81892,    // third column
	0,          0,          -0.1,       10
);*/


const glm::mat4 projectionMatrix = glm::mat4(
	1.30323,    0,          0,          0,
	0,          1.30323,    0,          0,
	0,          0,          -1.002,     -1,
	0,          0,          -0.2002,    0
);

const glm::mat4 viewMatrix = glm::mat4(
	1,      0,      -0,     0,
	-0,     1,      -0,     0,
	0,      0,      1,      0,
	-5,     -5,     -10,    1
);

const glm::mat4 frustumMatrix = projectionMatrix * viewMatrix;

const float FUDGED_RADIUS = RADIUS + (RADIUS * 0.275f) ;    // OK, so the xyzRadius vectors below are not orthogonal to the frustum planes...
															// ...so the test doesn't work perfectly and the radius needs to be fudged a little.
															// It's gross, but easy, and it should work well enough.
const glm::vec4 xRadius = glm::vec4(FUDGED_RADIUS, 0, 0, 0);
const glm::vec4 yRadius = glm::vec4(0, FUDGED_RADIUS, 0, 0);
const glm::vec4 zRadius = glm::vec4(0, 0, FUDGED_RADIUS, 0);


const float SCALE_FACTOR = 0.10f; // enlarge to get bigger spheres
const glm::vec3 SCALE_VECTOR = glm::vec3(SCALE_FACTOR, SCALE_FACTOR, SCALE_FACTOR);
const glm::mat4 ModelMatrix = glm::scale(glm::mat4(1.0), SCALE_VECTOR);

// Angle of each rotation, in PIs. I.e. if ANGLE == 0.65f, the angle will be between -0.65*PI and 0.65*PI
//const float ANGLE         = 0.30f;
float ANGLE;

// Number of changes in direction per second
//const float ROTATION_FREQUENCY    = 4.0f;
//float ROTATION_FREQUENCY;
float ROTATION_PERIOD;

// Used by MoveTargets2D for the version with a 2D grid of targets
//const float SPEED_2D          = 0.750f;
float SPEED_2D;

int NB_STATIC_CONDITIONS;

// Minimum and maximum angles for the rotation of directions.
float RANGE_MIN;
float RANGE_MAX;

// Number of trials per block
int NB_TRIALS;

vector<float> v_angles;
vector<float> v_periods;
vector<float> v_speeds;

std::vector<glm::vec3> v_conditions;

MessageWrite msgPosOut;
Message msgSpringIndex;
Message msgParams;

//StampInfo *PosMsgStamp = new flowvr::StampInfo("TypePositions",TypeInt::create());
bool threeDimensional = false;


void CreateConditionSet() {
	int nb_static_conditions = 5;

	for(unsigned int i=0; i<v_angles.size(); ++i)
		for(unsigned int j=0; j<v_periods.size(); ++j)
			for(unsigned int k=0; k<v_speeds.size(); ++k)
				v_conditions.push_back( glm::vec3(v_angles[i], v_periods[j], v_speeds[k]) );

	for(unsigned int i=0; i<v_conditions.size(); ++i)
		cout << "Condition " << i << ": " << v_conditions[i].x << ", " << v_conditions[i].y << ", " << v_conditions[i].z << endl;

	std::random_shuffle(v_conditions.begin(), v_conditions.end());

	cout << endl << endl;

	//for(int i=0; i<nb_static_conditions; ++i)
		v_conditions.insert(v_conditions.begin(), NB_STATIC_CONDITIONS, glm::vec3(1,1,0));

	for(unsigned int i=0; i<v_conditions.size(); ++i)
		cout << "Condition " << i << ": " << v_conditions[i].x << ", " << v_conditions[i].y << ", " << v_conditions[i].z << endl;
}

void ReadNBTargets(){ // And number!
	string line;
	const char* fname;
	if(threeDimensional)
		fname = "../../share/assets/config/parameters3d.txt";
	else
		fname = "../../share/assets/config/parameters.txt";
	ifstream myfile(fname);
	if (myfile.is_open()) {
		getline(myfile, line); // NB Targets
		getline(myfile, line);
		CHOSEN_NB_TARGETS_2D = atoi(line.c_str());
		myfile.close();
	}
	else {
		if(threeDimensional)
			cout << "Unable to open file parameters3d.txt." << endl;
		else
			cout << "Unable to open file parameters.txt." << endl;
		cout << "You might be in the wrong directory. Good luck figuring it out." << endl;
	}
}

float RandomFloat(float a, float b) {
	float random = ((float) rand()) / (float) RAND_MAX;
	float diff = b - a;
	float r = random * diff;
	return a + r;
}

void InitTargets2D(float* targets) {
	float x, y, z;
	//cout << "NB TARGETS:  " << CHOSEN_NB_TARGETS_2D << endl;
	for(int i=0; i<CHOSEN_NB_TARGETS_2D; ++i) {
		x = GRID_SIZE - MARGIN + static_cast <float> (rand()) / (static_cast <float> ( RAND_MAX/(XSTART + MARGIN - GRID_SIZE + MARGIN) ) );
		y = GRID_SIZE - MARGIN + static_cast <float> (rand()) / (static_cast <float> ( RAND_MAX/(YSTART + MARGIN - GRID_SIZE + MARGIN) ) );
		//z = GRID_SIZE - MARGIN + static_cast <float> (rand()) / (static_cast <float> ( RAND_MAX/(ZSTART + MARGIN - GRID_SIZE + MARGIN) ) );
		z = -20.0f;
		targets[3*i]        = x;
		targets[3*i + 1]    = y;
		targets[3*i + 2]    = z;
	}
}



void InitDirections2D(float* directions) {
	float vx, vy;
	for(int i=0; i<CHOSEN_NB_TARGETS_2D; ++i) {
		// vx and vy are random floats between -1.0f and +1.0f
		vx = -1.0f + static_cast <float> (rand()) / (static_cast <float> ( RAND_MAX/(2.0f) ) ); // 2.0f == MAX - MIN == 1.0f - (-1.0f)
		vy = -1.0f + static_cast <float> (rand()) / (static_cast <float> ( RAND_MAX/(2.0f) ) );
		directions[3*i]     = vx;
		directions[3*i + 1] = vy;
		// Nothing for 3*i + 2, since there's no movement on the Z axis
	}
}


void NormalizeDirections2D(float* directions) {
	float x, y, length;
	for(int i=0; i<CHOSEN_NB_TARGETS_2D; i++) {
		x = directions[3*i];
		y = directions[3*i + 1];
		length = sqrt(x*x + y*y);
		if(length != 0.0f) {
			directions[3*i]     = x/length;
			directions[3*i + 1] = y/length;
		} else {
			directions[3*i]     = 1.0f;
			directions[3*i + 1] = 1.0f;
			cerr << "Nil direction found; set to (1.0,1.0)." << endl;
		}
	}
}



void DisplayTargets2D(float* targets) {
	for(int i=0; i<CHOSEN_NB_TARGETS_2D; ++i) {
		cout << "Target " << i      << ": ";
		cout << targets[3*i]        << " | ";
		cout << targets[3*i + 1]    << endl;
		// No need to display z, but it's there.
	}
}

void DisplayDirections2D(float* directions){
	for(int i=0; i<CHOSEN_NB_TARGETS_2D; ++i) {
		cout << "Direction "<< i    << ": ";
		cout << directions[3*i]     << " | ";
		cout << directions[3*i + 1] << endl;
	}
}

void RecSpringToPosition(float* targets, double deltaT, int currentTarget, float expo) {
	if(expo < 1.0f)
		return;
	int nTarget = (currentTarget + 1) % CHOSEN_NB_TARGETS_2D;
	float xdist = targets[3*nTarget]        - targets[3*currentTarget];
	float ydist = targets[3*nTarget + 1]    - targets[3*currentTarget + 1];
	float currSqAmplitude = xdist*xdist + ydist*ydist;

	float gap = SQ_AMPLITUDE - currSqAmplitude;
	float expGap = std::pow(abs(gap), expo);
	if(gap < 0)
		expGap = -abs(expGap);
	else
		expGap = abs(expGap);

	if(abs(expGap) >= MAX_EXP_GAP) {
		//cout << "Big Gap: " << expGap << endl;
		if(expGap < 0)
			expGap = -MAX_EXP_GAP;
		else
			expGap = MAX_EXP_GAP;
	}
	expGap *= 4.0f*SPEED_2D;

	targets[3*nTarget]      += deltaT * SPRING_FACTOR * expGap * xdist;
	targets[3*nTarget + 1]  += deltaT * SPRING_FACTOR * expGap * ydist;

	// Next next target
	RecSpringToPosition(targets, deltaT, nTarget, expo*SHRINK_FACTOR);
}

void LogTimedConditions(ofstream* tclog, float angle, float period, float speed) {
	*tclog << glfwGetTime() << "    New Condition   " << angle*180.0f << "  " << 1.0f/period << "   " << speed << endl;
}

void LogTimedSelections(ofstream* tclog) {
	*tclog << glfwGetTime() << "    New Selection" << endl;
}

// Depends on direction vectors
void MoveTargets2D(float* targets, double deltaT, float* directions, int currentTarget) {
	float x, y;
	float bump = 0.001f; // To make sure targets stay in the window and don't just skirt long the edge, right outside.
	for(int i=0; i<CHOSEN_NB_TARGETS_2D; ++i) {
		targets[3*i]        += SPEED_2D * deltaT * directions[3*i];         // X
		targets[3*i + 1]    += SPEED_2D * deltaT * directions[3*i + 1];     // Y
		// z does not need to move, but it is still there

		RecSpringToPosition(targets, deltaT, currentTarget, 4.0f);

		// Bouncing is better
		bool xBelow = targets[3*i]      <= XSTART + RADIUS;
		bool xAbove = targets[3*i]      >= XSTART + GRID_SIZE - RADIUS;

		bool yBelow = targets[3*i + 1]  <= YSTART + RADIUS;
		bool yAbove = targets[3*i + 1]  >= YSTART + GRID_SIZE - RADIUS;

		if(xBelow) {
			targets[3*i] += bump;
			directions[3*i] = abs(directions[3*i]);
		} else if(xAbove) {
			targets[3*i] -= bump;
			directions[3*i] = -abs(directions[3*i]);
		}

		if(yBelow) {
			targets[3*i + 1] += bump;
			directions[3*i + 1] = abs(directions[3*i + 1]);
		} else if(yAbove) {
			targets[3*i + 1] -= bump;
			directions[3*i + 1] = -abs(directions[3*i + 1]);
		}
	}
}


glm::vec4 QuatMult(glm::vec4 q1, glm::vec4 q2) {
	glm::vec4 qr;
	qr.x = (q1.w * q2.x) + (q1.x * q2.w) + (q1.y * q2.z) - (q1.z * q2.y);
	qr.y = (q1.w * q2.y) - (q1.x * q2.z) + (q1.y * q2.w) + (q1.z * q2.x);
	qr.z = (q1.w * q2.z) + (q1.x * q2.y) - (q1.y * q2.x) + (q1.z * q2.w);
	qr.w = (q1.w * q2.w) - (q1.x * q2.x) - (q1.y * q2.y) - (q1.z * q2.z);

	return qr;
}

glm::vec4 QuatConjugate(glm::vec4 q) {
	return glm::vec4(-q.x, -q.y, -q.z, q.w);
}


// Rotates the directions of the targets every ROTATION_PERIOD
void PeriodicRotateDirections2D(float* directions) {
	float angle, dx, dy, cosa, sina;
	for(int i=0; i<CHOSEN_NB_TARGETS_2D; ++i) {
		// The angle should end up between range min and max.
		angle = RANGE_MIN + static_cast <float> (rand()) / (static_cast <float> ( RAND_MAX/(RANGE_MAX-RANGE_MIN) ) );
		dx = directions[3*i];
		dy = directions[3*i + 1];
		cosa = cos(angle);
		sina = sin(angle);
		directions[3*i]     = dx*cosa - dy*sina;
		directions[3*i + 1] = dx*sina + dy*cosa;
	}
}

inline bool IsTiny(float x) {
	return(abs(x) < 0.00001f);
}

// This function takes a vector v = (a,b,c) as three float values and finds a random, normalized vector orthogonal to v.
glm::vec3 FindOrthogonalVector(float a, float b, float c) {
	if( (a == 0) && (b == 0) & (c == 0) ) {
		cerr << "FindOrthogonalVector: nil vector as input, can't find orthogonal plane.";
		return(glm::vec3(0,0,0)); // exception would be much better
	}
	float x, y, z;
	// If I remember correctly, -100 and +100 are arbitrary and unimportant. What matters is that they have opposite signs and identical absolute values.
	// I should double-check this some time.
	if(IsTiny(c)) {
		if(IsTiny(b)) { //  # a,0,0
			x = 0.0f;
			y = RandomFloat(-100.0f, 100.0f);
			z = RandomFloat(-100.0f, 100.0f);
		} else {        //  # a,b,0
			x = RandomFloat(-100.0f, 100.0f);
			z = RandomFloat(-100.0f, 100.0f);
			y = -a*x/b;
		}
	} else {
		x = RandomFloat(-100.0f, 100.0f);
		y = RandomFloat(-100.0f, 100.0f);
		z = -(a*x + b*y)/c;
	}
	glm::vec3 res = glm::vec3(x,y,z);
	res = glm::normalize(res);
	return(res);
}

glm::vec4 AxisToQuat(glm::vec3 axis, float angle) {
	float half_angle = angle/2.0f;
	glm::vec4 q;
	q.x = axis.x * sin(half_angle);
	q.y = axis.y * sin(half_angle);
	q.z = axis.z * sin(half_angle);
	q.w = cos(half_angle);
	return(q);
}

void LogNewCondition(ofstream* condSwitches, float* parameters) {
	*condSwitches   << glfwGetTime() << "   New condition: " << parameters[0] << "  " << parameters[1] << " " << parameters[2] << endl;
}

void UpdateParams(float* oldPar, float* newPar) {
	oldPar[0] = newPar[0];
	oldPar[1] = newPar[1];
	oldPar[2] = newPar[2];
}


// Main loop. Moves the targets and sends them.
//void Communicate(   MessageWrite msgPosOut, float* targets, float* directions,
//					ModuleAPI* module, OutputPort* pOutPositions, InputPort* pInSpringIndex, InputPort* pInParams, ofstream* condSwitches) {
void Communicate(   MessageWrite msgPosOut, float* targets, float* directions,
					ModuleAPI* module, OutputPort* pOutPositions, InputPort* pInSpringIndex, ofstream* condSwitches) {
	double currentTime = glfwGetTime();
	double prevTime = currentTime;
	double prevLogTime = currentTime;
	double deltaT;
	double logDeltaT;
	double timeSinceLastRotation = 3600; // one hour, bogus value to trigger the first rotation
	int currentTarget = 0;
	float* parameters = new float[3];
	float oldParameters[3] = {-1.0f, -1.0f, -1.0f};
	int* springIndex;
	int nbIterations = 0;
	//int totalNbIt = 0;

	parameters[0] = 30.0f;
	parameters[1] = 1.0f;
	parameters[2] = 1.0f;
	// To launch the simulation, otherwise LeanRenderer just keeps waiting.
	module->put(pOutPositions,msgPosOut);

	timespec mySleep, remainder;
	mySleep.tv_sec = 0;
	mySleep.tv_nsec = 8000;
	int nanores = 42;

	while(module->wait()) {
		//nanores = nanosleep(&mySleep , &remainder);
		deltaT = currentTime - prevTime;
		logDeltaT = currentTime - prevLogTime;
		timeSinceLastRotation += deltaT;
		if(logDeltaT > 0) {
			if(nbIterations == 100) {
				//cout << float(nbIterations)/logDeltaT << " it/s" << endl;
				//totalNbIt += nbIterations;
				nbIterations = 0;
				prevLogTime = currentTime;
			}
		}
		++nbIterations;

		ANGLE				= parameters[0] / 180.0f;
		//cout << "Read angle:  " << ANGLE << endl;
		ROTATION_PERIOD		= parameters[1];
		SPEED_2D			= parameters[2];
		RANGE_MIN			= -M_PI * ANGLE;
		RANGE_MAX			= M_PI  * ANGLE;

		if(threeDimensional)
			MoveTargets3D(targets, deltaT, directions);
		else
			MoveTargets2D(targets, deltaT, directions, currentTarget);
		if(timeSinceLastRotation >= ROTATION_PERIOD) {
			if(threeDimensional)
				PeriodicRotateDirections3D(directions);
			else
				PeriodicRotateDirections2D(directions);
			timeSinceLastRotation = 0;
		}

		cout << targets[0] << endl;
		prevTime = currentTime;
		currentTime = glfwGetTime();
		module->put(pOutPositions, msgPosOut);
	}
}


int main(int argc, char *argv[]) {
	cout << "RandomMotion is ready to randomize motion.\n";


	if(argc >= 2) {
		if(atoi(argv[1]) == 3) {
			threeDimensional = true;
			cout << "3D mode, bitches!" << endl;
		}
	} else
		cout << "2D only, sorry." << endl;

	glfwInit();

	long int seed = (long int)glfwGetTime();
	long int* pSeed = &seed;
	srand (static_cast <unsigned> (time(pSeed))); // Seeds the random number generator

	ReadNBTargets();
	// Open in and out ports
	OutputPort* pOutPositions	= new OutputPort("atomPos");
	InputPort*  pInSpringIndex	= new InputPort("springIndex");
	BufferPool* pool			= new BufferPool();
	ModuleAPI* module = 0;

	vector<Port*> ports;
	ports.push_back(pOutPositions);


	// FlowVR initialization
	if (!(module = initModule(ports)))
		return -1;

	if(module == NULL){
		cerr << "Error: RandomMotion Module not initialized." << endl;
		return -1;
	}

	msgPosOut.data = pool->alloc(module->getAllocator(),CHOSEN_NB_TARGETS_2D * 3 * sizeof(float));
	float* positions = (float*) (msgPosOut.data.writeAccess());
	float* directions = new float[CHOSEN_NB_TARGETS_2D * 3 * sizeof(float)];
	if(threeDimensional)
		InitTargets3D(positions);
	else
		InitTargets2D(positions);

	if(threeDimensional) {
		InitDirections3D(directions);
		NormalizeDirections3D(directions);
	} else {
		InitDirections2D(directions);
		NormalizeDirections2D(directions);
	}

	ofstream conditionSwitches;
	conditionSwitches.open("condition_switches.txt");
	// Will loop until it's stopped.
	//Communicate(msgPosOut, positions, directions, module, pOutPositions, pInSpringIndex, pInParams, &conditionSwitches);
	Communicate(msgPosOut, positions, directions, module, pOutPositions, pInSpringIndex, &conditionSwitches);
	conditionSwitches.close();
	return 0;
}
