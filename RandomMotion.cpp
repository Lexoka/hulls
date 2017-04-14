#include <iostream>
#include <fstream> // files, reading and writing
#include <sstream>
#include <iterator>
#include <vector>
#include <stdlib.h> // rand
#define _USE_MATH_DEFINES // for pi
#include <math.h>   // sqrt
#include <stdio.h>
#include <time.h>
#include <flowvr/moduleapi.h>


#include </people/kouyoumdjian/local/glfw/include/GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

using namespace std;
using namespace flowvr;
using namespace glm;

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

const mat4 projectionMatrix = mat4(	1.30323,    0,          0,          0,
									0,          1.30323,    0,          0,
									0,          0,          -1.002,     -1,
									0,          0,          -0.2002,    0
);

const mat4 viewMatrix = mat4(	1,      0,      -0,     0,
								-0,     1,      -0,     0,
								0,      0,      1,      0,
								-5,     -5,     -10,    1
);

const mat4 frustumMatrix = projectionMatrix * viewMatrix;

const float FUDGED_RADIUS = RADIUS + (RADIUS * 0.275f) ;    // OK, so the xyzRadius vectors below are not orthogonal to the frustum planes...
															// ...so the test doesn't work perfectly and the radius needs to be fudged a little.
															// It's gross, but easy, and it should work well enough.
const vec4 xRadius = vec4(FUDGED_RADIUS, 0, 0, 0);
const vec4 yRadius = vec4(0, FUDGED_RADIUS, 0, 0);
const vec4 zRadius = vec4(0, 0, FUDGED_RADIUS, 0);


const float SCALE_FACTOR = 0.10f; // enlarge to get bigger spheres
const vec3 SCALE_VECTOR = vec3(SCALE_FACTOR, SCALE_FACTOR, SCALE_FACTOR);
const mat4 ModelMatrix = scale(mat4(1.0), SCALE_VECTOR);

// Angle of each rotation, in PIs. I.e. if ANGLE == 0.65f, the angle will be between -0.65*PI and 0.65*PI
float ANGLE;
float ROTATION_PERIOD;
float SPEED;

int CHOSEN_NB_TARGETS;

// Minimum and maximum angles for the rotation of directions.
float RANGE_MIN;
float RANGE_MAX;

void ReadNBTargets(){ // And number!
	string line;
	const char* fname;
	fname = "../../share/assets/config/parameters3d.txt";
	
	ifstream myfile(fname);
	if (myfile.is_open()) {
		getline(myfile, line); // NB Targets
		getline(myfile, line);
		CHOSEN_NB_TARGETS = atoi(line.c_str());
		myfile.close();
	}
	else {
		cout << "Unable to open file parameters3d.txt." << endl; 
		cout << "You might be in the wrong directory. Good luck figuring it out." << endl; 
	}       
}


float RandomFloat(float a, float b) {
	float random = ((float) rand()) / (float) RAND_MAX;
	float diff = b - a;
	float r = random * diff;
	return a + r;
}


void InitTargets3D(float* targets) {
	double x, y, z;
	for(int i=0; i<CHOSEN_NB_TARGETS; ++i) {
		x = RandomFloat(XSTART, XEND);
		y = RandomFloat(YSTART, YEND);
		z = RandomFloat(ZSTART3D, ZEND);
		targets[3*i]        = x;
		targets[3*i + 1]    = y;
		targets[3*i + 2]    = z;
		//cout << x << "    " << y << " " << z << endl;
		if(isnan(x) || isnan(y) || isnan(z))
			cout << "Init: NAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAN!:	" << x << " " << y << " " << z << endl;
	}
}


void InitDirections3D(float* directions) {
	float vx, vy, vz;
	for(int i=0; i<CHOSEN_NB_TARGETS; ++i) {
		// vx and vy are random floats between -1.0f and +1.0f
		vx = RandomFloat(-1.0f, 1.0f);
		vy = RandomFloat(-1.0f, 1.0f);
		vz = RandomFloat(-1.0f, 1.0f);
		directions[3*i]     = vx;
		directions[3*i + 1] = vy;
		directions[3*i + 2] = vz;
	}
}


void NormalizeDirections3D(float* directions) {
	float x, y, z, length;
	for(int i=0; i<CHOSEN_NB_TARGETS; i++) {
		x = directions[3*i];
		y = directions[3*i + 1];
		z = directions[3*i + 2];
		length = sqrt(x*x + y*y + z*z);
		if(length != 0.0f) {
			directions[3*i]     = x/length;
			directions[3*i + 1] = y/length;
			directions[3*i + 2] = z/length;
		} else {
			directions[3*i]     = 1.0f;
			directions[3*i + 1] = 1.0f;
			directions[3*i + 2] = 1.0f;
			cerr << "Nil direction found; set to (1.0,1.0)." << endl;
		}
	}
}


void Bounce3D(float* targets, double deltaT, float* directions, int i) {
	float bump = 0.001f; // To make sure targets stay in the window and don't just skirt long the edge, right outside.
	vec4 pos	= vec4(targets[3*i], targets[3*i + 1], targets[3*i + 2], 1.0f);

	// Positions on the surface of the sphere.
	vec4 leftEdge = pos - xRadius;
	vec4 rightEdge = pos + xRadius;

	vec4 bottomEdge = pos - yRadius;
	vec4 topEdge = pos + yRadius;

	vec4 frontEdge = pos + zRadius;
	vec4 backEdge = pos - zRadius;

	// Transforming them into clip space positions.
	leftEdge    = frustumMatrix * leftEdge;
	rightEdge   = frustumMatrix * rightEdge;

	bottomEdge  = frustumMatrix * bottomEdge;
	topEdge     = frustumMatrix * topEdge;

	frontEdge   = frustumMatrix * frontEdge;
	backEdge    = frustumMatrix * backEdge;

	// I may not need to add/substract the little 0.0001f.
	if(leftEdge.x < -leftEdge.w) {
		targets[3*i] += bump;
		directions[3*i] = abs(directions[3*i]);
	}
	else {
		if(rightEdge.x > rightEdge.w) {
			targets[3*i] -= bump;
			directions[3*i] = -abs(directions[3*i]);
		}
	}

	// Y
	if(bottomEdge.y < -bottomEdge.w) {
		targets[3*i + 1] += bump;
		directions[3*i + 1] = abs(directions[3*i + 1]);
	}
	else {
		if (topEdge.y > topEdge.w) {
			targets[3*i + 1] -= bump;
			directions[3*i + 1] = -abs(directions[3*i + 1]);
		}
	}

	// Z
	if(frontEdge.z <= 0.0f + CAMERA_MARGIN) {
		targets[3*i + 2] -= bump;
		directions[3*i + 2] = -abs(directions[3*i + 2]);
	} else {
		if (backEdge.z > backEdge.w) {
			targets[3*i + 2] += bump;
			directions[3*i + 2] = abs(directions[3*i + 2]);
		}
	}
}


void MoveTargets3D(float* targets, double deltaT, float* directions) {
	for(int i=0; i<CHOSEN_NB_TARGETS; ++i) {
		targets[3*i]        += SPEED * deltaT * directions[3*i];         // X
		targets[3*i + 1]    += SPEED * deltaT * directions[3*i + 1];     // Y
		targets[3*i + 2]    += SPEED * deltaT * directions[3*i + 2];     // Z

		//Bounce3D(targets, deltaT, directions, i);
	}
}


vec4 QuatMult(vec4 q1, vec4 q2) { 
	vec4 qr;
	qr.x = (q1.w * q2.x) + (q1.x * q2.w) + (q1.y * q2.z) - (q1.z * q2.y);
	qr.y = (q1.w * q2.y) - (q1.x * q2.z) + (q1.y * q2.w) + (q1.z * q2.x);
	qr.z = (q1.w * q2.z) + (q1.x * q2.y) - (q1.y * q2.x) + (q1.z * q2.w);
	qr.w = (q1.w * q2.w) - (q1.x * q2.x) - (q1.y * q2.y) - (q1.z * q2.z);
	return qr;
}


vec4 QuatConjugate(vec4 q) { 
	return vec4(-q.x, -q.y, -q.z, q.w); 
}


inline bool IsTiny(float x) {
	return(abs(x) < 0.00001f);
}


// This function takes a vector v = (a,b,c) as three float values and finds a random, normalized vector orthogonal to v.
vec3 FindOrthogonalVector(float a, float b, float c) {
	if( (a == 0) && (b == 0) & (c == 0) ) {
		cerr << "FindOrthogonalVector: nil vector as input, can't find orthogonal plane.";
		return(vec3(0,0,0)); // exception would be much better
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
	vec3 res = vec3(x,y,z);
	res = normalize(res);
	return(res);
}


vec4 AxisToQuat(vec3 axis, float angle) {
	float half_angle = angle/2.0f;
	vec4 q;
	q.x = axis.x * sin(half_angle);
	q.y = axis.y * sin(half_angle);
	q.z = axis.z * sin(half_angle);
	q.w = cos(half_angle);
	return(q);
}


void PeriodicRotateDirections3D(float* directions) {
	// The angle should end up between range min and max.
	double a, b, c;
	float angle, half;
	vec4 qResult, qDir, qRotArb, qConjArb;
	vec3 rotAxis;
	
	for(int i=0; i<CHOSEN_NB_TARGETS; ++i) {
		// The angle should end up between range min and max.
		a = directions[3*i];
		b = directions[3*i + 1];
		c = directions[3*i + 2];
		if(isnan(a) || isnan(b) || isnan(c)) {
			cout << "From directions float array: NAN!" << endl;
			cout << a << "  " << b << " " << c << endl;
		}
		angle = RANGE_MIN + static_cast <float> (rand()) / (static_cast <float> ( RAND_MAX/(RANGE_MAX-RANGE_MIN) ) );
		half = angle*0.5f;
		rotAxis = FindOrthogonalVector(a,b,c);      // Random, normalized vector orthogonal to current direction

		qDir		= vec4(a,b,c,0);                  // Current direction as quaternion
		qRotArb		= AxisToQuat(rotAxis, angle);       // Rotation quaternion from rotation axis and angle
		qConjArb	= QuatConjugate(qRotArb);          // Conjugate of rotation quaternion
		qResult		= QuatMult(qRotArb, qDir);          // Multiply rotation by direction
		qResult		= QuatMult(qResult, qConjArb);      // Multiply result by conjugate of rotation quaternion

		directions[3*i]     = float(qResult.x);
		directions[3*i + 1] = float(qResult.y);
		directions[3*i + 2] = float(qResult.z);
	}
}


void UpdateParams(float* oldPar, float* newPar) {
	oldPar[0] = newPar[0];
	oldPar[1] = newPar[1];
	oldPar[2] = newPar[2];
}


// Main loop. Moves the targets and sends them.
//void Communicate(   MessageWrite msgPosOut, float* targets, float* directions, 
//					ModuleAPI* module, OutputPort* pOutPositions, InputPort* pInParams) {
void Communicate(   MessageWrite msgPosOut, float* targets, float* directions, 
					ModuleAPI* module, OutputPort* pOutPositions) { //, BufferPool* pool, float* positions) {
	double currentTime = glfwGetTime();
	double prevTime = currentTime;
	double deltaT;
	double timeSinceLastRotation = 3600; // one hour, bogus value to trigger the first rotation
	float* parameters = new float[3];
	float oldParameters[3] = {-1.0f, -1.0f, -1.0f};

	parameters[0] = 30.0f;	// angle
	parameters[1] = 2.5f;	// frequency
	parameters[2] = 0.2f;	// speed

	// To launch the simulation, otherwise LeanRenderer just keeps waiting.
	timespec mySleep, remainder;
	mySleep.tv_sec = 0;
	mySleep.tv_nsec = 10000000;
	int nanores = 42;
	nanores = nanosleep(&mySleep , &remainder);
	module->put(pOutPositions,msgPosOut);
	
	/*
	msgPosOut.clear(); // fucks things up somehow; the next message becomes empty or something
	msgPosOut.data = pool->alloc(module->getAllocator(),CHOSEN_NB_TARGETS * 3 * sizeof(float));
	positions	= (float*) (msgPosOut.data.writeAccess());
	*/

	while(module->wait()) {
		//nanores = nanosleep(&mySleep , &remainder);
		deltaT = currentTime - prevTime;
		timeSinceLastRotation += deltaT;

		/*
		module->get(pInParams, msgParams);
		if(msgParams.data.getSize() > 0) {
			//parameters = (float*) msgParams.data.readAccess();
			memcpy(parameters, (float*) msgParams.data.readAccess(), msgParams.data.getSize());
			msgParams.clear();
			if(parameters != NULL) {
				ANGLE               = parameters[0] / 180.0f;
				//cout << "Read angle:  " << ANGLE << endl;
				ROTATION_PERIOD     = parameters[1];
				SPEED            = parameters[2];
				RANGE_MIN           = -M_PI * ANGLE;
				RANGE_MAX           = M_PI  * ANGLE;
				UpdateParams(oldParameters, parameters);
				}
			}
		}*/

		ANGLE				= parameters[0] / 180.0f;
		ROTATION_PERIOD		= parameters[1];
		SPEED				= parameters[2];
		RANGE_MIN			= -M_PI * ANGLE;
		RANGE_MAX			= M_PI  * ANGLE;
		
		MoveTargets3D(targets, deltaT, directions);
		
		if(timeSinceLastRotation >= ROTATION_PERIOD) {
			PeriodicRotateDirections3D(directions);
			timeSinceLastRotation = 0;
		}

		prevTime = currentTime;
		currentTime = glfwGetTime();
		module->put(pOutPositions, msgPosOut);
		/*
		msgPosOut.clear(); // fucks things up somehow; the next message becomes empty or something
		msgPosOut.data = pool->alloc(module->getAllocator(),CHOSEN_NB_TARGETS * 3 * sizeof(float));
		positions	= (float*) (msgPosOut.data.writeAccess());
		*/
	}
}


int main(int argc, char *argv[]) {
	glfwInit();

	long int seed = (long int)glfwGetTime();
	long int* pSeed = &seed;
	srand (static_cast <unsigned> (time(pSeed))); // Seeds the random number generator

	ReadNBTargets();

	// Open in and out ports
	OutputPort* pOutPositions	= new OutputPort("atomPos");
	//InputPort*  pInParams		= new InputPort("params");
	BufferPool* pool			= new BufferPool();
	ModuleAPI* module = 0;
	
	vector<Port*> ports;
	ports.push_back(pOutPositions);
	//pInParams->setNonBlockingFlag(true);
	//ports.push_back(pInParams);

	MessageWrite msgPosOut;
	Message msgParams;

	// FlowVR initialization
	if (!(module = initModule(ports)))
		return -1;

	if(module == NULL){
		cerr << "Error: RandomMotion Module not initialized." << endl;
		return -1;
	}   
	
	msgPosOut.data = pool->alloc(module->getAllocator(),CHOSEN_NB_TARGETS * 3 * sizeof(float));
	float* positions	= (float*) (msgPosOut.data.writeAccess());
	float* directions	= new float[CHOSEN_NB_TARGETS * 3 * sizeof(float)];
	InitTargets3D(positions);
	InitDirections3D(directions);
	NormalizeDirections3D(directions);
	
	//Communicate(msgPosOut, positions, directions, module, pOutPositions, pInParams, &conditionSwitches);
	Communicate(msgPosOut, positions, directions, module, pOutPositions); //, pool, positions);
	
	return 0;
}
