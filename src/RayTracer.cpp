// The main ray tracer.

#pragma warning (disable: 4786)

#include "RayTracer.h"
#include "scene/light.h"
#include "scene/material.h"
#include "scene/ray.h"

#include "parser/Tokenizer.h"
#include "parser/Parser.h"

#include "ui/TraceUI.h"
#include <cmath>
#include <algorithm>
#include <glm/glm.hpp>
#include <glm/gtx/io.hpp>
#include <string.h> // for memset

#include <iostream>
#include <fstream>
#include <thread>
#include <time.h>
#include <stdlib.h>

using namespace std;
extern TraceUI* traceUI;

// Use this variable to decide if you want to print out
// debugging messages.  Gets set in the "trace single ray" mode
// in TraceGLWindow, for example.
bool debugMode = false;

// Trace a top-level ray through pixel(i,j), i.e. normalized window coordinates (x,y),
// through the projection plane, and out into the scene.  All we do is
// enter the main ray-tracing method, getting things started by plugging
// in an initial ray weight of (0.0,0.0,0.0) and an initial recursion depth of 0.

glm::dvec3 RayTracer::trace(double x, double y)
{
	// Clear out the ray cache in the scene for debugging purposes,
	if (TraceUI::m_debug)
		scene->intersectCache.clear();

	ray r(glm::dvec3(0,0,0), glm::dvec3(0,0,0), glm::dvec3(1,1,1), ray::VISIBILITY);
	scene->getCamera().rayThrough(x,y,r);
	double dummy;
	glm::dvec3 ret = traceRay(r, glm::dvec3(1.0,1.0,1.0), traceUI->getDepth(), dummy);
	ret = glm::clamp(ret, 0.0, 1.0);
	return ret;
}

glm::dvec3 RayTracer::tracePixel(int i, int j)
{
	glm::dvec3 col(0,0,0);

	if( ! sceneLoaded() ) return col;

	double x = double(i)/double(buffer_width);
	double y = double(j)/double(buffer_height);

	unsigned char *pixel = buffer.data() + ( i + j * buffer_width ) * 3;
	col = trace(x, y);

	pixel[0] = (int)( 255.0 * col[0]);
	pixel[1] = (int)( 255.0 * col[1]);
	pixel[2] = (int)( 255.0 * col[2]);
	return col;
}

glm::dvec3 RayTracer::tracePixelAntiAlias(int i, int j)
{
	glm::dvec3 col(0,0,0);

	if( ! sceneLoaded() ) return col;

	double x = double(i)/double(buffer_width);
	double y = double(j)/double(buffer_height);

	//super sample
    int n_samples = traceUI->getSuperSamples();

	unsigned char *pixel = buffer.data() + ( i + j * buffer_width ) * 3;

    //take several pixels and obtain average
	if (traceUI->aaSwitch()) {
		for (int m = 0; m < n_samples; ++m) {
			for (int n = 0; n < n_samples; ++n) {
				x = (double(i) + double(m) / n_samples) / double(buffer_width);
				y = (double(j) + double(n) / n_samples) / double(buffer_height);
				col += trace(x,y);
			}
		}
		col /= n_samples * n_samples;
	} else {
		col = trace(x,y);
	}
    
	pixel[0] = (int)( 255.0 * col[0]);
	pixel[1] = (int)( 255.0 * col[1]);
	pixel[2] = (int)( 255.0 * col[2]);
	return col;
}


glm::dvec3 RayTracer::tracePixelJitter(int i, int j)
{
	glm::dvec3 col(0,0,0);
	srand((unsigned)time(NULL));

	if( ! sceneLoaded() ) return col;

	double x = (double(i)+0.5)/double(buffer_width);
	double y = (double(j)+0.5)/double(buffer_height);

	glm::dvec3 center_val = trace(x,y);

	double x_low = double(i);
	double y_low = double(j);

	//super sample
    int n_samples = traceUI->getSuperSamples();
	int threshold = traceUI->getAaThreshold();
	double x_new, y_new;
	int count = 0;
	if (traceUI->aaSwitch()){
		for (int m = 0; m < n_samples; m++){
			for (int n = 0; n < n_samples; n++){
				x_new = (x_low + (double(m) + rand()/double(RAND_MAX)) / n_samples) / double(buffer_width);
				y_new = (y_low + (double(n) + rand()/double(RAND_MAX)) / n_samples) / double(buffer_height);
				glm::dvec3 col_temp = trace(x_new, y_new);
				glm::dvec3 col_dif = (col_temp - center_val)*255.0;
				if (fabs(col_dif[0])<=threshold && fabs(col_dif[1])<=threshold && fabs(col_dif[2])<=threshold){
					col += col_temp;
					count ++;
				}
				// if (fabs(col_dif[0]) + fabs(col_dif[1]) + fabs(col_dif[2]) <= threshold){
				// 	col += col_temp;
				// 	count ++;
				// }
			}
		}
		col /= count;
	} else {
		col = center_val;
	}

	unsigned char *pixel = buffer.data() + ( i + j * buffer_width ) * 3;

	pixel[0] = (int)( 255.0 * col[0]);
	pixel[1] = (int)( 255.0 * col[1]);
	pixel[2] = (int)( 255.0 * col[2]);
	return col;
}

#define VERBOSE 0

// Do recursive ray tracing!  You'll want to insert a lot of code here
// (or places called from here) to handle reflection, refraction, etc etc.
glm::dvec3 RayTracer::traceRay(ray& r, const glm::dvec3& thresh, int depth, double& t )
{
	isect i;
	glm::dvec3 colorC;
#if VERBOSE
	std::cerr << "== current depth: " << depth << std::endl;
#endif

	if(scene->intersect(r, i)) {
		// YOUR CODE HERE

		// An intersection occurred!  We've got work to do.  For now,
		// this code gets the material for the surface that was intersected,
		// and asks that material to provide a color for the ray.

		// This is a great place to insert code for recursive ray tracing.
		// Instead of just returning the result of shade(), add some
		// more steps: add in the contributions from reflected and refracted
		// rays.

		const Material& m = i.getMaterial();
		//shade contribution
		colorC = m.shade(scene.get(), r, i);
		if (depth == 0){
			return colorC;
		}
		glm::dvec3 point_inter = r.at(i.getT());
		glm::dvec3 v_n = i.getN();
		glm::dvec3 v_d = r.getDirection();
		glm::dvec3 dir_reflect = glm::normalize(glm::reflect(v_d,v_n));
		glm::dvec3 sin_i = v_d - v_n * glm::dot(v_d , v_n);
		//Reflected ray
		if (m.kr(i) != glm::dvec3(0.0,0.0,0.0)){
		    ray ray_reflect = ray(point_inter, dir_reflect, r.getAtten(), ray::REFLECTION);
		    colorC += m.kr(i) * traceRay(ray_reflect, thresh, depth-1, t);
		}

		//Refracted ray
		if (m.kt(i) != glm::dvec3(0.0,0.0,0.0)){
			double face_dot = glm::dot(v_n , v_d);
			double n_i, n_r;
			glm::dvec3 i_normal;
			int input;
			if (face_dot < 0) { //coming into object
				n_i = 1;
				n_r = m.index(i);
				i_normal = v_n;
				input = 0;
			} else { //coming out of object
				n_i = m.index(i);
				n_r = 1;
				i_normal = -v_n;
				input = 1;
			}
			double n_transfer = n_r / n_i;
			if (glm::length(sin_i) < n_transfer) {
				glm::dvec3 sin_t = sin_i * n_i / n_r;
				double cos_val = glm::sqrt(1 - glm::dot(sin_t,sin_t));
				glm::dvec3 cos_t = -i_normal * cos_val;
				glm::dvec3 dir_refract = glm::normalize(sin_t + cos_t);
				ray ray_refract = ray(point_inter, dir_refract, r.getAtten(), ray::REFRACTION);
				isect j;
				glm::dvec3 kt_i = m.kt(i);
				// if (scene->intersect(ray_refract, j)) {
				// 	double dist = glm::distance(r.at(i.getT()) , ray_refract.at(j.getT()));
				//     glm::dvec3 kt_d = glm::dvec3(glm::pow(kt_i[0],dist) , glm::pow(kt_i[1],dist) , glm::pow(kt_i[2],dist));
                //     colorC += kt_d * traceRay(ray_refract, thresh, depth-1, t);
				// }	else {
                //     colorC +=  traceRay(ray_refract, thresh, depth-1, t);
				// }	
				colorC +=  traceRay(ray_refract, thresh, depth-1, t);		
			}
		}

	} else {
		// No intersection.  This ray travels to infinity, so we color
		// it according to the background color, which in this (simple) case
		// is just black.
		//
		// FIXME: Add CubeMap support here.
		// TIPS: CubeMap object can be fetched from traceUI->getCubeMap();
		//       Check traceUI->cubeMap() to see if cubeMap is loaded
		//       and enabled.

		if (traceUI->cubeMap()){
			colorC = traceUI->getCubeMap()->getColor(r);
		} else {
			colorC = glm::dvec3(0.0, 0.0, 0.0);
		}

		//colorC = glm::dvec3(0.0, 0.0, 0.0);
	}
#if VERBOSE
	std::cerr << "== depth: " << depth+1 << " done, returning: " << colorC << std::endl;
#endif
	return colorC;
}

RayTracer::RayTracer()
	: scene(nullptr), buffer(0), thresh(0), buffer_width(256), buffer_height(256), m_bBufferReady(false)
{
}

RayTracer::~RayTracer()
{
}

void RayTracer::getBuffer( unsigned char *&buf, int &w, int &h )
{
	buf = buffer.data();
	w = buffer_width;
	h = buffer_height;
}

double RayTracer::aspectRatio()
{
	return sceneLoaded() ? scene->getCamera().getAspectRatio() : 1;
}

bool RayTracer::loadScene(const char* fn)
{
	ifstream ifs(fn);
	if( !ifs ) {
		string msg( "Error: couldn't read scene file " );
		msg.append( fn );
		traceUI->alert( msg );
		return false;
	}

	// Strip off filename, leaving only the path:
	string path( fn );
	if (path.find_last_of( "\\/" ) == string::npos)
		path = ".";
	else
		path = path.substr(0, path.find_last_of( "\\/" ));

	// Call this with 'true' for debug output from the tokenizer
	Tokenizer tokenizer( ifs, false );
	Parser parser( tokenizer, path );
	try {
		scene.reset(parser.parseScene());
	}
	catch( SyntaxErrorException& pe ) {
		traceUI->alert( pe.formattedMessage() );
		return false;
	} catch( ParserException& pe ) {
		string msg( "Parser: fatal exception " );
		msg.append( pe.message() );
		traceUI->alert( msg );
		return false;
	} catch( TextureMapException e ) {
		string msg( "Texture mapping exception: " );
		msg.append( e.message() );
		traceUI->alert( msg );
		return false;
	}

	if (!sceneLoaded())
	return false;

	if (traceUI->kdSwitch()){
		//use kdTree
		scene->buildKdTree(traceUI->getMaxDepth(), traceUI->getLeafSize());
		scene->ifUseKD = true;
		//scene->buildKdTree(traceUI->getMaxDepth(), traceUI->getLeafSize());
		//scene->ifUseKD = true;
	} else {
		scene->ifUseKD - false;
	}

	return true;
}


void RayTracer::traceSetup(int w, int h)
{
	if (buffer_width != w || buffer_height != h)
	{
		buffer_width = w;
		buffer_height = h;
		bufferSize = buffer_width * buffer_height * 3;
		buffer.resize(bufferSize);
	}
	std::fill(buffer.begin(), buffer.end(), 0);
	m_bBufferReady = true;

	/*
	 * Sync with TraceUI
	 */

	threads = traceUI->getThreads();
	block_size = traceUI->getBlockSize();
	thresh = traceUI->getThreshold();
	samples = traceUI->getSuperSamples();
	aaThresh = traceUI->getAaThreshold();

	// YOUR CODE HERE
	// FIXME: Additional initializations
	threadDone.resize(traceUI->getThreads(), false);
}

/*
 * RayTracer::traceImage
 *
 *	Trace the image and store the pixel data in RayTracer::buffer.
 *
 *	Arguments:
 *		w:	width of the image buffer
 *		h:	height of the image buffer
 *
 */
void RayTracer::traceImage(int w, int h)
{
	// Always call traceSetup before rendering anything.
	traceSetup(w,h);

	// YOUR CODE HERE
	// FIXME: Start one or more threads for ray tracing
	//
	// TIPS: Ideally, the traceImage should be executed asynchronously,
	//       i.e. returns IMMEDIATELY after working threads are launched.
	//
	//       An asynchronous traceImage lets the GUI update your results
	//       while rendering.
	// for (int i=0; i<w; i++){
	// 	 for (int j=0; j<h; j++){
	// 		 glm::dvec3 prac = tracePixelAntiAlias(i,j);
	// 	 }
	// }

	for (int i=0; i<threads; i++){
		workThreads.push_back(std::thread(&RayTracer::renderThread,this, w, h, i));
	}
}

void RayTracer::renderThread(int width, int height, int index){
	for (int x = 0; x < width; x++){
		for (int y=0; y<height; y++){
			if(y % threads == index){
				//glm::dvec3 prac = this->tracePixelAntiAlias(x,y);
				glm::dvec3 pres = this->tracePixelJitter(x,y);
			}
		}
	}
	this->threadDone[index] = true;
}

int RayTracer::aaImage()
{
	// YOUR CODE HERE
	// FIXME: Implement Anti-aliasing here
	//
	// TIP: samples and aaThresh have been synchronized with TraceUI by
	//      RayTracer::traceSetup() function
	return 0;
}

bool RayTracer::checkRender()
{
	// YOUR CODE HERE
	// FIXME: Return true if tracing is done.
	//        This is a helper routine for GUI.
	//
	// TIPS: Introduce an array to track the status of each worker thread.
	//       This array is maintained by the worker threads.

	int count = 0;
	for (int i = 0; i < threads; i++){
		if (threadDone[i]){
			count++;
		}
	}
	if (count == threads){
		return true;
	}
	return false;
}

void RayTracer::waitRender()
{
	// YOUR CODE HERE
	// FIXME: Wait until the rendering  is done.
	//        This function is essential if you are using an asynchronous
	//        traceImage implementation.
	//
	// TIPS: Join all worker threads here.
	// for (int i=0; i < traceUI->getThreads()-1; i++){
	// 	workThreads[i].join();
	// }

}


glm::dvec3 RayTracer::getPixel(int i, int j)
{
	unsigned char *pixel = buffer.data() + ( i + j * buffer_width ) * 3;
	return glm::dvec3((double)pixel[0]/255.0, (double)pixel[1]/255.0, (double)pixel[2]/255.0);
}

void RayTracer::setPixel(int i, int j, glm::dvec3 color)
{
	unsigned char *pixel = buffer.data() + ( i + j * buffer_width ) * 3;

	pixel[0] = (int)( 255.0 * color[0]);
	pixel[1] = (int)( 255.0 * color[1]);
	pixel[2] = (int)( 255.0 * color[2]);
}


