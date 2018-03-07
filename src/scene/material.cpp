#include "material.h"
#include "../ui/TraceUI.h"
#include "light.h"
#include "ray.h"
extern TraceUI* traceUI;

#include <glm/gtx/io.hpp>
#include <iostream>
#include "../fileio/images.h"

using namespace std;
extern bool debugMode;

Material::~Material()
{
}

// Apply the phong model to this point on the surface of the object, returning
// the color of that point.
glm::dvec3 Material::shade(Scene* scene, const ray& r, const isect& i) const
{
	// YOUR CODE HERE

	// For now, this method just returns the diffuse color of the object.
	// This gives a single matte color for every distinct surface in the
	// scene, and that's it.  Simple, but enough to get you started.
	// (It's also inconsistent with the phong model...)

	// Your mission is to fill in this method with the rest of the phong
	// shading model, including the contributions of all the light sources.
	// You will need to call both distanceAttenuation() and
	// shadowAttenuation()
	// somewhere in your code in order to compute shadows and light falloff.
	//	if( debugMode )
	//		std::cout << "Debugging Phong code..." << std::endl;

	// When you're iterating through the lights,
	// you'll want to use code that looks something
	// like this:
	//
	// for ( const auto& pLight : scene->getAllLights() )
	// {
	//              // pLight has type unique_ptr<Light>
	// 		.
	// 		.
	// 		.
	// }
	glm::dvec3 res = ke(i) + ka(i) * scene->ambient();
	glm::dvec3 point_inters = r.at(i.getT());
	for (const auto& pLight : scene->getAllLights()){
		//Diffuse term
		glm::dvec3 vec_l = glm::normalize(pLight->getDirection(point_inters));
		glm::dvec3 temp = kd(i) * max(glm::dot(vec_l,i.getN()),0.0);
		glm::dvec3 temp_1 = pLight->shadowAttenuation(r, point_inters) * pLight->getColor() * pLight->distanceAttenuation(point_inters);
		glm::dvec3 add_diffuse = temp_1 * temp;
		//Specular term
		//glm::dvec3 vec_v = glm::normalize(scene->getCamera().getEye() - point_inters);
		glm::dvec3 vec_v = r.getDirection();
		glm::dvec3 vec_r0 = glm::reflect(vec_l,i.getN());
		glm::dvec3 vec_r = glm::normalize(vec_r0);
		glm::dvec3 temp1 = ks(i) * glm::pow(max(glm::dot(vec_v,vec_r) , 0.0),shininess(i));
		glm::dvec3 add_specular = temp_1 * temp1;
		res = res +  add_diffuse + add_specular;	
	}
	return res;
}

TextureMap::TextureMap(string filename)
{
	data = readImage(filename.c_str(), width, height);
	if (data.empty()) {
		width = 0;
		height = 0;
		string error("Unable to load texture map '");
		error.append(filename);
		error.append("'.");
		throw TextureMapException(error);
	}
}

glm::dvec3 TextureMap::getMappedValue(const glm::dvec2& coord) const
{
	// YOUR CODE HERE
	//
	// In order to add texture mapping support to the
	// raytracer, you need to implement this function.
	// What this function should do is convert from
	// parametric space which is the unit square
	// [0, 1] x [0, 1] in 2-space to bitmap coordinates,
	// and use these to perform bilinear interpolation
	// of the values.

	double x = coord[0] * width;
    double y = coord[1] * height;
	int x1 = (int)x;
	int y1 = (int)y;
	int x2 = x1 + 1;
    int y2 = y1 + 1;
	glm::dvec3 f11 = getPixelAt(x1, y1);
	glm::dvec3 f12 = getPixelAt(x1, y2);
	glm::dvec3 f21 = getPixelAt(x2, y1);
	glm::dvec3 f22 = getPixelAt(x2, y2);
	glm::dvec3 ff = (y2-y)*(x2-x)*f11 + (y2-y)*(x-x1)*f21 + (y-y1)*(x2-x)*f12 + (y-y1)*(x-x1)*f22;
	return ff;




	return glm::dvec3(1, 1, 1);
}

glm::dvec3 TextureMap::getPixelAt(int x, int y) const
{
	// YOUR CODE HERE
	//
	// In order to add texture mapping support to the
	// raytracer, you need to implement this function.

	if(data.size() == 0){
		return glm::dvec3(1.0,1.0,1.0);
	}
	if (x >= width){
		x = width - 1;
	} 
	if (y >= width){
		y = width - 1;
	}
	//Find the position in the big data array
	int position = (y * width + x)*3;
	return glm::dvec3(double(data[position])/255.0, double(data[position+1])/255.0, double(data[position+2])/255.0);

	//return glm::dvec3(1, 1, 1);
}

glm::dvec3 MaterialParameter::value(const isect& is) const
{
	if (0 != _textureMap)
		return _textureMap->getMappedValue(is.getUVCoordinates());
	else
		return _value;
}

double MaterialParameter::intensityValue(const isect& is) const
{
	if (0 != _textureMap) {
		glm::dvec3 value(
		        _textureMap->getMappedValue(is.getUVCoordinates()));
		return (0.299 * value[0]) + (0.587 * value[1]) +
		       (0.114 * value[2]);
	} else
		return (0.299 * _value[0]) + (0.587 * _value[1]) +
		       (0.114 * _value[2]);
}
