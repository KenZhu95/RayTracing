#include <cmath>
#include <iostream>

#include "light.h"
#include <glm/glm.hpp>
#include <glm/gtx/io.hpp>


using namespace std;

double DirectionalLight::distanceAttenuation(const glm::dvec3& P) const
{
	// distance to light is infinite, so f(di) goes to 0.  Return 1.
	return 1.0;
}


glm::dvec3 DirectionalLight::shadowAttenuation(const ray& r, const glm::dvec3& p) const
{
	// YOUR CODE HERE:
	// You should implement shadow-handling code here.
	isect i;
	ray shadow = ray(p , getDirection(p) , r.getAtten(), ray::SHADOW);
	if (getScene()->intersect(shadow , i)){
			const Material& material = i.getMaterial();
			if (material.Trans()) {
				//double dist = glm::distance(p,r.at(i.getT()));
				//glm::dvec3 kt_i = material.kt(i);
				//return glm::dvec3(glm::pow(kt_i[0],dist) , glm::pow(kt_i[1],dist) , glm::pow(kt_i[2],dist));
				return glm::dvec3(1.0,1.0,1.0);
			} else {
				return glm::dvec3(0.0,0.0,0.0);
			}
	} else {
	    return glm::dvec3(1.0, 1.0, 1.0);
	}
}

glm::dvec3 DirectionalLight::getColor() const
{
	return color;
}

glm::dvec3 DirectionalLight::getDirection(const glm::dvec3& P) const
{
	return -orientation;
}

double PointLight::distanceAttenuation(const glm::dvec3& P) const
{

	// YOUR CODE HERE

	// You'll need to modify this method to attenuate the intensity 
	// of the light based on the distance between the source and the 
	// point P.  For now, we assume no attenuation and just return 1.0
	double dist = glm::distance(position,P);
	return min(1.0 ,1/(constantTerm + linearTerm * dist + quadraticTerm * dist * dist));
}

glm::dvec3 PointLight::getColor() const
{
	return color;
}

glm::dvec3 PointLight::getDirection(const glm::dvec3& P) const
{
	return glm::normalize(position - P);
}


glm::dvec3 PointLight::shadowAttenuation(const ray& r, const glm::dvec3& p) const
{
	// YOUR CODE HERE:
	// You should implement shadow-handling code here.
	isect i;
	ray shadow = ray(p , this->getDirection(p) , r.getAtten(), ray::SHADOW);
	if (getScene()->intersect(shadow , i)){
		double dist = glm::distance(position,p);
		glm::dvec3 sect = shadow.at(i.getT());
		double disti = glm::distance(sect,p);
		if (disti < dist){
			const Material& material = i.getMaterial();
			//return material.kt(i);
			if (material.Trans()) {
				//return glm::dvec3(1.0,1.0,1.0);
				return material.kt(i);
			} else {
				return glm::dvec3(0.0,0.0,0.0);
			}
		} else {
			return glm::dvec3(1.0,1.0,1.0);
		}
	} else {
	    return glm::dvec3(1.0,1.0,1.0);
	}
}

#define VERBOSE 0

