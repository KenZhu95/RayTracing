#include <cmath>
#include<limits>
#include<stack>

#include "scene.h"
#include "light.h"
#include "kdTree.h"
#include "../ui/TraceUI.h"
#include "../SceneObjects/trimesh.h"
#include <glm/gtx/extended_min_max.hpp>
#include <iostream>
#include <glm/gtx/io.hpp>

using namespace std;

bool Geometry::intersect(ray& r, isect& i) const {
	double tmin, tmax;
	if (hasBoundingBoxCapability() && !(bounds.intersect(r, tmin, tmax))) return false;
	// Transform the ray into the object's local coordinate space
	glm::dvec3 pos = transform->globalToLocalCoords(r.getPosition());
	glm::dvec3 dir = transform->globalToLocalCoords(r.getPosition() + r.getDirection()) - pos;
	double length = glm::length(dir);
	dir = glm::normalize(dir);
	// Backup World pos/dir, and switch to local pos/dir
	glm::dvec3 Wpos = r.getPosition();
	glm::dvec3 Wdir = r.getDirection();
	r.setPosition(pos);
	r.setDirection(dir);
	bool rtrn = false;
	if (intersectLocal(r, i))
	{
		// Transform the intersection point & normal returned back into global space.
		i.setN(transform->localToGlobalCoordsNormal(i.getN()));
		i.setT(i.getT()/length);
		rtrn = true;
	}
	// Restore World pos/dir
	r.setPosition(Wpos);
	r.setDirection(Wdir);
	return rtrn;
}

bool Geometry::hasBoundingBoxCapability() const {
	// by default, primitives do not have to specify a bounding box.
	// If this method returns true for a primitive, then either the ComputeBoundingBox() or
    // the ComputeLocalBoundingBox() method must be implemented.

	// If no bounding box capability is supported for an object, that object will
	// be checked against every single ray drawn.  This should be avoided whenever possible,
	// but this possibility exists so that new primitives will not have to have bounding
	// boxes implemented for them.
	return false;
}

void Geometry::ComputeBoundingBox() {
    // take the object's local bounding box, transform all 8 points on it,
    // and use those to find a new bounding box.

    BoundingBox localBounds = ComputeLocalBoundingBox();
        
    glm::dvec3 min = localBounds.getMin();
    glm::dvec3 max = localBounds.getMax();

    glm::dvec4 v, newMax, newMin;

    v = transform->localToGlobalCoords( glm::dvec4(min[0], min[1], min[2], 1) );
    newMax = v;
    newMin = v;
    v = transform->localToGlobalCoords( glm::dvec4(max[0], min[1], min[2], 1) );
    newMax = glm::max(newMax, v);
    newMin = glm::min(newMin, v);
    v = transform->localToGlobalCoords( glm::dvec4(min[0], max[1], min[2], 1) );
    newMax = glm::max(newMax, v);
    newMin = glm::min(newMin, v);
    v = transform->localToGlobalCoords( glm::dvec4(max[0], max[1], min[2], 1) );
    newMax = glm::max(newMax, v);
    newMin = glm::min(newMin, v);
    v = transform->localToGlobalCoords( glm::dvec4(min[0], min[1], max[2], 1) );
    newMax = glm::max(newMax, v);
    newMin = glm::min(newMin, v);
    v = transform->localToGlobalCoords( glm::dvec4(max[0], min[1], max[2], 1) );
    newMax = glm::max(newMax, v);
    newMin = glm::min(newMin, v);
    v = transform->localToGlobalCoords( glm::dvec4(min[0], max[1], max[2], 1) );
    newMax = glm::max(newMax, v);
    newMin = glm::min(newMin, v);
    v = transform->localToGlobalCoords( glm::dvec4(max[0], max[1], max[2], 1) );
    newMax = glm::max(newMax, v);
    newMin = glm::min(newMin, v);
		
    bounds.setMax(glm::dvec3(newMax));
    bounds.setMin(glm::dvec3(newMin));
}

Scene::Scene()
{
	kdDepth = 0;
	kdLeafSize = 0;
	ifUseKD = false;
	kdRoot = nullptr;

}

Scene::~Scene()
{
}

void Scene::add(Geometry* obj) {
	obj->ComputeBoundingBox();
	sceneBounds.merge(obj->getBoundingBox());
	objects.emplace_back(obj);
	if (obj->hasBoundingBoxCapability()){
		boundedObjects.emplace_back(obj);
	}
}

void Scene::add(Light* light)
{
	lights.emplace_back(light);
}


// Get any intersection with an object.  Return information about the 
// intersection through the reference parameter.
bool Scene::intersect(ray& r, isect& i) const {
	double tmin = 0.0;
	double tmax = 0.0;
	bool have_one = false;
	if (this->ifUseKD && this->kdRoot != nullptr){
		if (this->kdRoot->bound.intersect(r, tmin, tmax)){
			intersectKdTree(r, i, this->kdRoot, have_one, tmin, tmax);
		}
		//have_one = intersectKdTree(r, i, this->kdRoot, tmin, tmax);
	} else {
		for(const auto& obj : objects) {
			isect cur;
			if( obj->intersect(r, cur) ) {
				if(!have_one || (cur.getT() < i.getT())) {
					i = cur;
					have_one = true;
				}
			}
		}
	}

	if(!have_one)
		i.setT(1000.0);
	// if debugging,
	if (TraceUI::m_debug)
		intersectCache.push_back(std::make_pair(new ray(r), new isect(i)));
	return have_one;
}

TextureMap* Scene::getTexture(string name) {
	auto itr = textureCache.find(name);
	if (itr == textureCache.end()) {
		textureCache[name].reset(new TextureMap(name));
		return textureCache[name].get();
	}
	return itr->second.get();
}


bool yzCompare(pair<Geometry*, int> first, pair<Geometry*, int> second)
{
	BoundingBox firstBB = first.first->getBoundingBox();
	BoundingBox secondBB = second.first->getBoundingBox();
	double firstX = firstBB.getMin()[0];
	double secondX = secondBB.getMin()[0];
	if (first.second == 1)
	{
		firstX = firstBB.getMax()[0];
	}
	if (second.second == 1)
	{
		secondX = secondBB.getMax()[0];
	}
	if (firstX == secondX)
	{
		if (first.second < second.second)
		{
			return true;
		}
		else if (first.second > second.second)
		{
			return false;
		}
		else
		{
			true;
		}
	}
	return (firstX < secondX);
}

bool xzCompare(pair<Geometry*, int> first, pair<Geometry*, int> second)
{
	BoundingBox firstBB = first.first->getBoundingBox();
	BoundingBox secondBB = second.first->getBoundingBox();
	double firstY = firstBB.getMin()[1];
	double secondY = secondBB.getMin()[1];
	if (first.second == 1)
	{
		firstY = firstBB.getMax()[1];
	}
	if (second.second == 1)
	{
		secondY = secondBB.getMax()[1];
	}
	if (firstY == secondY)
	{
		if (first.second < second.second)
		{
			return true;
		}
		else if (first.second > second.second)
		{
			return false;
		}
		else
		{
			true;
		}
	}
	return (firstY < secondY);
}

bool xyCompare(pair<Geometry*, int> first, pair<Geometry*, int> second)
{
	BoundingBox firstBB = first.first->getBoundingBox();
	BoundingBox secondBB = second.first->getBoundingBox();
	double firstZ = firstBB.getMin()[2];
	double secondZ = secondBB.getMin()[2];
	if (first.second == 1)
	{
		firstZ = firstBB.getMax()[2];
	}
	if (second.second == 1)
	{
		secondZ = secondBB.getMax()[2];
	}
	if (firstZ == secondZ)
	{
		if (first.second < second.second)
		{
			return true;
		}
		else if (first.second > second.second)
		{
			return false;
		}
		else
		{
			true;
		}
	}
	return (firstZ < secondZ);
}

void Scene::intersectKdTree(ray& r, isect& i, KdTree<Geometry>* current_node, bool& have_one, double t_min, double t_max) const{
    stack<KdStack> kdTreeStack;
	struct KdStack rootE = KdStack(current_node, t_min, t_max);
	kdTreeStack.push(rootE);
	while (!kdTreeStack.empty()){
		struct KdStack currentE = kdTreeStack.top();
		kdTreeStack.pop();
		if ((currentE.currentNode->getLeft() == nullptr) && (currentE.currentNode->getRight() == nullptr)) {
			//leaf node reached
			for (cgiter obj = currentE.currentNode->obj_vector.begin(); obj != currentE.currentNode->obj_vector.end(); obj++){
				isect current;
				if ((*obj)->isTrimesh()){
					Trimesh* trimesh = (Trimesh*)(*obj);
					double tri_min = 0.0;
					double tri_max = 0.0;
					bool trimeshHit = trimesh->triKdRoot->bound.intersect(r, tri_min, tri_max);
					intersectKdTree(r, i, trimesh->triKdRoot, have_one, tri_min, tri_max);
				} else {
					if ((*obj)->intersect(r,current)){
						if (!have_one || current.getT() < i.getT()){
							i = current;
							have_one = true;
						}
					}
				}
			}
		} else {
			KdTree<Geometry> *node_close, *node_far;
			if (r.getPosition()[currentE.currentNode->getIndex()] < currentE.currentNode->getSplitBoundingBox().getMin()[currentE.currentNode->getIndex()]){
				node_close = currentE.currentNode->getLeft();
				node_far = currentE.currentNode->getRight();
			} else {
				node_close = currentE.currentNode->getRight();
				node_far = currentE.currentNode->getLeft();
			}

			double t_star = ( currentE.currentNode->getSplitBoundingBox().getMin()[currentE.currentNode->getIndex()] - r.getPosition()[currentE.currentNode->getIndex()] )/ (r.getDirection()[currentE.currentNode->getIndex()]);

			if (t_star > currentE.tMax || t_star < 0){
				kdTreeStack.push(KdStack(node_close, t_min, t_max));
			} else if (t_star < currentE.tMin){
				kdTreeStack.push(KdStack(node_far, t_min, t_max));
			} else if (t_star >= currentE.tMin && t_star <= currentE.tMax){
				kdTreeStack.push(KdStack(node_close, t_min, t_star));
				kdTreeStack.push(KdStack(node_far, t_star, t_max));
			}
		}
	}
	return;
}



// bool Scene::intersectKdTree(ray& r, isect& i, KdTree<Geometry>* current_node, double t_min, double t_max) const{
//     stack<KdStack> kdTreeStack;
// 	struct KdStack rootE = KdStack(current_node, t_min, t_max);
// 	kdTreeStack.push(rootE);
// 	bool have_one = false;
// 	while (!kdTreeStack.empty()){
// 		struct KdStack currentE = kdTreeStack.top();
// 		kdTreeStack.pop();
// 		if ((currentE.currentNode->getLeft() == nullptr) && (currentE.currentNode->getRight() == nullptr)) {
// 			//leaf node reached
// 			for (int ii = 0; ii < currentE.currentNode->obj_vector.size(); ii++){
// 				isect current;
// 				if (currentE.currentNode->obj_vector[ii]->intersect(r, current)){
// 					if (!have_one || (current.getT()<i.getT())){
// 						i = current;
// 						have_one = true;
// 					}
// 				}
// 			}
// 		} else {
// 			//not a leaf node

// 			double taxis_max = r.at(t_max)[currentE.currentNode->getIndex()];
// 			double taxis_min = r.at(t_min)[currentE.currentNode->getIndex()];
// 			KdTree<Geometry>* left = currentE.currentNode->getLeft();
// 			KdTree<Geometry>* right = currentE.currentNode->getRight();

// 			double t_star = ( currentE.currentNode->getSplitBoundingBox().getMin()[currentE.currentNode->getIndex()] - r.getPosition()[currentE.currentNode->getIndex()] )/ (r.getDirection()[currentE.currentNode->getIndex()]);
// 			if (t_star > taxis_max && t_star > taxis_min){
// 				kdTreeStack.push(KdStack(left, t_min, t_max));
// 			} else if (t_star < taxis_max && t_star < taxis_min) {
// 				kdTreeStack.push(KdStack(right, t_min, t_max));
// 			} else {
// 				double rmin, rmax;
// 				if (currentE.currentNode->getRight()->bound.intersect(r,rmin,rmax)){
// 					kdTreeStack.push(KdStack(right, rmin, rmax));
// 				}
// 				double lmin, lmax;
// 				if (currentE.currentNode->getRight()->bound.intersect(r,rmin,rmax)){
// 					kdTreeStack.push(KdStack(left, lmin, lmax));
// 				}
// 			}

// 		}
// 	}
// 	return have_one;
// }


// bool Scene::intersectKdMain(ray& r, isect& i) const{
//     return false;
// }


void Scene::buildKdTree(int depth, int leaf_size){
	this->kdDepth = depth;
	this->kdLeafSize = leaf_size;
	this->kdRoot = new KdTree<Geometry>(true);
	for (const auto& pObj : this->getAllObjects()){
		if (pObj->hasBoundingBoxCapability()){
			kdRoot->obj_vector.push_back(&*pObj);
		}
	}
	kdRoot->setBoundingBox(this->bounds());
	std::vector<pair<Geometry*, int>> yzOrder;
	std::vector<pair<Geometry*, int>> xzOrder;
	std::vector<pair<Geometry*, int>> xyOrder;
	for (int i=0; i< kdRoot->getNoObjects(); i++){
		yzOrder.emplace_back(pair<Geometry*, int>(kdRoot->obj_vector[i], 0));
		yzOrder.emplace_back(pair<Geometry*, int>(kdRoot->obj_vector[i], 1));
		xzOrder.emplace_back(pair<Geometry*, int>(kdRoot->obj_vector[i], 0));
		xzOrder.emplace_back(pair<Geometry*, int>(kdRoot->obj_vector[i], 1));
		xyOrder.emplace_back(pair<Geometry*, int>(kdRoot->obj_vector[i], 0));
		xyOrder.emplace_back(pair<Geometry*, int>(kdRoot->obj_vector[i], 1));
	}

	vector<vector<pair<Geometry*, int>>> orderPlanes;
	orderPlanes.emplace_back(yzOrder);
	orderPlanes.emplace_back(xzOrder);
	orderPlanes.emplace_back(xyOrder);
	buildKdDetail(kdRoot, depth, leaf_size, orderPlanes);

}


void Scene::buildKdDetail(KdTree<Geometry>* kd_tree, int depth, int leaf_size, std::vector<std::vector<std::pair<Geometry*, int>>> orderedPlanes){
	sort(orderedPlanes[0].begin(), orderedPlanes[0].end(), yzCompare);
	sort(orderedPlanes[1].begin(), orderedPlanes[1].end(), xzCompare);
	sort(orderedPlanes[2].begin(), orderedPlanes[2].end(), xyCompare);
	if (depth == 0){
		return;
	} else {
		for (int i = 0; i < kd_tree->getNoObjects(); i++){
			if (kd_tree->obj_vector[i]->isTrimesh()){
				Trimesh *trimesh = (Trimesh*)(kd_tree->obj_vector[i]);
				if (!trimesh->kdTreeBuilt()){
					buildKdTrimesh(trimesh, this->kdDepth, this->kdLeafSize);
				}
			}
		}
	}
    
	double cost = numeric_limits<double>::max();
	double t_time = 15.0;
	double i_time = 80 * t_time;
	double total_area = kd_tree->bound.area();
	int final_index;
	glm::dvec3 final_point;
	for (int dim=0; dim<3; dim++){
		double area_left = 0;
		double area_right = 0;
		double num_left = 0;
		double num_right = 0;
		for (vector<pair<Geometry*,int>>::const_iterator it = orderedPlanes[dim].begin(); it != orderedPlanes[dim].end(); it++){
			BoundingBox box = (*it).first->getBoundingBox();
			if ((*it).second == 0){
				area_right += box.area();
				num_right += 1;
			}    	
		}

		for (vector<pair<Geometry*,int>>::const_iterator it = orderedPlanes[dim].begin(); it != orderedPlanes[dim].end(); it++){
			
			BoundingBox box = (*it).first->getBoundingBox();
			if ((*it).second == 0){
				area_left += box.area();
				num_left += 1;
			} else if ((*it).second == 1){
				num_right -= 1;
				area_right -= box.area();
			}
			double current_cost = t_time + (area_left/total_area)*num_left*i_time + (area_right/total_area)*num_right*i_time;
			if (current_cost < cost){
				cost = current_cost;
				final_index = dim;
				final_point = box.getMin();
				if ((*it).second == 1){
					final_point = box.getMax();
				}
			}
		}
	}

	KdTree<Geometry>* left_child = new KdTree<Geometry>();
	KdTree<Geometry>* right_child = new KdTree<Geometry>();
	vector<vector<pair<Geometry*, int>>> leftOrderPlanes(3);
	vector<vector<pair<Geometry*, int>>> rightOrderPlanes(3);
	leftOrderPlanes.emplace_back(vector<pair<Geometry*, int>>());
	leftOrderPlanes.emplace_back(vector<pair<Geometry*, int>>());
	leftOrderPlanes.emplace_back(vector<pair<Geometry*, int>>());
	rightOrderPlanes.emplace_back(vector<pair<Geometry*, int>>());
	rightOrderPlanes.emplace_back(vector<pair<Geometry*, int>>());
	rightOrderPlanes.emplace_back(vector<pair<Geometry*, int>>());
	kd_tree->planeNormal[final_index] = final_point[final_index];
	kd_tree->setIndex(final_index);
	glm::dvec3 min_point(-1.0e308, -1.0e308, -1.0e308);
    glm::dvec3 max_point(1.0e308, 1.0e308, 1.0e308);
	min_point[final_index] = final_point[final_index];
	max_point[final_index] = final_point[final_index];
	kd_tree->setSplitBoundingBox(BoundingBox(min_point, max_point));
	for (vector<Geometry*>::const_iterator it = kd_tree->obj_vector.begin(); it != kd_tree->obj_vector.end(); it++){
		BoundingBox box = (*it)->getBoundingBox();
		if (box.getMax()[final_index] < final_point[final_index]){
			left_child->obj_vector.emplace_back(*it);
			left_child->bound.merge(box);
			leftOrderPlanes[0].emplace_back(pair<Geometry*, int>(*it, 0));
			leftOrderPlanes[0].emplace_back(pair<Geometry*, int>(*it, 1));
			leftOrderPlanes[1].emplace_back(pair<Geometry*, int>(*it, 0));
			leftOrderPlanes[1].emplace_back(pair<Geometry*, int>(*it, 1));
			leftOrderPlanes[2].emplace_back(pair<Geometry*, int>(*it, 0));
			leftOrderPlanes[2].emplace_back(pair<Geometry*, int>(*it, 1));
		} else if (box.getMin()[final_index] > final_point[final_index]){
			right_child->obj_vector.emplace_back(*it);
			right_child->bound.merge(box);
			rightOrderPlanes[0].emplace_back(pair<Geometry*, int>(*it, 0));
			rightOrderPlanes[0].emplace_back(pair<Geometry*, int>(*it, 1));
			rightOrderPlanes[1].emplace_back(pair<Geometry*, int>(*it, 0));
			rightOrderPlanes[1].emplace_back(pair<Geometry*, int>(*it, 1));
			rightOrderPlanes[2].emplace_back(pair<Geometry*, int>(*it, 0));
			rightOrderPlanes[2].emplace_back(pair<Geometry*, int>(*it, 1));
		} else {
			left_child->obj_vector.emplace_back(*it);
			left_child->bound.merge(box);
			leftOrderPlanes[0].emplace_back(pair<Geometry*, int>(*it, 0));
			leftOrderPlanes[0].emplace_back(pair<Geometry*, int>(*it, 1));
			leftOrderPlanes[1].emplace_back(pair<Geometry*, int>(*it, 0));
			leftOrderPlanes[1].emplace_back(pair<Geometry*, int>(*it, 1));
			leftOrderPlanes[2].emplace_back(pair<Geometry*, int>(*it, 0));
			leftOrderPlanes[2].emplace_back(pair<Geometry*, int>(*it, 1));
			right_child->obj_vector.emplace_back(*it);
			right_child->bound.merge(box);
			rightOrderPlanes[0].push_back(pair<Geometry*, int>(*it, 0));
			rightOrderPlanes[0].push_back(pair<Geometry*, int>(*it, 1));
			rightOrderPlanes[1].push_back(pair<Geometry*, int>(*it, 0));
			rightOrderPlanes[1].push_back(pair<Geometry*, int>(*it, 1));
			rightOrderPlanes[2].push_back(pair<Geometry*, int>(*it, 0));
			rightOrderPlanes[2].push_back(pair<Geometry*, int>(*it, 1));
		}
	}

	kd_tree->setLeft(left_child);
	kd_tree->setRight(right_child);
	if (kd_tree->getLeft()->getNoObjects() > leaf_size){
		buildKdDetail(kd_tree->getLeft(), depth-1, leaf_size, leftOrderPlanes);
	} else {
		for (int i=0; i < kd_tree->getLeft()->getNoObjects(); i++){
			if (kd_tree->getLeft()->obj_vector[i]->isTrimesh()){
				Trimesh *trimesh = (Trimesh*)(kd_tree->getLeft()->obj_vector[i]);
				if (!trimesh->kdTreeBuilt()){
					buildKdTrimesh(trimesh, this->kdDepth, this->kdLeafSize);
				}
			}
		}
	}

	if (kd_tree->getRight()->getNoObjects() > leaf_size){
		buildKdDetail(kd_tree->getRight(), depth-1, leaf_size, rightOrderPlanes);
	} else {
		for (int i=0; i < kd_tree->getRight()->getNoObjects(); i++){
			if (kd_tree->getRight()->obj_vector[i]->isTrimesh()){
				Trimesh *trimesh = (Trimesh*)(kd_tree->getRight()->obj_vector[i]);
				if (!trimesh->kdTreeBuilt()){
					buildKdTrimesh(trimesh, this->kdDepth, this->kdLeafSize);
				}
			}
		}
	}
}

void Scene::buildKdTrimesh(Geometry* triMesh, int depth, int leaf_size){
	Trimesh *trimesh = (Trimesh*)(triMesh);
	trimesh->triKdRoot = new KdTree<Geometry>(true);
	for (int i=0; i < trimesh->getFaces().size(); i++){
		trimesh->triKdRoot->obj_vector.push_back(trimesh->getFaces()[i]);
	}

	trimesh->triKdRoot->setBoundingBox(trimesh->getBoundingBox());
	std::vector<pair<Geometry*, int>> yzOrder;
	std::vector<pair<Geometry*, int>> xzOrder;
	std::vector<pair<Geometry*, int>> xyOrder;
	for (int i=0; i<trimesh->triKdRoot->getNoObjects(); i++){
		yzOrder.push_back(pair<Geometry*, int>(trimesh->triKdRoot->obj_vector[i], 0));
		yzOrder.push_back(pair<Geometry*, int>(trimesh->triKdRoot->obj_vector[i], 1));
		xzOrder.push_back(pair<Geometry*, int>(trimesh->triKdRoot->obj_vector[i], 0));
		xzOrder.push_back(pair<Geometry*, int>(trimesh->triKdRoot->obj_vector[i], 1));
		xyOrder.push_back(pair<Geometry*, int>(trimesh->triKdRoot->obj_vector[i], 0));
		xyOrder.push_back(pair<Geometry*, int>(trimesh->triKdRoot->obj_vector[i], 1));
	}

	sort(yzOrder.begin(), yzOrder.end(), yzCompare);
	sort(xzOrder.begin(), xzOrder.end(), xzCompare);
	sort(xyOrder.begin(), xyOrder.end(), xyCompare);

	vector<vector<pair<Geometry*, int>>> orderPlanes;
	orderPlanes.push_back(yzOrder);
    orderPlanes.push_back(xzOrder);
    orderPlanes.push_back(xyOrder);
	buildKdDetail(trimesh->triKdRoot, depth, leaf_size, orderPlanes);
}



