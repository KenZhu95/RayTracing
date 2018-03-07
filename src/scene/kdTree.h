#pragma once

// Note: you can put kd-tree here

#ifndef __KDTREE_H__
#define __KDTREE_H__

#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <memory>

#include "ray.h"
#include "material.h"
#include "camera.h"
#include "bbox.h"


template <typename T>
class KdTree{
public:
    std::vector<T*> obj_vector;
    glm::dvec3 planeNormal;
    BoundingBox bound;
    BoundingBox splitBound;
    
    KdTree(){
        ifRoot = false;
        leftNode = nullptr;
        rightNode = nullptr;
        index = 4;
        planeNormal = glm::dvec3(0.0,0.0,0.0);
    }

    KdTree(bool isRoot){
        ifRoot = isRoot;
        leftNode = nullptr;
        rightNode = nullptr;
        index = 4;
        planeNormal = glm::dvec3(0.0,0.0,0.0);
    }

    // void setPlaneNormal(glm::dvec3 normal){
    //     planeNormal = normal;
    // }

    // glm::dvec3 getPlaneNormal(){
    //     return planeNormal;
    // }

    void setBoundingBox(BoundingBox box){
        bound = box;
    }

    BoundingBox getBoundingBox(){
        return bound;
    }

    void setSplitBoundingBox(BoundingBox box){
        splitBound = box;
    }

    BoundingBox getSplitBoundingBox(){
        return splitBound;
    }

    void setIfRoot(bool root){
        ifRoot = root;
    }

    bool getIfRoot(){
        return ifRoot;
    }

    void addObject(T* obj){
        obj_vector.push_back(obj);
    }

    T* getObject(int i){
        return obj_vector[i];
    }

    void setLeft(KdTree<T>* left){
        leftNode = left;
    }

    KdTree<T>* getLeft(){
        return leftNode;
    }

    void setRight(KdTree<T>* right){
        rightNode = right;
    }

    KdTree<T>* getRight(){
        return rightNode;
    }

    int getNoObjects(){
        return obj_vector.size();
    }

    void setIndex(int i){
        index = i;
    }

    int getIndex(){
        return index;
    }



private:
    KdTree<T>* leftNode;
    KdTree<T>* rightNode;
    bool ifRoot;


    int index;
};



#endif //__KDTREE_H__
