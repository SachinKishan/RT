#pragma once
#ifndef HITTABLE_H
#define HITTABLE_H

#include "ray.h"
#include "utility.h"
#include "Matrix.h"
#include "aabb.h"
class material;
struct hit_record {
    point3 p;//point of intersection
    vec3 normal;
    double t;

    double u;
    double v;//surface coordinates
    
    bool front_face;
    bool hit_one_point;
    shared_ptr<material> mat_ptr;
    inline void set_face_normal(const ray& r, const vec3& outward_normal) 
    {
        front_face = dot(r.direction(), outward_normal) < 0;
        normal = front_face ? outward_normal : -outward_normal;
    }

};

class hittable {
public:
    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
    virtual bool bounding_box(double time0, double time1, aabb& output_box) const = 0;
};



class translate : public hittable {
public:
    translate(shared_ptr<hittable> p, const vec3& displacement)
        : ptr(p), offset(displacement) {}

    virtual bool hit(
        const ray& r, double t_min, double t_max, hit_record& rec) const override;

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override;

public:
    shared_ptr<hittable> ptr;
    vec3 offset;
};

bool translate::hit(const ray& r, double t_min, double t_max, hit_record& rec) const
{
    ray moved_r(r.origin() - offset, r.direction());

    if (!ptr->hit(moved_r, t_min, t_max, rec))
        return false;

    //set new normal position
    rec.p += offset;
    rec.set_face_normal(moved_r, rec.normal);

    return true;
}

bool translate::bounding_box(double time0, double time1, aabb& output_box) const {
    if (!ptr->bounding_box(time0, time1, output_box))
        return false;
    output_box = aabb(
        output_box.min() + offset,
        output_box.max() + offset);
    return true;
}

class instance:public hittable
{
public:
    instance(shared_ptr<hittable> p,Matrix44<double> mat):ptr(p),transformationMatrix(mat),inverseTransformMatrix(mat.real_inverse()){}
    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override
    {
        vec3 rayOrigin = inverseTransformMatrix.multiplyVectorMatrix(r.orig);
        vec3 rayDir = inverseTransformMatrix.multiplyPointMatrix(r.dir);
        //rayDir = rayDir - rayOrigin;
       // rayDir = unit_vector(rayDir);
        ray inverseRay(rayOrigin, rayDir);

        if (ptr->hit(inverseRay, t_min, t_max, rec))
        {
            rec.normal = unit_vector(inverseTransformMatrix.returnTransposed().multiplyPointMatrix(rec.normal));
        	return true;
        }
    	return false;
    }

    virtual bool bounding_box(double time0, double time1, aabb& output_box) const override
    {
        if (!ptr->bounding_box(time0, time1, output_box))
            return false;
        output_box = aabb(
             inverseTransformMatrix.multiplyVectorMatrix(output_box.min()),
            inverseTransformMatrix.multiplyVectorMatrix(output_box.max()));
        return true;
    }
private:
    Matrix44<double> transformationMatrix;
    Matrix44<double> inverseTransformMatrix;
    shared_ptr<hittable> ptr;
};




#endif