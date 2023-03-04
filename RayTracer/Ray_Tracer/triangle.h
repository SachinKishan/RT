#pragma once
#include "hittable.h"
#include "VectorMath.h"

class triangle:public hittable
{
public:
	point3 v0, v1, v2;//counter clockwise
	shared_ptr<material> mat;
	
	triangle(point3 a):v0(a),v1(a),v2(a)
	{
	}

	triangle(point3 a,point3 b,point3 c, shared_ptr<material> m):v0(a),v1(b),v2(c),mat(m)
	{}

	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const override;


};

bool triangle::hit(const ray& r, double t_min, double t_max, hit_record& rec) const
{


	vec3 a = v1 - v0;
	vec3 b = v2 - v0;
	vec3 normal = unit_vector(cross(a, b));

	if (dot(r.dir, normal) < 0.001)return false;//ray and triangle are parallel- no intersection
	//if (dot(r.dir, normal) > 0)return false;

	float D = dot(normal, v0);
	float t = -(dot(normal, r.origin()) + D) / dot(normal, r.dir);
	if (t < 0)return false;//triangle is behind the camera
	vec3 hitPoint = r.at(t);


	// Step 2: inside-outside test
	vec3 C;  //vector perpendicular to triangle's plane 

	// edge 0
	vec3 edge0 = v1 - v0;
	vec3 vp0 = hitPoint - v0;
	C = cross(edge0,vp0);
	if (dot(C,normal) < 0) return false;  //P is on the right side 

	// edge 1
	vec3 edge1 = v2 - v1;
	vec3 vp1 = hitPoint - v1;
	C = cross(edge1, vp1);
	if (dot(C,normal) < 0)  return false;  //P is on the right side 

	// edge 2
	vec3 edge2 = v0 - v2;
	vec3 vp2 = hitPoint - v2;
	C = cross(vp2,edge2);
	if (dot(C,hitPoint) < 0) return false;  //P is on the right side; 

	return true;
}
