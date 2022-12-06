#pragma once
#ifndef MATERIAL_H
#define MATERIAL_H

#include "utility.h"
#include "VectorMath.h"
#include "texture.h"
#include "hittable.h"
struct hit_record;

class material {
public:
    virtual bool scatter(
        const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered
    ) const = 0;


    virtual color emitted(double u, double v, const point3& p) const {
        return color(0, 0, 0);
    }
};

class lambertian : public material {
public:
    lambertian(const color& a) : albedo(make_shared<solid_color>(a)) {}
    lambertian(shared_ptr<texture> a) : albedo(a) {}

    virtual bool scatter(
        const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered
    ) const override {
        auto scatter_direction = rec.normal + random_unit_vector();
        // Catch degenerate scatter direction
        if (scatter_direction.near_zero())
            scatter_direction = rec.normal;

        
        scattered = ray(rec.p, scatter_direction);
        attenuation = albedo->value(rec.u, rec.v, rec.p);        
        return true;
    }

public:
    shared_ptr<texture> albedo;
};

class metal : public material {
public:
    metal(const color& a, double f=0) : albedo(a), fuzz(f < 1 ? f : 1) {}

    

    virtual bool scatter(
        const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered
    ) const override {
        vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
        scattered = ray(rec.p, reflected + fuzz * random_in_unit_sphere());
        attenuation = albedo;
        return (dot(scattered.direction(), rec.normal) > 0);
    }

public:
    color albedo;
    double fuzz; 
};

class dielectric : public material {
public:
    dielectric(double index_of_refraction) : ir(index_of_refraction) {}

    virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered)const override 
    {
        attenuation = color(1.0, 1.0, 1.0);
        double refraction_ratio = rec.front_face ? (1.0 / ir) : ir;

        vec3 unit_direction = unit_vector(r_in.direction());
        double cos_theta = fmin(dot(-unit_direction, rec.normal), 1.0);
        double sin_theta = sqrt(1.0 - cos_theta * cos_theta);

        bool cannot_refract = refraction_ratio * sin_theta > 1.0;
        vec3 direction;

        if (cannot_refract || reflectance(cos_theta, refraction_ratio) > random_double())           
            direction = reflect(unit_direction, rec.normal);
        else
            direction = refract(unit_direction, rec.normal, refraction_ratio);

        scattered = ray(rec.p, direction);
        return true;
    }

public:
    double ir; // Index of Refraction

private:
    static double reflectance(double cosine, double ref_idx) {
        // Use Schlick's approximation for reflectance.
        auto r0 = (1 - ref_idx) / (1 + ref_idx);
        r0 = r0 * r0;
        return r0 + (1 - r0) * pow((1 - cosine), 5);
    }
};

class diffuse_light : public material {
public:
    diffuse_light(shared_ptr<texture> a) : emit(a) {}
    diffuse_light(color c) : emit(make_shared<solid_color>(c)) {}

    virtual bool scatter(
        const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered
    ) const override {
        return false;
    }

    virtual color emitted(double u, double v, const point3& p) const override {
        return emit->value(u, v, p);
    }

public:
    shared_ptr<texture> emit;
};

class basic_color :public material
{
public:
    basic_color(color a) :albedo(a)
    {}

    virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered) const override
    {
        attenuation = albedo;
        return true;
    }


    color albedo;

};

class lambert:public material
{
public:

    lambert(color c,float dc,float sc):diffuseCoefficient(dc),albedo(c),specularCoefficient(sc)
    {
	    
    }

	virtual bool scatter(const ray& r_in, const hit_record& rec, color& attenuation, ray& scattered) const override
	{
        vec3 directionOfLight1(1, 8, 7);
        vec3 directionOfLight2(9, 2, 9);
        color lightColor1(0, 0, 1);
        color lightColor2(1, 1, 1);
        float intensityOfLight1 = 0.1;
        float intensityOfLight2 = 0.4;

        vec3 normal = rec.normal;



        float l1 = diffuseCoefficient* intensityOfLight1 * std::max(0.0, dot(directionOfLight1, normal));
        float l2 = diffuseCoefficient * intensityOfLight2 * std::max(0.0, dot(directionOfLight2, normal));

        float p = 50;
        
        vec3 h1 = unit_vector(r_in.direction() + directionOfLight1);
        vec3 h2 = unit_vector(r_in.direction() + directionOfLight2);
        float blinn1 = intensityOfLight1 * std::pow((double)std::max(0.0, dot(normal, h1)),(double)p);
        float blinn2 = intensityOfLight2 * std::pow((double)std::max(0.0, dot(normal, h2)),(double)p);
        
    	attenuation = (blinn2*lightColor2)+albedo*l2+ (blinn1 * lightColor1) + albedo * l1;
        return true;
	}


    float diffuseCoefficient;
    float specularCoefficient;
    color albedo;

};

#endif