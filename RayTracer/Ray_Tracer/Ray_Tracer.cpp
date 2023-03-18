// Ray_Tracer.cpp : This file contains the 'main' function. Program execution begins and ends there.
//



///Todo
//////All original RT Features
///Camera
///Objects
///Scene manipulation
///
///
#include <algorithm>
#include <functional>
#include<stdlib.h>
#include "VectorMath.h"
#include "Matrix.h"
#include <iostream>
#include <random>

#include "aabb.h"
#include "box.h"
#include "bvh.h"
#include "camera.h"
#include "hittable.h"
#include "hittable_list.h"
#include "lodepng.h"
#include "material.h"
#include "perlin.h"
#include "ray.h"
#include "rectangle.h"
#include "sphere.h"
std::vector<unsigned char> backgroundimage; //the raw pixels
unsigned bgwidth, bgheight;

#pragma region PNG
void encodeOneStep(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height) {
	//Encode the image
	unsigned error = lodepng::encode(filename, image, width, height);

	//if there's an error, display it
	if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
}

void decodeOneStep(const char* filename) {

    //decode
    unsigned error = lodepng::decode(backgroundimage, bgwidth, bgheight, filename);

    //if there's an error, display it
    if (error) std::cout << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;

    //the pixels are now in the vector "image", 4 bytes per pixel, ordered RGBARGBA..., use it as texture, draw it, ...
}
#pragma endregion



color ray_color(const ray& r, const hittable_list& world, int depth, color defaultColor, float sx, float sy, bool scattering = true) {
    hit_record rec;
    std::vector<bool> shouldLight;
    defaultColor = Black;
    if (depth <= 0)
        return color(0, 0, 0);
    if (!world.hit(r, 0.001, infinity, rec))
    {
        return defaultColor;
    } 
    
        //ray scattered;
        //color attenuation;
        //color finalCol;
        //color shadowing;
        /*for (PointLight l : world.lights)
        {
            bool blocked_by_all = true;
            hit_record shadowRec;
            vec3 dir = vec3(
                l.lightDirection.x(),
                l.lightDirection.y(),
                l.lightDirection.z());

            Matrix44<double> a = build_local_coords(dir);

            ray shadow_ray(rec.p + (rec.normal * 1e-4), dir);

            //shadow_ray.dir = a.multiplyVectorMatrix(shadow_ray.dir);
            //shadow_ray.orig = a.multiplyVectorMatrix(shadow_ray.orig);

            if (world.hit(shadow_ray, 0.001, infinity, shadowRec))
            {
                shouldLight.push_back(false);
            }
            else { shouldLight.push_back(true); blocked_by_all = false; }
            shadowing = blocked_by_all ? Black : White;

        }*/
        ray scattered;
        color attenuation;
        color emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);
        if (!rec.mat_ptr->scatter(r, rec, attenuation, scattered,world.lights, shouldLight))
        {
            return emitted;
        }



        /*
        if (scattering)
            finalCol = attenuation * ray_color(scattered, world, depth - 1, defaultColor, sx, sy);
        else
            finalCol = attenuation;
        attenuation = finalCol;
        */

        return emitted + attenuation * ray_color(scattered, world, depth-1, defaultColor, sx, sy);
}

hittable_list earth() {
    auto earth_texture = make_shared<image_texture>("earthmap.jpg");
    auto light = make_shared<diffuse_light>(color(15, 15, 15));

    auto earth_surface = make_shared<lambertian>(Red);

    auto globe = make_shared<sphere>(point3(0, 0, 0), 4, earth_surface);
    auto l = make_shared<sphere>(point3(0, 4, 1), 0.5, light);
    Matrix44<double> f(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1); 
    ScaleMatrix scale(vec3(2, 1, 2));
    TranslationMatrix translation(0, 0, 0);
    f = f*scale* translation;
    for(int i=0;i<4;i++)
    {
	    for(int j=0;j<4;j++)
	    {
            std::cout << f[i][j] << " ";
	    }
        std::cout << std::endl;
    }
    auto spin = make_shared<instance>(globe, f);

    hittable_list world;
    world.add(l);
    world.add(spin);

    return world;
}
hittable_list cornell_box() {
    hittable_list objects;

    auto red = make_shared<lambertian>(color(.65, .05, .05));
    auto white = make_shared<lambertian>(color(.73, .73, .73));
    auto green = make_shared<lambertian>(color(.12, .45, .15));
    auto light = make_shared<diffuse_light>(color(15, 15, 15));

    objects.add(make_shared<yz_rect>(0, 555, 0, 555, 555, green));
     objects.add(make_shared<yz_rect>(0, 555, 0, 555, 0, red));
   objects.add(make_shared<xz_rect>(200, 343, 227, 320, 554, light));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 0, white));
    objects.add(make_shared<xz_rect>(0, 555, 0, 555, 555, white));
    objects.add(make_shared<xy_rect>(0, 555, 0, 555, 555, white));

    auto ball = make_shared<sphere>(point3(200, 200, 200), 50, red);

    auto cube=make_shared<box>(point3(130, 0, 65), point3(295, 165, 230), white);

    RotationMatrixZ zRot(10);
    RotationMatrixX XRot(1);
    Matrix44<double> f(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
    TranslationMatrix translation(0, 1, 1);

    f = f * XRot * zRot * translation;

	auto rotateCube = make_shared<instance>(cube, f);

    objects.add(cube);
    objects.add(rotateCube);
    objects.add(ball);

    return objects;
}


hittable_list cube_test_world()
{
    hittable_list world;

    auto ground_material = make_shared<lambert>(color(0.75, 0.75, 0.5), 0.5, 0);
    auto ground_material_l = make_shared<lambertian>(color(0.75, 0.75, 0.75));
    hittable_list boxes2;
    auto white = make_shared<lambertian>(.73);
    int ns = 100;
    //for (int j = 0; j < ns; j++) {
       // boxes2.add(make_shared<sphere>(point3::random(0, 5), 1, material3));
      //  boxes2.add(make_shared<sphere>(point3::random(0, 5), 1, material3));

    //world.add(make_shared<bvh_node>(boxes2, 0, 1));


    auto difflight = make_shared<diffuse_light>(color(4, 4, 4));
    auto cube = make_shared<box>(point3(4, 4, 4), point3(9, 9, 9), white);
    auto material2_lambert = make_shared<lambert>(color(0.4, 0.2, 0.1), 1.0, 0.2);
    world.add(make_shared<sphere>(point3(0, -999, 0), 1000, material2_lambert));


    auto l=(make_shared<sphere>(point3(7, 16, 7), 4, difflight));

    RotationMatrixZ zRot(10);
    RotationMatrixX XRot(1);
    ScaleMatrix scale(vec3(2,2,2));
    
    Matrix44<double> f(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
    TranslationMatrix translation(0, 0, 5);

    
    f = f * translation * scale;
    auto scaledLight = make_shared<instance>(l, f);
    auto rotateCube = make_shared<instance>(cube, RotationMatrixZ(45)*RotationMatrixX(45)*translation);
    

    world.add(l);
    world.add(scaledLight);
    world.add(cube);
    world.add(rotateCube);

    return world;
}


hittable_list random_scene() {
    hittable_list world;

    auto ground_material = make_shared<lambert>(color(0.75, 0.75, 0.5), 0.5, 0);
    auto ground_material_l = make_shared<lambertian>(color(0.75, 0.75, 0.75));
    
    vec3 directionOfLight1(7, 9, 9);
    vec3 directionOfLight2(-9, 10, 9);
    color lightColor1(1, 1, 1);
    color lightColor2(1, 1, 1);
    float intensityOfLight1 = 0.01;
    float intensityOfLight2 = 0.1;

    PointLight l1(intensityOfLight1, directionOfLight1, lightColor1);
    PointLight l2(intensityOfLight2, directionOfLight2, lightColor2);
    //world.lights.push_back(l1);
	//world.lights.push_back(l2);


	/*
    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            auto choose_mat = random_double();
            point3 center(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());

            if ((center - point3(4, 0.2, 0)).length() > 0.9) {
                shared_ptr<material> sphere_material;

                if (choose_mat < 0.8) {
                    // diffuse
                    auto albedo = color::random() * color::random();
                    sphere_material = make_shared<lambertian>(albedo);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
                else if (choose_mat < 0.95) {
                    // metal
                    auto albedo = color::random(0.5, 1);
                    auto fuzz = random_double(0, 0.5);
                    sphere_material = make_shared<metal>(albedo, fuzz);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
                else {
                    // glass
                    sphere_material = make_shared<dielectric>(1.5);
                    world.add(make_shared<sphere>(center, 0.2, sphere_material));
                }
            }
        }
    }*/

    //auto material1 = make_shared<dielectric>(1.5);
    //world.add(make_shared<sphere>(point3(0, 1, 0), 1, material1));

    auto material2 = make_shared<lambertian>(color(0.4, 0.2, 0.1));
    auto material2_lambert = make_shared<lambert>(color(0.4, 0.2, 0.1),1.0,0.2);
    world.add(make_shared<sphere>(point3(0, -999, 0), 1000, material2_lambert));

    //world.add(make_shared<sphere>(point3(-4, 1, 0), 1, material2));
    
    auto material4 = make_shared<lambert>(color(1, 0, 0),0.5,16);
    auto material5 = make_shared<lambert>(color(0, 0, 1),0.5,64);
    auto material6 = make_shared<lambert>(color(1, 1, 0),0.7,64);
    auto material3 = make_shared<metal>(color(1.0, 0.9, 0.5), 0.7,0.2,56);
    //world.add(make_shared<sphere>(point3(4, 1, 0), 1, material4));
    //world.add(make_shared<sphere>(point3(-4, 1, 0), 1, material5));
    //world.add(make_shared<sphere>(point3(2, 5, 0), 0.5, material3));
    

    double diff=0.35;
    double spec=128;
    auto john = make_shared<lambert>(Brown, diff, spec);
	auto mary = make_shared<lambert>(color(0.52, 0.8, 0.91), diff, spec);
	auto angel = make_shared<lambert>(color(0.7,0.7,0.7), diff, spec);
	auto jesus = make_shared<lambert>(White, diff, spec);
	auto king1 = make_shared<lambert>(Red, diff, spec);
	auto king2 = make_shared<lambert>(Green, diff , spec);
	auto king3 = make_shared<lambert>(Blue, diff, spec);
//    world.add(make_shared<sphere>(point3(0, 3, 0), 2, king3));

    //world.add(make_shared<sphere>(point3(-0.1, 2, 4.7), 2, jesus));
    //world.add(make_shared<sphere>(point3(-4, 4, 0), 4, john));
    //world.add(make_shared<sphere>(point3(3.2, 3, 0.5), 3, mary));
    //world.add(make_shared<sphere>(point3(13, 4, 3.4 - 4 * 2.4), 4, king1));
    //world.add(make_shared<sphere>(point3(13+4*2, 4, 3.4-4*1.2), 4, king2));
    //world.add(make_shared<sphere>(point3(13+4*4, 4, 3.4), 4, king3));
    //world.add(make_shared<sphere>(point3(-16, 2.5, 1.3), 2.5, angel));

    hittable_list boxes2;
    auto white = make_shared<lambertian>(.73);
    int ns = 100;
    //for (int j = 0; j < ns; j++) {
       // boxes2.add(make_shared<sphere>(point3::random(0, 5), 1, material3));
      //  boxes2.add(make_shared<sphere>(point3::random(0, 5), 1, material3));

	//world.add(make_shared<bvh_node>(boxes2, 0, 1));

    
        auto difflight = make_shared<diffuse_light>(color(4, 4, 4));
        world.add(make_shared<xy_rect>(3, 5, 1, 3, -2, difflight));
        world.add(make_shared<xy_rect>(-3, -5, 1, 3, -2, difflight));
        boxes2.add(make_shared<sphere>(point3(0, 2, 0), 2, material3));
        boxes2.add(make_shared<sphere>(point3(0, 5, 0), 0.5, difflight));

	world.add(make_shared<bvh_node>(boxes2, 0.0, 1.0));
    
    return world;
}

int main()
{
	//generate image
    const char* backgroundfilename = "FF13.png";

    decodeOneStep(backgroundfilename);
	const char* filename = "out1.png";
	std::vector<unsigned char> image;
	const auto aspect_ratio = 1;
	const int image_width = 512;
	const int image_height = static_cast<int>(image_width / aspect_ratio);
	int total = image_width * image_height;
    const int samples_per_pixel = 20;
    const int max_depth = 20;
	//image resizing
	image.resize(image_width * image_height * 4);

    // World

//    auto world = random_scene();
//    auto world = cornell_box();
    auto world = cube_test_world();


	// Camera
    
    //point3 lookfrom(6, 75, 150);
    vec3 lookfrom = point3(30, 10, -30);
    vec3 lookat = point3(7, 7, 7);
    vec3 vup(0, 1, 0);
    auto dist_to_focus = 100;
    auto aperture = 1;

    camera cam(lookfrom, lookat, vup, 40, aspect_ratio, aperture, dist_to_focus);

    color defaultColor = Blue;

    //soft shadows
    const size_t elements = samples_per_pixel;
    std::vector<float> px(elements);
    std::vector<float> py(elements);
    std::vector<float> sx(elements);
    std::vector<float> sy(elements);
    std::uniform_real_distribution<float> distribution(-1.0f, 1.0f); 
    std::mt19937 engine; // Mersenne twister MT19937
    auto generator = std::bind(distribution, engine);
    std::generate_n(px.begin(), elements, generator);
    std::generate_n(sx.begin(), elements, generator);
    std::generate_n(py.begin(), elements, generator);
    std::generate_n(sy.begin(), elements, generator);
    auto rng = std::default_random_engine{};
    std::shuffle(std::begin(sx), std::end(sx), rng);
    std::shuffle(std::begin(sy), std::end(sy), rng);
    //
    /*
	for (unsigned y = image_height - 1; y >0; y--)
	{
		for (unsigned x = 0; x < image_width; x++) 
        {
            color col(0, 0, 0);
            color defaultColor(0.3, 0.5, 0.7);
			int value= image_width * (image_height - y) + x;

			for (int s = 0; s < samples_per_pixel; ++s) {
				auto u = (x + px[s]) / (image_width - 1);
				auto v = (y + py[s]) / (image_height - 1);

				ray ra = cam.get_ray_perspective(u, v);
				col += ray_color(ra, world, max_depth,defaultColor,sx[s],sy[s]);

            }

            
            //std::cout << std::endl<<col;

#pragma region Image File Color Translation
			//conversion to unsigned char
            auto scale = 1.0 / samples_per_pixel;
            auto r = col.x();
            auto g = col.y();
            auto b = col.z();
            r = sqrt(scale * r);
            g = sqrt(scale * g);
            b = sqrt(scale * b);
            
			const unsigned char rf = static_cast<unsigned char>(std::min(1., r) * 255);
			const unsigned char gf = static_cast<unsigned char>(std::min(1., g) * 255);
			const unsigned char bf = static_cast<unsigned char>(std::min(1., b) * 255);
			int val = 4 * image_width * (image_height - y) + 4 * x;
            //allot color
			//int val = 4 * image_width * (image_height- y) + 4 * x;
			image[val + 0] = rf;
			image[val + 1] = gf;
			image[val + 2] = bf;
			image[val + 3] = 255;//alpha
#pragma endregion


        }
        if (y % 10 == 0)
        {
            system("cls");
            std::cout <<"Progress: " << int(100 * (double(image_height - y) / image_height)) << "%" << std::endl;
        }
	}
	encodeOneStep(filename, image, image_width, image_height);
    */
    Matrix44<double> m = ScaleMatrix(vec3(2, 4, 8)) * TranslationMatrix(1, 2, 3);

    for(int i=0;i<4;i++)
    {
	    for(int j=0;j<4;j++)
	    {
            std::cout << m[i][j] << " ";

	    }    std::cout << std::endl;

    }


	return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
