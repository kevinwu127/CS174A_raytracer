//
// template-rt.cpp
//

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include "matm.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits>
#include <cstring>
using namespace std;

#define MAX_DEPTH 3

int g_width;
int g_height;

struct Ray
{
    vec4 origin;
    vec4 dir;
};

struct Sphere
{
    string name;

    // Position of Sphere
    float pos_x;
    float pos_y;
    float pos_z;

    // Scaling of Sphere
    float scl_x;
    float scl_y;
    float scl_z;

    // Color
    float color_r;
    float color_g;
    float color_b;

    // K values
    float k_a;
    float k_d;
    float k_s;
    float k_r;

    // Specular Exponent
    float spec_n;
};

struct Light
{
    string name;

    // Position
    float pos_x;
    float pos_y;
    float pos_z;

    // Intensity
    float i_r;
    float i_g;
    float i_b;
};

struct Back
{
    // Background Color
    float color_r;
    float color_g;
    float color_b;
};

struct Ambient
{
    // Ambient Intensity
    float i_r;
    float i_g;
    float i_b;
};

vector<Sphere> spheres;
vector<Light> lights;
Back back;
Ambient ambient;

string output_name;

vector<vec4> g_colors;

float g_left;
float g_right;
float g_top;
float g_bottom;
float g_near;


// -------------------------------------------------------------------
// Input file parsing

vec4 toVec4(const string& s1, const string& s2, const string& s3)
{
    stringstream ss(s1 + " " + s2 + " " + s3);
    vec4 result;
    ss >> result.x >> result.y >> result.z;
    result.w = 1.0f;
    return result;
}

float toFloat(const string& s)
{
    stringstream ss(s);
    float f;
    ss >> f;
    return f;
}

void parseLine(const vector<string>& vs)
{
    int i = 0;

    if (vs[i] == "NEAR") { g_near = toFloat(vs[i + 1]); }

    else if (vs[i] == "LEFT") { g_left = toFloat(vs[i + 1]); }

    else if (vs[i] == "RIGHT") { g_right = toFloat(vs[i + 1]); }

    else if (vs[i] == "BOTTOM") { g_bottom = toFloat(vs[i + 1]); }

    else if (vs[i] == "TOP") { g_top = toFloat(vs[i + 1]); }

    else if (vs[i] == "RES")
    {
        g_width = (int)toFloat(vs[i + 1]);
        g_height = (int)toFloat(vs[i + 2]);
        g_colors.resize(g_width * g_height);
    }

    else if (vs[i] == "SPHERE")
    {
        Sphere sphere;
        sphere.name = vs[i + 1];
        sphere.pos_x = toFloat(vs[i + 2]);
        sphere.pos_y = toFloat(vs[i + 3]);
        sphere.pos_z = toFloat(vs[i + 4]);
        sphere.scl_x = toFloat(vs[i + 5]);
        sphere.scl_y = toFloat(vs[i + 6]);
        sphere.scl_z = toFloat(vs[i + 7]);
        sphere.color_r = toFloat(vs[i + 8]);
        sphere.color_g = toFloat(vs[i + 9]);
        sphere.color_b = toFloat(vs[i + 10]);
        sphere.k_a = toFloat(vs[i + 11]);
        sphere.k_d = toFloat(vs[i + 12]);
        sphere.k_s = toFloat(vs[i + 13]);
        sphere.k_r = toFloat(vs[i + 14]);
        sphere.spec_n = toFloat(vs[i + 15]);
        spheres.push_back(sphere);
    }

    else if (vs[i] == "LIGHT")
    {
        Light light;
        light.name = vs[i + 1];
        light.pos_x = toFloat(vs[i + 2]);
        light.pos_y = toFloat(vs[i + 3]);
        light.pos_z = toFloat(vs[i + 4]);
        light.i_r = toFloat(vs[i + 5]);
        light.i_g = toFloat(vs[i + 6]);
        light.i_b = toFloat(vs[i + 7]);
        lights.push_back(light);
    }

    else if (vs[i] == "BACK")
    {
        back.color_r = toFloat(vs[i + 1]);
        back.color_g = toFloat(vs[i + 2]);
        back.color_b = toFloat(vs[i + 3]);
    }

    else if (vs[i] == "AMBIENT")
    {
        ambient.i_r = toFloat(vs[i + 1]);
        ambient.i_g = toFloat(vs[i + 2]);
        ambient.i_b = toFloat(vs[i + 3]);
    }

    else if (vs[i] == "OUTPUT") { output_name = vs[i + 1]; }   
}

void loadFile(const char* filename)
{
    ifstream is(filename);
    if (is.fail())
    {
        cout << "Could not open file " << filename << endl;
        exit(1);
    }
    string s;
    vector<string> vs;
    while(!is.eof())
    {
        vs.clear();
        getline(is, s);
        istringstream iss(s);
        while (!iss.eof())
        {
            string sub;
            iss >> sub;
            vs.push_back(sub);
        }
        parseLine(vs);
    }
}


// -------------------------------------------------------------------
// Utilities

void setColor(int ix, int iy, const vec4& color)
{
    int iy2 = g_height - iy - 1; // Invert iy coordinate.
    g_colors[iy2 * g_width + ix] = color;
}

void dump() // Debugging purposes only
{
    int spheres_size = spheres.size();
    int lights_size = lights.size();

    cout << "NEAR " << g_near << endl;
    cout << "LEFT " << g_left << endl;
    cout << "RIGHT " << g_right << endl;
    cout << "BOTTOM " << g_bottom << endl;
    cout << "TOP " << g_top << endl;
    cout << "RES " << g_width << " " << g_height << endl;
    for (int i = 0; i < spheres_size; i++)
    {
        cout << "SPHERE " << spheres[i].name << " "
                          << spheres[i].pos_x << " "
                          << spheres[i].pos_y << " "
                          << spheres[i].pos_z << " "
                          << spheres[i].scl_x << " "
                          << spheres[i].scl_y << " "
                          << spheres[i].scl_z << " "
                          << spheres[i].color_r << " "
                          << spheres[i].color_g << " "
                          << spheres[i].color_b << " "
                          << spheres[i].k_a << " "
                          << spheres[i].k_d << " "
                          << spheres[i].k_s << " "
                          << spheres[i].k_r << " "
                          << spheres[i].spec_n << endl;
    }
    for (int i = 0; i < lights_size; i++)
    {
        cout << "LIGHT " << lights[i].name << " "
                         << lights[i].pos_x << " "
                         << lights[i].pos_y << " "
                         << lights[i].pos_z << " "
                         << lights[i].i_r << " "
                         << lights[i].i_g << " "
                         << lights[i].i_b << endl;
    }
    cout << "BACK " << back.color_r << " "
                    << back.color_g << " "
                    << back.color_b << endl;
    cout << "AMBIENT " << ambient.i_r << " "
                       << ambient.i_g << " "
                       << ambient.i_b << endl;
    cout << "OUTPUT " << output_name << endl;
}

bool lightIntersection(const vec4& light, const vec4& closest_hit_point, bool inside);
vec4 trace(const Ray& ray, int depth, bool inside);

vec4 light(const Ray& ray, int& sphere_num, const vec4& closest_hit_point, vec4& normal, const int depth, bool inside)
{
    int size_l = lights.size();
    vec4 diffuse_factor(0.0f,0.0f,0.0f,1.0f);
    vec4 specular_factor(0.0f,0.0f,0.0f,1.0f);
    vec4 reflected_color(0.0f,0.0f,0.0f,1.0f);
    vec4 ambient_factor(ambient.i_r * spheres[sphere_num].k_a * spheres[sphere_num].color_r,
                        ambient.i_g * spheres[sphere_num].k_a * spheres[sphere_num].color_g,
                        ambient.i_b * spheres[sphere_num].k_a * spheres[sphere_num].color_b, 1.0f);
    vec4 light;
    vec4 light_norm;
    vec4 reflection_norm;
    vec4 final_light;
    vec4 ray_norm = normalize(ray.origin - closest_hit_point);
    vec4 normal_norm;

    bool light_inside = false;

    if (inside)
        normal_norm = -normalize(vec4(normal.x, normal.y, normal.z, 0.0f));
    else
        normal_norm = normalize(vec4(normal.x, normal.y, normal.z, 0.0f));

    for (int i = 0; i < size_l; i++)
    {
        // Get the normalized light vector from the hit point to the light source
        light = vec4(lights[i].pos_x, lights[i].pos_y, lights[i].pos_z, 1.0f);
        light_norm = normalize(light - closest_hit_point);

        // Get the Reflection vector at the hit point
        // R = 2*N(N . L) - L
        reflection_norm = 2 * normal_norm * dot(normal_norm, light_norm) - light_norm;


        // Check if the light intersects any other objects
        // If true, then set shadow by not setting any illumination
        if (lightIntersection(light, closest_hit_point, inside)) { continue; }

        // Set the diffuse factor for illumination
        float n_dot_l = dot(normal_norm, light_norm);
        if (n_dot_l >= 0)
        {
            diffuse_factor = vec4(diffuse_factor.x + (lights[i].i_r * spheres[sphere_num].k_d * n_dot_l * spheres[sphere_num].color_r),
                                  diffuse_factor.y + (lights[i].i_g * spheres[sphere_num].k_d * n_dot_l * spheres[sphere_num].color_g),
                                  diffuse_factor.z + (lights[i].i_b * spheres[sphere_num].k_d * n_dot_l * spheres[sphere_num].color_b), 1.0f);

            // Set the specular factor for illumination
            float r_dot_v = dot(reflection_norm, ray_norm);
            if (r_dot_v >= 0 && !inside)
            {
                specular_factor = vec4(specular_factor.x + (lights[i].i_r * spheres[sphere_num].k_s * pow(r_dot_v, spheres[sphere_num].spec_n)),
                                       specular_factor.y + (lights[i].i_g * spheres[sphere_num].k_s * pow(r_dot_v, spheres[sphere_num].spec_n)),
                                       specular_factor.z + (lights[i].i_b * spheres[sphere_num].k_s * pow(r_dot_v, spheres[sphere_num].spec_n)), 1.0f);
            }
        }
        else { continue; }

        // Recursive reflected ray
        if (depth > 0)
        {
            // Find the reflected ray
            Ray reflected_ray;
            reflected_ray.origin = closest_hit_point;
            reflected_ray.dir = -2 * normal_norm * dot(normal_norm, ray.dir) + ray.dir;

            // Recursively find the color returned by trace
            reflected_color = trace(reflected_ray, depth, inside) * spheres[sphere_num].k_r;      
        }
        
    }

    // Add together all the values for illumination
    final_light = vec4(diffuse_factor.x + specular_factor.x + ambient_factor.x + reflected_color.x, 
                       diffuse_factor.y + specular_factor.y + ambient_factor.y + reflected_color.y, 
                       diffuse_factor.z + specular_factor.z + ambient_factor.z + reflected_color.z, 1.0f);
    return final_light;
}


// -------------------------------------------------------------------
// Intersection routine

bool solveQuadratic(const float& a, const float& b, const float& c, 
                    float &t0, float &t1)
{
    // Discriminant
    float discr = b * b - (a * c);

    if (discr < 0) { return false; }             // No Solution

    else if (discr == 0) {  t0 = t1 = -b / a; }  // One Solution

    else
    {
                                                 // Two Solutions    
        float q = (b > 0) ?
            -(b + sqrt(discr)) :
            -(b - sqrt(discr));
        t0 = q / a;
        t1 = c / q; 
    }
    
    if (t0 > t1) { swap(t0, t1); } // Make sure t0 is where the ray enters sphere

    return true;
}

vec3 toVec3(const vec4& v4)
{
    vec3 v3;
    v3.x = v4.x;
    v3.y = v4.y;
    v3.z = v4.z;
    return v3;
}

void invert(const Ray& ray, const mat4& matrix, Ray& inverse_ray, mat4& inverse_matrix)
{
     // Find the inverse transformed ray & the inverse matrix
    if (InvertMatrix(matrix, inverse_matrix))
    {
        inverse_ray.origin = inverse_matrix * ray.origin;
        inverse_ray.dir = inverse_matrix * ray.dir;
    }
    else
    {
        cout << "Matrix not invertible" << endl;
        exit(1);
    }
}


float intersection(const Ray& ray, const mat4& matrix, const int depth)
{
    Ray inverse_ray;
    mat4 inverse_matrix;

    invert(ray, matrix, inverse_ray, inverse_matrix);

    // Convert to vec3 to use in quadratic equation
    vec3 S_prime = toVec3(inverse_ray.origin);
    vec3 c_prime = toVec3(inverse_ray.dir);

    // Get A, B, and C values for quadratic
    float a = dot(c_prime, c_prime);
    float b = dot(S_prime, c_prime);
    float c = dot(S_prime, S_prime) - 1;

    // Intersection Points
    float t0;
    float t1;

    // Value to check depending on depth of recursion
    float epsilon;
    if (depth == MAX_DEPTH - 1) { epsilon = 1.0f; }
    else { epsilon = 0.001f; }

    // Quadratic Equation
    if (solveQuadratic(a, b, c, t0, t1)) { return t0 < t1 && t0 > epsilon ? t0 : t1; }
        
    else { return -1; } // No Intersections
}

bool lightIntersection(const vec4& light, const vec4& closest_hit_point, bool inside)
{
    vec4 inverse_origin;
    vec4 inverse_ray;
    mat4 matrix;
    mat4 inverse_matrix;

    int size = spheres.size();

    // Find the vector starting from the hit point to the light source
    vec4 light_vector = light - closest_hit_point;


    // Iterate through all the spheres
    for (int i = 0; i < size; i++)
    {
        // Generate Transformation Matrix
        matrix = Translate(spheres[i].pos_x, spheres[i].pos_y, spheres[i].pos_z);
        matrix *= Scale(spheres[i].scl_x, spheres[i].scl_y, spheres[i].scl_z);

        // Invert the matrix 
        if (!InvertMatrix(matrix, inverse_matrix))
        {
            cout << "Matrix not invertible" << endl;
            exit(1);
        }

        // Find the start point and vector for the unit sphere
        inverse_origin = inverse_matrix * closest_hit_point;
        inverse_ray    = inverse_matrix * light_vector;

        // Convert to vec3 to use in quadratic equation
        vec3 S_prime = toVec3(inverse_origin);
        vec3 c_prime = toVec3(inverse_ray);

        // Get A, B, and C values for quadratic
        float a = dot(c_prime, c_prime);
        float b = dot(S_prime, c_prime);
        float c = dot(S_prime, S_prime) - 1;

        // Intersection Points
        float t0;
        float t1;

        // The bottom limit for checking for intersection
        float epsilon = 0.0001f;

        if (solveQuadratic(a, b, c, t0, t1))
        {
            // returns true if there IS an intersection
            if( (t0 < 1 && t0 > epsilon) || (t1 < 1 && t1 > epsilon) )
                 return true;
        }
    }

    // if no intersections with ANY sphere, then return false
    return false;
}

vec4 getNormal(const Ray& ray, const mat4& matrix, const float& time)
{
    Ray inverse_ray;
    mat4 inverse_matrix;
    vec4 inverse_hit_point;
    vec4 inverse_normal;

    // Find the inverse ray and matrix
    invert(ray, matrix, inverse_ray, inverse_matrix);

    // Get the hit point on the unit sphere
    inverse_hit_point = inverse_ray.origin + inverse_ray.dir * time;

    // Acquire the normal vector on the unit sphere
    inverse_normal = vec4(inverse_hit_point.x, inverse_hit_point.y, inverse_hit_point.z, 0.0f);

    // Return the transformed normal vector
    return transpose(inverse_matrix) * inverse_normal;
}


// -------------------------------------------------------------------
// Ray tracing

vec4 trace(const Ray& ray, int depth, bool inside)
{
    mat4 matrix; 
    vec4 hit_point;
    vec4 closest_hit_point;
    vec4 normal;
    vec4 color;

    int size_s = spheres.size();
    float closest_z = 5000;
    bool intersect_exists = false;
    int sphere_num;
    float time;


    // Decrement depth for recursion
    depth--;

    // Iterate through every sphere
    for (int i = 0; i < size_s; i++)
    {
        // Generate Transformation Matrix
        matrix = Translate(spheres[i].pos_x, spheres[i].pos_y, spheres[i].pos_z);
        matrix *= Scale(spheres[i].scl_x, spheres[i].scl_y, spheres[i].scl_z);


            if ((spheres[i].pos_z + spheres[i].scl_z) > 0)
                inside = true;
            else
                inside = false;
        

        // Find first intersection
        time = intersection(ray, matrix, depth);

        // Calculate the hit point on the transformed sphere
        hit_point = ray.origin + (time * ray.dir);

        // If the hit point is the closest hit point in the z direction
        // And the time value is greater than 0.001f
        if (abs(hit_point.z - ray.origin.z) < closest_z && time > 0.0001f)
        {
            // If at the top layer of recursion, don't acquire any color for
            // spheres in front of the near plane
            if (depth == MAX_DEPTH - 1 && time < 1) { continue; }

            // Set the new closest hit point;
            closest_z = abs(hit_point.z - ray.origin.z);
            closest_hit_point = hit_point;

            // Get the normal of this hit point
            normal = getNormal(ray, matrix, time);

            // Set the sphere number for use in lighting function
            sphere_num = i;
            intersect_exists = true;

            // Find the color to be used on this pixel
            color = light(ray, sphere_num, closest_hit_point, normal, depth, inside);
        
        }
    }
    // If there is no intersect, then set the background color to the input value
    if (!intersect_exists && depth == MAX_DEPTH - 1) { return vec4(back.color_r, back.color_g, back.color_b, 1.0f); }

    // Otherwise, return the color vector found in the previous for loop
    else { return color; }
    
}

vec4 getDir(int ix, int iy)
{
    float x = g_left + ((1.*ix) / (1.*g_width)) * (g_right - g_left);
    float y = g_bottom + ((1.*iy) / (1.*g_height)) * (g_top - g_bottom);
    float z = -g_near;
    return vec4(x, y, z, 0.0f);
}

void renderPixel(int ix, int iy)
{
    Ray ray;
    ray.origin = vec4(0.0f, 0.0f, 0.0f, 1.0f);
    ray.dir = getDir(ix, iy);
    bool inside = false;
    vec4 color = trace(ray, MAX_DEPTH, inside);
    setColor(ix, iy, color);
}

void render()
{
    for (int iy = 0; iy < g_height; iy++)
        for (int ix = 0; ix < g_width; ix++)
            renderPixel(ix, iy);
}


// -------------------------------------------------------------------
// PPM saving

void savePPM(int Width, int Height, char* fname, unsigned char* pixels) 
{
    FILE *fp;
    const int maxVal=255;

    printf("Saving image %s: %d x %d\n", fname, Width, Height);
    fp = fopen(fname,"wb");
    if (!fp) {
        printf("Unable to open file '%s'\n", fname);
        return;
    }
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", Width, Height);
    fprintf(fp, "%d\n", maxVal);

    for(int j = 0; j < Height; j++) {
        fwrite(&pixels[j*Width*3], 3, Width, fp);
    }

    fclose(fp);
}

void saveFile()
{
    // Convert color components from floats to unsigned chars.
    unsigned char* buf = new unsigned char[g_width * g_height * 3];
    for (int y = 0; y < g_height; y++)
        for (int x = 0; x < g_width; x++)
            for (int i = 0; i < 3; i++)
            {
                if (g_colors[y*g_width + x][i] > 1.0) { g_colors[y*g_width + x][i] = 1.0f; }
                if (g_colors[y*g_width + x][i] < 0.0) { g_colors[y*g_width + x][i] = 0.0f; }
                buf[y*g_width * 3 + x * 3 + i] = (unsigned char)(((float*)g_colors[y*g_width + x])[i] * 255.9f);
            }
    

    char * name = new char[output_name.length() + 1];
    strcpy(name, output_name.c_str());
    savePPM(g_width, g_height, name, buf);
    delete[] buf;
    delete[] name;
}


// -------------------------------------------------------------------
// Main

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        cout << "Usage: template-rt <input_file.txt>" << endl;
        exit(1);
    }
    loadFile(argv[1]);
    render();
    saveFile();
	return 0;
}

