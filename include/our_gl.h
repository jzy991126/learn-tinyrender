#pragma once

#include "geometry.h"
#include "tgaimage.h"

extern Matrix ModelView;
extern Matrix Viewport;
extern Matrix Projection;

struct FragInfo {
  vec3 world_pos_;
};

void viewport(int x, int y, int w, int h);
void projection(float coeff = 0.f); // coeff = -1/c
void lookat(Vec3f eye, Vec3f center, Vec3f up);

struct IShader {
  virtual ~IShader();
  virtual Vec4f vertex(int iface, int nthvert, FragInfo &info) = 0;
  virtual bool fragment(Vec3f bar, TGAColor &color, const FragInfo &info) = 0;
};

void triangle(Vec4f *pts, IShader &shader, TGAImage &image, float *zbuffer);
