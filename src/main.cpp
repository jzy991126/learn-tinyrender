// #include "geometry.h"
// #include "model.h"
// #include "tgaimage.h"
// #include <algorithm>
// #include <fstream>
// #include <limits>

// const TGAColor white = TGAColor(255, 255, 255, 255);
// const TGAColor red = TGAColor(255, 0, 0, 255);

// const int width = 800;
// const int height = 800;
// const int depth = 255;

// const vec3 eye(1, 1, 3);
// const vec3 center(0, 0, 0);

// // struct Ivec2 {
// //   Ivec2() = default;
// //   Ivec2(int x, int y) : x(x), y(y) {}
// //   int x, y;
// //   Ivec2 operator-(const Ivec2 &b) const { return Ivec2{x - b.x, y - b.y};
// }
// // };

// Vec3f barycentric(const Vec3f &a, const Vec3f &b, const Vec3f &c,
//                   const Vec3f &p) {
//   Vec3f ac = c - a, ab = b - a, pa = a - p;
//   Vec3f x(ab.x, ac.x, pa.x), y(ab.y, ac.y, pa.y);
//   Vec3f uv = cross(x, y);
//   if (abs(uv.z) < 1e-5) {
//     return Vec3f(-1, -1, -1);
//   } else {
//     uv = uv / uv.z;
//     return Vec3f(1 - uv.x - uv.y, uv.x, uv.y);
//   }
// }

// // void line(Ivec2 start, Ivec2 end, TGAImage &image, const TGAColor &color)
// {

// //   Ivec2 vec = end - start;
// //   bool change = false;

// //   if (abs(vec.y) > abs(vec.x)) {
// //     std::swap(start.x, start.y);
// //     std::swap(end.x, end.y);
// //     change = true;
// //   }
// //   if (start.x > end.x)
// //     std::swap(start, end);
// //   vec = end - start;
// //   int flag = vec.y > 0 ? 1 : -1;

// //   int dx = vec.x, dy = vec.y;
// //   int error = abs(dy) * 2;
// //   int y = start.y;
// //   int sume = 0;
// //   for (int x = start.x; x <= end.x; x++) {
// //     if (change) {
// //       image.set(y, x, color);
// //     } else {
// //       image.set(x, y, color);
// //     }
// //     sume += error;
// //     if (sume > dx) {
// //       y += flag;
// //       sume -= 2 * dx;
// //     }
// //   }
// // }

// // void mtriangle(Ivec2 &a, Ivec2 &b, Ivec2 &c, TGAImage &image,
// //                const TGAColor &color) {
// //   if (a.y > b.y)
// //     std::swap(a, b);
// //   if (a.y > c.y)
// //     std::swap(a, c);
// //   if (b.y > c.y)
// //     std::swap(b, c);
// //   Ivec2 vec1 = c - a, vec2 = b - a, vec3 = c - b;
// //   float k1 = (float)vec1.x / (vec1.y), k2 = (float)vec2.x / (vec2.y),
// //         k3 = (float)vec3.x / (vec3.y);
// //   float sx = a.x, ex = a.x;
// //   for (int y = a.y; y < b.y; y++) {
// //     line({int(sx), y}, {int(ex), y}, image, color);
// //     sx += k1, ex += k2;
// //   }
// //   for (int y = b.y; y < c.y; y++) {
// //     line({int(sx), y}, {int(ex), y}, image, color);
// //     sx += k1, ex += k3;
// //   }
// // }
// void triangle(const Vec3f &a, const Vec3f &b, const Vec3f &c, TGAImage
// &image,
//               float *zbuffer, Vec2f uva, Vec2f uvb, Vec2f uvc,
//               const TGAImage &tex, float intense) {
//   int minx = std::min(std::min(a.x, b.x), c.x);
//   minx = std::max(0, minx);
//   int miny = std::min(std::min(a.y, b.y), c.y);
//   miny = std::max(0, miny);
//   int maxx = std::max(std::max(a.x, b.x), c.x);
//   maxx = std::min(width - 1, maxx);
//   int maxy = std::max(std::max(a.y, b.y), c.y);
//   maxy = std::min(height - 1, maxy);

//   for (int x = minx; x <= maxx; x++) {
//     for (int y = miny; y <= maxy; y++) {
//       Vec3f p(x, y, 0);
//       Vec3f baycenter = barycentric(a, b, c, p);
//       if (baycenter.x < 0 || baycenter.y < 0 || baycenter.z < 0) {
//         continue;
//       }
//       p.z = baycenter.x * a.z + baycenter.y * b.z + baycenter.z * c.z;
//       Vec2f uv = uva * baycenter.x + uvb * baycenter.y + uvc * baycenter.z;
//       if (zbuffer[y * width + x] < p.z) {
//         zbuffer[y * width + x] = p.z;
//         TGAColor color = tex.get(uv.x * tex.width(), uv.y * tex.height());
//         image.set(x, y,
//                   TGAColor(intense * color[2], intense * color[1],
//                            intense * color[0]));
//       }
//     }
//   }
// }

// Matrix v2m(Vec3f v) {
//   Matrix m;
//   m[0][0] = v.x;
//   m[1][0] = v.y;
//   m[2][0] = v.z;
//   m[3][0] = 1.f;
//   return m;
// }
// Vec3f m2v(Matrix m) {
//   return Vec3f(m[0][0] / m[3][0], m[1][0] / m[3][0], m[2][0] / m[3][0]);
// }

// Matrix viewport(int x, int y, int w, int h) {
//   Matrix m = Matrix::identity();
//   m[0][3] = x + w / 2.f;
//   m[1][3] = y + h / 2.f;
//   m[2][3] = depth / 2.f;

//   m[0][0] = w / 2.f;
//   m[1][1] = h / 2.f;
//   m[2][2] = depth / 2.f;
//   return m;
// }

// Matrix look_at(Vec3f eye, Vec3f center, Vec3f up) {
//   Vec3f z = (eye - center).normalize();
//   vec3 x = cross(up, z).normalize();
//   vec3 y = cross(z, x).normalize();
//   Matrix rot = Matrix::identity();
//   Matrix trans = Matrix::identity();
//   for (int i = 0; i < 3; i++) {
//     rot[0][i] = x[i];
//     rot[1][i] = y[i];
//     rot[2][i] = z[i];
//     trans[i][3] = -eye[i];
//   }
//   return rot * trans;
// }

// int main(int argc, char **argv) {

//   Model *model;
//   if (2 == argc) {
//     model = new Model(argv[1]);
//   } else {
//     model = new Model("obj/african_head/african_head.obj");
//   }
//   TGAImage image(width, height, TGAImage::RGB);
//   TGAImage texture;
//   texture.read_tga_file("objs/african_head_diffuse.tga");

//   //   for (int i = 0; i < model->nfaces(); i++) {
//   //     std::vector<int> face = model->face(i);
//   //     for (int j = 0; j < 3; j++) {
//   //       Vec3f v0 = model->vert(face[j]);
//   //       Vec3f v1 = model->vert(face[(j + 1) % 3]);
//   //       int x0 = (v0.x + 1.) * width / 2.;
//   //       int y0 = (v0.y + 1.) * height / 2.;
//   //       int x1 = (v1.x + 1.) * width / 2.;
//   //       int y1 = (v1.y + 1.) * height / 2.;
//   //       line({x0, y0}, {x1, y1}, image, white);
//   //     }
//   //   }
//   float *zbuffer = new float[width * height];

//   std::fill(zbuffer, zbuffer + width * height,
//             -std::numeric_limits<float>::max());
//   Vec3f light_dir(0, 0, -1);
//   for (int i = 0; i < model->nfaces(); i++) {
//     Vec3f screen_coords[3];
//     Vec3f world_coords[3];
//     Vec2f uvs[3];
//     Matrix Projection = Matrix::identity();
//     Matrix ViewPort = viewport(0, 0, width, height);
//     Matrix ModelView = look_at(eye, center, Vec3f(0, 1, 0));
//     Projection[3][2] = -1.f / (eye - center).norm();
//     for (int j = 0; j < 3; j++) {
//       Vec3f world_coord = model->vert(i, j);
//       world_coords[j] = world_coord;
//       uvs[j] = model->uv(i, j);
//       //   screen_coords[j] =
//       //       Vec3f((world_coord.x + 1.) * width / 2.,
//       //             (world_coord.y + 1.) * height / 2., world_coord.z);
//       auto s = ViewPort * Projection;
//       screen_coords[j] =
//           m2v(ViewPort * Projection * ModelView * v2m(world_coord));
//     }
//     Vec3f normal = cross((world_coords[2] - world_coords[0]),
//                          (world_coords[1] - world_coords[0]));
//     normal.normalize();
//     float intens = normal * light_dir;
//     if (intens > 0) {
//       triangle(screen_coords[0], screen_coords[1], screen_coords[2], image,
//                zbuffer, uvs[0], uvs[1], uvs[2], texture, intens);
//     }
//   }

//   // triangle(a, b, c, image, white);
//   image.write_tga_file("output.tga");
//   return 0;
// }

#include <iostream>
#include <vector>

#include "geometry.h"
#include "model.h"
#include "our_gl.h"
#include "tgaimage.h"

Model *model = NULL;
const int width = 800;
const int height = 800;

Vec3f light_dir(1, 1, 1);
Vec3f eye(1, 1, 3);
Vec3f center(0, 0, 0);
Vec3f up(0, 1, 0);

struct GouraudShader : public IShader {
  Vec3f varying_intensity; // written by vertex shader, read by fragment shader
  Vec3f world_pos_[3];
  mat<2, 3, float> varying_uv;

  Vec4f vertex(int iface, int nthvert, FragInfo &info) {
    varying_intensity[nthvert] =
        std::max(0.f, model->normal(iface, nthvert) *
                          light_dir); // get diffuse lighting intensity
    varying_uv.set_col(nthvert, model->uv(iface, nthvert));
    Vec4f gl_Vertex =
        embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
    world_pos_[nthvert] = model->vert(iface, nthvert);
    return Viewport * Projection * ModelView *
           gl_Vertex; // transform it to screen coordinates
  }

  bool fragment(Vec3f bar, TGAColor &color, const FragInfo &info) {

    float intensity =
        varying_intensity * bar; // interpolate intensity for the current pixel
    float inten = 1.0;
    if (intensity > .85f) {
      inten = 1.f;
    } else if (intensity > 0.65f) {
      inten = 0.8;
    } else if (intensity > 0.45f) {
      inten = 0.6;
    } else if (intensity > 0.25f) {
      inten = 0.3;
    } else
      inten = 0.0f;

    vec2 uv = varying_uv * bar;
    vec3 normal = model->normal(uv).normalize();
    vec3 r = ((normal * (light_dir * normal) * 2.f) - light_dir);
    intensity = std::max(0.f, normal * light_dir);
    float spec = pow(std::max(r.z, 0.0f), model->specular(uv));
    float diff = std::max(0.f, normal * light_dir);
    TGAColor c = model->diffuse(uv);
    color = c;
    for (int i = 0; i < 3; i++)
      color[i] = std::min<float>(5 + c[i] * (diff + .6 * spec), 255);
    return false; // no, we do not discard this pixel
  }
};

int main(int argc, char **argv) {
  if (2 == argc) {
    model = new Model(argv[1]);
  } else {
    model = new Model("obj/african_head/african_head.obj");
  }

  lookat(eye, center, up);
  viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
  projection(-1.f / (eye - center).norm());
  light_dir.normalize();

  TGAImage image(width, height, TGAImage::RGB);
  TGAImage zbuffer(width, height, TGAImage::GRAYSCALE);

  GouraudShader shader;
  for (int i = 0; i < model->nfaces(); i++) {
    Vec4f screen_coords[3];
    FragInfo info;
    for (int j = 0; j < 3; j++) {
      screen_coords[j] = shader.vertex(i, j, info);
    }
    triangle(screen_coords, shader, image, zbuffer);
  }

  image.flip_vertically(); // to place the origin in the bottom left corner
  // of
  //  the image
  zbuffer.flip_vertically();
  image.write_tga_file("output.tga");
  zbuffer.write_tga_file("zbuffer.tga");

  delete model;
  return 0;
}