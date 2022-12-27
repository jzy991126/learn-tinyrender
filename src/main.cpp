#include "model.h"
#include "tgaimage.h"
#include <fstream>

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);

const int width = 800;
const int height = 800;

struct Ivec {
  int x, y;
  Ivec operator-(const Ivec &b) { return Ivec{x - b.x, y - b.y}; }
};

void line(Ivec start, Ivec end, TGAImage &image, const TGAColor &color) {

  Ivec vec = end - start;
  bool change = false;

  if (abs(vec.y) > abs(vec.x)) {
    std::swap(start.x, start.y);
    std::swap(end.x, end.y);
    change = true;
  }
  if (start.x > end.x)
    std::swap(start, end);
  vec = end - start;
  int flag = vec.y > 0 ? 1 : -1;

  int dx = vec.x, dy = vec.y;
  int error = abs(dy) * 2;
  int y = start.y;
  int sume = 0;
  for (int x = start.x; x <= end.x; x++) {
    if (change) {
      image.set(y, x, color);
    } else {
      image.set(x, y, color);
    }
    sume += error;
    if (sume > dx) {
      y += flag;
      sume -= 2 * dx;
    }
  }
}

int main(int argc, char **argv) {

  Model *model;
  if (2 == argc) {
    model = new Model(argv[1]);
  } else {
    model = new Model("objs/african_head.obj");
  }
  TGAImage image(width, height, TGAImage::RGB);

  for (int i = 0; i < model->nfaces(); i++) {
    std::vector<int> face = model->face(i);
    for (int j = 0; j < 3; j++) {
      Vec3f v0 = model->vert(face[j]);
      Vec3f v1 = model->vert(face[(j + 1) % 3]);
      int x0 = (v0.x + 1.) * width / 2.;
      int y0 = (v0.y + 1.) * height / 2.;
      int x1 = (v1.x + 1.) * width / 2.;
      int y1 = (v1.y + 1.) * height / 2.;
      line({x0, y0}, {x1, y1}, image, white);
    }
  }
  // line({0,0},{0,50},image,red);
  // image.flip_vertically(); // i want to have the origin at the left bottom
  //  corner of the image
  image.write_tga_file("output.tga");
  return 0;
}