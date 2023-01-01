

#include <iostream>
#include <vector>

#include "geometry.h"
#include "model.h"
#include "our_gl.h"
#include "tgaimage.h"

Model *model = NULL;

const int width = 800;
const int height = 800;
const int depth = 255;

TGAImage shadow_buffer(width, height, TGAImage::RGB);
TGAImage use_buffer(width, height, TGAImage::RGB);

Vec3f light_dir(1, 1, 1);
Vec3f eye(0, 0, 3);
Vec3f center(0, 0, 0);
Vec3f up(0, 1, 0);

struct DepthShader : public IShader {
  mat<3, 3, float> varing_tri;
  Vec4f vertex(int iface, int nthvert, FragInfo &info) {
    Vec4f vertex = embed<4>(model->vert(iface, nthvert));
    vertex = Viewport * Projection * ModelView * vertex;
    varing_tri.set_col(nthvert, proj<3>(vertex / vertex[3]));
    return vertex;
  }
  bool fragment(Vec3f bar, TGAColor &color, const FragInfo &info) {
    Vec3f p = varing_tri * bar;
    color = TGAColor(255, 255, 255) * (p.z / depth);
    return false;
  }
};

struct GouraudShader : public IShader {
  Vec3f varying_intensity; // written by vertex shader, read by fragment shader
  Vec3f world_pos_[3];
  vec2 uvs[3];
  vec3 ndc[3];
  mat<4, 4, float> shadow_mat;
  mat<3, 3, float> varing_nrm;
  mat<3, 3, float> varying_tri;
  mat<2, 3, float> varying_uv;

  Vec4f vertex(int iface, int nthvert, FragInfo &info) {
    varying_intensity[nthvert] =
        std::max(0.f, model->normal(iface, nthvert) *
                          light_dir); // get diffuse lighting intensity
    varying_uv.set_col(nthvert, model->uv(iface, nthvert));

    varing_nrm.set_col(nthvert, model->normal(iface, nthvert));
    uvs[nthvert] = model->uv(iface, nthvert);
    Vec4f gl_Vertex =
        embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
    world_pos_[nthvert] = model->vert(iface, nthvert);
    auto t_pos = Projection * ModelView * gl_Vertex;
    // ndc[nthvert] =
    //     vec3(t_pos[0] / t_pos[3], t_pos[1] / t_pos[3], t_pos[2] /
    //     t_pos[3]);
    ndc[nthvert] = model->vert(iface, nthvert);
    gl_Vertex = Viewport * Projection * ModelView * gl_Vertex;
    varying_tri.set_col(nthvert, proj<3>(gl_Vertex / gl_Vertex[3]));
    return gl_Vertex; // transform it to screen coordinates
  }

  bool fragment(Vec3f bar, TGAColor &color, const FragInfo &info) {

    float intensity =
        varying_intensity * bar; // interpolate intensity for the current pixel

    vec3 bn = (varing_nrm * bar).normalize();
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
    //--------------------------------------------
    mat<2, 2, float> a;
    vec2 uv1 = uvs[1] - uvs[0];
    vec2 uv2 = uvs[2] - uvs[0];

    a[0][0] = uv1.x;
    a[0][1] = uv1.y;
    a[1][0] = uv2.x;
    a[1][1] = uv2.y;

    auto ai = a.invert();
    mat<2, 3, float> e;
    e[0] = ndc[1] - ndc[0];
    e[1] = ndc[2] - ndc[0];

    mat<2, 3, float> t = ai * e;
    vec3 t1 = t[0].normalize();
    vec3 t2 = t[1].normalize();
    vec3 norm = cross(t1, t2).normalize();

    mat<3, 3, float> tbn;
    tbn.set_col(0, t1);
    tbn.set_col(1, t2);
    tbn.set_col(2, bn);
    //------------------------------------------------------

    //----------------------------------------------------
    // mat<3, 3, float> A;
    // A[0] = ndc[1] - ndc[0];
    // A[1] = ndc[2] - ndc[0];
    // A[3] = bn;
    // auto AI = A.invert();
    // vec3 i = AI * vec3(uv1[0], uv2[0], 0);
    // vec3 j = AI * vec3(uv1[1], uv2[1], 0);

    //--------------------------------------------

    Vec4f sb_p = shadow_mat * embed<4>(varying_tri * bar);
    sb_p = sb_p / sb_p[3];
    float shadow_z = use_buffer.get(sb_p[0], sb_p[1])[0];
    float shadow = (sb_p[2] - shadow_z) > 0 ? 1 : 0;
    shadow = 0.3 + 0.7 * shadow;
    vec2 uv = varying_uv * bar;
    vec3 normal = bn.normalize();
    vec3 r = ((normal * (light_dir * normal) * 2.f) - light_dir);
    intensity = std::max(0.f, normal * light_dir);
    float spec = pow(std::max(r.z, 0.0f), model->specular(uv));
    float diff = std::max(0.f, normal * light_dir);
    // TGAColor c = model->diffuse(uv);
    // color = c;
    // for (int i = 0; i < 3; i++)
    //   color[i] = std::min<float>(5 + c[i] * (diff + .6 * spec), 255);
    color = TGAColor(255, 255, 255, 255) * shadow * intensity;
    return false; // no, we do not discard this pixel
  }
};

int main(int argc, char **argv) {
  if (2 == argc) {
    model = new Model(argv[1]);
  } else {
    model = new Model("obj/african_head/african_head.obj");
  }
  light_dir.normalize();
  {

    TGAImage depth(width, height, TGAImage::RGB);
    lookat(light_dir, center, up);
    viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
    projection(0);

    DepthShader depthshader;
    Vec4f screen_coords[3];
    FragInfo info;
    for (int i = 0; i < model->nfaces(); i++) {
      for (int j = 0; j < 3; j++) {
        screen_coords[j] = depthshader.vertex(i, j, info);
      }
      triangle(screen_coords, depthshader, use_buffer, shadow_buffer);
    }
    // use_buffer.flip_vertically(); // to place the origin in the bottom left
    //  corner of the image
    use_buffer.write_tga_file("use.tga");

    shadow_buffer.flip_horizontally();

    shadow_buffer.write_tga_file("shadow.tga");
  }

  Matrix M = Viewport * Projection * ModelView;
  lookat(eye, center, up);
  viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
  projection(-1.f / (eye - center).norm());

  TGAImage image(width, height, TGAImage::RGB);
  TGAImage zbuffer(width, height, TGAImage::GRAYSCALE);

  GouraudShader shader;
  shader.shadow_mat = M * (Viewport * Projection * ModelView).invert();
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