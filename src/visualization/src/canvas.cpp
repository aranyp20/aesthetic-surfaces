#include "canvas.h"
#include "curvaturecalculator.h"
#include "organize.hpp"
#include "settings.h"
#include "type_converter.h"
#include <OpenGL/OpenGL.h>
#include <OpenGL/gl.h>
#include <QPainter>
#include <QtCore/qnamespace.h>


//#define SPLINE_RENDER

Canvas::Canvas(QWidget *parent) : QOpenGLWidget(parent)
{
   

}

Canvas::~Canvas()
{}

//TODO Fix it, its not yaw-pitch-roll
Eigen::Matrix4d Canvas::modelRotMatrix() const
{

  Eigen::Matrix4d r_z;
  r_z << std::cos(model_yaw), -std::sin(model_yaw), 0, 0,
          std::sin(model_yaw), std::cos(model_yaw), 0, 0,
          0,0,1,0,
          0,0,0,1;

  Eigen::Matrix4d r_x;
  r_x << 1,0,0,0,
        0, std::cos(model_roll), -std::sin(model_roll), 0,
        0, std::sin(model_roll), std::cos(model_roll), 0,
        0,0,0,1;

  Eigen::Matrix4d r_y;
  r_y << std::cos(model_pitch), 0, std::sin(model_pitch), 0,
          0,1,0,0,
          -std::sin(model_pitch),0, std::cos(model_pitch),0,
          0,0,0,1;

  return r_x * r_y * r_z;
}

void Canvas::changeYaw(double diff)
{
  model_yaw += diff;
}

void Canvas::changePitch(double diff)
{
  model_pitch += diff;
}

void Canvas::changeRoll(double diff)
{
  model_roll += diff;
}

void Canvas::setHighlightEdges(bool status)
{
  highlight_edges = status;
}

void Canvas::initializeGL()
{
    sp = new QOpenGLShaderProgram();
    sp->addShaderFromSourceCode(QOpenGLShader::Vertex,
                                "#version 450\n"

                                "uniform mat4  M, V, P;\n"
                                "uniform float offset;\n"


                                "in vec3 position;\n"
                                "in vec3 color;\n"
                                "out vec4 fragColor;\n"

                                "vec4 pos = vec4(position,1);\n"

                                "void main(){\n"
                                "fragColor = vec4(color,1.0);\n"
                                "pos = P * V * M * pos;\n"
                                "gl_Position = vec4(pos.x, pos.y, pos.z - offset, pos.w);\n"
                                "}");
    sp->addShaderFromSourceCode(QOpenGLShader::Fragment,
                                "#version 450\n"
                                "in vec4 fragColor;\n"
                                "out vec4 finalColor;\n"

                                "void main(){\n"
                                "finalColor = fragColor;\n"
                                "}");
    sp->link();



    vao.create();
    vao.bind();
    vbo.create();
    vbo.setUsagePattern(QOpenGLBuffer::DynamicDraw);

    vao.release();
    vbo.release();

    glEnable(GL_DEPTH_TEST);
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
}



void Canvas::setPrintable(const std::shared_ptr<common::BaseMesh> _printable_mesh)
{
  printable_mesh = _printable_mesh;
}


void Canvas::setCurvaturToHueAttributes(const common::BaseMesh& mesh, double outlier)
{
  max_color_val = 0.0;
  std::vector<double> curvatures;
  for (common::MyMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
  {
    auto& vh = *v_it;
    auto& vertexPos = mesh.point(vh);

    OpenMesh::VPropHandleT<double> myprop;
    if(!printable_mesh->get_property_handle(myprop, "demo_color")){
      //std::cout<<"Prop not found."<<std::endl;
    }
    else {
      if(printable_mesh->property(myprop, vh) > max_color_val) {
	max_color_val = printable_mesh->property(myprop, vh);
      }
      curvatures.push_back(printable_mesh->property(myprop, vh));
    }
    
  }

  if(curvatures.empty())
    return;
  
  std::sort(curvatures.begin(), curvatures.end());

  const auto start_rate = (1 - outlier) / 2;
  const auto end_rate = outlier + start_rate;
  const size_t start_index = std::floor(curvatures.size() * start_rate);
  const size_t end_index = std::ceil(curvatures.size() * end_rate);
  const auto& start_value = curvatures[start_index];
  const auto& end_value = curvatures[end_index];

  hue_divider = end_value - start_value;
  hue_offset = -(start_value / hue_divider);

  
}

std::vector<Canvas::qGlVertex> Canvas::printableMeshToLines() const //edges
{
  std::vector<qGlVertex> retval;

  if (!printable_mesh) {
    std::cout<< "Canvas::printableMeshToTriangles called while printable_mesh was nullptr."<<std::endl;
    return retval;
  }

  auto& mesh = *printable_mesh;  


  std::ofstream f("aaaa.obj");

  for(auto it : mesh.vertices()) {
    auto p = mesh.point(it);
    f << "v "<< p[0]<<" "<<p[1]<<" "<<p[2]<<std::endl;
  }
  
  for (common::MyMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it) {
    common::MyMesh::EdgeHandle edge = *e_it;

    // Get the vertices of the edge
    common::MyMesh::VertexHandle v0 = mesh.from_vertex_handle(mesh.halfedge_handle(edge, 0));
    common::MyMesh::VertexHandle v1 = mesh.to_vertex_handle(mesh.halfedge_handle(edge, 0));

    #ifdef SPLINE_RENDER
    if(v0.idx()==0 || v1.idx()==0){
      continue;
    }
    #endif
    
    common::MyMesh::Point vertex_position1 = printable_mesh->point(v0);
    common::MyMesh::Point vertex_position2 = printable_mesh->point(v1);



    retval.push_back({{vertex_position1[0], vertex_position1[1], vertex_position1[2]}, {0.0, 0.0, 1.0}});
    retval.push_back({{vertex_position2[0], vertex_position2[1], vertex_position2[2]}, {0.0, 0.0, 1.0}});


    f << "l "<< v0.idx() +1 <<" "<< v1.idx()+1 <<std::endl;



  }

  return retval;
}


std::vector<Canvas::qGlVertex> Canvas::printableFaceToTriangles(const common::MyMesh::FaceHandle& fh) const
{
  std::vector<Canvas::qGlVertex> retval;
  for (common::MyMesh::ConstFaceHalfedgeIter fhe_it = printable_mesh->cfh_iter(fh); fhe_it.is_valid(); ++fhe_it) {
      common::MyMesh::VertexHandle vh = printable_mesh->to_vertex_handle(*fhe_it);
      common::MyMesh::Point vertex_position = printable_mesh->point(vh);



      //double color = 0;
      //bool has_curvature = false;

      Eigen::Vector3d rgb_curvature(1,0,0);
      
      OpenMesh::VPropHandleT<double> myprop;
      if(printable_mesh->get_property_handle(myprop, "demo_color")){
	const auto color = printable_mesh->property(myprop, vh);
	if(color < -1e-10) {
	  rgb_curvature=Eigen::Vector3d(0,1,0);
	}
	else if(color > 1e-10) {
	  rgb_curvature=Eigen::Vector3d(0,0,1);
	}
	else {
	  rgb_curvature=Eigen::Vector3d(0,0,0);
	  //const auto val = (color + 10.0) / 20.0;
	  //rgb_curvature=Eigen::Vector3d(val,val,val);
	}
	
	//has_curvature = true;
      }
	 

          
      //const auto rgb_curvature = has_curvature ? common::color::hsvToRgb({color / hue_divider + hue_offset, 1.0, 1.0}) : Eigen::Vector3d(1.0, 0.0, 0.0);
      //const auto rgb_curvature = has_curvature ? common::color::hsvToRgb({color / max_color_val, 1.0, 1.0}) : Eigen::Vector3d(1.0, 0.0, 0.0);

      retval.push_back({{vertex_position[0], vertex_position[1], vertex_position[2]}, {static_cast<float>(rgb_curvature[0]), static_cast<float>(rgb_curvature[1]), static_cast<float>(rgb_curvature[2])}});
  }

  if (retval.size() == 4) {
    retval.push_back(retval[1]);
    retval.push_back(retval[2]);
    retval.push_back(retval[3]);
    retval.erase(retval.begin() + 2);
  }


  return retval;
}

//TODO optimize it to use TRIANGLE_STRIP
std::vector<Canvas::qGlVertex> Canvas::printableMeshToTriangles() const
{
  std::vector<Canvas::qGlVertex> retval;

  if (!printable_mesh) {
    std::cout<< "Canvas::printableMeshToTriangles called while printable_mesh was nullptr."<<std::endl;
    return retval;
  }


  common::MyMesh::FaceIter face_it, face_end(printable_mesh->faces_end());

  for (face_it = printable_mesh->faces_begin(); face_it != face_end; ++face_it)
  {
      common::MyMesh::FaceHandle fh = *face_it;

      const auto face_data = printableFaceToTriangles(fh);

      retval.insert(retval.end(), face_data.begin(), face_data.end());
  
  }

  return retval;
}


void Canvas::paintGL()
{
  if (!printable_mesh) {
    return;
  }

    const auto& v = camera.V();
    const auto& p = camera.P();
    QMatrix4x4 q_v;
    QMatrix4x4 q_p;


    const auto m = modelRotMatrix();
    QMatrix4x4 q_m;


    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        q_v(i, j) = static_cast<float>(v(i, j));
      }
    }

    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        q_p(i, j) = static_cast<float>(p(i, j));
      }
    }

    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        q_m(i, j) = static_cast<float>(m(i, j));
      }
    }
    //setCurvaturToHueAttributes(*printable_mesh, 1.0);
    std::vector<qGlVertex> pp = printableMeshToTriangles();
/*
    sp->bind();
    sp->setUniformValue("V", q_v);
    sp->setUniformValue("P", q_p);
    sp->setUniformValue("M", q_m);

    sp->setUniformValue("offset", 0.0f);


    vao.bind();
    vbo.bind();

    const void *printable_data = pp.data();

    vbo.allocate(printable_data, sizeof(qGlVertex) * pp.size());

    sp->enableAttributeArray("position");
    sp->enableAttributeArray("color");

    sp->setAttributeBuffer(0, GL_FLOAT, offsetof(qGlVertex, position), 3, sizeof(qGlVertex));
    sp->setAttributeBuffer(1, GL_FLOAT, offsetof(qGlVertex, color), 3, sizeof(qGlVertex));

*/
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    std::vector<qGlVertex> edgpoints = printableMeshToLines();
/*
    glDrawArrays(GL_TRIANGLES, 0, pp.size());

  ////////////////

    const void *printable_data2 = edgpoints.data();

    vbo.allocate(printable_data2, sizeof(qGlVertex) * edgpoints.size());

    sp->enableAttributeArray("position");
    sp->enableAttributeArray("color");

    sp->setUniformValue("offset", 0.01f);


    sp->setAttributeBuffer(0, GL_FLOAT, offsetof(qGlVertex, position), 3, sizeof(qGlVertex));
    sp->setAttributeBuffer(1, GL_FLOAT, offsetof(qGlVertex, color), 3, sizeof(qGlVertex));
    glDrawArrays(GL_LINES, 0, edgpoints.size());


*/

#ifdef SPLINE_RENDER

    glLineWidth(0.2f);
    glBegin(GL_LINES);

    glColor3f(1.0, 0.0, 0.0);
    glVertex2f(0, 0);
    glVertex2f(1, 0);

    glColor3f(0.0, 0.0, 1.0);
    glVertex2f(0, 0);
    glVertex2f(0, 1);
    
    glEnd();

#endif
    
    
    if (highlight_edges) {      
      glLineWidth(0.1f);
      glBegin(GL_LINES);
      for(const auto& line : edgpoints) {
	Eigen::Vector4d c_pos = p*v*m *  Eigen::Vector4d(line.position[0], line.position[1], line.position[2], 1.0);
	c_pos = c_pos/ c_pos[3];
	glColor3f(1.0, 1.0, 1.0);
	glVertex3f(c_pos[0], c_pos[1], c_pos[2] - 0.01f);
      }
      glEnd();
    }
    
    glBegin(GL_TRIANGLES);
    for(const auto& side : pp) {
#ifdef SPLINE_RENDER
      continue;
#endif
      Eigen::Vector4d c_pos = p*v*m *  Eigen::Vector4d(side.position[0], side.position[1], side.position[2], 1.0);
      c_pos = c_pos/ c_pos[3];
      glColor3f(side.color[0], side.color[1], side.color[2]);
      glVertex3f(c_pos[0], c_pos[1], c_pos[2]);
    }
    glEnd();


    if(common::settings::show_vertex_ids) {
      for(const auto& vertex : printable_mesh->vertices()) {
	const auto pos = common::converter::meshPointToEigen(printable_mesh->point(vertex));
	Eigen::Vector4d c_pos = p*v*m *  Eigen::Vector4d(pos[0], pos[1], pos[2], 1.0);
	c_pos = c_pos/ c_pos[3];

	int x_screen = (int)((c_pos[0] + 1.0f) * 0.5f * width());
	int y_screen = (int)((1.0f - c_pos[1]) * 0.5f * height());

	QPainter painter(this);
	glDisable(GL_DEPTH_TEST);
	painter.setPen(Qt::blue);
	painter.setFont(QFont("Arial", 10));
	painter.drawText(x_screen, y_screen, std::to_string(vertex.idx()).c_str());
	glEnable(GL_DEPTH_TEST);
	painter.end();
      }

    }
}

void Canvas::resizeGL(int w, int h)
{
}

void Canvas::mousePressEvent(QMouseEvent *event)
{

}

void Canvas::mouseMoveEvent(QMouseEvent *event)
{

}


