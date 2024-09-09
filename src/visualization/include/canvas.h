#pragma once

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QMouseEvent>


#include <iostream>
#include <vector>

#include "common_defines.h"
#include "camera.h"
#include "mesh.h"
#include "objectloader.h"





class Canvas : public QOpenGLWidget
{
  Q_OBJECT

  QOpenGLBuffer vbo{QOpenGLBuffer::VertexBuffer};

  QOpenGLVertexArrayObject vao;
  QOpenGLShaderProgram *sp;



  Camera camera;

  struct qGlVertex
  {
    QVector3D position;
    QVector3D color;
  };

public:

  Canvas(QWidget *parent);

  virtual ~Canvas();

  void setPrintable(const std::shared_ptr<common::BaseMesh> _printable_mesh);

  void changeYaw(double diff);
  void changePitch(double diff);
  void changeRoll(double diff);

  void setHighlightEdges(bool state);

protected:

  void initializeGL()override;

  void resizeGL(int w, int h)override;

  void paintGL()override;

  void mousePressEvent(QMouseEvent *event)override;

  void mouseMoveEvent(QMouseEvent *event)override;

private:

  std::shared_ptr<const common::BaseMesh> printable_mesh = nullptr;
  double model_yaw = 0;
  double model_pitch = 0;
  double model_roll = 0;

  Eigen::Matrix4d modelRotMatrix() const;


  double max_color_val = 0.0;
  // divides THEN offsets
  double hue_divider = 1;
  double hue_offset = 0;
  void setCurvaturToHueAttributes(const common::BaseMesh& mesh, double outlier = 0.95);
  std::vector<qGlVertex> printableFaceToTriangles(const common::MyMesh::FaceHandle& fh) const;
  std::vector<qGlVertex> printableMeshToTriangles() const; //faces
  std::vector<qGlVertex> printableMeshToLines() const; //edges


  bool highlight_edges = false;
};




