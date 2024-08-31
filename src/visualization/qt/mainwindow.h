#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QtWidgets/qprogressbar.h>
#include <memory>
#include <shared_mutex>
#include "canvas.h"
#include "common_defines.h"
#include "mesh.h"
#include "objectloader.h"
#include "discretefairer.h"
#include "subdivider.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();


  void cbSliceBarProgressed(int val);
  static QProgressBar* mainProgressBar;

private:
  Ui::MainWindow *ui;


  std::shared_ptr<Canvas> active_canvas;
  
  std::shared_ptr<common::MyMesh>* m_mesh = &m_mesh_a;

  std::shared_ptr<common::MyMesh> m_mesh_a;

  std::shared_ptr<common::MyMesh> m_mesh_b;

  

  framework::ObjectLoader object_loader;
  core::DiscreteFairer discrete_fairer;
  core::Subdivider subdivider;

  
  int df_subdivision_count = 0;
  int df_iteration_count = 0;


  void connectSignalsAndSlots();
  void initWidgets();

public slots:

  void yawPlus();
  void yawMinus();
  void pitchPlus();
  void pitchMinus();
  void rollPlus();
  void rollMinus();
  void setHighlightEdges(int status);
  void setShowVertexIds(int stataus);
  void loadModel();
  void changeLoadedModel(int index);
  void resetModel();
  void changeAlgorithm(int index);
  void changeMeshSlot(bool index);
  void performMethod();
  void performSubdivision();
  void setDFSubdivisionCount(int n);
  void setDFIterationCount(int n);
  void setLogAestheticAlpha(double);
  void exportModel();
  void changeCurvatureType(int index);

};
#endif // MAINWINDOW_H
