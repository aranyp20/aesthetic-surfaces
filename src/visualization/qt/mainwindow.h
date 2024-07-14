#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "objectloader.h"
#include "discretefairer.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:
    Ui::MainWindow *ui;

  std::shared_ptr<common::MyMesh> m_mesh;

    framework::ObjectLoader object_loader;
  core::DiscreteFairer discrete_fairer;

  int df_subdivision_count = 0;
  int df_iteration_count = 0;
public slots:

    void yawPlus();
    void yawMinus();
    void pitchPlus();
    void pitchMinus();
    void rollPlus();
    void rollMinus();
  void setHighlightEdges(int status);
  void changeLoadedModel(int index);

  void performMethod();
  void setDFSubdivisionCount(int n);
  void setDFIterationCount(int n);

};
#endif // MAINWINDOW_H
