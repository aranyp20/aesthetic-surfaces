#include "mainwindow.h"
#include "common_defines.h"
#include "mesh.h"
#include "settings.h"
#include "ui_mainwindow.h"

#include "canvas.h"
#include <QtCore/qglobal.h>
#include <QtCore/qobject.h>
#include <QtCore/qstring.h>
#include <QtCore/qvariant.h>
#include <QtWidgets/qcombobox.h>
#include <QtWidgets/qprogressbar.h>
#include <QtWidgets/qpushbutton.h>
#include <QtWidgets/qradiobutton.h>
#include <QtWidgets/qspinbox.h>
#include <memory>
#include "parametriclogaesthetic1.h"
#include "logaesthetic_spline.h"

QProgressBar *MainWindow::mainProgressBar = nullptr;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    //ui->main_openGL_widget = ui->openGLWidget;
    ui->openGLWidget = ui->main_openGL_widget;
    MainWindow::mainProgressBar = ui->algProgressBar;

    connectSignalsAndSlots();


    initWidgets();


    //loadModel();
}

MainWindow::~MainWindow()
{
  delete ui;
}

void MainWindow::connectSignalsAndSlots()
{
  QObject::connect(ui->modelSelectionComboBox, qOverload<int>(&QComboBox::currentIndexChanged), this, &MainWindow::changeLoadedModel);
  QObject::connect(ui->algSelectorComboBox, qOverload<int>(&QComboBox::currentIndexChanged), this, &MainWindow::changeAlgorithm);
  QObject::connect(ui->curvatureSelectorComboBox, qOverload<int>(&QComboBox::currentIndexChanged), this, &MainWindow::changeCurvatureType);
  QObject::connect(ui->exportModelPushButton, &QPushButton::pressed, this, &MainWindow::exportModel);
  QObject::connect(ui->loadModelPushButton, &QPushButton::pressed, this, &MainWindow::loadModel);
  QObject::connect(ui->logAlphaSpinBox, qOverload<double>(&QDoubleSpinBox::valueChanged), this, &MainWindow::setLogAestheticAlpha);
  QObject::connect(ui->radioButton, qOverload<bool>(&QRadioButton::toggled), this, &MainWindow::changeMeshSlot);
  QObject::connect(ui->yawplus, &QPushButton::pressed, this, &MainWindow::yawPlus);
  QObject::connect(ui->yawminus, &QPushButton::pressed, this, &MainWindow::yawMinus);
  QObject::connect(ui->pitchplus, &QPushButton::pressed, this, &MainWindow::pitchPlus);
  QObject::connect(ui->pitchminus, &QPushButton::pressed, this, &MainWindow::pitchMinus);
  QObject::connect(ui->rollplus, &QPushButton::pressed, this, &MainWindow::rollPlus);
  QObject::connect(ui->rollminus, &QPushButton::pressed, this, &MainWindow::rollMinus);
  QObject::connect(ui->highlightEdgesCheckBox, qOverload<int>(&QCheckBox::stateChanged), this, &MainWindow::setHighlightEdges);
  QObject::connect(ui->subdivisionCountBox, qOverload<int>(&QSpinBox::valueChanged), this, &MainWindow::setDFSubdivisionCount);
  QObject::connect(ui->iterationCountBox, qOverload<int>(&QSpinBox::valueChanged), this, &MainWindow::setDFIterationCount);
  QObject::connect(ui->executeButton, &QPushButton::pressed, this, &MainWindow::performMethod);
  QObject::connect(ui->subdivideButton, &QPushButton::pressed, this, &MainWindow::performSubdivision);
  QObject::connect(ui->resetModelButton, &QPushButton::pressed, this, &MainWindow::resetModel);
  QObject::connect(ui->vertexIdsCheckBox, qOverload<int>(&QCheckBox::stateChanged), this, &MainWindow::setShowVertexIds);
  QObject::connect(ui->logGenPushButton, &QPushButton::pressed, this, &MainWindow::generateLogAesthetic);
}

void MainWindow::initWidgets()
{
  ui->logAlphaSpinBox->setMinimum(-5000);
  ui->logAlphaSpinBox->setMaximum(5000);
  ui->logAlphaSpinBox->setSingleStep(0.1);
  ui->logAlphaSpinBox->setValue(common::settings::log_aesthetic_alpha);

  const auto file_options = object_loader.loadFileOptions();

  size_t tmpi = 0;
  for (const auto& file_option : file_options) {
    ui->modelSelectionComboBox->addItem(QString::fromStdString(file_option), QVariant((int)tmpi++));
    if (file_option == common::settings::default_model) {
      ui->modelSelectionComboBox->setCurrentIndex(tmpi-1);
    }
  }

  ui->algSelectorComboBox->addItem(QString::fromStdString(common::settings::alg_name.at(common::settings::BASIC)), QVariant((int)common::settings::BASIC));
  ui->algSelectorComboBox->addItem(QString::fromStdString(common::settings::alg_name.at(common::settings::BEZIER)), QVariant((int)common::settings::BEZIER));
  ui->algSelectorComboBox->addItem(QString::fromStdString(common::settings::alg_name.at(common::settings::LOG_AESTHETIC)), QVariant((int)common::settings::LOG_AESTHETIC));

  ui->curvatureSelectorComboBox->addItem(QString::fromStdString(common::settings::curvature_name.at(common::settings::MEAN)), QVariant((int)common::settings::MEAN));
  ui->curvatureSelectorComboBox->addItem(QString::fromStdString(common::settings::curvature_name.at(common::settings::GAUSSIAN)), QVariant((int)common::settings::GAUSSIAN));

  ui->curvatureSelectorComboBox->setCurrentIndex(common::settings::selected_curvature);

  ui->algSelectorComboBox->setCurrentIndex(2);
  
  ui->radioButton->setChecked(true);
  ui->subdivisionCountBox->setValue(1);
  ui->iterationCountBox->setValue(1);
  ui->highlightEdgesCheckBox->setChecked(true);
  ui->algProgressBar->setValue(0);
  ui->vertexIdsCheckBox->setChecked(common::settings::show_vertex_ids);



  //ui->exportModelNameLineEdit->setPlaceholderText("filename.txt");

}

void MainWindow::yawPlus()
{
    ui->main_openGL_widget->changeYaw(0.2);
    ui->main_openGL_widget->update();
}
void MainWindow::yawMinus()
{
    ui->main_openGL_widget->changeYaw(-0.2);
    ui->main_openGL_widget->update();
}

void MainWindow::pitchPlus()
{
    ui->main_openGL_widget->changePitch(0.2);
    ui->main_openGL_widget->update();
}
void MainWindow::pitchMinus()
{
    ui->main_openGL_widget->changePitch(-0.2);
    ui->main_openGL_widget->update();
}

void MainWindow::rollPlus()
{
    ui->main_openGL_widget->changeRoll(0.2);
    ui->main_openGL_widget->update();
}
void MainWindow::rollMinus()
{
    ui->main_openGL_widget->changeRoll(-0.2);
    ui->main_openGL_widget->update();
}

void MainWindow::setHighlightEdges(int status)
{
  if(status) {
    ui->main_openGL_widget->setHighlightEdges(true);
  }
  else {
    ui->main_openGL_widget->setHighlightEdges(false);    
  }
  ui->main_openGL_widget->update();

}

void MainWindow::loadModel()
{
  *m_mesh = object_loader.loadFromFile( ui->modelSelectionComboBox->currentIndex());
  ui->main_openGL_widget->setPrintable(*m_mesh);
  ui->main_openGL_widget->update();
}

void MainWindow::resetModel()
{
  if(!*m_mesh) {
    std::cout << "No loaded mesh."<< std::endl;
    return;
  }
  *m_mesh = std::make_shared<common::MyMesh>((**m_mesh).original_state);
  ui->main_openGL_widget->setPrintable(*m_mesh);
  ui->main_openGL_widget->update();
}

void MainWindow::performMethod()
{
  discrete_fairer.execute(**m_mesh, df_iteration_count, [this](int v){this->cbSliceBarProgressed(v);});
  ui->main_openGL_widget->update();
}

void MainWindow::performSubdivision()
{
  subdivider.execute(**m_mesh, df_subdivision_count);
  ui->main_openGL_widget->update();
}

void MainWindow::setDFSubdivisionCount(int n)
{
  df_subdivision_count = n;
}

void MainWindow::setDFIterationCount(int n)
{
  df_iteration_count = n;
}

void MainWindow::changeAlgorithm(int n)
{
  common::settings::selected_alg = common::settings::Algorithm(n);
}

void MainWindow::changeCurvatureType(int n)
{
  common::settings::selected_curvature = common::settings::CurvatureType(n);
}

void MainWindow::changeMeshSlot(bool index)
{
  if (index) {
    m_mesh = &m_mesh_a;
  }
  else {
    m_mesh = &m_mesh_b;
  }

  ui->main_openGL_widget->setPrintable(*m_mesh);
  ui->main_openGL_widget->update();
}

void MainWindow::setLogAestheticAlpha(double v)
{
  common::settings::log_aesthetic_alpha = v;
}

void MainWindow::changeLoadedModel(int index)
{
  *m_mesh = object_loader.loadFromFile(index);
  ui->main_openGL_widget->setPrintable(*m_mesh);
  ui->main_openGL_widget->update();
}

void MainWindow::cbSliceBarProgressed(int val)
{
  (MainWindow::mainProgressBar)->setValue(val);
}

void MainWindow::setShowVertexIds(int status)
{
  common::settings::show_vertex_ids = static_cast<bool>(status);
  ui->main_openGL_widget->update();
}

void MainWindow::exportModel()
{
  if(!*m_mesh) {
    return;
  }

  QString fileName = ui->exportModelNameLineEdit->text();

  if(fileName.isEmpty()) {
    fileName = "result_mesh";
  }

  OpenMesh::IO::write_mesh(**m_mesh, "./out/" + fileName.toStdString() + ".obj");

  
}

void MainWindow::generateLogAesthetic()
{
  //core::ParametricLogAesthetic1 surf;
  core::LogAestheticSpline surf;
  //m_log_gen = std::make_shared<common::BaseMesh>(surf.tessellate(50));
  m_log_gen = std::make_shared<common::BaseMesh>(surf.execute());
  ui->openGLWidget->setPrintable(m_log_gen);
  ui->openGLWidget->update();
}

