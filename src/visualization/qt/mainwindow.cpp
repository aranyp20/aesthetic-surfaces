#include "mainwindow.h"
#include "common_defines.h"
#include "settings.h"
#include "ui_mainwindow.h"

#include "canvas.h"
#include <QtCore/qobject.h>
#include <QtCore/qstring.h>
#include <QtCore/qvariant.h>
#include <QtWidgets/qpushbutton.h>
#include <QtWidgets/qradiobutton.h>
#include <memory>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
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


    ui->logAlphaSpinBox->setMinimum(-5);
    ui->logAlphaSpinBox->setMaximum(5);
    ui->logAlphaSpinBox->setSingleStep(0.5);
    ui->logAlphaSpinBox->setValue(common::settings::log_aesthetic_alpha);

    const auto file_options = object_loader.loadFileOptions();

    size_t tmpi = 0;
    for (const auto& file_option : file_options) {
      ui->modelSelectionComboBox->addItem(QString::fromStdString(file_option), QVariant((int)tmpi++));
    }

    ui->algSelectorComboBox->addItem(QString::fromStdString(common::settings::alg_name.at(common::settings::BASIC)), QVariant((int)common::settings::BASIC));
    ui->algSelectorComboBox->addItem(QString::fromStdString(common::settings::alg_name.at(common::settings::BEZIER)), QVariant((int)common::settings::BEZIER));
    ui->algSelectorComboBox->addItem(QString::fromStdString(common::settings::alg_name.at(common::settings::LOG_AESTHETIC)), QVariant((int)common::settings::LOG_AESTHETIC));

    QObject::connect(ui->modelSelectionComboBox, qOverload<int>(&QComboBox::currentIndexChanged), this, &MainWindow::changeLoadedModel);
    QObject::connect(ui->algSelectorComboBox, qOverload<int>(&QComboBox::currentIndexChanged), this, &MainWindow::changeAlgorithm);

    
    QObject::connect(ui->radioButton, qOverload<bool>(&QRadioButton::toggled), this, &MainWindow::changeMeshSlot);
    
    changeAlgorithm(0);


    ui->radioButton->setChecked(true);

}

MainWindow::~MainWindow()
{
    delete ui;
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

void MainWindow::changeLoadedModel(int index)
{
  *m_mesh = object_loader.loadFromFile(index);
  *m_loaded_mesh = std::make_shared<common::MyMesh>(**m_mesh);
  ui->main_openGL_widget->setPrintable(*m_mesh);
  ui->main_openGL_widget->update();
}

void MainWindow::resetModel()
{
  if(!m_loaded_mesh) {
    std::cout << "No loaded mesh."<< std::endl;
    return;
  }
  *m_mesh = std::make_shared<common::MyMesh>(**m_loaded_mesh);
  ui->main_openGL_widget->setPrintable(*m_mesh);
  ui->main_openGL_widget->update();
}

void MainWindow::performMethod()
{
  discrete_fairer.execute(**m_mesh, df_iteration_count);
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

void MainWindow::changeMeshSlot(bool index)
{
  if (index) {
    m_mesh = &m_mesh_a;
    m_loaded_mesh = &m_loaded_mesh_a;
  }
  else {
    m_mesh = &m_mesh_b;
    m_loaded_mesh = &m_loaded_mesh_b;
  }

  ui->main_openGL_widget->setPrintable(*m_mesh);
  ui->main_openGL_widget->update();
}
