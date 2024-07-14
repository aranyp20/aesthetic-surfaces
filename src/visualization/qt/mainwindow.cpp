#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "canvas.h"

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



    const auto file_options = object_loader.loadFileOptions();

    size_t tmpi = 0;
    for (const auto& file_option : file_options) {
      ui->modelSelectionComboBox->addItem(QString::fromStdString(file_option), QVariant((int)tmpi++));
    }
    

    QObject::connect(ui->modelSelectionComboBox, qOverload<int>(&QComboBox::currentIndexChanged), this, &MainWindow::changeLoadedModel);


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
  m_mesh = object_loader.loadFromFile(index);
  ui->main_openGL_widget->setPrintable(m_mesh);
  ui->main_openGL_widget->update();
}

void MainWindow::performMethod()
{
  discrete_fairer.execute(*m_mesh, df_subdivision_count, df_iteration_count);
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
