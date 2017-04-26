#ifndef MAIN_WINDOW_H
#define MAIN_WINDOW_H

#include <QtGui>
#include "ui_main_window.h"

/**
 * Main window class for the application
 */
class MainWindow : public QMainWindow
{
  Q_OBJECT

  public:
    MainWindow(QMainWindow* parent = 0);

  private:
    Ui::MainWindow ui;
};

#endif
