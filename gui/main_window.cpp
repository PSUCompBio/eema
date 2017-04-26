#include "main_window.h"
#include "moc_main_window.cpp"

/**
 * Constructor
 */
MainWindow::MainWindow(QMainWindow* parent)
  : QMainWindow(parent)
{
  this->ui.setupUi(this);
}
